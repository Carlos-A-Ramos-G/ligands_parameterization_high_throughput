#!/usr/bin/env python3
"""
Ligand parameterization runner.

Commands
--------
  python run_param.py run     [--ligands L ...] [--jobs N] [--force]
  python run_param.py setup   [--ligands L ...] [--force] [--submit]
  python run_param.py submit  [--ligands L ...]
  python run_param.py status  [--ligands L ...]

  run     Run parameterization locally without SLURM.  Each ligand runs
          antechamber -> parmchk2 -> tleap directly via subprocess.
          --jobs N  runs N ligands in parallel (default: 1).
          Requires antechamber, parmchk2, and tleap in PATH.
  setup   Write tleap input files and one SLURM job script per ligand.
          Pass --submit to also launch the jobs immediately after generation.
  submit  Submit existing SLURM scripts without regenerating any files.
  status  Show per-step timing and SLURM state for each ligand.

Each ligand is processed in a separate CPU-only SLURM job (jobs run in parallel):

    antechamber -> parmchk2 -> tleap

The residue name is read directly from the RESNAME column of each ligand PDB —
no hardcoded name is needed. This makes the tool universal for docking outputs
(UNK), co-crystallized ligands (their own 3-letter code), or anything else.

Outputs per ligand (written to output_dir/<ligand>/):
    <ligand>.mol2    – GAFF atom types + AM1-BCC charges
    <ligand>.frcmod  – missing GAFF parameters
    <ligand>.lib     – AMBER library file

These three files are what the MMGBSA protocol needs to build the complex topology.

Charge detection requires RDKit (primary).
Install: conda install -c conda-forge rdkit
"""

import argparse
import concurrent.futures
import logging
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

import yaml

log = logging.getLogger(__name__)


# =============================================================================
# _tick timing helper — embedded verbatim in every SLURM script
# =============================================================================

_TICK_FN = r"""_tick() {
    local label=$1
    local step_sec=$(( SECONDS - STEP_START ))
    local total_sec=$(( SECONDS - JOB_START ))
    printf "%s\t%d\t%d\t%s\n" \
        "$label" "$step_sec" "$total_sec" "$(date '+%Y-%m-%dT%H:%M:%S')" >> "$TIMING_LOG"
    printf "[%s] %-20s  step %02d:%02d:%02d  total %02d:%02d:%02d\n" \
        "$(date '+%H:%M:%S')" "$label" \
        $(( step_sec/3600 )) $(( step_sec%3600/60 )) $(( step_sec%60 )) \
        $(( total_sec/3600 )) $(( total_sec%3600/60 )) $(( total_sec%60 ))
    STEP_START=$SECONDS
}"""


# =============================================================================
# SLURM script helpers
# =============================================================================

def _sbatch_header(job_name: str, log_prefix: Path, resources: dict) -> str:
    lines = [
        "#!/bin/bash",
        f"#SBATCH --job-name={job_name}",
        f"#SBATCH --output={log_prefix}.out",
        f"#SBATCH --error={log_prefix}.err",
    ]
    for key, value in resources.items():
        lines.append(f"#SBATCH --{key}={value}")
    return "\n".join(lines)


def _module_block(cfg: dict) -> str:
    lines = [""]
    cpu_module = cfg.get("amber", {}).get("cpu_module", "").strip()
    if cpu_module:
        lines.append(f"module load {cpu_module}")
    lines.append("")
    return "\n".join(lines)


def _write_exe(path: Path, content: str) -> None:
    path.write_text(content)
    path.chmod(0o755)


# =============================================================================
# Utilities
# =============================================================================

def _setup_logging(base: Path) -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler(base / "param_run.log"),
        ],
    )


def _load_config(path: str) -> dict:
    with open(path) as fh:
        return yaml.safe_load(fh)


def _discover_ligands(ligands_dir: Path) -> list:
    pdbs = sorted(ligands_dir.glob("*.pdb"), key=lambda p: p.stem)
    if not pdbs:
        raise FileNotFoundError(f"No .pdb files found in {ligands_dir}")
    return [p.stem for p in pdbs]


def _tail(path: Path, n: int = 30) -> str:
    lines = path.read_text().splitlines()
    return "\n".join(lines[-n:])


# =============================================================================
# PDB residue name detection
# =============================================================================

_SOLVENT_RESNAMES = {"HOH", "WAT", "SOL", "TIP", "T3P", "Na+", "Cl-", "NA", "CL"}


def _read_resname(pdb_path: Path) -> str:
    """Read the residue name from the first non-solvent HETATM or ATOM record.

    PDB format columns (1-indexed): 18-20 = residue name.
    Python slice (0-indexed): line[17:20].
    """
    with open(pdb_path) as fh:
        for line in fh:
            if line.startswith(("HETATM", "ATOM  ")):
                resname = line[17:20].strip()
                if resname and resname not in _SOLVENT_RESNAMES:
                    return resname
    raise RuntimeError(
        f"No valid residue name found in {pdb_path}. "
        "Ensure the file contains HETATM or ATOM records with a 3-letter residue name."
    )


# =============================================================================
# Charge detection
# =============================================================================

def _infer_net_charge(pdb_path: Path) -> int:
    """Infer net formal charge from a ligand PDB using RDKit.

    Tries candidate total charges until DetermineBonds produces a self-consistent
    assignment. Raises RuntimeError if none converge (bad geometry or missing H).
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import rdDetermineBonds
    except ImportError:
        raise RuntimeError(
            "RDKit is required for charge detection. "
            "Install with: conda install -c conda-forge rdkit"
        )

    mol = Chem.MolFromPDBFile(str(pdb_path), removeHs=False, sanitize=False)
    if mol is None:
        raise RuntimeError(f"RDKit could not parse {pdb_path}")

    for charge in [0, -1, 1, -2, 2, -3, 3]:
        try:
            mol_try = Chem.RWMol(Chem.Mol(mol))
            rdDetermineBonds.DetermineBonds(mol_try, charge=charge)
            Chem.SanitizeMol(mol_try)
            if Chem.GetFormalCharge(mol_try) == charge:
                return charge
        except Exception:
            continue

    raise RuntimeError(
        f"Could not determine net charge for {pdb_path.name}. "
        "Check that the PDB contains all hydrogens and has correct geometry, "
        "or set net_charge to an explicit integer in config.yaml."
    )


def _resolve_charge(pdb_path: Path, lcfg: dict) -> int:
    raw = lcfg.get("net_charge", "auto")
    if str(raw).lower() == "auto":
        return _infer_net_charge(pdb_path)
    return int(raw)


# =============================================================================
# Input file writers
# =============================================================================

def _write_tleap_params(lig: str, cfg: dict, work: Path, resname: str) -> None:
    ff     = cfg.get("forcefield", {})
    mol2   = work / f"{lig}.mol2"
    frcmod = work / f"{lig}.frcmod"
    lib    = work / f"{lig}.lib"
    (work / "tleap_params.in").write_text(
        f"source {ff.get('ligand', 'leaprc.gaff')}\n"
        f"{resname} = loadmol2 {mol2.resolve()}\n"
        f"check {resname}\n"
        f"loadamberparams {frcmod.resolve()}\n"
        f"saveoff {resname} {lib.resolve()}\n"
        "quit\n"
    )


# =============================================================================
# SLURM script generator
# =============================================================================

def _gen_param_script(lig: str, cfg: dict, base: Path,
                      charge: int, resname: str) -> str:
    """Generate the full SLURM bash script for one ligand (CPU-only)."""
    lcfg   = cfg.get("ligand", {})
    mult   = lcfg.get("multiplicity", 1)
    method = lcfg.get("charge_method", "bcc")

    ligands_dir = (base / cfg.get("ligands_dir", "ligands")).resolve()
    output_dir  = (base / cfg.get("output_dir",  "output")).resolve()
    ligand_pdb  = (ligands_dir / f"{lig}.pdb").resolve()
    work        = (output_dir / lig).resolve()
    logs_dir    = (base / "logs").resolve()

    timing_log = work / "timing.log"
    mol2       = work / f"{lig}.mol2"
    frcmod     = work / f"{lig}.frcmod"
    tleap_in   = work / "tleap_params.in"

    p = str

    lines = [
        _sbatch_header(
            f"{lig}_param",
            logs_dir / lig,
            cfg["slurm"]["cpu"],
        ),
        _module_block(cfg),
        "set -euo pipefail",
        "",
        "# ---- Timing",
        f"TIMING_LOG={p(timing_log)}",
        "JOB_START=$SECONDS",
        "STEP_START=$SECONDS",
        "",
        _TICK_FN,
        "",
        "_tick start",
        "",
        f"cd {p(work)}",
        "",
        "# ---- Step 1: antechamber  (GAFF atom types + AM1-BCC charges)",
        f"antechamber -i {p(ligand_pdb)} -fi pdb \\",
        f"    -o {p(mol2)} -fo mol2 \\",
        f"    -c {method} -s 2 -nc {charge} -m {mult} -rn {resname}",
        "_tick antechamber",
        "",
        "# ---- Step 2: parmchk2  (check for missing GAFF parameters)",
        f"parmchk2 -i {p(mol2)} -f mol2 -o {p(frcmod)}",
        "_tick parmchk2",
        "",
        "# ---- Step 3: tleap  (generate AMBER library file)",
        f"tleap -f {p(tleap_in)}",
        "_tick tleap",
        "",
    ]
    return "\n".join(lines)


# =============================================================================
# Status command
# =============================================================================

_STEP_ORDER  = ["start", "antechamber", "parmchk2", "tleap"]
_STEP_LABELS = {
    "start":       "job start",
    "antechamber": "antechamber",
    "parmchk2":    "parmchk2",
    "tleap":       "tleap",
}


def _hms(secs: int) -> str:
    return f"{secs // 3600:02d}:{secs % 3600 // 60:02d}:{secs % 60:02d}"


def status(cfg: dict, base: Path, ligs: list) -> None:
    running = {}
    try:
        sq = subprocess.run(
            ["squeue", "--me", "-o", "%j %T", "--noheader"],
            capture_output=True, text=True, check=False,
        )
        for line in sq.stdout.splitlines():
            parts = line.split()
            if len(parts) >= 2:
                running[parts[0]] = parts[1]
    except FileNotFoundError:
        pass  # squeue not available on this machine

    output_dir = base / cfg.get("output_dir", "output")

    for lig in ligs:
        slurm_state = running.get(f"{lig}_param", "not submitted")
        work        = output_dir / lig
        timing_log  = work / "timing.log"

        outputs = {
            "mol2":   (work / f"{lig}.mol2").exists(),
            "frcmod": (work / f"{lig}.frcmod").exists(),
            "lib":    (work / f"{lig}.lib").exists(),
        }
        done_mark = " [DONE]" if all(outputs.values()) else ""
        print(f"\n{lig}  [{slurm_state}]{done_mark}")
        print(f"  mol2={outputs['mol2']}  frcmod={outputs['frcmod']}  lib={outputs['lib']}")

        if not timing_log.exists():
            print("  (no timing log yet)")
            continue

        entries = {}
        for line in timing_log.read_text().splitlines():
            parts = line.strip().split("\t")
            if len(parts) == 4:
                label, step_s, total_s, ts = parts
                entries[label] = (int(step_s), int(total_s), ts)

        print(f"  {'Step':<20} {'Step time':>10} {'Cumulative':>10}  Timestamp")
        print(f"  {'-'*20} {'-'*10} {'-'*10}  {'-'*19}")
        for key in _STEP_ORDER:
            if key not in entries:
                continue
            step_sec, total_sec, ts = entries[key]
            label = _STEP_LABELS.get(key, key)
            if key == "start":
                print(f"  {label:<20} {'--':>10} {'--':>10}  {ts}")
            else:
                print(f"  {label:<20} {_hms(step_sec):>10} {_hms(total_sec):>10}  {ts}")


# =============================================================================
# Setup command
# =============================================================================

def setup(cfg: dict, base: Path, ligs: list,
          force: bool = False, submit_after: bool = False) -> None:
    ligands_dir = base / cfg.get("ligands_dir", "ligands")
    output_dir  = base / cfg.get("output_dir",  "output")
    logs_dir    = base / "logs"
    scripts_dir = base / "scripts"

    for d in (output_dir, logs_dir, scripts_dir):
        d.mkdir(parents=True, exist_ok=True)

    print(f"Base dir   : {base}")
    print(f"Ligands    : {len(ligs)} -- {ligs}")
    print(f"Output dir : {output_dir}")
    print()

    failed = []
    for lig in ligs:
        script = scripts_dir / f"{lig}_param.sh"
        if script.exists() and not force:
            log.info("[%s] script exists -- skipping (use --force to regenerate)", lig)
            continue

        pdb_src = ligands_dir / f"{lig}.pdb"
        if not pdb_src.exists():
            log.error("[%s] ligand PDB not found: %s", lig, pdb_src)
            failed.append(lig)
            continue

        try:
            resname = _read_resname(pdb_src)
            log.info("[%s] residue name = %s", lig, resname)

            lcfg   = cfg.get("ligand", {})
            charge = _resolve_charge(pdb_src, lcfg)
            log.info("[%s] net charge = %+d", lig, charge)

            work = output_dir / lig
            work.mkdir(parents=True, exist_ok=True)

            _write_tleap_params(lig, cfg, work, resname)
            _write_exe(script, _gen_param_script(lig, cfg, base, charge, resname))
            log.info("[%s] script written: %s", lig, script)

        except Exception as exc:
            log.error("[%s] setup failed: %s", lig, exc)
            failed.append(lig)

    if failed:
        log.error("Setup failed for: %s", failed)
        sys.exit(1)

    print("\nSetup complete.")
    if submit_after:
        submit(cfg, base, ligs)
    else:
        print("Run 'python run_param.py submit' to launch the jobs.")


# =============================================================================
# Submit command
# =============================================================================

def submit(cfg: dict, base: Path, ligs: list) -> None:
    scripts_dir = base / "scripts"
    for lig in ligs:
        script = scripts_dir / f"{lig}_param.sh"
        if not script.exists():
            sys.exit(
                f"Script not found: {script}\n"
                "Run 'python run_param.py setup' first."
            )
        result = subprocess.run(
            ["sbatch", str(script)], capture_output=True, text=True
        )
        if result.returncode == 0:
            print(f"  [{lig}] {result.stdout.strip()}")
        else:
            print(f"  [{lig}] sbatch failed: {result.stderr.strip()}",
                  file=sys.stderr)


# =============================================================================
# Local run helpers
# =============================================================================

def _run_single(lig: str, cfg: dict, base: Path,
                charge: int, resname: str, force: bool) -> tuple:
    """Run antechamber -> parmchk2 -> tleap locally for one ligand.

    Returns (success: bool, message: str).
    stdout/stderr are captured to logs/<lig>_local.log.
    """
    lcfg   = cfg.get("ligand", {})
    mult   = lcfg.get("multiplicity", 1)
    method = lcfg.get("charge_method", "bcc")

    ligands_dir = (base / cfg.get("ligands_dir", "ligands")).resolve()
    output_dir  = (base / cfg.get("output_dir",  "output")).resolve()
    logs_dir    = (base / "logs").resolve()

    ligand_pdb = ligands_dir / f"{lig}.pdb"
    work       = output_dir / lig
    mol2       = work / f"{lig}.mol2"
    frcmod     = work / f"{lig}.frcmod"
    lib        = work / f"{lig}.lib"
    timing_log = work / "timing.log"
    local_log  = logs_dir / f"{lig}_local.log"

    if not force and mol2.exists() and frcmod.exists() and lib.exists():
        return True, f"[{lig}] already done -- skipping (use --force to redo)"

    work.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)
    timing_log.unlink(missing_ok=True)

    _write_tleap_params(lig, cfg, work, resname)

    t = [time.monotonic(), time.monotonic()]  # [job_start, step_start]

    def tick(label: str) -> None:
        now       = time.monotonic()
        step_sec  = int(now - t[1])
        total_sec = int(now - t[0])
        ts = datetime.now().strftime("%Y-%m-%dT%H:%M:%S")
        with open(timing_log, "a") as fh:
            fh.write(f"{label}\t{step_sec}\t{total_sec}\t{ts}\n")
        t[1] = now

    tick("start")

    steps = [
        ("antechamber", [
            "antechamber",
            "-i", str(ligand_pdb), "-fi", "pdb",
            "-o", str(mol2), "-fo", "mol2",
            "-c", method, "-s", "2",
            "-nc", str(charge), "-m", str(mult), "-rn", resname,
        ]),
        ("parmchk2", [
            "parmchk2",
            "-i", str(mol2), "-f", "mol2",
            "-o", str(frcmod),
        ]),
        ("tleap", [
            "tleap", "-f", str(work / "tleap_params.in"),
        ]),
    ]

    with open(local_log, "w") as log_fh:
        for step_name, cmd in steps:
            log.info("[%s] running %s ...", lig, step_name)
            result = subprocess.run(
                cmd, cwd=str(work),
                stdout=log_fh, stderr=subprocess.STDOUT,
            )
            tick(step_name)
            if result.returncode != 0:
                tail = _tail(local_log)
                return False, (
                    f"[{lig}] {step_name} failed (exit {result.returncode}).\n"
                    f"Log: {local_log}\n"
                    f"--- last lines ---\n{tail}"
                )

    return True, f"[{lig}] done  (log: {local_log})"


# =============================================================================
# Local run command
# =============================================================================

def run(cfg: dict, base: Path, ligs: list,
        jobs: int = 1, force: bool = False) -> None:
    ligands_dir = base / cfg.get("ligands_dir", "ligands")

    print(f"Base dir   : {base}")
    print(f"Ligands    : {len(ligs)} -- {ligs}")
    print(f"Parallel   : {jobs} job(s)")
    print()

    # Pre-flight: resolve resname and charge for all ligands before starting any job.
    tasks = []
    preflight_failed = []
    for lig in ligs:
        pdb_src = ligands_dir / f"{lig}.pdb"
        if not pdb_src.exists():
            log.error("[%s] ligand PDB not found: %s", lig, pdb_src)
            preflight_failed.append(lig)
            continue
        try:
            resname = _read_resname(pdb_src)
            charge  = _resolve_charge(pdb_src, cfg.get("ligand", {}))
            log.info("[%s] residue name = %s, net charge = %+d", lig, resname, charge)
            tasks.append((lig, charge, resname))
        except Exception as exc:
            log.error("[%s] pre-flight failed: %s", lig, exc)
            preflight_failed.append(lig)

    if preflight_failed:
        log.error("Pre-flight failed for: %s", preflight_failed)
        sys.exit(1)

    run_failed = []

    if jobs == 1:
        for lig, charge, resname in tasks:
            ok, msg = _run_single(lig, cfg, base, charge, resname, force)
            print(msg)
            if not ok:
                run_failed.append(lig)
    else:
        def _worker(task):
            return _run_single(task[0], cfg, base, task[1], task[2], force)

        with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as pool:
            futures = {pool.submit(_worker, t): t[0] for t in tasks}
            for fut in concurrent.futures.as_completed(futures):
                ok, msg = fut.result()
                print(msg)
                if not ok:
                    run_failed.append(futures[fut])

    if run_failed:
        log.error("Parameterization failed for: %s", run_failed)
        sys.exit(1)

    output_dir = base / cfg.get("output_dir", "output")
    print(f"\nAll ligands parameterized successfully.")
    print(f"Outputs in: {output_dir}")


# =============================================================================
# Entry point
# =============================================================================

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Ligand parameterization runner",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--config", default="config.yaml",
        help="YAML config file (default: config.yaml)",
    )

    sub = parser.add_subparsers(dest="command", required=True)

    p_run = sub.add_parser("run", help="Run parameterization locally (no SLURM)")
    p_run.add_argument("--ligands", nargs="+",
        help="Ligand IDs to process (default: all PDBs in ligands_dir)")
    p_run.add_argument("--jobs", type=int, default=1, metavar="N",
        help="Number of ligands to run in parallel (default: 1)")
    p_run.add_argument("--force", action="store_true",
        help="Redo parameterization even if outputs already exist")

    p_setup = sub.add_parser("setup", help="Write input files and SLURM scripts")
    p_setup.add_argument("--ligands", nargs="+",
        help="Ligand IDs to process (default: all PDBs in ligands_dir)")
    p_setup.add_argument("--force", action="store_true",
        help="Overwrite existing scripts and input files")
    p_setup.add_argument("--submit", action="store_true",
        help="Submit SLURM jobs immediately after setup")

    p_submit = sub.add_parser("submit", help="Submit existing SLURM scripts")
    p_submit.add_argument("--ligands", nargs="+",
        help="Ligand IDs to submit (default: all PDBs in ligands_dir)")

    p_status = sub.add_parser("status", help="Print step timings and SLURM state")
    p_status.add_argument("--ligands", nargs="+",
        help="Ligand IDs to query (default: all PDBs in ligands_dir)")

    args = parser.parse_args()

    config_path = Path(args.config)
    if not config_path.exists():
        sys.exit(f"Config file not found: {config_path}")

    cfg  = _load_config(str(config_path))
    base = Path(cfg.get("base_dir", ".")).resolve()
    _setup_logging(base)

    ligands_dir = base / cfg.get("ligands_dir", "ligands")
    ligs = args.ligands if args.ligands else _discover_ligands(ligands_dir)

    if args.command == "run":
        run(cfg, base, ligs, jobs=args.jobs, force=args.force)
    elif args.command == "setup":
        setup(cfg, base, ligs, force=args.force, submit_after=args.submit)
    elif args.command == "submit":
        submit(cfg, base, ligs)
    elif args.command == "status":
        status(cfg, base, ligs)


if __name__ == "__main__":
    main()
