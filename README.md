Ligand Parameterization
=======================

High-throughput AMBER/GAFF parameterization for organic protein-ligand binders.
Generates .mol2/.frcmod/.lib files for entire compound libraries — in parallel
locally or via SLURM on HPC clusters.

This protocol describes the parameterization of a ligand set used to study 
the binding free energy of these compounds to the caspase-1 enzyme, a key 
inflammatory target implicated in cellular inflammation-related diseases such 
as rheumatoid arthritis, type 2 diabetes, atherosclerosis and Alzheimer's.
Accurate parameterization is essential to ensure reliable molecular dynamics 
simulations and free energy calculations, enabling a rigorous characterization 
of ligand–enzyme interactions.

The complete study can be found here:
https://www.sciencedirect.com/science/article/abs/pii/S0223523419301886

PREREQUISITES
-------------

  Tool              Purpose
  ----------------  ----------------------------------------
  AmberTools >= 22  antechamber, parmchk2, tleap
  RDKit             Automatic net charge detection
  PyYAML            Config file parsing
  SLURM             Job scheduler (cluster only)

Install Python dependencies using the provided requirements.txt:

    conda install -c conda-forge --file requirements.txt

  Or with pip:

    pip install -r requirements.txt

  Note: if you hit binary compatibility issues with rdkit (e.g. PythonB
  version mismatches), prefer the conda-forge install over pip.


WORKFLOW
--------

  ligands/<lig>.pdb
      |
      +-- residue name  ->  read from RESNAME column (automatic)
      +-- net charge    ->  RDKit auto-detect (or set in config.yaml)
      |
      antechamber       ->  output/<lig>/<lig>.mol2
      |
      parmchk2          ->  output/<lig>/<lig>.frcmod
      |
      tleap             ->  output/<lig>/<lig>.lib

Each ligand runs as an independent CPU-only SLURM job (all in parallel),
or locally via subprocess with optional parallelism (--jobs N).


QUICK START -- LOCAL MACHINE
-----------------------------

1. Drop ligand PDB files into ligands/

       ligands/
           compoundA.pdb
           compoundB.pdb
           compoundC.pdb

   PDB files must contain all hydrogens and correct 3D geometry.
   The residue name is read from the RESNAME column (cols 18-20)
   automatically -- no configuration needed.

2. Activate your conda environment

       conda activate AmberTools25   # or your env name

3. Run

       # One ligand at a time (safe on a laptop)
       python run_param.py run

       # 4 ligands in parallel (good for a workstation)
       python run_param.py run --jobs 4

       # Only specific ligands
       python run_param.py run --ligands compoundA compoundB

       # Redo already-completed ligands
       python run_param.py run --force

   Progress is printed per ligand. Full stdout/stderr for each tool is
   saved to logs/<lig>_local.log -- check there if something fails.

4. Check results

       python run_param.py status


QUICK START -- HPC CLUSTER (SLURM)
------------------------------------

1. Drop ligand PDB files into ligands/ (same as above)

2. Edit config.yaml -- minimally required:
     - slurm.cpu.partition  : your cluster CPU partition name
     - amber.cpu_module     : AmberTools module (or leave empty for conda)
     - forcefield.ligand    : must match the force field used in your MD protocol

3. Generate SLURM scripts

       python run_param.py setup

   Reads all .pdb files in ligands/, detects residue names and charges,
   and writes:
     - output/<lig>/tleap_params.in
     - scripts/<lig>_param.sh

4. Submit jobs

       python run_param.py submit

   Or combine steps 3 and 4:

       python run_param.py setup --submit

5. Check progress

       python run_param.py status

   Example output:

       compoundA  [RUNNING]
         mol2=True  frcmod=False  lib=False
         Step                 Step time  Cumulative  Timestamp
         -------------------  ---------  ----------  -------------------
         job start                   --          --  2025-05-20T10:01:03
         antechamber           00:34:21    00:34:21  2025-05-20T10:35:24

       compoundB  [PENDING]
         mol2=False  frcmod=False  lib=False
         (no timing log yet)

       compoundC  [DONE]
         mol2=True  frcmod=True  lib=True


CLI REFERENCE
-------------

  python run_param.py [--config CONFIG] <command> [options]

  Commands:
    run     [--ligands L ...]  [--jobs N]  [--force]
    setup   [--ligands L ...]  [--force]   [--submit]
    submit  [--ligands L ...]
    status  [--ligands L ...]

  Flag              Applies to     Description
  ----------------  -------------  ----------------------------------------
  --ligands L ...   all            Ligand IDs to process (default: all PDBs)
  --jobs N          run            Ligands to run in parallel (default: 1)
  --force           run, setup     Redo even if outputs/scripts exist
  --submit          setup          Submit SLURM jobs right after setup


OUTPUT STRUCTURE
----------------

  ligand_parameterization/
      ligands/                   <- place input .pdb files here
      output/
          compoundA/
              compoundA.mol2     GAFF atom types + AM1-BCC charges
              compoundA.frcmod   missing GAFF parameters
              compoundA.lib      AMBER library file
              tleap_params.in
              timing.log
          compoundB/
              ...
      scripts/                   generated SLURM job scripts
      logs/                      SLURM .out/.err  or  local run logs


USING THE OUTPUTS
-----------------

Each output/<lig>/ directory contains the two files that tleap needs
to build any AMBER topology:

    loadamberparams  /path/to/output/<lig>/<lig>.frcmod
    loadoff          /path/to/output/<lig>/<lig>.lib

Point your downstream protocol (MMGBSA, FEP, plain MD, ...) at the
output/ folder and load these files before assembling the complex.


CONFIGURATION REFERENCE
------------------------

  Key                   Default         Description
  --------------------  --------------  ------------------------------------
  base_dir              .               Root of the parameterization repo
  ligands_dir           ligands         Folder with input .pdb files
  output_dir            output          Folder for parameterization outputs
  amber.cpu_module      ""              HPC module (e.g. ambertools/24)
  slurm.cpu.*           see config      Any #SBATCH option, emitted verbatim
  ligand.net_charge     auto            auto (RDKit) or explicit integer
  ligand.multiplicity   1               Spin multiplicity
  ligand.charge_method  bcc             bcc, resp, or gas
  forcefield.ligand     leaprc.gaff     Must match your downstream MD protocol


NOTES
-----

  GAFF vs GAFF2
    leaprc.gaff is the default. If you switch to leaprc.gaff2 here,
    use the same force field in your downstream MD protocol as well.
    Mixing GAFF and GAFF2 parameters in the same system is incorrect.

  Residue name
    Automatically read from the RESNAME column of each PDB. Works with
    docking outputs (UNK), co-crystallized ligands, or any naming convention.

  AM1-BCC wall time
    Budget 30-90 minutes per ligand depending on molecular size. Use at
    least time: "04:00:00" for typical drug-like molecules, "24:00:00"
    for larger peptide-like compounds.

  antechamber intermediates
    ANTECHAMBER_*, ATOMTYPE.INF, sqm.* are written to output/<lig>/ and
    can be deleted after the job succeeds (.gitignore excludes them).
