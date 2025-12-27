# Running Molecular Dynamics

StreaMD runs preparation, equilibration, production, continuation, and analysis. Activate your environment first (`conda activate md`).

## Quick Start Examples

### CPU: Protein in Water (end-to-end)
```bash
run_md -p protein_H_HIS.pdb --md_time 1 --nvt_time 1000 --npt_time 1000

# Extend to 2 ns using the same working directory
run_md --wdir_to_continue md_files/md_run/protein_H_HIS_ligand1_replica1 --md_time 2
```
Outputs: see {doc}`outputs` for the `md_files/` structure with `md_preparation/` and `md_run/<system>/` plus `md_analysis/` summaries.

### GPU: Protein-Ligand with Continuation
```bash
# Initial 1 ns run on a GPU
run_md -p protein_H_HIS.pdb -l ligand.mol --md_time 1 --device gpu
```
Outputs: `md_files/md_run/protein_H_HIS_ligand1_replica1/` contains GROMACS outputs and `md_analysis/`; see {doc}`outputs` for file details.

## Command-Line Reference (`run_md -h`)
```
run_md -h

usage: run_md [-h] [--config FILENAME] [-p FILENAME] [-d WDIR] [-l FILENAME] [--cofactor FILENAME] [--clean_previous_md] [--hostfile FILENAME] [-c INTEGER]
              [--mdrun_per_node INTEGER] [--device cpu] [--gpu_ids GPU ID [GPU ID ...]] [--ntmpi_per_gpu int] [--topol topol.top]
              [--topol_itp topol_chainA.itp topol_chainB.itp [topol_chainA.itp topol_chainB.itp ...]] [--posre posre.itp [posre.itp ...]]
              [--protein_forcefield amber99sb-ildn] [--noignh] [--md_time ns] [--npt_time ps] [--nvt_time ps] [--seed int] [--replicas INTEGER] [--no_dr]
              [--not_clean_backup_files] [--steps [STEPS ...]] [--mdp_dir Path to a directory with specific MDP files] [--save_traj_without_water]
              [--wdir_to_continue DIRNAME [DIRNAME ...]] [-o OUT_SUFFIX] [--deffnm prefix for MD files] [--tpr FILENAME] [--cpt FILENAME] [--xtc FILENAME]
              [--ligand_list_file all_ligand_resid.txt] [--ligand_id UNL] [--activate_gaussian module load Gaussian/09-d01]
              [--gaussian_exe g09 or /apps/all/Gaussian/09-d01/g09/g09] [--gaussian_basis B3LYP/6-31G*] [--gaussian_memory 120GB] [--metal_resnames [MN ...]]
              [--metal_cutoff 2.8] [--metal_charges {MN:2, ZN:2, CA:2}]

Run or continue MD simulation. Allowed systems: Protein, Protein-Ligand, Protein-Cofactors (multiple), Protein-Ligand-Cofactors (multiple)

options:
  -h, --help            show this help message and exit
  --config FILENAME     Path to YAML configuration file with default arguments
  -o OUT_SUFFIX, --out_suffix OUT_SUFFIX
                        User-defined unique suffix for output files

Standard Molecular Dynamics Simulation Run:
  -p FILENAME, --protein FILENAME
                        Input file of protein. Supported formats: *.pdb or gro
  -d WDIR, --wdir WDIR  Working directory. If not set, the current directory will be used.
  -l FILENAME, --ligand FILENAME
                        Input file with compound(s). Supported formats: *.mol or sdf or mol2
  --cofactor FILENAME   Input file with compound(s). Supported formats: *.mol or sdf or mol2
  --clean_previous_md   Remove a production MD simulation directory if it exists to re-initialize production MD setup
  --hostfile FILENAME   Text file with addresses of nodes of a Dask SSH cluster. Typically passed as the $PBS_NODEFILE variable from inside a PBS script.
                        The first line in this file will be the address of the scheduler running on the standard port 8786. If omitted, calculations will run on a
                        single machine as usual.
  -c INTEGER, --ncpu INTEGER
                        Number of CPU cores per server. Use all available CPUs by default.
  --mdrun_per_node INTEGER
                        Number of simulations to run per node. If multiple simulations per node are requested, the available CPU cores will be split evenly across them.
                        By default, run only one simulation per node and use all available CPUs.
  --device cpu          Calculate bonded and non-bonded interactions on: auto, cpu, or gpu.
  --gpu_ids GPU ID [GPU ID ...]
                        List of unique GPU device IDs available to use. Specify when using multiple GPUs or targeting exact GPU devices.
  --ntmpi_per_gpu int   The number of thread-MPI ranks to start per GPU. The default, 1, will start one rank per GPU.
  --topol topol.top     Topology file (required if a gro-file is provided for the protein). All output files obtained from gmx2pdb should preserve the original names
  --topol_itp topol_chainA.itp topol_chainB.itp [topol_chainA.itp topol_chainB.itp ...]
                        itp files for individual protein chains (required if a gro-file is provided for the protein). All output files obtained from gmx2pdb should
                        preserve the original names
  --posre posre.itp [posre.itp ...]
                        posre file(s) (required if a gro-file is provided for the protein). All output files obtained from gmx2pdb should preserve the original names
  --protein_forcefield amber99sb-ildn
                        Force field for protein preparation. Available FF can be found at Miniconda3/envs/md/share/gromacs/top.
  --noignh              By default, StreaMD uses gmx pdb2gmx -ignh, which re-adds hydrogens using residue names (the correct protonation states must be provided by
                        the user) and ignores the original hydrogens. If the --noignh argument is used, the original hydrogen atoms will be preserved during the
                        preparation, although there may be problems with recognition of atom names by GROMACS.
  --md_time ns          Time of MD simulation in ns. Default: 1 ns.
  --npt_time ps         Time of NPT equilibration in ps. Default: 1000 ps.
  --nvt_time ps         Time of NVT equilibration in ps. Default: 1000 ps.
  --seed int            Random seed.
  --replicas INTEGER    Number of replicate simulations to run per complex
  --no_dr               Turn off the acdoctor mode and do not check/diagnose problems in the input ligand file in the next attempt if the regular antechamber run for
                        ligand preparation fails (ligand_mol2prep.sh script related issues). Use this argument carefully and ensure that you provide valid structures.
  --not_clean_backup_files
                        Do not remove backups of MD files.
  --steps [STEPS ...]   Run particular step(s) of the StreaMD run. Options: 1 - preparation (protein, ligand, cofactor); 2 - MD equilibration (minimization, NVT,
                        NPT); 3 - MD simulation; 4 - MD analysis. Example: 3 4. If steps 2/3/4 are used, --wdir_to_continue should be used to provide directories
                        with files obtained during step 1.
  --mdp_dir Path to a directory with specific MDP files
                        To use specific MD settings, the user can provide a path to a directory that contains any of the following .mdp files: ions.mdp, minim.mdp,
                        nvt.mdp, npt.mdp, md.mdp. Missing .mdp files will be replaced by default StreaMD files. Provided .mdp files will be used as templates,
                        although the system StreaMD parameters (seed, nvt_time, npt_time, md_time, and tc-grps (cannot be changed by the user)) will override the
                        ones provided. Warning: The names of the files must be strictly preserved.
  --save_traj_without_water
                        Save additional md_out_nowater.tpr and md_fit_nowater.xtc files for more memory-efficient analysis.
  --wdir_to_continue DIRNAME [DIRNAME ...]
                        Single or multiple directories containing simulations created by StreaMD. Use with steps 2, 3, or 4 to continue a run. Directories should
                        contain tpr, cpt, xtc, and optional all_ligand_resid.txt files. If you want to continue your own simulation not created by the tool, use
                        --tpr, --cpt, --xtc, and --wdir arguments (--ligand_list_file is optional and required to run MD analysis after the simulation).

Continue or Extend Molecular Dynamics Simulation:
  --deffnm prefix for MD files
                        Prefix for the MD files. Used to run, extend, or continue the simulation. If --wdir_to_continue is used, files named deffnm.tpr, deffnm.cpt,
                        and deffnm.xtc will be read from those directories.
  --tpr FILENAME        Use an explicit tpr file to continue a non-StreaMD simulation
  --cpt FILENAME        Use an explicit cpt file to continue a non-StreaMD simulation
  --xtc FILENAME        Use an explicit xtc file to continue a non-StreaMD simulation
  --ligand_list_file all_ligand_resid.txt
                        To run automatic MD analysis for ligands after continuing a non-StreaMD simulation, set a ligand_list file. Format (no headers): user_ligand_id
                        gromacs_ligand_id. Example: my_ligand UNL. Can be set via CLI or placed in --wdir_to_continue directories.
  --ligand_id UNL       Set ligand_id (if not default UNL) to run automatic MD analysis for the ligand after continuing a simulation.

Boron-containing molecules or MCPBPY usage (use together with Standard Molecular Dynamics Simulation Run arguments group):
  --activate_gaussian module load Gaussian/09-d01
                        String to load the Gaussian module if necessary.
  --gaussian_exe g09 or /apps/all/Gaussian/09-d01/g09/g09
                        Path to Gaussian executable or alias. Required to run preparation of boron-containing compounds.
  --gaussian_basis B3LYP/6-31G*
                        Gaussian basis.
  --gaussian_memory 120GB
                        Gaussian memory usage.

MCPBPY usage (Use together with Standard Molecular Dynamics Simulation Run and Boron-containing molecules arguments group):
  --metal_resnames [MN ...]
                        Metal residue names to run the MCPB.py procedure. Start MCPB.py procedure only if gaussian_exe and activate_gaussian arguments are set up; otherwise
                        standard gmx2pdb procedure will be run.
  --metal_cutoff 2.8    Metal residue cutoff to run MCPB.py procedure
  --metal_charges {MN:2, ZN:2, CA:2}
                        Metal residue charges in dictionary format. Start the MCPB.py procedure only if metal_resnames and gaussian_exe and activate_gaussian arguments are
                        set up; otherwise standard gmx2pdb procedure will be run.

```

## Standard Workflows

### Protein in Water (use all available CPUs, by default)
```bash
run_md -p protein_H_HIS.pdb --md_time 1 --nvt_time 1000 --npt_time 1000
```

### Protein-Ligand
```bash
# Single ligand
run_md -p protein_H_HIS.pdb -l ligand.mol --md_time 1 --ncpu 128
```

StreaMD can run multiple simulations for a set of ligands bound to the same protein. All ligands should have valid, aligned 3D coordinates.
If any ligand preparation fails, those systems are skipped and only successfully prepared systems continue to the next steps.
```bash
run_md -p protein_H_HIS.pdb -l ligands.sdf
```

### Protein-Cofactor
All cofactor molecules must be present in a simulated system; if any cofactor preparation fails, StreaMD stops and does not continue to the next step.
```bash
run_md -p protein_H_HIS.pdb --cofactor cofactors.sdf --md_time 1
```

### Using a Configuration File
Arguments passed via CLI take precedence over `config.yml`.
```bash
run_md --config config.yml
```
See `configuration.md` for details.

### Boron-Containing Molecules (Gaussian)
See {doc}`advanced_features` for Gaussian setup, options, and an example command.

### Ligand-Binding Metalloproteins with MCPB.py
See {doc}`advanced_features` for MCPB.py requirements and command examples.

## Replicas
Run several independent simulations of one prepared system with `--replicas`. The system is prepared once and copied for each replica under `md_files/md_run/<complex>_replicaN`. Replica seeds increment from the value passed via `--seed`; if `--seed -1` (default) is provided, all replicas keep `-1`.

If a replica directory already exists, StreaMD reuses existing files (with warnings). Running later with a higher `--replicas` count adds only new replicas.
```bash
run_md -p protein_H_HIS.pdb -l ligand.mol --md_time 1 --replicas 3 --seed 1024

# Add 1 extra replica to 3 already existing
run_md -p protein_H_HIS.pdb -l ligand.mol --md_time 1 --replicas 4 --seed 1024

```

## Multi-Server Execution
Provide a hostfile (first line is the Dask scheduler).
```bash
# PBS
run_md -p protein_H_HIS.pdb -l molecules.sdf --cofactor cofactors.sdf --md_time 1 \
  --hostfile "$PBS_NODEFILE" --ncpu 128

# SLURM
srun hostname | sort | uniq > hostfile
run_md -p protein_H_HIS.pdb -l molecules.sdf --cofactor cofactors.sdf --md_time 1 \
  --hostfile hostfile --ncpu 128
```

## Continue or Extend Simulations
StreaMD resumes automatically when checkpoint files exist.
```bash
# Continue or extend a StreaMD run
run_md --wdir_to_continue md_files/md_run/protein_H_HIS_ligand_*/ --md_time 2

# Continue a non-StreaMD run with explicit files
run_md --wdir mdrun --md_time 3 \
  --tpr protein_H_HIS_ligand_1/md_out.tpr \
  --cpt protein_H_HIS_ligand_1/md_out.cpt \
  --xtc protein_H_HIS_ligand_1/md_out.xtc \
  --ligand_list_file protein_H_HIS_ligand_1/all_ligand_resid.txt

# Skip preparation and run specific steps
run_md --wdir_to_continue md_files/md_run/protein_H_HIS_ligand_1/ --md_time 3 --steps 3 4
```

You can rerun the same command to resume an interrupted run; StreaMD detects checkpoint files and continues from the stop point.
```
run_md -p protein_H_HIS.pdb -l ligand.mol --md_time 10
```

## GPU Usage
StreaMD supports GPU acceleration for minimization, NVT/NPT equilibration, and production. Performance depends on hardware and system size; monitor GPU load (e.g., `nvidia-smi`). 
For tuning guidance, see the GROMACS GPU documentation: https://manual.gromacs.org/documentation/current/user-guide/mdrun-performance.html

- Select device: `--device gpu`, `--device cpu`, or `--device auto` (GROMACS decides).
  - With `--device gpu`, StreaMD offloads nb, update, PME, bonded, and PME-FFT calculations (all GPU-capable tasks) to the GPU.
  - Letting GROMACS auto-offload with `--device auto` may be optimal when CPUs are relatively powerful compared to GPUs.

To improve performance, one can use multiple GPUs or start multiple ranks per GPU.
> [!WARNING]
> Increasing the number of GPUs does not always improve performance. Each simulation will use all provided GPUs.

- Increase thread-MPI ranks per GPU: `--ntmpi_per_gpu 2`.
  - The `ntmpi_per_gpu` argument determines the total thread-MPI ranks (`gmx mdrun -ntmpi X`, where **X** = ntmpi_per_gpu * number of GPUs). The default is 1; 2 ranks per GPU may perform better in some cases.
- Use specific multiple GPUs: `--gpu_ids 0 1 2 3`.

For better GPU usage (e.g., multiple small systems), you can start multiple simulations on one or more nodes; the tool automatically splits available CPU cores across simulations:
- Run multiple simulations per node (splits CPUs): `--mdrun_per_node 2`.

> [!WARNING] 
> All simulations will still utilize all provided GPUs, which can lead to suboptimal GPU load. This feature is still under development.
- Multiple tasks on the same node using multiple GPUs (experimental): `--mdrun_per_node 2 --gpu_ids 0 1`.

Examples:
```bash
# Single GPU
run_md -p protein_HIS.pdb -l ligand.mol --md_time 1 --device gpu

# Multiple GPUs per simulation
run_md -p protein_HIS.pdb -l ligand.mol --md_time 1 --device gpu --gpu_ids 0 1 2 3

# Increase thread-MPI ranks per GPU
run_md -p protein_HIS.pdb -l ligand.mol --md_time 1 --device gpu --ntmpi_per_gpu 2

# Multiple runs per node with GPU
run_md -p protein_HIS.pdb -l ligands.sdf --md_time 1 --device gpu --mdrun_per_node 2

# Multiple tasks on same node using multiple GPUs (experimental)
run_md -p protein_HIS.pdb -l ligands.sdf --md_time 1 --device gpu --mdrun_per_node 2 --gpu_ids 0 1

# CPU-only even if GPUs are present
run_md -p protein_HIS.pdb -l ligand.mol --md_time 1 --device cpu

# Let GROMACS auto-offload
run_md -p protein_HIS.pdb -l ligand.mol --md_time 1 --device auto --gpu_ids 0
```

## Custom MDP Files and Step Control
See {doc}`advanced_features` for custom `.mdp` usage, step selection, trajectory size options, and other advanced flags.

## Additional Analysis Tools
- Binding free energies: see {doc}`mm_pbsa`
- Interaction fingerprints: see {doc}`plif`
- Trajectory convergence: see {doc}`trajectory_convergence`
