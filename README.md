# StreaMD: a tool to perform high-throughput automated molecular dynamics simulations

## Table Of Contents
- [Installation](#installation)
- [Description](#description)
- [Features](#features)
- [Molecular Dynamics](#run-molecular-dynamics-simulations)
  - [Data Preparation](#data-preparation)
    - [Target preparation](#1-target-preparation)
    - [Docking](#2-docking-procedure)
  - [MD simulations run](#usage)
    - [Different systems](#run-simulation-for-different-systems)
      - [Protein in water](#protein-in-water)
      - [Protein-Ligand](#protein---ligand)
      - [Protein-cofactors](#protein---cofactors)
      - [Boron-containing compounds](#simulations-with-boron-containing-compounds)
      - [Ligand Binding Metalloprotein with MCPB.py](#simulations-of-ligand-binding-metalloprotein-with-mcpbpy)
    - [Multiple servers](#simulations-using-multiple-servers)
    - [Continue the interrupted simulations](#continue-the-interrupted-simulations) 
    - [Extend the simulation](#extend-the-simulation)
    - [GPU usage](#gpu-usage)
      - [Single GPU](#run-using-single-gpu)
      - [Multiple GPUs](#run-each-simulation-using-multiple-gpus)
      - [The number of thread-MPI ranks per GPU](#increase-the-number-of-thread-mpi-ranks-per-gpu)
      - [Multiple runs per node](#multiple-runs-per-node)
      - [Multiple tasks on the same node while using multiple GPUs](#run-multiple-tasks-on-the-same-node-while-using-multiple-gpus)
      - [Use CPUs only](#run-simulations-only-using-cpus-on-a-server-where-gpus-are-available)
      - [Automatically distribute calculations between the CPU and GPU](#automatically-offload-calculations-across-cpu-and-gpu)
    - [StreaMD output files](#output)
      - [Simulations output](#--md-output-files)
      - [Analysis files](#--analysis-output-files)
- [MM-PBSA/MM-GBSA energy calculation](#mm-pbsamm-gbsa-energy-calculation)
  - [Usage](#usage-1) 
  - [Examples](#examples)
    - [Protein-ligand system](#protein-ligand-system)
    - [Protein-ligand-cofactors system](#protein-ligand-cofactors-system)
  - [Output](#output-1)
- [ProLIF Protein-Ligand Interaction Fingerprints](#prolif-protein-ligand-interaction-fingerprints)
  - [Usage](#usage-2)
  - [Examples](#examples-1)
    - [Protein-ligand system](#protein-ligand-system-1)
    - [Protein-ligand-cofactors system](#protein-ligand-cofactors-system-1)
    - [Effective parallel processing](#effective-parallel-processing)
  - [Output](#output-2)
  - [Supplementary run_prolif analysis scripts](#supplementary-run_prolif-scripts)
- [Trajectory convergence analysis](#trajectory-convergence-analysis)
  - [Usage](#usage-3)
  - [Examples](#examples-2)
  - [Output](#output-3)
- [Logging](#logging)
- [License](#license)
- [Citation](#citation)
    

## installation
Please choose the yaml file for cpu only (**env.yml**) or gpu/cpu (**env_gpu.yml**) usage.  

The **env_gpu.yml** file can only be installed on a machine with an available GPU.
```
conda env create --file env.yml
conda activate md
```
```
pip install streamd
```
or the latest version from github
```
pip install git+https://github.com/ci-lab-cz/streamd.git
```

[Return to the Table Of Contents](#table-of-contents)  

## **Description**
#### **Fully automated pipeline for molecular dynamics**

## Features:  
- supports run of multiple simultaneous molecular dynamics simulations
- supports simulation for different systems:  
    - Protein in Water;  
    - Protein - Ligand;  
    - Protein - Cofactor (multiple);  
    - Protein - Ligand - Cofactor (multiple);  
  
- supports of simulations of boron-containing molecules using Gaussian Software
- supports of simulations of ligand binding metalloproteins with MCPB.py
- supports distributed computing using dask library
- supports of running of parallel simulations on multiple servers
- supports to extend time of MD simulations 
- supports to continue of interrupted MD simulation
- interrupted MD preparation can be restarted by invoking the same command
- implemented tools for end-state free energy calculations (gmx_MMPBSA) and Protein-Ligand Interaction analysis (ProLIF)
- supports customized .mdp files
- interactive trajectory convergence analysis for multiple complexes
- GPU support

[Return to the Table Of Contents](#table-of-contents)  

## **Data preparation**
Before run MD simulation it is important to prepare protein by yourself to make sure you simulate correct system.
> [!CAUTION]  
> To prevent problems with recognition of different atom names by Gromacs, we use the `-ignh` argument, which re-adds the hydrogen atoms, ignoring the original one,
> while converting the protein to a gro file (`gmx pdb2gmx -water tip3p -ignh`). Make sure to rename residues (CYS-CYX, HIS-HIE-HIP) by proper protonation states (_more details in Target Preparation section_).  
> As an option, the user can use the `--noignh` argument when running StreaMD, although there may be problems with atom name recognition that users must resolve themselves or, as another option, provide to Streamd the already prepared protein.gro, topol.top and posre.itp files.   
#### Example of preparation steps before MD: 
#### 1) Target Preparation:  
*Manual preparation:*
- **Fill missing residues and loops**  

*Using Chimera:*   

``Tools -> Sequence -> Structure -> Modeller (loops/refinement)``  
``Tools -> Structure Editing -> Dock Prep ``
* **Explicit water molecules as well as cofactors from a crystal structure can be removed, or if necessarily retained manually;**

* **Remove co-crystallizated ligands;**

* **Add hydrogens based on protonation states.**
  * Check states of histidines and put proper aliases HIE, HID or HIP instead of HIS (otherwise protonation can be distorted during MD preparation stage)  

*type into Chimera cmd:*
```
setattr r type HID :HIS@HD1,DD1,TD1,HND
setattr r type HIP :HID@HE2,DE2,TE2
setattr r type HIE :HIS@HE2
```

#### 2) Docking procedure
Required to obtain relevant poses of the ligand(s) if needed
* **Perform docking procedure by [EasyDock - automate molecular docking tool](https://github.com/ci-lab-cz/easydock)**<br> 

[Return to the Table Of Contents](#table-of-contents)

### Run molecular dynamics simulations   
``` source activate md ```  
### **Usage**
```
run_md -h
usage: run_md [-h] [-p FILENAME] [-d WDIR] [-l FILENAME] [--cofactor FILENAME] [--clean_previous_md] [--hostfile FILENAME] [-c INTEGER] [--mdrun_per_node INTEGER]
              [--device cpu] [--gpu_ids GPU ID [GPU ID ...]] [--ntmpi_per_gpu int] [--topol topol.top]
              [--topol_itp topol_chainA.itp topol_chainB.itp [topol_chainA.itp topol_chainB.itp ...]] [--posre posre.itp [posre.itp ...]]
              [--protein_forcefield amber99sb-ildn] [--noignh] [--md_time ns] [--npt_time ps] [--nvt_time ps] [--seed int] [--no_dr] [--not_clean_backup_files]
              [--steps [STEPS ...]] [--mdp_dir Path to a directory with specific mdp files] [--wdir_to_continue DIRNAME [DIRNAME ...]] [-o OUT_SUFFIX]
              [--save_traj_without_water] [--deffnm preffix for md files] [--tpr FILENAME] [--cpt FILENAME] [--xtc FILENAME] [--ligand_list_file all_ligand_resid.txt]
              [--ligand_id UNL] [--activate_gaussian module load Gaussian/09-d01] [--gaussian_exe g09 or /apps/all/Gaussian/09-d01/g09/g09]
              [--gaussian_basis B3LYP/6-31G*] [--gaussian_memory 120GB] [--metal_resnames [MN ...]] [--metal_cutoff 2.8] [--metal_charges {MN:2, ZN:2, CA:2}]

Run or continue MD simulation. Allowed systems: Protein, Protein-Ligand, Protein-Cofactors(multiple), Protein-Ligand-Cofactors(multiple)

options:
  -h, --help            show this help message and exit
  -o OUT_SUFFIX, --out_suffix OUT_SUFFIX
                        Suffix for output files

Standard Molecular Dynamics Simulation Run:
  -p FILENAME, --protein FILENAME
                        Input file of protein. Supported formats: *.pdb or gro
  -d WDIR, --wdir WDIR  Working directory. If not set the current directory will be used.
  -l FILENAME, --ligand FILENAME
                        Input file with compound(s). Supported formats: *.mol or sdf or mol2
  --cofactor FILENAME   Input file with compound(s). Supported formats: *.mol or sdf or mol2
  --clean_previous_md   Remove a production MD simulation directory if it exists to re-initialize production MD setup
  --hostfile FILENAME   Text file with addresses of nodes of dask SSH cluster. The most typical, it can be passed as $PBS_NODEFILE variable from inside a PBS script.
                        The first line in this file will be the address of the scheduler running on the standard port 8786. If omitted, calculations will run on a
                        single machine as usual.
  -c INTEGER, --ncpu INTEGER
                        Number of CPU per server. Use all available cpus by default.
  --mdrun_per_node INTEGER
                        Number of simulations to run per 1 server.In case of multiple simulations per node, the available CPU cores will be split evenly across these
                        simulations.By default, run only 1 simulation per node and use all available cpus
  --device cpu          Calculate bonded and non-bonded interactions on: auto, cpu, gpu
  --gpu_ids GPU ID [GPU ID ...]
                        List of unique GPU device IDs available to use. Use in case of multiple GPUs usage or using exact GPU devices.
  --ntmpi_per_gpu int   The number of thread-MPI ranks to start per GPU. The default, 1, will start one rank per GPU
  --topol topol.top     Topology file (required if a gro-file is provided for the protein).All output files obtained from gmx2pdb should preserve the original names
  --topol_itp topol_chainA.itp topol_chainB.itp [topol_chainA.itp topol_chainB.itp ...]
                        itp files for individual protein chains (required if a gro-file is provided for the protein).All output files obtained from gmx2pdb should
                        preserve the original names
  --posre posre.itp [posre.itp ...]
                        posre file(s) (required if a gro-file is provided for the protein).All output files obtained from gmx2pdb should preserve the original names
  --protein_forcefield amber99sb-ildn
                        Force Field for protein preparation.Available FF can be found at Miniconda3/envs/md/share/gromacs/top
  --noignh              By default, Streamd uses gmx pdb2gmx -ignh, which re-adds hydrogens using residue names (the correct protonation states must be provided by
                        user) and ignores the original hydrogens. If the --noignh argument is used, the original hydrogen atoms will be preserved during the
                        preparation, although there may be problems with recognition of atom names by GROMACS.
  --md_time ns          Time of MD simulation in ns
  --npt_time ps         Time of NPT equilibration in ps
  --nvt_time ps         Time of NVT equilibration in ps
  --seed int            seed
  --no_dr               Turn off the acdoctor mode and do not check/diagnose problems in the input ligand file in the next attempt if the regular antechamber run for
                        ligand preparation fails (ligand_mol2prep.sh script related issues). Use this argument carefully and ensure that you provide valid structures
  --not_clean_backup_files
                        Not to remove all backups of md files
  --steps [STEPS ...]   Run a particular step(s) of the StreaMD run. Options:1 - run preparation step (protein, ligand, cofactor preparation) 2 - run MD equilibration
                        step (minimization, NVT, NPT) 3 - run MD simulation 4 - run MD analysis. Ex: 3 4 If 2/3/4 step(s) are used --wdir_to_continue argument should
                        be used to provide directories with files obtained during the step 1
  --mdp_dir Path to a directory with specific mdp files
                        To use specific MD settings, the user can provide a path to a directory that contains any of the following .mdp files: ions.mdp, minim.mdp,
                        nvt.mdp, npt.mdp, md.mdp. Missing .mdp files will be replaced by default StreaMD files. Provided .mdp files will be used as templates,
                        although the system StreaMD parameters (seed, nvt_time, npt_time, md_time, and tc-grps (can not be changed by user)) will override the ones
                        provided. Warning: The names of the files must be strictly preserved.
  --save_traj_without_water
                        Save additional md_out_nowater.tpr and md_fit_nowater.xtc files for more memory efficient analysis.
  --wdir_to_continue DIRNAME [DIRNAME ...]
                        Single or multiple directories contain simulations created by the tool. Use with steps 2,3,4 to continue run. ' Should consist of: tpr, cpt,
                        xtc and all_ligand_resid.txt files. File all_ligand_resid.txt is optional and used to run md analysis for the ligands. If you want to continue
                        your own simulation not created by the tool use --tpr, --cpt, --xtc and --wdir or arguments (--ligand_list_file is optional and required to
                        run md analysis after simulation )

Continue or Extend Molecular Dynamics Simulation:
  --deffnm preffix for md files
                        Preffix for the md files. Used to run, extend or continue the simulation. If --wdir_to_continue is used files as deffnm.tpr, deffnm.cpt,
                        deffnm.xtc will be searched from --wdir_to_continue directories
  --tpr FILENAME        Use explicit tpr arguments to continue a non-StreaMD simulation
  --cpt FILENAME        Use explicit cpt arguments to continue a non-StreaMD simulation
  --xtc FILENAME        Use explicit xtc arguments to continue a non-StreaMD simulation
  --ligand_list_file all_ligand_resid.txt
                        If you want automatic md analysis for ligands was run after continue of non-StreaMD simulation you should set ligand_list file. Format of the
                        file (no headers): user_ligand_id gromacs_ligand_id. Example: my_ligand UNL. Can be set up or placed into --wdir_to_continue directory(ies)
  --ligand_id UNL       If you want to run an automatic md analysis for the ligand after continue of simulation you can set ligand_id if it is not UNL default value

Boron-containing molecules or MCPBPY usage (use together with Standard Molecular Dynamics Simulation Run arguments group):
  --activate_gaussian module load Gaussian/09-d01
                        String to load gaussian module if necessary
  --gaussian_exe g09 or /apps/all/Gaussian/09-d01/g09/g09
                        Path to gaussian executable or alias. Required to run preparation of boron-containing compounds.
  --gaussian_basis B3LYP/6-31G*
                        Gaussian Basis
  --gaussian_memory 120GB
                        Gaussian Memory Usage

MCPBPY usage (Use together with Standard Molecular Dynamics Simulation Run and Boron-containing molecules arguments group):
  --metal_resnames [MN ...]
                        Metal residue names to run MCPB.py procedure. Start MCPBPY procedure only if gaussian_exe and activate_gaussian arguments are set up,Otherwise
                        standard gmx2pdb procedure will be run.
  --metal_cutoff 2.8    Metal residue cutoff to run MCPB.py procedure
  --metal_charges {MN:2, ZN:2, CA:2}
                        Metal residue charges in dictionary formatStart MCPBPY procedure only if metal_resnames and gaussian_exe and activate_gaussian arguments are
                        set up, otherwise standard gmx2pdb procedure will be run

``` 

#### **Run simulation for different systems:**
##### Protein in Water
```
run_md -p protein_H_HIS.pdb --md_time 0.1 --nvt_time 100 --npt_time 100 --ncpu 128 
```

##### Protein - Ligand
```
run_md -p protein_H_HIS.pdb -l ligand.mol --md_time 0.1 --nvt_time 100 --npt_time 100 --ncpu 128 
```

##### Protein - Cofactors
All molecules should present in simulated system, so any problem with preparation of cofactors will interrupt the program. 
```
run_md -p protein_H_HIS.pdb --cofactor cofactors.sdf --md_time 0.1 --nvt_time 100 --npt_time 100 --ncpu 128 

```

##### **Simulations with boron-containing compounds**  
*Gaussian Software* should be available.  
Gaussian optimization and charge calculation will be run only for molecules with boron atoms, other molecules will be processed by regular procedure by Antechamber.
If Gaussian cannot be load boron-containing molecules will be skipped.  
Any `--ligand` or `--cofactor` files can consist of boron-containing compounds
```
run_md -p protein_H_HIS.pdb -l molecules.sdf --cofactor cofactors.sdf --md_time 0.1 --npt_time 10 --nvt_time 10 --activate_gaussian "module load Gaussian/09-d01" --gaussian_exe g09 --ncpu 128

```

##### **Simulations of Ligand Binding Metalloprotein with MCPB.py**  
*Gaussian Software* should be available.
```
run_md -p protein_H_HIS.pdb -l molecules.sdf --cofactor cofactors.sdf --md_time 0.1 --npt_time 10 --nvt_time 10 --activate_gaussian "module load Gaussian/09-d01" --gaussian_exe g09 --ncpu 128 --metal_resnames ZN

```
[Return to the Table Of Contents](#table-of-contents)<br>  

#### **Simulations using multiple servers**
```
PBS:
run_md -p protein_H_HIS.pdb -l molecules.sdf --cofactor cofactors.sdf --md_time 0.1 --npt_time 10 --nvt_time 10 --hostfile $PBS_NODEFILE --ncpu 128

SLURM:
srun hostname | sort | uniq > hostfile  
run_md -p protein_H_HIS.pdb -l molecules.sdf --cofactor cofactors.sdf --md_time 0.1 --npt_time 10 --nvt_time 10 --hostfile hostfile --ncpu 128

```
[Return to the Table Of Contents](#table-of-contents)<br>  


#### **Continue the interrupted simulations**  
You can continue the interrupted run by re-executing the previous command. The tool will recognize the checkpoint files and continue the run from the unfinished step.  

[Return to the Table Of Contents](#table-of-contents)<br>  


#### **Extend the simulation**  
you can continue your simulation unlimited times. As the `--md_time` argument user should set up the overall time of the simulation
```
run_md --wdir_to_continue md_files/md_run/protein_H_HIS_ligand_*/ --md_time 0.2
```
or use explicit `--tpr`, `--cpt` and `--xtc` arguments to continue a non-StreaMD simulation
```
run_md --wdir_to_continue md_files/md_run/protein_H_HIS_ligand_1/  --md_time 0.3 --tpr protein_H_HIS_ligand_1/md_out.tpr --cpt protein_H_HIS_ligand_1/md_out.cpt --xtc protein_H_HIS_ligand_1/md_out.xtc
```
in case you don't want to check/run all preparation steps with using non-StreaMD simulations you can use `--steps` argument 
```
run_md --wdir_to_continue md_files/md_run/protein_H_HIS_ligand_1/ --md_time 0.3 --steps 3 4
```

[Return to the Table Of Contents](#table-of-contents)<br>  


#### GPU usage
The StreaMD tool supports running of energy minimization, NVT, and NPT equilibration steps, as well as production simulations on GPU(s).   
When run with `--device gpu` argument the StreaMD offloads nb, update, pme, bonded, pmefft (all which can be run on GPU) computations to GPU.  
More details on: https://manual.gromacs.org/documentation/current/user-guide/mdrun-performance.html

The performance and gpu usage strongly depend on the type of hardware and size of the system.  
It is always good to check the GPU usage (for example, by `nvidia-smi` command).

##### Run using single GPU
```
run_md -p protein_HIS.pdb -l ligand.mol --md_time 1 --device gpu
```

To improve the performance one can use multiple GPUs or start multiple ranks per GPU.

##### Run each simulation using multiple GPUs 
> [!WARNING]
> Increasing the number of GPUs does not always improve performance.
The each single simulation will use all provided GPUs. 
```
run_md -p protein_HIS.pdb -l ligand.mol --md_time 1 --device gpu --gpu_ids 0 1 2 3
```
> [!WARNING]
> If you want to split the simulations across multiple GPUs but still run the task on the same node use the --mdrun_per_node argument (see below).

##### Increase the number of thread-MPI ranks per GPU  
```-ntmpi_per_gpu_``` argument is used to calculate _total number of thread-MPI ranks_ (```gmx mdrun -ntmpi_ X```, where **X**=ntmpi_per_gpu*number of GPUs to use). By default, ```ntmpi_per_gpu``` equals 1, although usage of 2 thread-MPI ranks per GPU may return better performance.
```
run_md -p protein_HIS.pdb -l ligand.mol --md_time 1 --device gpu -ntmpi_per_gpu 2
```

##### Multiple runs per node
Sometimes for more full GPU usage user can start multiple simulations on a single/multiple nodes and the tool automatically splits the available CPU cores across these simulations:
```
run_md -p protein_HIS.pdb -l ligands.sdf --md_time 1 --device gpu --mdrun_per_node 2
```
##### Run multiple tasks on the same node while using multiple GPUs
> [!WARNING]
> All simulations will still be utilizing all provided GPUs which can lead to suboptimal GPU load. This feature is still under development.
```
run_md -p protein_HIS.pdb -l ligands.sdf --md_time 1 --device gpu --mdrun_per_node 2 --gpu_ids 0 1
```
##### Run simulations only using CPUs on a server where GPUs are available
```
run_md -p protein_HIS.pdb -l ligand.mol --md_time 1 --device cpu
```
##### Automatically offload calculations across CPU and GPU
To let GROMACS automatically offload calculations between CPU and GPU may be optimal on hardware where the CPUs are relatively powerful compared to the GPUs.
```
run_md -p protein_HIS.pdb -l ligand.mol --md_time 1 --device auto --gpu_ids 0
```
   
[Return to the Table Of Contents](#table-of-contents)<br>  


#### **Output**   
*each run creates in the working directory (or in the current directory if wdir argument was not set up):*
 1) a unique streaMD log file which name contains name of the protein, ligand file, cofactor file and time of run.  
 log_*protein-fname*\_*ligand-fname*\_*cofactor-fname*\_*start-time*.log  
 Contains important information/warnings/errors about the main program run.
 2) a unique bash log file.  
 streamd_bash_*protein-fname*\_*ligand-fname*\_*cofactor-fname*\_*start-time*.log  
 Contains stdout from Gromacs and Antechamber.  
 
will be created the next folders:
```
md_files/
- md_preparation/
 -- protein/
 -- ligands/
 -- cofactors/
- md_run/
 -- protein-id_ligand-id
 --- md_analysis
```

```
md_files/md_preparation/protein/:
protein.gro  posre.itp  topol.top

OR for multiple chain protein:
md_files/md_preparation/protein/:
protein.gro 
topol.top
posre_Protein_chain_A.itp    
posre_Protein_chain_B.itp  
topol_Protein_chain_A.itp
topol_Protein_chain_B.itp
```

```
md_files/md_preparation/ligands/:
all_resid.txt

ligand_1/
ligand_1.frcmod  ligand_1.lib     ligand_1.top        sqm.in
ligand_1.gro     ligand_1.mol     leap.log            sqm.out
ligand_1.inpcrd  ligand_1.mol2    posre_ligand_1.itp  sqm.pdb
ligand_1.itp     ligand_1.prmtop  resid.txt           tleap.in

ligand_2/
..
```

```
md_files/md_preparation/cofactors/:
all_resid.txt

cofactor_1/
cofactor_1.frcmod  cofactor_1.lib     cofactor_1.top        sqm.in
cofactor_1.gro     cofactor_1.mol     leap.log              sqm.out
cofactor_1.inpcrd  cofactor_1.mol2    posre_cofactor_1.itp  sqm.pdb
cofactor_1.itp     cofactor_1.prmtop  resid.txt             tleap.in

cofactor_2/
```

```
md_files/md_run/protein_H_HIS_ligand_1/
ligand_1.itp          em.trr       ions.tpr    md_out.edr         md_out.tpr             npt.cpt   npt.tpr  nvt.log    
cofactor_1.itp        em.edr       frame.pdb   md_out.gro         md_out.xtc             npt.edr   npt.trr  nvt.mdp  topol.top
all.itp               em.gro       md_fit.xtc  md_out.log         md_short_forcheck.xtc  npt.gro   nvt.cpt  nvt.tpr    
all_ligand_resid.txt  em.log       index.ndx   md.mdp             mdout.mdp              minim.mdp npt.log  nvt.edr  solv.gro
complex.gro           em.tpr       ions.mdp    md_out.cpt         newbox.gro             npt.mdp   nvt.gro  posre.itp solv_ions.gro

```
##### - **Analysis output files**
```
md_files/md_run/protein_H_HIS_ligand_1/md_analysis

potential_protein_HIS_ligand_1.{csv, png, xtc} - potential energy of Energy Minimization step calculated by gmx energy 
temperature_protein_HIS_ligand_1.{csv, png, xtc} - system temperature of NVT simulation calculated by gmx energy
density_protein_HIS_ligand_1.{csv, png, xtc}  - total density of NPT simulations calculated by gmx energy
pressure_protein_HIS_ligand_1.{csv, png, xtc} - system pressure of NPT simulations calculated by gmx energy
rmsd_protein_HIS_ligand_1.{csv, png} - Root mean square deviation of atomic positions for backbone and ligand and Active Site (default 5A) if Protein-Ligand simulation was performed
rmsf_protein_HIS_ligand_1.{csv, png, xtc, pdb} - root mean square fluctuation (RMSF, i.e. standard deviation) of atomic positions in the trajectory
gyrate_protein_HIS_ligand_1.{csv, png, xtc} - radius of gyration
```
##### - **MD output files**
```
md_fit.xtc - MD trajectory with removed PBC and fitted into Protein or Protein-Ligand group
md_short_forcheck.xtc - short trajectory to check if simulation was valid
frame.pdb - a frame for topology

```

[Return to the Table Of Contents](#table-of-contents)  


## Supplementary tools
### MM-PBSA/MM-GBSA energy calculation
*The tool is based on [gmx_MMPBSA](https://valdes-tresanco-ms.github.io/gmx_MMPBSA/dev/) tool.*  
Calculation parameters can be changed/added by customized [mmpbsa.in](https://valdes-tresanco-ms.github.io/gmx_MMPBSA/dev/input_file/) file 
> [!NOTE] 
> The user can control the length of the trajectory for analysis by editing 
> the startframe, endframe, and interval arguments (in mmpbsa.in file), where every 100 frames equals 1 ns.

#### **Usage**
```
run_gbsa -h
usage: run_gbsa [-h] [-i DIRNAME [DIRNAME ...]] [--topol topol.top] [--tpr md_out.tpr] [--xtc md_fit.xtc] [--index index.ndx] [-m mmpbsa.in] [-d WDIR]
                [--out_files OUT_FILES [OUT_FILES ...]] [--hostfile FILENAME] [-c INTEGER] [--ligand_id UNL] [-a [STRING ...]] [--clean_previous]

Run MM-GBSA/MM-PBSA calculation using gmx_MMPBSA tool

options:
  -h, --help            show this help message and exit
  -i DIRNAME [DIRNAME ...], --wdir_to_run DIRNAME [DIRNAME ...]
                        single or multiple directories for simulations. Should consist of: tpr, xtc, ndx files
  --topol topol.top     topol file from the the MD simulation. Will be ignored if --wdir_to_run is used
  --tpr md_out.tpr      tpr file from the the MD simulation. Will be ignored if --wdir_to_run is used
  --xtc md_fit.xtc      xtc file of the simulation. Trajectory should have no PBC and be fitted on the Protein_Ligand group. Will be ignored if --wdir_to_run is used
  --index index.ndx     Gromacs index file from the simulation. Will be ignored if --wdir_to_run is used
  -m mmpbsa.in, --mmpbsa mmpbsa.in
                        MMPBSA input file. If not set up default template will be used.
  -d WDIR, --wdir WDIR  Working directory for program output. If not set the current directory will be used.
  --out_files OUT_FILES [OUT_FILES ...]
                        gmxMMPBSA out files (FINAL*.dat) to parse. If set will be used over other variables.
  --hostfile FILENAME   text file with addresses of nodes of dask SSH cluster. The most typical, it can be passed as $PBS_NODEFILE variable from inside a PBS script.
                        The first line in this file will be the address of the scheduler running on the standard port 8786. If omitted, calculations will run on a
                        single machine as usual.
  -c INTEGER, --ncpu INTEGER
                        number of CPU per server. Use all cpus by default.
  --ligand_id UNL       Ligand residue ID
  -a [STRING ...], --append_protein_selection [STRING ...]
                        residue IDs whuch will be included in the protein system (cofactors).Example: ZN MG
  --clean_previous      Clean previous temporary gmxMMPBSA files
  -o string, --out_suffix string
                        Unique suffix for output files. By default, start-time_unique-id.
                        Unique suffix is used to separate outputs from different runs.

  
```

#### **Examples**
#### Protein-ligand system
```
run_gbsa  --wdir_to_run md_files/md_run/protein_H_HIS_ligand_1 md_files/md_run/protein_H_HIS_ligand_2  -c 128 -m mmpbsa.in
```

#### Protein-ligand-cofactors system

In case, you have a cofactor-protein system, the ```--ligand_id``` and ```--append_protein_selection``` arguments can be used
in this scenario. The system residue names of both the ligand and cofactors can be found in the _md_files/md_run/protein-ligand/all_ligand_resid.txt_ file. 

To calculate the affinity between a protein-cofactor system and a ligand, use `--ligand_id 'UNL'` and `--append_protein_selection 'CFR'`. 
By default, StreaMD uses 'UNL' as the ligand system residue name, but it is recommended to verify the exact residue name in the _all_ligand_resid.txt_ file.
Additionally, specify `--append_protein_selection 'CFR'` to include cofactor into the protein system for the calculations (you can find the exact cofactor residue name also in the _all_ligand_resid.txt_ file).
```
run_gbsa  --wdir_to_run md_files/md_run/protein_H_HIS_ligand_*  --append_protein_selection MG GTP
```
To calculate binding free energy between the protein system and the cofactor, use `--ligand_id 'CFR'` instead.
```
run_gbsa  --wdir_to_run md_files/md_run/protein_H_HIS_ligand_*  --append_protein_selection MG --ligand_id GTP
```


#### **Output**   
*each run creates in the working directory (or in the current directory if wdir argument was not set up):*
Unique suffix is used to separate outputs from different runs.
 1) a unique streaMD log file  
 log_mmpbsa_*unique-suffix*.log 
 Contains important information/warnings/errors about the main run_gbsa program run.
 2) a unique bash log file. 
 log_mmpbsa_bash_*unique-suffix*.log   
 Contains stdout from gmx_MMPBSA 
 3) GBSA_output_*unique-suffix*.csv with summary csv if MMGBSA method was run
 4) PBSA_output_*unique-suffix*.csv with summary csv if MMPBSA method was run
 
 each wdir_to_run has FINAL_RESULTS_MMPBSA_*unique-suffix*.csv with GBSA/PBSA output. 
 
[Return to the Table Of Contents](#table-of-contents)  


### ProLIF Protein-Ligand Interaction Fingerprints
#### **Usage**
```
usage: run_prolif [-h] [-i DIRNAME [DIRNAME ...]] [--xtc FILENAME] [--tpr FILENAME] [-l STRING] [-s INTEGER] [--protein_selection STRING] [-a STRING] [-d WDIR] [-v]
                  [--hostfile FILENAME] [-c INTEGER] [--n_jobs INTEGER] [--width FILENAME] [--height FILENAME] [--occupancy float] [--not_save_pics] [-o string]

Get protein-ligand interactions from MD trajectories using ProLIF module.

options:
  -h, --help            show this help message and exit
  -i DIRNAME [DIRNAME ...], --wdir_to_run DIRNAME [DIRNAME ...]
                        single or multiple directories for simulations.
                                                     Should consist of: md_out.tpr and md_fit.xtc files (default: None)
  --xtc FILENAME        input trajectory file (XTC). Will be ignored if --wdir_to_run is used (default: None)
  --tpr FILENAME        input topology file (TPR). Will be ignored if --wdir_to_run is used (default: None)
  -l STRING, --ligand STRING
                        residue name of a ligand in the input trajectory. (default: UNL)
  -s INTEGER, --step INTEGER
                        step to take every n-th frame. ps (default: 1)
  --protein_selection STRING
                        The protein selection atoms. Example: "protein" or "protein and byres around 20.0 resname UNL" (default: protein)
  -a STRING, --append_protein_selection STRING
                        the string which will be concatenated to the protein selection atoms. Example: "resname ZN or resname MG". (default: None)
  -d WDIR, --wdir WDIR  Working directory for program output. If not set the current directory will be used. (default: None)
  -v, --verbose         print progress. (default: False)
  --hostfile FILENAME   text file with addresses of nodes of dask SSH cluster. The most typical, it can be passed as $PBS_NODEFILE variable from inside a PBS script. The first line in this file will be the address of the scheduler running on the standard port 8786. If omitted, calculations will run on a single machine as usual. (default: None)
  -c INTEGER, --ncpu INTEGER
                        number of CPU per server. Use all available cpus by default.
  --n_jobs INTEGER      Number of processes to run per each trajectory. Provided CPUs (--ncpu arg) will be distributed between number of trajectories and number of processes per each trajectory (--n_jobs arg). (default: 1)
  --width FILENAME      width of the output pictures (default: 15)
  --height FILENAME     height of the output pictures (default: 10)
  --occupancy float     occupancy of the unique contacts to show. Applied for plifs_occupancyX.html (for each complex) and prolif_output_occupancyX.png (all systems aggregated plot) (default: 0.6)
  --not_save_pics       not create html and png files (by frames) for each unique trajectory. Only overall prolif png file will be created. (default: False)
  -o string, --out_suffix string
                        Unique suffix for output files. By default, start-time_unique-id.Unique suffix is used to separate outputs from different runs.

```

#### **Examples**
##### Protein-ligand system
```
run_prolif  --wdir_to_run md_files/md_run/protein_H_HIS_ligand_1 md_files/md_run/protein_H_HIS_ligand_2  -c 128 -v -s 5
```
##### Protein-ligand-cofactors system
In case, you have a cofactor-protein system, the ```--ligand_id``` and ```--append_protein_selection``` arguments can be used in this scenario. 
The system residue names of both the ligand and cofactors can be found in the _md_files/md_run/protein-ligand/all_ligand_resid.txt_ file. 

To obtain interaction fingerprints between a protein-cofactor system and a ligand, use `--ligand_id 'UNL'` and `--append_protein_selection 'CFR'` arguments. 
By default, StreaMD uses 'UNL' as the ligand system residue name, but it is recommended to verify the exact residue name in the _all_ligand_resid.txt_ file.
Additionally, specify `--append_protein_selection 'CFR'` to include cofactor into the protein system for the calculations (you can find the exact cofactor residue name also in the _all_ligand_resid.txt_ file).
```
run_prolif  --wdir_to_run md_files/md_run/protein_H_HIS_ligand_* --append_protein_selection MG GTP
```
To calculate interaction fingerprints between protein system and cofactor, use `--ligand_id 'CFR'` instead.
```
run_prolif  --wdir_to_run md_files/md_run/protein_H_HIS_ligand_*  --ligand_id 'GTP'
```
#### Effective parallel processing
There are 2 arguments which users can control to achieve more effective parallelization:
1) `--ncpu`
The overall maximum number of cores available for usage. By default, StreaMD utilizes all available CPUs.
2) `--n_jobs`
Number of processes to run per each interaction analysis tasks. Equivalent of _n_jobs_ for parallel processing in ProLIF.
By default, StreaMD distributes the specified number of cores (`--ncpu`) evenly
between the available CPUs and the number of tasks to execute (e.g., multiple directories provided via `--wdir_to_run`).
However, by default, the `--n_jobs` value is limited to 12 to avoid the [bottleneck issue](https://github.com/chemosim-lab/ProLIF/issues/110) described by the ProLIF authors. 
Users can override this limitation by explicitly specifying the `--n_jobs` argument value.

#### **Output**  
1) in each directory where xtc file is located  *plifs.csv*, *plifs.png*,*plifs_map.png*, *plifs.html* file for each simulation will be created
2) *prolif_output_*unique-suffix*.csv/png* - aggregated csv/png output file for all analyzed simulations. Unique suffix is used to separate outputs from different runs.

#### Supplementary run_prolif scripts
_run_prolif applies all this scripts automatically. Use it if you want more detailed analysis or to change the picture/fonts sizes._  
**prolif_drawmap**  
Draw prolif plot for analysis binding mode of multiple ligands
````
prolif_drawmap  -h
usage: prolif_drawmap [-h] -i FILENAME [FILENAME ...] [-o FILENAME] [--width FILENAME] [--height FILENAME] [--base_size FILENAME]

Draw prolif plot for analysis binding mode of multiple ligands

options:
  -h, --help            show this help message and exit
  -i FILENAME [FILENAME ...], --input FILENAME [FILENAME ...]
                        input file with prolif output for the set of molecules. Supported formats: *.csv
                        Ex: prolif_output.csv
  --occupancy float
                        minimum occupancy of the unique contacts to show
  --width int      width of the output picture
  --height int     height of the output picture
  --base_size int  base size of the output picture

````

**prolif_draw_by_frame**

```
prolif_draw_by_frame -h
usage: prolif_draw_by_frame [-h] -i [FILENAME ...] [-o FILENAME] [--filt_only_H] [--width FILENAME] [--height FILENAME] [--base_size FILENAME]

options:
  -h, --help            show this help message and exit
  -i [FILENAME ...], --input [FILENAME ...]
                        input file with prolif output for the unique molecule. Supported formats: *.csv
                        Ex: plifs.csv
  --occupancy float
                        minimum occupancy of the unique contacts to show. Show all contacts by default.
  --filt_only_H         filt residues where only hydrophobic contacts occur
  --width int      width of the output picture
  --height int     height of the output picture
  --base_size int  base size of the output picture
```

[Return to the Table Of Contents](#table-of-contents)  


### Trajectory convergence analysis 

To identify converged segments of molecular dynamics trajectories _run_analysis_ module calculates
the average root-mean-square deviation (RMSD) of the ligand, protein, and active site residues within 5Å of the ligand, as well as the standard deviation of RMSD for the same trajectory segment.
The average RMSD should provide an insight into the ligand movement or rotation relative to its initial pose, while the standard deviation reflects the stability of the ligand pose within the selected trajectory segment. 
The conclusions can be valuable for subsequent MM-GBSA/PBSA calculations.

> [!NOTE]
> Such analysis is automatically performed as a part of the run_md default run. 
> For the separate analysis the user can use either run_module with the --wdir_to_continue _directory_path_ and --steps 4 arguments to perform full analysis of already obtained simulations or the run_analysis script itself (the rmsd file in csv format from run_md will be still needed).

#### **USAGE**
````
usage: run_analysis.py [-h] [-i FILENAME [FILENAME ...]] [--rmsd_type backbone [backbone ...]] [--time_ranges 0-1 5-10 9-10 [0-1 5-10 9-10 ...]] [-d dirname]
                       [--paint_by PAINT_BY] [-o OUT_SUFFIX] [--title RMSD Mean vs RMSD Std]

Run rmsd analysis for StreaMD output files

options:
  -h, --help            show this help message and exit
  -i FILENAME [FILENAME ...], --input FILENAME [FILENAME ...]
                        input file(s) with rmsd. Supported formats: *.csv. Required columns: time(ns) ligand_name system. Sep: /\t.
  --rmsd_type backbone [backbone ligand ActiveSite]
                        Column names in the input file to use for the analysis 
  --time_ranges 0-1 5-10 9-10 [0-1 5-10 9-10 ...]
                        Time ranges in nanoseconds. Default: Default: start-end, middle-end, end-1 (in nanoseconds)
  -d dirname, --wdir dirname
                        Output files directory
  --paint_by csv_file   File to paint by additional column. Required columns: - Protein-ligand simulations: system\tligand_name. - Protein only in water simulations: system. Sep: /\t. The plot will be painted by any other than system and
                        ligand_name column.
  -o OUT_SUFFIX, --out_suffix OUT_SUFFIX
                        Suffix for output files
  --title RMSD Mean vs RMSD Std
                        Title for html plot. Default: RMSD Mean vs RMSD Std
````
#### **Examples**
preferred way:
```
run_md  --wdir_to_continue md_files/md_run/protein_H_HIS_ligand_1 md_files/md_run/protein_H_HIS_ligand_2  --steps 4 -d md_files
```
by the script itself:
```
run_rmsd_analysis -i rmsd_all_systems.csv --rmsd_type "backbone" "ligand" "ActiveSite5.0A" --paint_by exp_data.csv -o protein --title "Protein. RMSD Mean vs RMSD Std" --time_ranges 0-1 0-2 0-10 5-10 9-10

```
#### **Output**
1) csv output file containing output data  
rmsd_mean_std_time-ranges_*start-time*.csv
2) html file with interactive visualization  
rmsd_mean_std_time-ranges_*start-time*.html

[Return to the Table Of Contents](#table-of-contents)  

## Logging  
All system information or errors are saved into logging files which would be placed into your main working directory (the current working directory or the path which was passed through --wdir argument):  
**run_md:**
```
log_protein-fname_ligand-fname_cofactor-fname_*start-time*.log - StreaMD logging user info (status of the )
streamd_bash_protein-fname_ligand-fname_cofactor-fname_*start-time*.log - StreaMD bash system logging info
```

**run_gbsa:**
```
log_mmpbsa_*unique-suffix*.log - StreaMD logging user info
log_mmpbsa_bash_*unique-suffix*.log - StreaMD bash system logging info
```

**run_prolif:**
```
log_prolif_*unique-suffix*.log - StreaMD logging user info
```

[Return to the Table Of Contents](#table-of-contents)   


## License
MIT
## Citation
Ivanova A, Mokshyna O, Polishchuk P.  
StreaMD: the toolkit for high-throughput molecular dynamics simulations.  
*J. Cheminf.* **2024**, 16 (1), 123.  
https://doi.org/10.1186/s13321-024-00918-w

[Return to the Table Of Contents](#table-of-contents)  