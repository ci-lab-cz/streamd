# StreaMD: a tool to perform high-throughput automated molecular dynamics simulations

## installation
*Source: https://valdes-tresanco-ms.github.io/gmx_MMPBSA/installation/*

[env.yml](https://github.com/ci-lab-cz/streamd/blob/master/env.yml)
```
conda env create -n md --file env.yml

```

## **Description**
#### **Fully automatic pipeline for molecular dynamics**

#### Features:  
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

### **USAGE**
```
run_md -h
usage: run_md [-h] [-p FILENAME] [-d WDIR] [-l FILENAME] [--cofactor FILENAME] [--clean_previous_md] [--hostfile FILENAME] [-c INTEGER] [--topol topol.top]
              [--topol_itp topol_chainA.itp topol_chainB.itp [topol_chainA.itp topol_chainB.itp ...]] [--posre posre.itp [posre.itp ...]]
              [--protein_forcefield amber99sb-ildn] [--md_time ns] [--npt_time ps] [--nvt_time ps] [--seed int] [--not_clean_log_files] [--steps [STEPS ...]]
              [--wdir_to_continue DIRNAME [DIRNAME ...]] [--deffnm preffix for md files] [--tpr FILENAME] [--cpt FILENAME] [--xtc FILENAME]
              [--ligand_list_file all_ligand_resid.txt] [--ligand_id UNL] [--activate_gaussian module load Gaussian/09-d01]
              [--gaussian_exe g09 or /apps/all/Gaussian/09-d01/g09/g09] [--gaussian_basis B3LYP/6-31G*] [--gaussian_memory 120GB] [--metal_resnames [MN ...]]
              [--metal_cutoff 2.8] [--metal_charges {MN:2, ZN:2, CA:2}]

Run or continue MD simulation. Allowed systems: Protein, Protein-Ligand, Protein-Cofactors(multiple), Protein-Ligand-Cofactors(multiple)

options:
  -h, --help            show this help message and exit

Standard Molecular Dynamics Simulation Run:
  -p FILENAME, --protein FILENAME
                        input file of protein. Supported formats: *.pdb or gro
  -d WDIR, --wdir WDIR  Working directory. If not set the current directory will be used.
  -l FILENAME, --ligand FILENAME
                        input file with compound. Supported formats: *.mol or sdf
  --cofactor FILENAME   input file with compound. Supported formats: *.mol or sdf
  --clean_previous_md   remove a production MD simulation directory if it exists to re-initialize production MD setup
  --hostfile FILENAME   text file with addresses of nodes of dask SSH cluster. The most typical, it can be passed as $PBS_NODEFILE variable from inside a PBS script.
                        The first line in this file will be the address of the scheduler running on the standard port 8786. If omitted, calculations will run on a
                        single machine as usual.
  -c INTEGER, --ncpu INTEGER
                        number of CPU per server. Use all cpus by default.
  --topol topol.top     topology file (required if a gro-file is provided for the protein).All output files obtained from gmx2pdb should preserve the original names
  --topol_itp topol_chainA.itp topol_chainB.itp [topol_chainA.itp topol_chainB.itp ...]
                        Itp files for individual protein chains (required if a gro-file is provided for the protein).All output files obtained from gmx2pdb should
                        preserve the original names
  --posre posre.itp [posre.itp ...]
                        posre file(s) (required if a gro-file is provided for the protein).All output files obtained from gmx2pdb should preserve the original names
  --protein_forcefield amber99sb-ildn
                        Force Field for protein preparation
  --md_time ns          time of MD simulation in ns
  --npt_time ps         time of NPT equilibration in ps
  --nvt_time ps         time of NVT equilibration in ps
  --seed int            seed
  --not_clean_log_files
                        Not to remove all backups of md files
  --steps [STEPS ...]   Run a particular step(s) of the StreaMD run. Options:
                        1 - run preparation step (protein, ligand, cofactor preparation)
                        2 - run MD equilibration step (minimization, NVT, NPT)
                        3 - run MD simulation
                        4 - run MD analysis
                        Ex: 3 4. 
                        If 2 or 3 or 4 step(s) are used --wdir_to_continue argument should be
                        used to provide directories with files obtained during the previous steps
  --wdir_to_continue DIRNAME [DIRNAME ...]
                        single or multiple directories contain simulations created by the tool. 
                        Use with steps 2,3,4 to continue the run. Should consist of: tpr, cpt,
                        xtc and all_ligand_resid.txt files. File all_ligand_resid.txt is optional and used to run md analysis for the ligands. If you want to continue
                        your own simulation not created by the tool use --tpr, --cpt, --xtc and --wdir or arguments (--ligand_list_file is optional and required to
                        run md analysis after simulation )

Continue or Extend Molecular Dynamics Simulation:
  --deffnm preffix for md files
                        Used to run, extend or continue the simulation.
                        If --wdir_to_continue is used files as deffnm.tpr, deffnm.cpt, deffnm.xtc will be searched from --wdir_to_continue directories
  --tpr FILENAME        use explicit tpr to continue a non-StreaMD simulation
  --cpt FILENAME        use explicit cpt to continue a non-StreaMD simulation
  --xtc FILENAME        use explicit xtc to continue a non-StreaMD simulation
  --ligand_list_file all_ligand_resid.txt
                        If you want automatic md analysis for ligands was run after continue of non-StreaMD simulation you should set ligand_list file. Format of the file (no
                        headers): user_ligand_id gromacs_ligand_id. Example: my_ligand UNL. Can be set up or placed into --wdir_to_continue directory(ies)
  --ligand_id UNL       If you want to run an automatic md analysis for a ligand after continue of simulation you can set ligand_id if it is not UNL as a default value

Boron-containing molecules or MCPBPY usage (use together with Standard Molecular Dynamics Simulation Run arguments group):
  --activate_gaussian module load Gaussian/09-d01
                        string that load gaussian module if necessary
  --gaussian_exe g09 or /apps/all/Gaussian/09-d01/g09/g09
                        path to gaussian executable or alias. Requred to run preparation of boron-containing compounds.
  --gaussian_basis B3LYP/6-31G*
                        Gaussian Basis
  --gaussian_memory 120GB
                        Gaussian Memory Usage

MCPBPY usage (use together with Standard Molecular Dynamics Simulation Run and Boron-containing molecules arguments group):
  --metal_resnames [MN ...]
                        Metal residue names to run MCPB.py procedure. Start MCPBPY procedure only if gaussian_exe and activate_gaussian arguments are set up,Otherwise
                        standard gmx2pdb procedure will be run.
  --metal_cutoff 2.8    Metal residue cutoff to run MCPB.py procedure
  --metal_charges {MN:2, ZN:2, CA:2}
                        Metal residue charges in dictionary formatStart MCPBPY procedure only if metal_resnames and gaussian_exe and activate_gaussian arguments are
                        set up,Otherwise standard gmx2pdb procedure will be run.

```

### **Examples**
! Before run MD simulation it is important to prepare protein by yourself to make sure you simulate correct system.
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
Required to obtain relevant poses of the ligand if need
* **Perform docking procedure**  
https://github.com/ci-lab-cz/easydock

#### Run molecular dynamics simulation
``` source activate md ```  
 
**Run simulation for different sytems:**
- Protein in Water
```
run_md -p protein_H_HIS.pdb --md_time 0.1 --nvt_time 100 --npt_time 100 --ncpu 128 
```

- Protein - Ligand
```
run_md -p protein_H_HIS.pdb -l ligand.mol --md_time 0.1 --nvt_time 100 --npt_time 100 --ncpu 128 
```

- Protein - Cofactor
All molecules should present in simulated system, so any problem with preparation of cofactors will interrupt the program. 
```
run_md -p protein_H_HIS.pdb --cofactor cofactors.sdf --md_time 0.1 --nvt_time 100 --npt_time 100 --ncpu 128 

```

**To run simulations with boron-containing compounds**  
*Gaussian Software* should be available.  
Gaussian optimization and charge calculation will be run only for molecules with boron atoms, other molecules will be processed by regular procedure by Antechamber.
If Gaussian cannot be load boron-containing molecules will be skipped.  
Any --ligand or --cofactor files can consist of boron-containing compounds
```
run_md -p protein_H_HIS.pdb -l molecules.sdf --cofactor cofactors.sdf --md_time 0.1 --npt_time 10 --nvt_time 10 --activate_gaussian "module load Gaussian/09-d01" --gaussian_exe g09 --ncpu 128

```

**To run simulations with MCPBPY using parametrization procedure**  
*Gaussian Software* should be available.
```
run_md -p protein_H_HIS.pdb -l molecules.sdf --cofactor cofactors.sdf --md_time 0.1 --npt_time 10 --nvt_time 10 --activate_gaussian "module load Gaussian/09-d01" --gaussian_exe g09 --ncpu 128 --metal_resnames ZN

```
**To run simulations using multiple servers**
```
PBS:
run_md -p protein_H_HIS.pdb -l molecules.sdf --cofactor cofactors.sdf --md_time 0.1 --npt_time 10 --nvt_time 10 --hostfile $PBS_NODEFILE --ncpu 128

SLURM:
srun hostname | sort | uniq > hostfile  
run_md -p protein_H_HIS.pdb -l molecules.sdf --cofactor cofactors.sdf --md_time 0.1 --npt_time 10 --nvt_time 10 --hostfile hostfile --ncpu 128

```

**To extend the simulation**
you can continue your simulation unlimited times. As the --md_time argument user should set up the overall time of the simulation
```
run_md --wdir_to_continue md_files/md_run/protein_H_HIS_ligand_*/ --md_time 0.2
```
or use explicit tpr, cpt and xtc arguments to continue a non-StreaMD simulation
```
run_md --wdir_to_continue md_files/md_run/protein_H_HIS_ligand_1/  --md_time 0.3 --tpr protein_H_HIS_ligand_1/md_out.tpr --cpt protein_H_HIS_ligand_1/md_out.cpt --xtc protein_H_HIS_ligand_1/md_out.xtc
```
in case you don't want to check/run all preparation steps with using non-StreaMD simulations you can use --steps argument 
```
run_md --wdir_to_continue md_files/md_run/protein_H_HIS_ligand_1/ --md_time 0.3 --steps 3 4
```
  
**Output**   
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
md_files/md_run/

protein_H_HIS_ligand_1/
ligand_1.itp          density.xvg  em.trr      ions.tpr                md_out.edr         md_out.tpr             npt.cpt  npt.tpr  nvt.log    potential.xvg         rmsd.xvg       temperature.xvg
cofactor_1.itp        em.edr       frame.pdb   md_centermolsnoPBC.xtc  md_out.gro         md_out.xtc             npt.edr  npt.trr  nvt.mdp    pressure.xvg          rmsf.pdb       topol.top
all.itp               em.gro       gyrate.xvg  md_fit.xtc              md_out.log         md_short_forcheck.xtc  npt.gro  nvt.cpt  nvt.tpr    rmsd_cofactor_1.xvg   rmsf.xvg
all_ligand_resid.txt  em.log       index.ndx   md.mdp                  mdout.mdp          minim.mdp              npt.log  nvt.edr  nvt.trr    rmsd_ligand_1.xvg     solv.gro
complex.gro           em.tpr       ions.mdp    md_out.cpt              md_out_noj_noPBC.xtc  newbox.gro          npt.mdp  nvt.gro  posre.itp  rmsd_xtal.xvg         solv_ions.gro

protein_H_HIS_ligand_2/
```
- **MD output files**
```
md_fit.xtc - MD trajectory with removed PBC and fitted into Protein or Protein-Ligand group
md_short_forcheck.xtc - short trajectory to check if simulation was valid
frame.pdb - a frame for topology

```
- **Analysis data**  
```
potential.png 
temperature.png
pressure.png
density.png
rmsd.png - rmsd of the protein against minimized structure
rmsd_xtal.png - rmsd of the protein against crystal structure
rmsd_cofactor_1.png - rmsd of cofactor against minimized structure
rmsd_cofactor_1_xtal.png - rmsd of the ligand against crystal structure
rmsd_ligand_1.png - rmsd of the ligand against minimized structure
rmsd_ligand_1_xtal.png - rmsd of the ligand against crystal structure
rmsf.png - root mean square fluctuation (RMSF, i.e. standard deviation) of atomic positions in the trajectory
gyrate.png - radius of gyration
```

## Supplementary tools
### MM-PBSA/MM-GBSA energy calculation
#### The tool is based on [gmx_MMPBSA](https://valdes-tresanco-ms.github.io/gmx_MMPBSA/dev/)
Calculation arguments can be changed/added by customized [mmpbsa.in](https://valdes-tresanco-ms.github.io/gmx_MMPBSA/dev/input_file/) file 
#### **USAGE**
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
```

### **Examples**
```
run_gbsa  --wdir_to_run md_files/md_run/protein_H_HIS_ligand_1 md_files/md_run/protein_H_HIS_ligand_2  -c 128 -m mmpbsa.in
```
**Output**   
*each run creates in the working directory (or in the current directory if wdir argument was not set up):*
 1) a unique streaMD log file  
 log_mmpbsa_*start-time*.log 
 Contains important information/warnings/errors about the main run_gbsa program run.
 2) a unique bash log file. 
 log_mmpbsa_bash_*start-time*.log   
 Contains stdout from gmx_MMPBSA 
 3) GBSA_output_*start-time*.csv with summary csv if MMGBSA method was run
 4) PBSA_output_*start-time*.csv with summary csv if MMPBSA method was run
 
 each wdir_to_run has FINAL_RESULTS_MMPBSA_*start-time*.csv with GBSA/PBSA output. 
 
### ProLIF Protein-Ligand Interaction Fingerprints
#### **USAGE**
```
run_prolif -h
usage: run_prolif [-h] [-i DIRNAME [DIRNAME ...]] [--xtc FILENAME] [--tpr FILENAME] [-l STRING] [-s INTEGER] [-a STRING] [-d WDIR] [-v] [--hostfile FILENAME]
                  [-c INTEGER] [--width FILENAME] [--height FILENAME] [-o FILENAME] [--not_save_pics]

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
  -a STRING, --append_protein_selection STRING
                        the string which will be concatenated to the protein selection atoms. Example: "resname ZN or resname MG". (default: None)
  -d WDIR, --wdir WDIR  Working directory for program output. If not set the current directory will be used. (default: None)
  -v, --verbose         print progress. (default: False)
  --hostfile FILENAME   text file with addresses of nodes of dask SSH cluster. The most typical, it can be passed as $PBS_NODEFILE variable from inside a PBS script. The first line in this file will be the address of the scheduler running on the standard port 8786. If omitted, calculations will run on a single machine as usual. (default: None)
  -c INTEGER, --ncpu INTEGER
                        number of CPU per server. Use all cpus by default. (default: 32)
  --width FILENAME      width of the output pictures (default: 15)
  --height FILENAME     height of the output pictures (default: 10)
  -o FILENAME, --occupancy FILENAME
                        occupancy of the unique contacts to show (default: 0.6)
  --not_save_pics       not create html and png files (by frames) for each unique trajectory. Only overall prolif png file will be created. (default: False)

```

### **Examples**
```
run_prolif  --wdir_to_run md_files/md_run/protein_H_HIS_ligand_1 md_files/md_run/protein_H_HIS_ligand_2  -c 128 -v -s 5
```
**Output**  
1) in each directory where xtc file is located  *plifs.csv*, *plifs.png*,*plifs_map.png*, *plifs.html' file for each simulation will be created
2) *prolif_output.csv/png* - aggregated csv/png output file for all analyzed simulations

### Logging
all system info or errors are saved into logging files which would be placed into your main working directory (the current working directory or the path which was passed through --wdir argument):
```
run_md:
log_protein-fname_ligand-fname_cofactor-fname_current-date.log - StreaMD logging user info (status of the )
streamd_bash_protein-fname_ligand-fname_cofactor-fname_current-date.log - StreaMD bash system logging info

run_gbsa:
log_mmpbsa_current-date.log - StreaMD logging user info
log_mmpbsa_bash_current-date.log - StreaMD bash system logging info

```
### Licence
MIT