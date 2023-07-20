# Python module to automate molecular dynamics

## Dependency
* **Python**
* **AmberTools**
* **RDKit**
* **gmx_MMPBSA**
* **mpi4py**
* **ParmEd**
* **Dask**

## installation
*Source: https://valdes-tresanco-ms.github.io/gmx_MMPBSA/installation/*

```
will be fixed

conda update conda
conda create -n gmxMMPBSA python=3.8 -y -q
conda activate gmxMMPBSA
conda install mpi4py
conda install -c conda-forge ambertools=21 compilers
conda install -c conda-forge rdkit

python -m pip install git+https://github.com/ParmEd/ParmEd.git@16fb236    
python -m pip install gmx_MMPBSA

```
## **Description**
#### **Fully automatic pipeline for molecular dynamics**

#### Features:  
- supports simulation for different systems:  
    - Protein in Water;  
    - Protein - Ligand;  
    - Protein - Multiple coenzymes;  
    - Protein - Ligand - Multiple coenzymes;  
    
- supports distributed computing using dask library
- supports of running of parallel simulations on multiple servers
- supports to extend time of MD simulations 
- supports to continue of interrupted MD simulation
- interrupted MD preparation can be restarted by invoking the same command


### **USAGE**
```
$ python run_md.py -h
usage: run_md.py [-h] [-p FILENAME] [-d WDIR] [-l FILENAME]
                 [--cofactor FILENAME] [--not_add_H] [--clean_previous_md]
                 [--hostfile FILENAME] [-c INTEGER] [--topol topol.top]
                 [--topol_itp topol.top [topol.top ...]]
                 [--posre posre.itp [posre.itp ...]] [--md_time ns]
                 [--npt_time ps] [--nvt_time ps] [--not_clean_log_files]
                 [--tpr FILENAME] [--cpt FILENAME] [--xtc FILENAME]
                 [--wdir_to_continue DIRNAME [DIRNAME ...]]
                 [--deffnm preffix for md files]

Run or continue MD simulation. Allowed systems: Protein in water, Protein - Ligand, 
Protein - multiple Coenzymes, Protein - Ligand - multiple Coenzymes, 

options:
  -h, --help            show this help message and exit
  -p FILENAME, --protein FILENAME
                        input file of the protein. Supported formats: *.pdb or gro
  -d WDIR, --wdir WDIR  working directory. If not set the current directory will be used.
  -l FILENAME, --ligand FILENAME
                        input file with compounds. Supported formats: *.mol or sdf.
                        If molecules have names it will be used for file names.  
  --cofactor FILENAME   input file with compound. Supported formats: *.mol or sdf.
                        If molecules have names it will be used for file names. 
  --gromacs_version     Default: GROMACS/2021.4-foss-2020b-PLUMED-2.7.3 .
  --not_add_H           disable to add hydrogens to molecules for simulation. By default always add
  --clean_previous_md   remove all previous prepared for the current system MD files. 
                        Prepared ligand and protein files will be intact and reused if it exists.
  --hostfile FILENAME   text file with addresses of nodes of dask SSH cluster.
                        The most typical, it can be passed as $PBS_NODEFILE
                        variable from inside a PBS script. The first line in
                        this file will be the address of the scheduler running
                        on the standard port 8786. If omitted, calculations
                        will run on a single machine as usual.
  -c INTEGER, --ncpu INTEGER
                        number of CPU per server. Use all cpus by default.
  --topol topol.top     required if gro file of the protein is provided
  --topol_itp topol.top [topol_chainA.itp topol_chainB.itp]
                        if the protein has multiple chain there are multiple topologies for each of them.
                        Itp files for each chain are required
  --posre posre.itp [posre.itp ...]
                        required if gro file of the protein is provided. Can
                        be multiple files for each chain.
  --md_time ns          time of MD simulation in ns
  --npt_time ps         time of NPT equilibration in ps
  --nvt_time ps         time of NVT equilibration in ps
  --not_clean_log_files Not to remove all backups of md files
  --tpr FILENAME        tpr file from the previous MD simulation. Use to extend or continue the simulation
  --cpt FILENAME        cpt file from previous simulation. Use to extend or continue the simulation
  --xtc FILENAME        xtc file from previous simulation. Use to extend or continue the simulation
  --wdir_to_continue DIRNAME [DIRNAME ...]
                        directories for the previous simulation. Use to extend or continue the simulation
                        Should consist of: tpr, cpt, xtc files
  --deffnm preffix for md files
                        preffix for the previous md files. Use to extend or continue the simulation.
                        Only if wdir_to_continue is used. Use if each --tpr, --cpt, --xtc arguments are not set up. 
                        Files deffnm.tpr, deffnm.cpt, deffnm.xtc will be used from wdir_to_continue.
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
``` source activate gmxMMPBSA ```  
 
**Run simulation for different sytems:**
- Protein in Water
```
python run_md.py -p protein_H_HIS.pdb --md_time 0.1 --nvt_time 100 --npt_time 100 --c 128 
```

- Protein - Ligand
```
python run_md.py -p protein_H_HIS.pdb -l ligand.mol --md_time 0.1 --nvt_time 100 --npt_time 100 --c 128 
```

- Protein - Coenzyme

```
python run_md.py -p protein_H_HIS.pdb --cofactor cofactors.sdf --md_time 0.1 --nvt_time 100 --npt_time 100 --c 128 

```

**To extend the simulation**
```
python run_md.py --wdir_to_continue md_preparation/md_files/protein_H_HIS_ligand_1/ md_preparation/md_files/protein_H_HIS_ligand_1/ --md_time 0.2 --deffnm md_out
```
you can continue your simulation unlimited times just use previous extended deffnm OR tpr, cpt and xtc arguments 
```
python run_md.py --wdir_to_continue md_preparation/md_files/protein_H_HIS_ligand_1/  --md_time 0.3 --deffnm md_out_0.2
```
or 
```
python run_md.py --wdir_to_continue md_preparation/md_files/protein_H_HIS_ligand_1/  --md_time 0.3 --tpr md_preparation/md_files/protein_H_HIS_ligand_1/md_out_0.2.tpr --cpt md_preparation/md_files/protein_H_HIS_ligand_1/md_out_0.2.cpt --xtc md_preparation/md_files/protein_H_HIS_ligand_1/md_out_0.2.xtc
```
  
**Output**   
In the working directory (or in the current directory if wdir argument was not set up) will be created the next folders:
```
md_preparation/
md_preparation/protein/
md_preparation/var_lig/
md_preparation/system_lig/
md_preparation/md_files/
```

```
md_preparation/protein/:
protein.gro  posre.itp  topol.top

OR for multiple chain protein:
md_preparation/protein/:
protein.gro 
topol.top
posre_Protein_chain_A.itp    
posre_Protein_chain_B.itp  
topol_Protein_chain_A.itp
topol_Protein_chain_B.itp
```

```
md_preparation/var_lig/:
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
md_preparation/system_lig/:
all_resid.txt

cofactor_1/
cofactor_1.frcmod  cofactor_1.lib     cofactor_1.top        sqm.in
cofactor_1.gro     cofactor_1.mol     leap.log              sqm.out
cofactor_1.inpcrd  cofactor_1.mol2    posre_cofactor_1.itp  sqm.pdb
cofactor_1.itp     cofactor_1.prmtop  resid.txt             tleap.in

cofactor_2/
...
```

```
md_preparation/md_files/

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
potential.xvg
temperature.xvg
pressure.xvg
density.xvg
rmsd.xvg - rmsd of the protein against minimized structure
rmsd_xtal.xvg - rmsd of the protein against crystal structure
rmsd_cofactor_1.xvg - rmsd of each cofactors
rmsd_ligand_1.xvg - rmsd of the ligand
rmsf.xvg - root mean square fluctuation (RMSF, i.e. standard deviation) of atomic positions in the trajectory
gyrate.xvg - radius of gyration
```
_All analysis results you can visualize with xmgrace:_  
`` for i in *.xvg; do gracebat $i;done `` 


###Licence
BSD-3