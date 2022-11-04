# Protein-ligand molecular dynamics simulation + PB(GB)SA calculation
## Scripts:  
* **01_complex_preparation_md.pbs**  
* **02_pbsa.pbs**
* **03_sum_result_pbsa.pbs**

## Dependency
* **Python (3.8)**
* **AmberTools**
* **RDKit**
* **gmx_MMPBSA (1.4.3)**
* **mpi4py**
* **ParmEd**

## installation
*Source: https://valdes-tresanco-ms.github.io/gmx_MMPBSA/installation/*

```

conda update conda
conda create -n gmxMMPBSA python=3.8 -y -q
conda activate gmxMMPBSA
conda install mpi4py
conda install -c conda-forge ambertools=21 compilers
conda install -c conda-forge rdkit

python -m pip install git+https://github.com/ParmEd/ParmEd.git@16fb236    
python -m pip install gmx_MMPBSA==1.4.3

```

### **Practical Example**

#### 1) Target Preparation:  
*Manual preparation:*
* **Fill missing residues and loops;**
  * Using Chimera: ```Tools -> Sequence -> Structure -> Modeller (loops/refinement)```  
  
*All next steps can be perfomed by Chimera command:*  ```Tools -> Structure Editing -> Dock Prep ```
* **Explicit water molecules as well as cofactors from a crystal structure can be removed, or if necessarily retained manually;**

* **Remove co-crystallizated ligands;**

* **Add hydrogens based on protonation states.**
  * For Docking procedure better to remove non-polar hydrogens:
  ```
  Select -> Chemistry -> IDATM type -> HC
  Actions -> Atoms/Bonds -> delete
  ```
  * Check states of histidines and put proper aliases HIE, HID or HIP instead of HIS (otherwise protonation can be distorted during MD preparation stage)
  * For Molecular dynamics and PBSA calculation - all hydrogens should be preserved. 
  
*Target preparation can be automatically performed by Schrodinger*
https://www.schrodinger.com/kb/29

#### 2) Docking procedure
Required to obtain relevant pose
* **Perform docking procedure**  
Using script from https://github.com/DrrDom/rdkit-scripts/blob/master/vina_dock.py   

#### 3) Ligand preparation 
* **Prepare docked poses for molecular dynamics procedure**  
Using scripts from https://github.com/DrrDom/rdkit-scripts/ and chemaxon
```
# extract doking poses from db
python get_sdf_from_db.py -i grow_procedure.db -o exmp.sdf

# add All Hydrogens
standardize -c "addexplicitH" exmp.sdf -f sdf -o exmp_H.sdf

# Convert sdf to mol format (Amber uses only pdb or mol)
molconvert mol exmp_H.sdf -o ligand.mol -m

# Just for convenience rename ligandN.mol to Mol.mol
rename_func() { molid=$(head -1 $1); cp $1 $molid.mol; }
export -f rename_func
parallel rename_func ::: *.mol
```

#### 4) Protein-Ligand molecular dynamics simulation
*01_complex_preparation_md.pbs*
1) Put all prepared above *.mol files into *ligands* directory
2) Here we use GNU parallel module but it can be replaced by for or xargs commands
Run each ligand-protein systems separately on each cluster node  

**! Important**  
Before to run do not forget to set appropriate number of cpu (Default:128) in each pbs scripts:  
```
#PBS -l select=1:ncpus=set_your_number:ompthreads=2 
```
**To run:**  
``module load parallel``
Absolute path:
```
parallel -j1  "test -d {1/.} || mkdir {1/.}; qsub -v lfile= {1},pfile= path/protein_H_full.pdb,script_path=path/scripts,wdir=/workdir/{1/.},mdtime=1  /home/md-scripts/Protein-Ligand_MD_PBGBSA/01_complex_preparation_md.pbs" ::: path/ligands/*.mol
OR
for i in mols/*.mol; do echo $i; molfile=${i##*/}; molname=${molfile%.*}; echo $molname; qsub -A OPEN-25-38  -l walltime=48:00:00 -v lfile=$(pwd)/$i,pfile=$(pwd)/4o2b_full_clean_ready_H_HIS.pdb,script_path=/home/md-scripts/scripts/,wdir=$molname,mdtime=10,gromacs_version='GROMACS/2021.4-foss-2020b-PLUMED-2.7.3'  /home/md-scripts/Protein-Ligand_MD_PBGBSA/01-1_complex_preparation_md.pbs ;done

```
Relative path can be set as well, but remember that script will be running from working directory and all relative paths should consider it.
Like:
```
parallel -j1  "test -d {1/.} || mkdir {1/.}; qsub -v lfile=../{1},pfile=../protein_H_full.pdb,script_path=../../scripts,wdir=/workdir/{1/.},mdtime=1 ../../scripts/01_complex_preparation_md.pbs" ::: ligands/*.mol
```

**Description:**  
Command creates separate folders named as each molname.mol if there were not created previously, and then runs all simulations simultaneous. 
- _lfile_ - ligand file (use directory path upper because script would be run in the created molname directory).
- _pfile_ - protein file prepared for MD (Missing residues and loops problems were resolved).  
If pdb protein file is used all hydrogens will be ignored and re-added by gmx (-ignh) 
due to of extreme intolerance of the gmx pdb2gmx command to the variability in the numeration and typing of hydrogens.  
Usually it works well but if some unusual case of your target protonation is observed you can set previously prepared gro file
 (gro file can be obtained by command`` gmx pdb2gmx -f $pfile -o $PNAME.gro -water tip3p -ignh ``)   
- _mdtime_ - time (in _ns_ units) for molecular dynamics simulation.
- _script_path_ - path to _scripts_ dir from this repository. Dir consists of all molecular dynamics parameters (*.mdp files), 
additional scripts, PBSA(GBSA) parameters (mmpbsa.in).
- _wdir_ - working directory. If skipped, current directory will be used. Should be different from _script_path_ argument.
- _gromacs_version_ - GROMACS version to use (Default: _GROMACS/2018.1-intel-2017c-hybrid-single-PLUMED_)

**Output** 
- **MD trajectories**  
md_out.[gro,tpr,cpt,xtc]
- **MD trajectory with removed PBC.**  
md_out_noPBC.xtc
- **Intermediate analysis data**  
*.xvg files. 
- **Final analysis data**  
rmsd for both protein and ligand, rmsf, radius of gyration.
  
_All analysis results you can visualize with xmgrace:_  
`` for i in *.xvg; do gracebat $i;done `` 

#### 5) gmx_MMPBSA energy evaluation of the protein - ligand system
``parallel -j1  "test -d {1/.} || mkdir {1/.}; qsub -v tpr=md_out.tpr,xtc=md_fit.xtc,script_path=path/scripts,wdir=wdir/{1/.},NP=128 02_pbsa.pbs" ::: ligands/*.mol ``

**Description:**  
Command creates separate folders named as each molname.mol if there were not created previously, and then runs all calculations simultaneous. 
- _tpr_ - md_out.tpr
- _xtc_ - md_fit.xtc Complex trajectories. Make sure the trajectory is fitted and PBC have been removed.
- _script_path_ - path to _scripts_ dir from this repository. Dir consists of PBSA(GBSA) parameters (mmpbsa.in).
- _wdir_ - working directory. If skipped, current directory will be used. Should be different from _script_path_ argument.
- _NP_ - number of processors. Should be less than total number of using frames. Can be skipped. Default: 128  
- _topol_ - topol.top. Can be skipped. Default: $wdir/topol.top.
- _index_ - index.ndx. Can be skipped. Default: $wdir/index.top
- _LNAME_ - Ligand ID. Previous script (_01_complex_preparation_md.pbs_) has prepared files with fixed names, so by default it can be skipped. Default: UNL
- _gromacs_version_ - GROMACS version to use. Can be skipped. Default: _GROMACS/2018.1-intel-2017c-hybrid-single-PLUMED_.

**Output** 

#### 6) MMPBSA result collecting
``qsub -v wdir=$(pwd) 03_sum_result_pbsa.pbs``

**Description:**  
Search all dirs in the working directory for FINAL_RESULTS_MMPBSA.dat files. Save necessary information (name of directory, PBSA, GBSA and Interaction energy) to result files.

**Output**    
- _Sum_Result_POISSON_BOLTZMANN.csv_
- _Sum_Result_GENERALIZED_BORN.csv_
- _Sum_Result_IE.csv_

### Calculation parameters
All molecular dynamics (*.mdp files) and PBSA (*.in) parameters files are located in the _scripts_ directory. 
All (*.mdp) parameters are ready for running in real calculations.  
The most variable part is PB(GB)SA calculation. Here we set the most optimal parameter values.  
**But it is highly recommended to set the most critical parameters according to what is better for your system:**
````
# General namelist variables
&general
  startframe           = 1                                              # First frame to analyze
  endframe             = 9999999                                        # Last frame to analyze
  interval             = 1                                              # Number of frames between adjacent frames analyzed
# optimal number of frames better be from 200 up to 1000, more frames gives more accurate calculation but become more time-consuming
  
  ie_segment           = 25                                             # Trajectory segment to calculate interaction entropy
# last 25% of frames (set above) will be used for calculation of interaction entropy 

&pb
  indi                 = 1.0                                            # Internal dielectric constant
# Value depends on total charge of the system
# indi = 1  - default
# indi = 2 - 1-2 charged amino acids or ion-ion interaction(s)
# indi = 4 - highly charged system (several charged amino acids) or strong ion-ion interaction(s) 
# !Important. if you get predominantly positive energies try to raise up indi constant  

  fillratio            = 4.0                                            # See "fillratio" in AmberTools/PBSA manual
# value should be smaller for large proteins, such as membrane protein. More details: https://valdes-tresanco-ms.github.io/gmx_MMPBSA/examples/Protein_membrane/?h=fillratio
````
More details about PB(GB)SA input variables you can find there:
- https://valdes-tresanco-ms.github.io/gmx_MMPBSA/input_file/
- https://ambermd.org/doc12/Amber21.pdf
    
### License
BSD-3

### References:
* Valdés-Tresanco, M.S., Valdés-Tresanco, M.E., Valiente, P.A. and Moreno E. gmx_MMPBSA: A New Tool to Perform End-State Free Energy Calculations with GROMACS. Journal of Chemical Theory and Computation, 2021 17 (10), 6281-6291. https://pubs.acs.org/doi/10.1021/acs.jctc.1c00645.   
* Bill R. Miller, T. Dwight McGee, Jason M. Swails, Nadine Homeyer, Holger Gohlke, and Adrian E. Roitberg. MMPBSA.py: An Efficient Program for End-State Free Energy Calculations. Journal of Chemical Theory and Computation, 2012 8 (9), 3314-3321. https://pubs.acs.org/doi/10.1021/ct300418h. 
