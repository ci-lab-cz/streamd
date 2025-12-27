# Outputs and Files

StreaMD organizes each run into predictable folders with logs, intermediates, and analysis outputs.

## Logs
- `log_*protein-fname*_*ligand-fname*_*cofactor-fname*_*start-time*.log`: StreaMD status, warnings, and errors
- `streamd_bash_*protein-fname*_*ligand-fname*_*cofactor-fname*_*start-time*.log`: stdout from GROMACS and Antechamber

Detailed logs for each ligand/cofactor preparation runs
- `md_preparation/ligands/ligandN/streamd_bash_*protein-fname*_*ligand-fname*_*cofactor-fname*_*start-time*.log`: stdout for Antechamber or Gaussian ligand preparation run 
- `md_preparation/cofactors/cofactorN/streamd_bash_*protein-fname*_*ligand-fname*_*cofactor-fname*_*start-time*.log`: stdout for Antechamber or Gaussian cofactor preparation run 

## Directory Layout
Each run creates in the working directory (or current directory if `--wdir` is not set):
```
md_files/
- md_preparation/
 -- protein/
 -- ligands/
 -- cofactors/
- md_run/
 -- protein-id_ligand-id_replicaN
 --- md_analysis
```

### `md_preparation/protein/`
```
protein.gro  posre.itp  topol.top
```
For multiple chains:
```
protein.gro
topol.top
posre_Protein_chain_A.itp
posre_Protein_chain_B.itp
topol_Protein_chain_A.itp
topol_Protein_chain_B.itp
```

### `md_preparation/ligands/`
```
all_resid.txt

ligand_1/
ligand_1.frcmod  ligand_1.lib     ligand_1.top        sqm.in
ligand_1.gro     ligand_1.mol     leap.log            sqm.out
ligand_1.inpcrd  ligand_1.mol2    posre_ligand_1.itp  sqm.pdb
ligand_1.itp     ligand_1.prmtop  resid.txt           tleap.in
streamd_bash_1o2i_H_HIS_ligand_1_date.log

ligand_2/
...
```

### `md_preparation/cofactors/`
```
all_resid.txt

cofactor_1/
cofactor_1.frcmod  cofactor_1.lib     cofactor_1.top        sqm.in
cofactor_1.gro     cofactor_1.mol     leap.log              sqm.out
cofactor_1.inpcrd  cofactor_1.mol2    posre_cofactor_1.itp  sqm.pdb
cofactor_1.itp     cofactor_1.prmtop  resid.txt             tleap.in
streamd_bash_1o2i_H_HIS_cofactor_1_date.log

cofactor_2/
...
```

### `md_preparation/systems/`
```
protein_H_HIS_ligand_1.itp  complex.gro  ions.tpr   minim.mdp   nvt.mdp                     solv.gro                         topol.top
all.itp                     index.ndx    md.mdp     newbox.gro  posre_protein_ligand_1.itp  solv_ions.gro
all_ligand_resid.txt        ions.mdp     mdout.mdp  npt.mdp     posre.itp                   streamd_bash_protein_HIS_ligand_1_date.log
```

### `md_run/<system>[_replicaN]/`
```
ligand_1.itp          em.trr       ions.tpr    md_out.edr         md_out.tpr             npt.cpt   npt.tpr  nvt.log
cofactor_1.itp        em.edr       frame.pdb   md_out.gro         md_out.xtc             npt.edr   npt.trr  nvt.mdp  topol.top
all.itp               em.gro       md_fit.xtc  md_out.log         md_short_forcheck.xtc  npt.gro   nvt.cpt  nvt.tpr  nvt.edr  solv.gro
complex.gro           em.tpr       ions.mdp    md_out.cpt         newbox.gro             npt.mdp   npt.log  index.ndx
all_ligand_resid.txt  em.log       md_fit.xtc  md_out.xtc         mdout.mdp              minim.mdp npt.tpr  posre.itp solv_ions.gro
```

## Analysis Outputs   
`md_analysis/`
- `potential_*.{csv,png,xtc}`: potential energy from energy minimization (gmx energy)
- `temperature_*.{csv,png,xtc}`: system temperature during NVT
- `density_*.{csv,png,xtc}`: total density during NPT
- `pressure_*.{csv,png,xtc}`: system pressure during NPT
- `rmsd_*.{csv,png}`: RMSD for backbone, ligand, and active site (5 A default)
- `rmsf_*.{csv,png,xtc,pdb}`: RMSF traces and structures
- `gyrate_*.{csv,png,xtc}`: radius of gyration

See {doc}`mm_pbsa` for MM-PBSA/MM-GBSA outputs, {doc}`plif` for interaction fingerprints, and {doc}`trajectory_convergence` for RMSD convergence analysis.

## MD Outputs    
`md_run/<system>/`
- `md_fit.xtc`: PBC-removed, protein- or protein-ligand-fitted trajectory
- `md_short_forcheck.xtc`: short trajectory for sanity checks
- `frame.pdb`: starting frame (0.1 ns)
- `last_frame.pdb`: last frame of the trajectory
- Standard GROMACS files (`md_out.*`, `em.*`, `nvt.*`, `npt.*`, `index.ndx`, `topol.top`, etc.)
