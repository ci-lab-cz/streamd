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
all_ligand_resid.txt  em.log       mdout.mdp   minim.mdp          posre.itp              solv_ions.gro
```

## Analysis Outputs   
`md_analysis/`
- `potential_*.{xvg,png}`: potential energy from energy minimization (gmx energy)
- `temperature_*.{xvg,png}`: system temperature during NVT
- `density_*.{xvg,png}`: total density during NPT
- `pressure_*.{xvg,png}`: system pressure during NPT
- `rmsd_*.{csv,png}`: RMSD for backbone and Cα, ligand, active site (5 A default), and optional local-pocket ligand RMSD
- `rmsf_*.{xvg,png,pdb}`: per-residue Cα RMSF traces and structures (fitted on and computed over Cα atoms)
- `gyrate_*.{xvg,png}`: radius of gyration

RMSD columns are reported after global protein-backbone alignment unless noted otherwise:
- `backbone`: protein backbone (N, CA, C, O) RMSD; the backbone is the superposition (fit) group.
- `CA`: protein Cα RMSD (the Cα subset, reported in the same backbone-aligned frame).
- `ligand`: ligand heavy-atom RMSD after global protein-backbone alignment. This measures ligand motion relative to the globally aligned protein.
- `ActiveSite5.0A`: RMSD of reference-defined binding-site backbone atoms. The pocket is defined from protein residues within the active-site distance of the ligand in the reference frame, then kept fixed for all trajectory frames.
- `ligand_local`: ligand heavy-atom RMSD after local alignment on the same reference-defined pocket backbone. This measures ligand pose stability relative to its local binding pocket and is useful for flexible or multidomain proteins where global protein motion can inflate the ordinary `ligand` RMSD.

`ActiveSite5.0A` and `ligand_local` are only defined for a single binding site. If more than one copy of the ligand residue is present (e.g. one per chain in a multimer), both are skipped for that system with a warning, while `CA`, `backbone` and `ligand` are still reported.

See {doc}`mm_pbsa` for MM-PBSA/MM-GBSA outputs, {doc}`plif` for interaction fingerprints, and {doc}`trajectory_convergence` for RMSD convergence analysis.

## MD Outputs    
`md_run/<system>/`
- `md_fit.xtc`: PBC-corrected trajectory centered using the complete molecular complex and fitted using Protein + primary ligand for systems containing a ligand. Protein-only systems use `Protein`; cofactor-only systems use `Protein` plus the cofactor groups.
- `md_short_forcheck.xtc`: short trajectory for sanity checks
- `frame.pdb`: representative start frame around 10 ps (0.01 ns)
- `last_frame.pdb`: last frame of the trajectory
- `md_fit_nowater.xtc`, `md_out_nowater.tpr`, `md_out_nowater.gro`: optional no-water analysis files when `--save_traj_without_water` is enabled
- Standard GROMACS files (`md_out.*`, `em.*`, `nvt.*`, `npt.*`, `index.ndx`, `topol.top`, etc.)
