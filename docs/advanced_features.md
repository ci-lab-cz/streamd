# Advanced Features

Specialized options for complex systems and performance tuning.

## Boron-Containing Molecules (Gaussian)
- Requires Gaussian. Only molecules with boron trigger Gaussian optimization; others follow the standard Antechamber workflow.
- If Gaussian cannot be loaded, boron-containing molecules are skipped while other ligands/cofactors proceed.
```bash
run_md -p protein_H_HIS.pdb -l molecules.sdf --cofactor cofactors.sdf --md_time 1 \
  --activate_gaussian "module load Gaussian/09-d01" --gaussian_exe g09 --ncpu 128
```

## Metalloproteins with MCPB.py
- Requires Gaussian plus metal residue names/charges.
- Automatically starts MCPB.py when `--metal_resnames` and Gaussian options are present; otherwise falls back to standard preparation.
```bash
run_md -p protein_H_HIS.pdb -l molecules.sdf --md_time 1 \
  --activate_gaussian "module load Gaussian/09-d01" --gaussian_exe g09 \
  --metal_resnames ZN --metal_cutoff 2.8 --metal_charges "{ZN:2}"
```

## Custom `.mdp` Files
Use `--mdp_dir` to supply any of `ions.mdp`, `minim.mdp`, `nvt.mdp`, `npt.mdp`, or `md.mdp`. Missing files fall back to StreaMD defaults.
```bash
run_md -p protein.pdb -l ligand.mol --mdp_dir custom_mdp/ --md_time 5
```
StreaMD preserves user parameters but overrides seed and timing options if provided via CLI. Keep filenames unchanged to ensure they are detected.

## Run Specific Steps
Use `--steps` to operate on existing outputs without repeating the full pipeline.
```
1 - preparation (protein/ligand/cofactor)
2 - equilibration (minimization, NVT, NPT)
3 - production MD
4 - analysis
```
```bash
run_md --wdir_to_continue md_files/md_run/protein_H_HIS_ligand_1 --steps 3 4 --md_time 3
```

## Additional Controls
- `--save_traj_without_water` saves water-stripped trajectories for lighter analysis.
- `--no_dr` skips acdoctor checks for problematic ligands (use with care).
- `--deffnm` sets custom prefixes when continuing/extending runs.
