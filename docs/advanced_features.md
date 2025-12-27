# Advanced Features

Specialized options for complex systems and performance tuning. Use this page as a checklist; full CLI details stay in {doc}`running_md`.

## Gaussian for Boron-Containing Molecules
- Requires Gaussian. Only ligands/cofactors with boron trigger Gaussian optimization and charge calculation; others use the standard Antechamber path.
- If Gaussian cannot be loaded, boron-containing molecules are skipped while other inputs continue.
- Key flags: `--activate_gaussian`, `--gaussian_exe`, `--gaussian_basis`, `--gaussian_memory`.
- Omit `--activate_gaussian` if Gaussian is already available on the node; use it only to run a load command such as `module load Gaussian/09-d01`.
```bash
run_md -p protein_H_HIS.pdb -l molecules.sdf --cofactor cofactors.sdf --md_time 1 \
  --activate_gaussian "module load Gaussian/09-d01" --gaussian_exe g09 --ncpu 128
```

## MCPB.py for Metalloproteins
- Also requires Gaussian. MCPB.py runs when metal residue names and Gaussian options are provided; otherwise StreaMD falls back to standard `gmx2pdb`.
- Key flags: `--metal_resnames`, `--metal_cutoff`, `--metal_charges` plus Gaussian flags above.
```bash
run_md -p protein_H_HIS.pdb -l molecules.sdf --cofactor cofactors.sdf --md_time 1 \
  --activate_gaussian "module load Gaussian/09-d01" --gaussian_exe g09 --metal_resnames ZN
```

## Custom MDP Files
- Supply `ions.mdp`, `minim.mdp`, `nvt.mdp`, `npt.mdp`, or `md.mdp` via `--mdp_dir`. Missing files fall back to defaults.
- StreaMD may override seed and timing if set via CLI; keep filenames unchanged for detection.
```bash
run_md -p protein.pdb -l ligand.mol --mdp_dir custom_mdp/ --md_time 5
```

## StreaMD Step Control
- `--steps` runs selected pipeline stages when continuing existing runs (`1` preparation, `2` equilibration, `3` production, `4` analysis).
```bash
run_md --wdir_to_continue md_files/md_run/protein_H_HIS_ligand_1 --md_time 3 --steps 3 4
```

## Trajectory Size and Diagnostics
- `--save_traj_without_water` writes water-stripped trajectories for lighter analysis.
- `--no_dr` skips acdoctor ligand diagnostics (use with care).
- `--deffnm` sets a custom prefix for continuing/extending runs.

## Parallelism and GPUs
- Combine `--device`, `--gpu_ids`, `--ntmpi_per_gpu`, and `--mdrun_per_node` to balance CPU/GPU use; see examples in {doc}`running_md`. 
- Set the `--ncpu` option to limit CPU usage; by default, StreaMD uses all available CPU cores.
- Distributed runs use `--hostfile` for Dask servers (PBS/SLURM examples in {doc}`running_md`).
