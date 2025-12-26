# PLIF Analysis

StreaMD bundles ProLIF-based tools to extract protein-ligand interaction fingerprints (PLIFs) from trajectories.

## Run ProLIF
Provide one or more simulation directories (each containing `md_out.tpr` and `md_fit.xtc`) or explicit files.
```bash
# Multiple protein-ligand systems
run_prolif --wdir_to_run md_files/md_run/protein_H_HIS_ligand_1 md_files/md_run/protein_H_HIS_ligand_2 \
  --ncpu 128 --step 5 --ligand UNL --verbose
```

## Cofactors and Custom Selections
Use `--append_protein_selection` to include cofactors in the protein group or swap ligand/cofactor roles by changing `--ligand`.
```bash
run_prolif --wdir_to_run md_files/md_run/protein_H_HIS_ligand_* \
  --append_protein_selection MG GTP   # include cofactors in the protein selection
run_prolif --wdir_to_run md_files/md_run/protein_H_HIS_ligand_* \
  --ligand GTP                        # treat a cofactor as the ligand of interest
```

## Parallelism
- `--ncpu`: maximum CPUs available; distributed across trajectories
- `--n_jobs`: ProLIF processes per trajectory (defaults to 1; capped to avoid bottlenecks unless explicitly set)

## Outputs
- Per-trajectory: `plifs.csv`, `plifs.png`, `plifs_map.png`, `plifs.html` (optional HTML/PNG controlled by `--not_save_pics`)
- Aggregated: `prolif_output_<suffix>.{csv,png}` summarizing all processed systems

## Supplementary Plotting
- `prolif_drawmap`: build summary plots across many ligands from one or more `prolif_output.csv` files
- `prolif_draw_by_frame`: visualize frame-level contacts from a single `plifs.csv`

Use `-h` with any command to see the full option set.
