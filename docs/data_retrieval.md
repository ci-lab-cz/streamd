# Data Retrieval

StreaMD organizes every run into predictable folders with logs, intermediates, and analysis outputs.

## Directory Layout
```
md_files/
  md_preparation/
    protein/        # prepared protein topology/posre files
    ligands/        # per-ligand parameters and residue IDs
    cofactors/      # per-cofactor parameters and residue IDs
    systems/        # combined system (complex.gro, mdp files, index.ndx, topologies)
  md_run/
    <system>[_replicaN]/
      md_analysis/  # analysis outputs for the run
```

## Logs
- `log_*protein*_*ligand*_*cofactor*_*start*.log`: StreaMD status and warnings
- `streamd_bash_*protein*_*ligand*_*cofactor*_*start*.log`: raw output from GROMACS/Antechamber

## Analysis Outputs (`md_analysis/`)
- `potential_*.{csv,png}`, `temperature_*.{csv,png}`, `density_*.{csv,png}`, `pressure_*.{csv,png}`: energy and thermodynamic profiles
- `rmsd_*.{csv,png}`: RMSD for backbone, ligand, and active site (5 A default)
- `rmsf_*.{csv,png,xtc,pdb}`: RMSF traces and structures
- `gyrate_*.{csv,png,xtc}`: radius of gyration summaries

## MD Outputs (`md_run/<system>/`)
- `md_fit.xtc`: PBC-removed, protein- or protein-ligand-fitted trajectory
- `md_short_forcheck.xtc`: short trajectory for sanity checks
- `frame.pdb`: starting frame (0.1 ns)
- `last_frame.pdb`: final frame from the trajectory
- Standard GROMACS files (`md_out.*`, `em.*`, `nvt.*`, `npt.*`, `index.ndx`, `topol.top`, etc.)

Use these artifacts as inputs for downstream analysis tools such as MM/PB(GB)SA or ProLIF.

## MM-PBSA/MM-GBSA (`run_gbsa`)
Run binding free-energy calculations with gmx_MMPBSA using StreaMD trajectories.
```bash
# Multiple protein-ligand systems
run_gbsa --wdir_to_run md_files/md_run/protein_H_HIS_ligand_1 md_files/md_run/protein_H_HIS_ligand_2 \
  -c 128 -m mmpbsa.in

# Include cofactors in the protein selection
run_gbsa --wdir_to_run md_files/md_run/protein_H_HIS_ligand_* \
  --append_protein_selection MG GTP

# Treat a cofactor as the ligand of interest
run_gbsa --wdir_to_run md_files/md_run/protein_H_HIS_ligand_* \
  --ligand_id GTP
```
Outputs include `GBSA_output_<suffix>.csv`, `PBSA_output_<suffix>.csv`, and `FINAL_RESULTS_MMPBSA_<suffix>.csv` per system, plus decomposition tables when `mmpbsa.in` contains a `&decomp` block.

## Trajectory Convergence (`run_rmsd_analysis`)
Identify stable segments by comparing RMSD means and standard deviations.
```bash
# Analyze existing runs (preferred)
run_md --wdir_to_continue md_files/md_run/protein_H_HIS_ligand_1 md_files/md_run/protein_H_HIS_ligand_2 --steps 4 -d md_files

# Directly on RMSD CSV files
run_rmsd_analysis -i rmsd_all_systems.csv --rmsd_type backbone ligand ActiveSite5.0A \
  --paint_by exp_data.csv --time_ranges 0-1 0-2 0-10 5-10 9-10 --out_suffix protein
```
Results: `rmsd_mean_std_time-ranges_<suffix>.csv` and an interactive HTML plot.
