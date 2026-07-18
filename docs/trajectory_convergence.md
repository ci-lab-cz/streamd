# Trajectory Convergence Analysis     
`run_rmsd_analysis`

`run_rmsd_analysis` identifies converged trajectory segments by comparing RMSD means and standard deviations for ligand, protein, active site residues (5 Angstrom default), and other RMSD columns present in the input files. The same analysis runs automatically as part of `run_md` step 4 from {doc}`running_md`.

StreaMD RMSD files can contain two ligand-pose metrics. `ligand` is ligand heavy-atom RMSD after global protein-backbone alignment, so it measures ligand motion relative to the globally aligned protein. `ligand_local` is ligand heavy-atom RMSD after local alignment on the reference-defined binding-site backbone, so it measures ligand pose stability relative to its pocket. `ligand_local` is most useful for flexible or multidomain proteins where movement of the binding-site domain relative to the rest of the protein can make the global `ligand` RMSD look high even when the ligand remains stable in the pocket.

Each input file must contain the required base columns `time(ns)`, `system`, and `ligand_name`; a missing base column raises a clear per-file error. Requested RMSD metrics are handled leniently so one incomplete system cannot abort a whole batch: a metric (such as `ligand_local` or `ActiveSite5.0A`) missing from only some input files is still aggregated for the files that have it, with those systems left empty and a warning logged; a metric missing from every input file is dropped with a warning. An error is raised only when none of the requested metrics are present anywhere.

For protein-only or legacy RMSD files, keep the `ligand_name` column present with empty values.

By default, StreaMD analyses:

1. the complete available trajectory interval;
2. the final half of the available interval, beginning at the first whole-nanosecond boundary at or after the exact midpoint. If that boundary would fall at or beyond the trajectory end, the preceding whole-nanosecond boundary is used;
3. the final 1 ns, when the trajectory duration exceeds 1 ns.

The paint-by file must contain at least one value column in addition to `protein_name` and, for ligand systems, `ligand_name`.

## Usage   
`run_rmsd_analysis -h`
```
usage: run_rmsd_analysis [-h] -i FILENAME [FILENAME ...] [--rmsd_type COLUMN [COLUMN ...]] [--time_ranges 0-1 5-10 9-10 [0-1 5-10 9-10 ...]] [-d dirname]
                       [--paint_by PAINT_BY] [-o OUT_SUFFIX] [--title RMSD Mean vs RMSD Std]

Run rmsd analysis for StreaMD output files

options:
  -h, --help            show this help message and exit
  -i FILENAME [FILENAME ...], --input FILENAME [FILENAME ...]
                        Input RMSD TSV file(s). Required base columns: time(ns), system, ligand_name. A column requested through --rmsd_type that is missing from some files is skipped for those systems (with a warning); one missing from all files is dropped.
  --rmsd_type COLUMN [COLUMN ...]
                        RMSD column names to aggregate, for example CA, backbone, ligand, ligand_local, or ActiveSite5.0A
  --time_ranges 0-1 5-10 9-10 [0-1 5-10 9-10 ...]
                        Time ranges in nanoseconds. Default: start-end, middle-end, and the final 1 ns.
  -d dirname, --wdir dirname
                        Output files directory
  --paint_by PAINT_BY   File to paint by additional column. Required columns: protein-ligand simulations: protein_name, ligand_name; protein-only simulations: protein_name. The plot is coloured using the first column other than
                        protein_name and ligand_name.
  -o OUT_SUFFIX, --out_suffix OUT_SUFFIX
                        Suffix for output files
  --title RMSD Mean vs RMSD Std
                        Title for html plot. Default: RMSD Mean vs RMSD Std
```

## Examples
Preferred (uses `run_md` analysis step):
```bash
run_md --wdir_to_continue md_files/md_run/protein_H_HIS_ligand_1 md_files/md_run/protein_H_HIS_ligand_2 --steps 4 -d md_files
```

Direct usage:
```bash
run_rmsd_analysis -i rmsd_all_systems.csv --rmsd_type "CA" "backbone" "ligand" "ligand_local" "ActiveSite5.0A" \
  --paint_by exp_data.csv -o protein --title "Protein. RMSD Mean vs RMSD Std" --time_ranges 0-1 0-2 0-10 5-10 9-10
```

## Outputs
1) `rmsd_mean_std_time-ranges_*start-time*.csv` - output data
2) `rmsd_mean_std_time-ranges_*start-time*.html` - interactive visualization

See {doc}`outputs` for where RMSD analysis artifacts are stored within each run.
