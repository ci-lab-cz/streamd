# Trajectory Convergence Analysis     
`run_rmsd_analysis`

`run_rmsd_analysis` identifies converged trajectory segments by comparing RMSD means and standard deviations for ligand, protein, and active site residues (5 Angstrom default). The same analysis runs automatically as part of `run_md` step 4 from {doc}`running_md`.

## Usage (`run_rmsd_analysis -h`)
```
usage: run_rmsd_analysis [-h] [-i FILENAME [FILENAME ...]] [--rmsd_type backbone [backbone ...]] [--time_ranges 0-1 5-10 9-10 [0-1 5-10 9-10 ...]] [-d dirname]
                       [--paint_by PAINT_BY] [-o OUT_SUFFIX] [--title RMSD Mean vs RMSD Std]

Run rmsd analysis for StreaMD output files

options:
  -h, --help            show this help message and exit
  -i FILENAME [FILENAME ...], --input FILENAME [FILENAME ...]
                        input file(s) with rmsd. Supported formats: *.csv. Required columns: time(ns) ligand_name system. Sep: /\t.
  --rmsd_type backbone [backbone ligand ActiveSite]
                        Column names in the input file to use for the analysis 
  --time_ranges 0-1 5-10 9-10 [0-1 5-10 9-10 ...]
                        Time ranges in nanoseconds. Default: Default: start-end, middle-end, end-1 (in nanoseconds)
  -d dirname, --wdir dirname
                        Output files directory
  --paint_by csv_file   File to paint by additional column. Required columns: - Protein-ligand simulations: system\tligand_name. - Protein only in water simulations: system. Sep: /\t. The plot will be painted by any other than system and
                        ligand_name column.
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
run_rmsd_analysis -i rmsd_all_systems.csv --rmsd_type "backbone" "ligand" "ActiveSite5.0A" \
  --paint_by exp_data.csv -o protein --title "Protein. RMSD Mean vs RMSD Std" --time_ranges 0-1 0-2 0-10 5-10 9-10
```

## Outputs
1) `rmsd_mean_std_time-ranges_*start-time*.csv` - output data
2) `rmsd_mean_std_time-ranges_*start-time*.html` - interactive visualization

See {doc}`outputs` for where RMSD analysis artifacts are stored within each run.
