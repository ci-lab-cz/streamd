# MM-PBSA/MM-GBSA (`run_gbsa`)

StreaMD wraps gmx_MMPBSA to compute binding free energies from trajectories.

## Usage (`run_gbsa -h`)
```
run_gbsa -h
usage: run_gbsa [-h] [-i DIRNAME [DIRNAME ...]] [--topol topol.top] [--tpr md_out.tpr] [--xtc md_fit.xtc] [--index index.ndx] [-m mmpbsa.in] [-d WDIR]
                [--out_files OUT_FILES [OUT_FILES ...]] [--hostfile FILENAME] [-c INTEGER] [--ligand_id UNL] [-a [STRING ...]] [--clean_previous]

Run MM-GBSA/MM-PBSA calculation using gmx_MMPBSA tool

options:
  -h, --help            show this help message and exit
  -i DIRNAME [DIRNAME ...], --wdir_to_run DIRNAME [DIRNAME ...]
                        single or multiple directories for simulations. Should consist of: tpr, xtc, ndx files
  --topol topol.top     topol file from the the MD simulation. Will be ignored if --wdir_to_run is used
  --tpr md_out.tpr      tpr file from the the MD simulation. Will be ignored if --wdir_to_run is used
  --xtc md_fit.xtc      xtc file of the simulation. Trajectory should have no PBC and be fitted on the Protein_Ligand group. Will be ignored if --wdir_to_run is used
  --index index.ndx     Gromacs index file from the simulation. Will be ignored if --wdir_to_run is used
  -m mmpbsa.in, --mmpbsa mmpbsa.in
                        MMPBSA input file. If not set up default template will be used.
  -d WDIR, --wdir WDIR  Working directory for program output. If not set the current directory will be used.
  --out_files OUT_FILES [OUT_FILES ...]
                        gmxMMPBSA out files (FINAL*.dat) to parse. If set will be used over other variables.
  --hostfile FILENAME   text file with addresses of nodes of dask SSH cluster. The most typical, it can be passed as $PBS_NODEFILE variable from inside a PBS script.
                        The first line in this file will be the address of the scheduler running on the standard port 8786. If omitted, calculations will run on a
                        single machine as usual.
  -c INTEGER, --ncpu INTEGER
                        number of CPU per server. Use all cpus by default.
  --ligand_id UNL       Ligand residue ID
  -a [STRING ...], --append_protein_selection [STRING ...]
                        residue IDs whuch will be included in the protein system (cofactors).Example: ZN MG
  --clean_previous      Clean previous temporary gmxMMPBSA files
  -o string, --out_suffix string
                        Unique suffix for output files. By default, start-time_unique-id.
                        Unique suffix is used to separate outputs from different runs.
```

Calculation parameters can be changed/added by customized `mmpbsa.in` (see https://valdes-tresanco-ms.github.io/gmx_MMPBSA/dev/input_file/). Control trajectory length via `startframe`, `endframe`, and `interval`. Required inputs (tpr/xtc/ndx) are produced by {doc}`running_md`.

## Examples

### Protein-ligand system
```bash
run_gbsa --wdir_to_run md_files/md_run/protein_H_HIS_ligand_1 md_files/md_run/protein_H_HIS_ligand_2 -c 128 -m mmpbsa.in
```

### Protein-ligand-cofactor system
Residue names for ligands/cofactors are listed in `md_files/md_run/protein-ligand/all_ligand_resid.txt`.
```bash
# Include a cofactor in the protein selection
run_gbsa --wdir_to_run md_files/md_run/protein_H_HIS_ligand_* --append_protein_selection MG GTP

# Treat a cofactor as the ligand
run_gbsa --wdir_to_run md_files/md_run/protein_H_HIS_ligand_* --append_protein_selection MG --ligand_id GTP
```

### Per-residue decomposition
Add a `&decomp` section to `mmpbsa.in` with the desired decomposition parameters. `run_gbsa` detects this and produces `FINAL_DECOMP_MMPBSA_<suffix>.dat` along with aggregated `GBSA_decomp_avg_<suffix>.csv` and `PBSA_decomp_avg_<suffix>.csv` when PBSA is used.

## Outputs
Each run creates (in `--wdir` or current directory):
1) `log_mmpbsa_*unique-suffix*.log` - StreaMD logging
2) `log_mmpbsa_bash_*unique-suffix*.log` - gmx_MMPBSA stdout
3) `GBSA_output_*unique-suffix*.csv` and `PBSA_output_*unique-suffix*.csv` (when applicable) - summary CSV
4) Per `--wdir_to_run` directory: `FINAL_RESULTS_MMPBSA_*unique-suffix*.csv` (and decomposition files when enabled)

See {doc}`outputs` for folder layout and other analysis artifacts.
