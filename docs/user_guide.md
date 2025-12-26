# Usage Guide

Use StreaMD to prepare biomolecular systems, launch simulations, and recover from interruptions. The examples below mirror the most common workflows.

## Prepare Inputs
- Clean proteins before running: fill missing loops, remove crystallographic ligands, and assign correct protonation (e.g., rename histidines to HID/HIE/HIP).
- StreaMD re-adds hydrogens by default (`pdb2gmx -ignh`). Use `--noignh` only if you are sure atom names match the chosen force field.
- Align ligand/cofactor coordinates to the protein frame; multi-ligand runs expect all poses to be aligned.

## Standard Runs
Activate your environment (`conda activate md`) before running the commands.

### Protein in Water
```bash
run_md -p protein_H_HIS.pdb --md_time 1 --nvt_time 1000 --npt_time 1000 --ncpu 32
```

### Protein-Ligand
```bash
# Single ligand
run_md -p protein_H_HIS.pdb -l ligand.mol --md_time 1 --ncpu 32

# Multiple aligned ligands in one file
run_md -p protein_H_HIS.pdb -l ligands.sdf --md_time 1 --ncpu 32
```

### Protein-Cofactor
All cofactors must be present; preparation errors stop the run.
```bash
run_md -p protein_H_HIS.pdb --cofactor cofactors.sdf --md_time 1 --ncpu 32
```

### Using a Configuration File
CLI flags override values from YAML.
```bash
run_md --config config.yml --ncpu 32
```

## Replicas
Prepare once, then launch independent replicas with incremented seeds.
```bash
run_md -p protein_H_HIS.pdb -l ligand.mol --md_time 1 --replicas 3 --seed 1024
```
Re-running with a larger `--replicas` count adds only missing replicas.

## Distributed Runs
Provide a hostfile (first line is the Dask scheduler).
```bash
# PBS
run_md -p protein_H_HIS.pdb -l molecules.sdf --cofactor cofactors.sdf --md_time 1 --hostfile "$PBS_NODEFILE" --ncpu 128

# SLURM
srun hostname | sort | uniq > hostfile
run_md -p protein_H_HIS.pdb -l molecules.sdf --md_time 1 --hostfile hostfile --ncpu 128
```

## Continue or Extend
StreaMD resumes automatically if checkpoint files exist. To extend a finished run, point to the prior directory or explicit files.
```bash
# Continue or extend a StreaMD run
run_md --wdir_to_continue md_files/md_run/protein_H_HIS_ligand_1 --md_time 2

# Continue a non-StreaMD run
run_md --wdir mdrun --md_time 3 --tpr md_out.tpr --cpt md_out.cpt --xtc md_out.xtc --ligand_list_file all_ligand_resid.txt

# Skip preparation/analysis for an existing run
run_md --wdir_to_continue md_files/md_run/protein_H_HIS_ligand_1 --md_time 3 --steps 3 4
```

## GPU Usage
- Select device: `--device gpu` (GPU), `--device cpu` (force CPU), or `--device auto` (GROMACS decides).
- Use specific GPUs: `--gpu_ids 0 1`.
- Increase ranks per GPU: `--ntmpi_per_gpu 2`.
- Run multiple simulations per node (splits CPUs): `--mdrun_per_node 2`.
- Combine multiple GPUs with multiple tasks (experimental): `--mdrun_per_node 2 --gpu_ids 0 1`.
```bash
run_md -p protein_HIS.pdb -l ligand.mol --md_time 1 --device gpu --gpu_ids 0 1 --ntmpi_per_gpu 2
```

## Output at a Glance
- Logs: `log_*` (StreaMD status) and `streamd_bash_*` (GROMACS/Antechamber output)
- Working directories: `md_files/md_preparation/` and `md_files/md_run/<system>[_replicaN]/`
- Analysis folder per run: `md_analysis/` with energy, RMSD, RMSF, and radius of gyration summaries

See **Data Retrieval** for detailed file listings.
