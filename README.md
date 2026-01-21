![StreaMD Logo](./streamd_logo.png)
# StreaMD:  a tool to perform high-throughput automated molecular dynamics simulations

StreaMD provides an end-to-end molecular dynamics workflow that takes a PDB structure as input and automatically performs system preparation, equilibration, production and continuation runs, and analysis,
producing XTC trajectories together with ready-to-use plots and CSV outputs.

## Features:  
- Run multiple simultaneous molecular dynamics simulations
- Run multiple replicas of the same system for multiple complexes in a single command
- Simulation for different systems:  
    - Protein in Water;  
    - Protein - Ligand;  
    - Protein - Cofactor (multiple);  
    - Protein - Ligand - Cofactor (multiple);  

- Simulations of boron-containing molecules using Gaussian software
- Simulations of ligand-binding metalloproteins with MCPB.py
- Distributed computing using dask library
- Running parallel simulations on multiple servers
- Extending the time of MD simulations 
- Continuing interrupted MD simulations
- Restarting interrupted MD preparation by invoking the same command
- Implemented tools for end-state free energy calculations (gmx_MMPBSA) and protein–ligand interaction analysis (ProLIF)
- Support for customized .mdp files
- Interactive trajectory convergence analysis for multiple complexes
- GPU support


## Quick start
```bash
# Create environment (choose CPU-only or GPU)
conda env create --file env.yml -n md          # or env_gpu.yml on GPU-capable hosts
conda activate md

# Install
pip install streamd
# or latest main branch
pip install git+https://github.com/ci-lab-cz/streamd.git

```
   
**Minimal protein-ligand run (1 ns)**
```
run_md -p protein.pdb -l ligand.mol --md_time 1
```
   
**Protein - multiple ligands multiple replicas runs (1 ns)**
```
run_md -p protein.pdb -l ligands.sdf --md_time 1 --replicas 3 --seed 1024
```
  
**Extend successfully finished simulations**
```
run_md --wdir_to_continue md_files/md_run/protein_H_HIS_ligand_*/ --md_time 10
```
    
**GPU-accelerated simulations**
```
run_md -p protein_HIS.pdb -l ligand.mol --md_time 1 --device gpu --ncpu 32
```

**MM-PBSA/MM-GBSA calculation support**
```
run_gbsa --wdir_to_run md_files/md_run/protein_H_HIS_ligand_1 md_files/md_run/protein_H_HIS_ligand_2 -c 128 -m mmpbsa.in
```
This functionality is based on the [gmx_MMPBSA tool](https://valdes-tresanco-ms.github.io/gmx_MMPBSA/dev)

**ProLIF (Protein-Ligand Interaction Fingerprints) Analysis**
```
run_prolif --wdir_to_run md_files/md_run/protein_H_HIS_ligand_1 md_files/md_run/protein_H_HIS_ligand_2 -c 128 -s 5
```
This functionality is based on the [ProLIF tool](https://github.com/chemosim-lab/ProLIF)


More examples can be found in the [documentation](https://streamd.readthedocs.io/)


## Documentation
https://streamd.readthedocs.io/

## License
MIT

## Ready-to-use containers (Apptainer)
Pre-built `.sif` images are available (CPU and GPU) in the [Zenodo record](https://zenodo.org/records/18176058)  
**CPU:** `apptainer run --cleanenv streamd_cpu.sif run_md --help`  
**GPU:** `apptainer run --nv --cleanenv streamd_gpu.sif run_md --help`

The provided `.sif` images are intended for Apptainer on Linux/HPC systems.
GPU usage requires an NVIDIA GPU node and launching with `--nv` and `run_md --device gpu`.

## Citation
Ivanova A, Mokshyna O, Polishchuk P.  
StreaMD: the toolkit for high-throughput molecular dynamics simulations.  
*J. Cheminf.* **2024**, 16 (1), 123.  
https://doi.org/10.1186/s13321-024-00918-w