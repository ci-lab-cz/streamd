# StreaMD Documentation

Welcome to the documentation for **StreaMD**, a toolkit for high-throughput molecular dynamics simulations built to take you from system preparation to automated analysis.

## Overview
StreaMD automates preparation, equilibration, production, and analysis of biomolecular systems. It can handle proteins in water or complex protein-ligand/cofactor systems, manage replicas, and resume interrupted jobs while keeping results organized.

## Key Features
- High-throughput orchestration for protein-only, protein-ligand, and protein-cofactor simulations
- Gaussian/MCPB.py support for boron-containing ligands and metalloproteins
- Distributed execution via Dask with GPU/CPU balancing and per-node throttling
- Built-in analysis (RMSD, RMSF, radius of gyration) plus MM/PB(GB)SA and ProLIF integrations
- Restartable workflows that safely extend or continue prior runs

## Quick Start
```bash
conda env create -f env.yml -n md        # or use env_gpu.yml for GPU builds
conda activate md
pip install streamd                      # or pip install git+https://github.com/ci-lab-cz/streamd.git

# Minimal protein-ligand run (1 ns)
run_md -p protein.pdb -l ligand.mol --md_time 1 --ncpu 32
```

## Getting Help
- Repository: https://github.com/ci-lab-cz/streamd
- Issues: report problems or feature requests on GitHub
- Citation: see the **Citation** section for publication details

## Documentation Structure
- **Installation**: environment setup and package installation
- **Usage Guide**: common workflows for preparing and running simulations
- **Advanced Features**: Gaussian/MCPB.py, distributed runs, and GPU tuning
- **Data Retrieval**: folder layout, logs, and built-in analysis outputs
- **PLIF Analysis**: ProLIF-based interaction fingerprints and plotting helpers
- **Python API**: module-level reference for programmatic use
- **Customization**: configuration files and custom `.mdp` settings
- **Citation** and **Changelog**

```{toctree}
:maxdepth: 2
:caption: Contents

installation
user_guide
advanced_features
data_retrieval
plif
api
customization
citation
changelog
```
