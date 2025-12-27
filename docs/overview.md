# Overview

StreaMD is a toolkit to run high-throughput molecular dynamics simulations end-to-end: preparation, equilibration, production, continuation, and analysis.

## Description
- Fully automated pipeline for molecular dynamics
- Supports proteins in water, protein-ligand, protein-cofactor, and combined protein-ligand-cofactor systems
- Gaussian and MCPB.py integration for boron-containing compounds and metalloproteins
- Distributed computing via Dask with CPU/GPU balancing
- Multiple replicas in a single command
- Restartable: extend or continue interrupted runs safely
- Built-in analysis (RMSD, RMSF, radius of gyration) plus MM/PB(GB)SA and ProLIF utilities
- GPU support and custom `.mdp` files

## Features
- Multiple simultaneous simulations
- Launch multiple replicas of the same system through a single command
- Simulation for different systems:  
    - Protein in Water;  
    - Protein - Ligand;  
    - Protein - Cofactor (multiple);  
    - Protein - Ligand - Cofactor (multiple);

- Support for custom `.mdp` files
- Boron-containing molecule support through Gaussian
- MCPB.py workflow for ligand-binding metalloproteins
- Distributed computing with Dask and multi-server support
- Extend or continue MD simulations and preparation stages
- End-state free energy (gmx_MMPBSA) and ProLIF interaction analysis
- Interactive trajectory convergence analysis
- GPU support with control over GPUs and thread-MPI ranks

    
## Next Steps
- Set up the environment: see {doc}`installation`
- Prepare inputs: see {doc}`data_preparation`
- Run simulations: see {doc}`running_md`
- Explore outputs and analysis: see {doc}`outputs`, {doc}`mm_pbsa`, {doc}`plif`, and {doc}`trajectory_convergence`
