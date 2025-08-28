**0.2.6**
- Added GPU support
- Use only available CPUs
- Fixed logging
- Added mdrun_per_node argument

**0.2.7**
- Fixed an issue with ligand restraints during NPT and NVT
- RMSD of ligand, backbone, Active Site groups are calculated by MDanalysis and saved in 1 common output file
- Use xtc and tpr files containing no water for RMSD analysis (to prevent memory issues)
- Added additional trajectory convergence analysis (html interactive files requires plotly)
- Added plotly dependency in env.yml
- Transform xvg analysis files to csv
- Save successfully finished names of the systems to the text output file
- Calculate MMPBSA intermediate files in tmp directories for multiple run in one directory
- Add unique suffix for run_md and run_gbsa output files for simultaneous runs in the same working directory

**0.2.8**
- Added n_jobs and protein_selection arguments to run_prolif
- Added save_traj_without_water argument to run_md

**0.2.9**
- Prolif n_jobs automatic calculation
- Fixed bug with interrupted continuation runs
- Fixed bug with protein only in water simulations analysis
- Add directory information into rmsd output files for replicate runs
- Separate analysis output subdirectory
- Tests added

**0.3**
- Set default temperature to 300K
- Add Position Restraints for MCPB.py pipline
- Logo
- Prolif draw_map point_size argument
- Prolif and GBSA added directory column
- Change nvt/npt/md time and seed in user-defined mdp files only if the `--seed`, `--nvt_time`, `--npt_time` and `--md_time` arguments were explicitly set up in the StreaMD arguments list otherwise user mdp parameters will be used

**future version 0.3.1 (the latest version from github)**
- allow to keep all gmx_MMPBSA temporary files for debugging or for `gmx_MMPBSA_ana` usage (`--debug argument`)
- allow to provide arguments through config.yml for run_md, run_prolif, run_gbsa
- fix minor bug with steps - allow to run run_md steps 2,3,4 without --wdir_to_continue if step 1 is also provided
- added gmx_MMPBSA residue decomposition analysis feature
- allow specifying solvation box shape and edge distance from the CLI
