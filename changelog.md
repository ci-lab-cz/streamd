**0.6**
- add argument `ligand_forcefield` - AmberTools parameter set (`gaff` or `gaff2`) for standard organic ligand parameterization (default `gaff`). 
- add argument `salt_concentration` - salt concentration in mol/L passed to `gmx genion -conc` (default: not set, only charge neutralization is performed)
- add argument `ion_pname` - positive ion name passed to `gmx genion -pname` (default `NA`)
- add argument `ion_nname` - negative ion name passed to `gmx genion -nname` (default `CL`)
- add argument `water_model` - water model passed to `gmx pdb2gmx -water` (default `tip3p`)
- update GROMACS to 2025.2 and MDAnalysis to 2.10.0
- change `ActiveSite` RMSD definition: the pocket is now the complete backbone of protein residues within `active_site_dist` of the ligand heavy atoms (previously only backbone atoms lying within the cutoff), so values are not directly comparable to earlier versions
- RMSF is now computed on Cα atoms (`gmx rmsf` selects `C-alpha`, used for both the fit and the per-residue output) instead of `Protein`
- protein RMSD is reported as both `backbone` (N/CA/C/O, the backbone superposition/fit group) and `CA` (the protein Cα subset)
- fix `last_frame.pdb`: the final frame was selected by passing the number of frames to `gmx trjconv -dump`, which interprets the value as a time in ps and extracted an early frame instead of the last one; the actual final trajectory time is now parsed from `gmx check` and passed to `-dump` (a clear error is raised if it cannot be determined)
- `last_frame.pdb` is now written as a visualization-ready snapshot: molecules are made whole across periodic boundaries and the protein–ligand–cofactor complex is centered in a compact unit cell (`gmx trjconv -pbc mol -center -ur compact`), with the full `System` (water and ions included) written
- add `start_frame.pdb` - the first production frame written with the same visualization-ready treatment as `last_frame.pdb` (`-pbc mol -center -ur compact`, full `System`), for visual comparison of the start and end of the simulation
- improve RMSD plot readability: move the legend outside the axes (below, left-aligned) so it no longer overlaps the curves, and draw lines with transparency (new `alpha` argument to `plot_rmsd`, default `0.7`) so overlapping RMSD curves remain visible
- `run_prolif` major speed-up: restrict the ProLIF protein to binding-site residues (union of residues within the cutoff of the ligand over the trajectory) instead of converting the whole protein to RDKit every frame; results are identical (~15–20× faster for large solvated systems). New argument `binding_site_cutoff` (default `12` Å; set `0` to use the full protein)
- `run_prolif` add argument `parallel_strategy` (default `chunk`): ProLIF auto-selects `queue` for solvated systems
- `run_prolif` make `--water_bridge` usable: restrict waters to an updating per-frame selection near the ligand instead of converting every water every frame. New argument `water_cutoff` (default `8` Å, automatically widened for higher `water_bridge_order`; set `0` to use all waters)
- `run_prolif` add argument `ligand_sdf`: optional reference `.sdf`/`.mol` bond-order template for the ligand (usually unnecessary — ligand chemistry is inferred from the topology and the interpreted SMILES is logged for verification)
- `run_prolif`/`prolif_drawmap` show the occupancy percentage above each dot in `prolif_output_occupancyX.png` by default (use `--no-show_percentage` to disable)

**0.5**
- add argument `box_type` - simulation box type (`triclinic`, `cubic`, `dodecahedron`, `octahedron`) (default `cubic`) defined using `gmx editconf -bt`
- add argument `box_padding_nm` - minimum solute-to-box edge distance (default `1 nm = 10 A`) defined using `gmx editconf -d`

**0.4.1** 
- added documentation
- fixed bug with --no_dr and gaussian run

**0.4**
- allow to keep all gmx_MMPBSA temporary files for debugging or for `gmx_MMPBSA_ana` usage (`--debug argument`)
- allow to provide arguments through config.yml for run_md, run_prolif, run_gbsa
- fix minor bug with steps - allow to run run_md steps 2,3,4 without --wdir_to_continue if step 1 is also provided
- added gmx_MMPBSA residue decomposition analysis feature
- add support for running replicate simulations of the same complexes
- save the last frame into `last_frame.pdb`

**0.3**
- Set default temperature to 300K
- Add Position Restraints for MCPB.py pipline
- Logo
- Prolif draw_map point_size argument
- Prolif and GBSA added directory column
- Change nvt/npt/md time and seed in user-defined mdp files only if the `--seed`, `--nvt_time`, `--npt_time` and `--md_time` arguments were explicitly set up in the StreaMD arguments list otherwise user mdp parameters will be used

**0.2.9**
- Prolif n_jobs automatic calculation
- Fixed bug with interrupted continuation runs
- Fixed bug with protein only in water simulations analysis
- Add directory information into rmsd output files for replicate runs
- Separate analysis output subdirectory
- Tests added

**0.2.8**
- Added n_jobs and protein_selection arguments to run_prolif
- Added save_traj_without_water argument to run_md

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

**0.2.6**
- Added GPU support
- Use only available CPUs
- Fixed logging
- Added mdrun_per_node argument
