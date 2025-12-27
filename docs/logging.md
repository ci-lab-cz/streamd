# Logging

StreaMD writes user-facing logs and raw tool output alongside your working directory (current directory or `--wdir`).

**run_md**
```
log_protein-fname_ligand-fname_cofactor-fname_*start-time*.log     # StreaMD status
streamd_bash_protein-fname_ligand-fname_cofactor-fname_*start-time*.log  # GROMACS/Antechamber output

Detailed logs for each ligand/cofactor preparation runs
- md_preparation/ligands/ligandN/streamd_bash_*protein-fname*_*ligand-fname*_*cofactor-fname*_*start-time*.log`: stdout for Antechamber or Gaussian ligand preparation run 
- md_preparation/cofactors/cofactorN/streamd_bash_*protein-fname*_*ligand-fname*_*cofactor-fname*_*start-time*.log`: stdout for Antechamber or Gaussian cofactor preparation run 

```

**run_gbsa**
```
log_mmpbsa_*unique-suffix*.log          # StreaMD status
log_mmpbsa_bash_*unique-suffix*.log     # gmx_MMPBSA output
```

**run_prolif**
```
log_prolif_*unique-suffix*.log          # StreaMD status
```

More detailed information about the logs can be found in the {doc}`outputs` file.
