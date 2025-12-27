# Logging

StreaMD writes user-facing logs and raw tool output alongside your working directory (current directory or `--wdir`).

**run_md**
```
log_protein-fname_ligand-fname_cofactor-fname_*start-time*.log     # StreaMD status
streamd_bash_protein-fname_ligand-fname_cofactor-fname_*start-time*.log  # GROMACS/Antechamber output
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
