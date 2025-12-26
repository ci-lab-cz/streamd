# Customization

Tune StreaMD behavior with configuration files and custom MD settings.

## Configuration Files
`run_md`, `run_prolif`, and `run_gbsa` accept `--config` with a YAML file. Command-line flags take precedence.

Example `config.yml` for `run_md`:
```yaml
protein: protein.pdb
ligand: ligand.mol
steps: 1 2 3
md_time: 2
seed: 1024
```
Run with:
```bash
run_md --config config.yml --ncpu 32
```

## Custom MD Parameters
- Provide `.mdp` templates via `--mdp_dir` (see **Advanced Features**).
- Override trajectory handling with `--save_traj_without_water` for lighter analysis.

## Logging
Each command writes a user-facing log and a bash log with underlying tool output. Logs live alongside your working directory and are useful when resuming runs or filing bug reports.
