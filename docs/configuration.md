# Configuration

`run_md`, `run_prolif`, and `run_gbsa` accept a `--config` option that supplies default arguments from a YAML file. Command-line arguments override YAML values.

Example `config.yml` for `run_md`:
```yaml
protein: protein.pdb
ligand: ligand.mol
steps: 3 4
md_time: 2
seed: 1024
```

Run using the configuration:
```bash
run_md --config config.yml --ncpu 32
```

Notes:
- Use CLI flags to override specific settings.
- Combine with `--steps` to skip preparation or analysis when continuing an existing run (see {doc}`running_md`).
