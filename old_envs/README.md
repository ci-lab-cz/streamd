# Archived conda environment files

Previous versions of `env.yml` and `env_gpu.yml`, kept for reproducing older
StreaMD setups. The **current** environment files live at the repository root
(`../env.yml`, `../env_gpu.yml`). Files here are named
`<name>_<date>_<commit>.yml` and sort chronologically.

To use one, e.g.:

```bash
conda env create -f old_envs/env_2024-10-14_698933c.yml
```

## `env.yml` history (oldest → newest)

| File | Commit | Date | Change |
|------|--------|------|--------|
| | `env_2024-10-14_698933c.yml` | 698933c | 2024-10-14 | add plotly for rmsd analysis |
| *(current)* `../env.yml` | e9d703b | 2026-07-08 | update packages version (newer MDAnalysis support) #85 |

## `env_gpu.yml` history (oldest → newest)

| File | Commit | Date | Change |
|------|--------|------|--------|
| `env_gpu_2024-10-14_698933c.yml` | 698933c | 2024-10-14 | add plotly for rmsd analysis |
| *(current)* `../env_gpu.yml` | e9d703b | 2026-07-08 | update packages version (newer MDAnalysis support) #85 |
