# Installation

Set up StreaMD with the provided conda environments and install the Python package.

## Requirements
- Conda (or Mamba) to create the supplied environments
- GPU-enabled machines should use `env_gpu.yml`; CPU-only machines should use `env.yml`
- Access to Gaussian if you plan to run Gaussian-dependent workflows

## Install
```bash
# Create environment (choose one)
conda env create --file env.yml -n md
# OR
conda env create --file env_gpu.yml -n md

conda activate md

# Install StreaMD
pip install streamd
# or the latest from GitHub
pip install git+https://github.com/ci-lab-cz/streamd.git
```

## Verify
```bash
run_md -h
run_gbsa -h
run_prolif -h
```
If the commands print their help text, your installation is ready.

Next steps: see {doc}`data_preparation` for input prep and {doc}`running_md` for running simulations.
