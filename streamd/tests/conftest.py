"""Shared pytest fixtures for StreaMD tests."""

import logging
from glob import glob
import os
import shutil
# ruff: noqa: E402
#, TemporaryDirectory

import sys
import types
import importlib


def _stub_if_missing(name: str, attrs: dict | None = None) -> None:
    """Register a lightweight stub for *name* if importing it fails.

    Parameters
    ----------
    name:
        Fully qualified module name to check.
    attrs:
        Optional mapping of attributes to set on the created stub module.
    """

    if name in sys.modules:
        return
    try:
        importlib.import_module(name)
        return
    except Exception:
        mod = types.ModuleType(name)
        if attrs:
            for key, value in attrs.items():
                setattr(mod, key, value)
        sys.modules[name] = mod


# Stub heavy optional dependencies for tests so that modules importing them
# during test collection do not fail. Only lightweight placeholders are
# provided – real modules, when installed, are left untouched.
_stub_if_missing("MDAnalysis")

# ProLIF-related stubs
_stub_if_missing("prolif")
_stub_if_missing("prolif.plotting")
_stub_if_missing("prolif.plotting.barcode", {"Barcode": type("Barcode", (), {})})
_stub_if_missing("prolif.plotting.network", {"LigNetwork": type("LigNetwork", (), {})})

_stub_if_missing("matplotlib", {"rcParams": {}})
_stub_if_missing("matplotlib.pyplot", {"ioff": lambda: None})

_stub_if_missing(
    "streamd.prolif.prolif2png", {"convertprolif2png": lambda *a, **k: None}
)
_stub_if_missing(
    "streamd.prolif.prolif_frame_map",
    {
        "create_framemap": lambda *a, **k: None,
        "convertplifbyframe2png": lambda *a, **k: None,
    },
)

# Other optional libraries
_stub_if_missing("dask")
class _DummyClient:
    def __init__(self, *a, **k):
        pass

    def forward_logging(self, *a, **k):
        pass

    def retire_workers(self, *a, **k):
        pass

    def shutdown(self):
        pass

    def scheduler_info(self):
        return {'workers': {}}

    def run(self, *a, **k):
        pass


_stub_if_missing(
    "dask.distributed",
    {
        "Client": _DummyClient,
        "SSHCluster": type("SSHCluster", (), {}),
    },
)

_stub_if_missing("MDAnalysis.analysis", {"rms": object()})

_stub_if_missing("seaborn", {"set_context": lambda *a, **k: None})
_stub_if_missing("plotly", {"__path__": []})
_stub_if_missing("plotly.offline")
_stub_if_missing("plotly.express")

_stub_if_missing("parmed", {"Structure": type("Structure", (), {})})

Mol = type(
    "Mol",
    (),
    {
        "__init__": lambda self: setattr(self, "props", {}),
        "HasProp": lambda self, name: name in self.props,
        "GetProp": lambda self, name: self.props.get(name, ""),
        "SetProp": lambda self, name, value: self.props.__setitem__(name, value),
    },
)

_stub_if_missing("rdkit")
_stub_if_missing(
    "rdkit.Chem",
    {
        "Mol": Mol,
        "MolFromMolFile": lambda *a, **k: Mol(),
        "SanitizeMol": lambda *a, **k: None,
        "SDMolSupplier": lambda *a, **k: [],
        "PropertyPickleOptions": type("PropertyPickleOptions", (), {"AllProps": 0}),
        "SetDefaultPickleProperties": lambda *a, **k: None,
    },
)
_stub_if_missing("rdkit.Chem.rdmolops")

_stub_if_missing("pandas", {"DataFrame": type("DataFrame", (), {})})
_stub_if_missing("numpy")

import pytest
from streamd.utils.utils import temporary_directory_debug
# from _pytest.mark import Mark

import streamd


logging.basicConfig(level=logging.ERROR)


def pytest_addoption(parser):
    parser.addoption("--run-preparation", action="store_true", default=False,
        help="Run slow ligand and protein preparation tests")
    parser.addoption("--run-md", action="store_true", default=False,
        help="Run md tests")
    parser.addoption("--run-md-full", action="store_true", default=False,
        help="Run detailed md tests")
    parser.addoption("--run-gbsa", action="store_true", default=False, help="Run gbsa tests")
    parser.addoption("--run-gbsa-full", action="store_true", default=False,
        help="Run detailed gbsa tests")
    parser.addoption("--run-prolif", action="store_true", default=False,
        help="Run prolif tests")
    parser.addoption("--run-prolif-full", action="store_true", default=False,
        help="Run detailed prolif tests")
    parser.addoption("--run-analysis", action="store_true", default=False,
        help="Run full analysis tests")
    parser.addoption("--not-cleanup", action="store_true", default=False,
        help="Do not clean tmp directories with test output for more detailed investigation")


# For the debug purposes tests temporary directories can support optional removal (remove=True by default, use --not_clean otherwise)
# The current version of Streamd is limited to python 3.10 and TemporaryDirectory(remove=False) supports optional remove only from python 3.12



# @pytest.fixture(scope="session")
def pytest_configure(config):
    """Configure pytest globals for test suite."""
    pytest.streamd_directory = os.path.dirname(streamd.__file__)
    pytest.script_directory = os.path.join(pytest.streamd_directory, 'scripts')
    pytest.data_directory = os.path.join(pytest.streamd_directory, 'tests', 'data')
    pytest.mdp_directory = os.path.join(pytest.script_directory, 'mdp')

    pytest.system_name = 'protein_HIS_1ke7_LS3'
    pytest.ff = 'amber99sb-ildn'
    pytest.cleanup = not config.getoption("--not-cleanup")

# provide arguments to fixture dir_with_input_for_protein_preparation
# def get_marker(request, name) -> Mark:
#     markers =
#     if len(markers) > 1:
#         pytest.fail(f"Found multiple markers for {name}")
#     return markers[0]

@pytest.fixture
def dir_with_input_for_preparation() -> str:
    """Create temporary directory with input files for preparation tests."""
    with temporary_directory_debug(remove=pytest.cleanup, suffix='_preparation') as dirname:
        shutil.copyfile(os.path.join(pytest.data_directory, 'protein_HIS.pdb'),
                        os.path.join(dirname, 'protein_HIS.pdb'))
        shutil.copyfile(os.path.join(pytest.data_directory, 'ligand.mol'),
                        os.path.join(dirname, 'ligand.mol'))
        yield dirname


@pytest.fixture
def dir_with_control_files_for_preparation() -> str:
    """Provide control files for preparation workflow."""
    with temporary_directory_debug(remove=pytest.cleanup, suffix='_control_files') as dirname:
        prot_files = ['protein_HIS.gro', 'topol.top', 'posre.itp']
        for f in prot_files:
            shutil.copyfile(os.path.join(pytest.data_directory, 'mdrun_test',
                                     'md_preparation', 'protein', 'protein_HIS', f),
                        os.path.join(dirname, f))
        ligand_files = ['1ke7_LS3.gro', '1ke7_LS3.mol2', '1ke7_LS3.top','posre_1ke7_LS3.itp', 'resid.txt']
        for f in ligand_files:
            shutil.copyfile(os.path.join(pytest.data_directory, 'mdrun_test',
                                         'md_preparation', 'ligands', '1ke7_LS3', f),
                            os.path.join(dirname, f))
        yield dirname


@pytest.fixture
def dir_with_streamd_output_for_analysis() -> str:
    """Temporary directory mimicking MD run output for analysis tests."""
    with temporary_directory_debug(remove=pytest.cleanup, suffix='_analysis') as dirname:
        for file_name in ['md_out.xtc', 'md_out.tpr', 'index.ndx', 'solv_ions.gro']:
            shutil.copyfile(os.path.join(pytest.data_directory, 'mdrun_test',
                                         'md_run', pytest.system_name, file_name),
                            os.path.join(dirname, file_name))
        yield dirname


@pytest.fixture
def dir_and_rmsd_files() -> (str, str):
    """Prepare directory with RMSD CSV files for analysis tests."""
    rmsd_file_list = glob(os.path.join(pytest.data_directory, 'rmsd_files', '*'))
    with temporary_directory_debug(remove=pytest.cleanup, suffix='_rmsd-analysis') as dirname:
        tmp_rmsd_file_list = []
        for rmsd_file in rmsd_file_list:
            tmp_rmsd_file = os.path.join(dirname, os.path.basename(rmsd_file))
            shutil.copyfile(rmsd_file, tmp_rmsd_file)
            tmp_rmsd_file_list.append(tmp_rmsd_file)
        yield dirname, tmp_rmsd_file_list


@pytest.fixture
def tmp_experimental_file_for_html_paintby(dir_and_rmsd_files) -> str:
    """Return a temporary experimental data file for HTML painting tests."""
    dirname, _ = dir_and_rmsd_files
    paint_file = os.path.join(pytest.data_directory, 'paint_by_file_exp.csv')
    tmp_paint_file = os.path.join(dirname, 'paint_by_file_exp.csv')
    shutil.copyfile(paint_file, tmp_paint_file)
    return tmp_paint_file


@pytest.fixture
def dir_with_streamd_output_for_gbsa() -> str:
    """Create directory with prepared files for GBSA tests."""
    with temporary_directory_debug(remove=pytest.cleanup, suffix='_gbsa') as dirname:
        for file_name in ['md_fit.xtc', 'md_out.tpr', 'index.ndx', 'topol.top',
                          'all.itp', '1ke7_LS3.itp']:
            shutil.copyfile(os.path.join(pytest.data_directory, 'mdrun_test',
                                         'md_run', pytest.system_name, file_name),
                            os.path.join(dirname, file_name))
        shutil.copyfile(os.path.join(pytest.data_directory, 'mmpbsa.in'),
                        os.path.join(dirname, 'mmpbsa.in'))
        yield dirname


@pytest.fixture
def dir_with_streamd_output_for_prolif() -> str:
    """Create directory with prepared files for ProLIF tests."""
    with temporary_directory_debug(remove=pytest.cleanup, suffix='_prolif') as dirname:
        for file_name in ['md_fit.xtc', 'md_out.tpr']:
            shutil.copyfile(os.path.join(pytest.data_directory, 'mdrun_test',
                                         'md_run', pytest.system_name, file_name),
                            os.path.join(dirname, file_name))
        yield dirname

