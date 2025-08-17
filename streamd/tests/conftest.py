"""Shared pytest fixtures for StreaMD tests."""

import logging
from glob import glob
import os
import shutil

import pytest
from streamd.utils.utils import temporary_directory_debug

import streamd


logging.basicConfig(level=logging.ERROR)


def pytest_addoption(parser):
    """Add custom command-line options for selective test runs."""
    parser.addoption("--run-preparation", action="store_true", default=False,
                     help="Run slow ligand and protein preparation tests")
    parser.addoption("--run-md", action="store_true", default=False,
                     help="Run md tests")
    parser.addoption("--run-md-full", action="store_true", default=False,
                     help="Run detailed md tests")
    parser.addoption("--run-gbsa", action="store_true", default=False,
                     help="Run gbsa tests")
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


def pytest_configure(config):
    """Configure pytest globals for test suite."""
    pytest.streamd_directory = os.path.dirname(streamd.__file__)
    pytest.script_directory = os.path.join(pytest.streamd_directory, 'scripts')
    pytest.data_directory = os.path.join(pytest.streamd_directory, 'tests', 'data')
    pytest.mdp_directory = os.path.join(pytest.script_directory, 'mdp')

    pytest.system_name = 'protein_HIS_1ke7_LS3'
    pytest.ff = 'amber99sb-ildn'
    pytest.cleanup = not config.getoption("--not-cleanup")


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

