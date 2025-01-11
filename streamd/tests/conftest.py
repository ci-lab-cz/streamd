import os
from tempfile import TemporaryDirectory

import pytest

import streamd

def pytest_addoption(parser):
    parser.addoption(
        "--run-preparation",
        action="store_true",
        default=False,
        help="Run slow ligand and protein preparation tests",
    )


#@pytest.fixture(scope="session")
def pytest_configure():
    pytest.streamd_directory = os.path.dirname(streamd.__file__)
    pytest.script_directory = os.path.join(pytest.streamd_directory, 'scripts')
    pytest.mdp_directory = os.path.join(pytest.script_directory, 'mdp')

    pytest.protein_file = os.path.join(pytest.streamd_directory, 'tests', 'data', 'protein_HIS.pdb')
    pytest.ligand_mol_file = os.path.join(pytest.streamd_directory, 'tests', 'data', 'ligand.mol')

    pytest.ff = 'amber99sb-ildn'
#@pytest.fixture(scope="session")
#def
