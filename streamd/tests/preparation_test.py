import os

import pytest

from streamd.utils.utils import run_check_subprocess, get_mol_resid_pair
from streamd.preparation.ligand_preparation import prepare_input_ligands

preparation_test = pytest.mark.skipif(
    "not config.getoption('--run-preparation')",
    reason="Only run when --run-preparation is given",
)


@preparation_test
def test_prepare_protein_gro(dir_with_input_for_preparation, dir_with_control_files_for_preparation):
    dirname = dir_with_input_for_preparation

    assert not os.path.isfile(f'{os.path.join(dirname, "protein.gro")}')
    assert not os.path.isfile(f'{os.path.join(dirname, "topol.top")}')

    cmd = f'gmx pdb2gmx -f {os.path.join(dirname, "protein_HIS.pdb")} -o {os.path.join(dirname, "protein.gro")} ' \
          f'-water tip3p -ignh ' \
          f'-i {os.path.join(dirname, "posre.itp")} ' \
          f'-p {os.path.join(dirname, "topol.top")} -ff {pytest.ff}'

    assert run_check_subprocess(cmd, env=os.environ.copy(), ignore_error=True) is not None
    
    assert os.path.isfile(f'{os.path.join(dirname, "protein.gro")}')
    assert os.path.isfile(f'{os.path.join(dirname, "topol.top")}')

    # check if output is identical




# @pytest.mark.skip(reason="Test directly with -k test_prepare_ligand. Takes some time")
@preparation_test
def test_prepare_ligand(dir_with_input_for_preparation):
    expected_mol_id = '1ke7_LS3'
    dirname = dir_with_input_for_preparation

    assert not os.path.isfile(os.path.join(dirname, expected_mol_id, "resid.txt"))
    assert not os.path.isfile(os.path.join(dirname, expected_mol_id, f"{expected_mol_id}.itp"))
    assert not os.path.isfile(os.path.join(dirname, expected_mol_id, f"{expected_mol_id}.mol2"))
    assert not os.path.isfile(os.path.join(dirname, expected_mol_id, f"{expected_mol_id}.gro"))

    # list with prepared directory or empty list
    result = prepare_input_ligands(ligand_fname=os.path.join(dirname, 'ligand.mol'),
                                   preset_resid='UNL',
                                   protein_resid_set={},
                                   script_path=pytest.script_directory,
                                   project_dir=pytest.streamd_directory,
                                   wdir_ligand=dirname,
                                   no_dr=False,
                                   gaussian_exe=None,
                                   activate_gaussian=None,
                                   gaussian_basis=None,
                                   gaussian_memory=None,
                                   hostfile=None, ncpu=1,
                                   bash_log='test')

    assert result
    assert result[0] == os.path.join(dirname, expected_mol_id)

    assert os.path.isfile(f'{os.path.join(dirname, expected_mol_id, "resid.txt")}')
    molid, resid = next(get_mol_resid_pair(f'{os.path.join(dirname, expected_mol_id, "resid.txt")}'))
    assert molid == expected_mol_id

    assert os.path.isfile(f'{os.path.join(dirname, expected_mol_id, f"{molid}.itp")}')
    assert os.path.isfile(f'{os.path.join(dirname, expected_mol_id, f"{molid}.mol2")}')
