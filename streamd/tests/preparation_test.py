import os
from tempfile import TemporaryDirectory
import pytest

from streamd.utils.utils import run_check_subprocess, get_mol_resid_pair
from streamd.preparation.ligand_preparation import prepare_input_ligands

preparation_test = pytest.mark.skipif(
    "not config.getoption('--run-preparation')",
    reason="Only run when --run-preparation is given",
)

@preparation_test
def test_prepare_protein_gro():
    with TemporaryDirectory(dir=os.path.curdir) as dirname:
        cmd = f'gmx pdb2gmx -f {pytest.protein_file} -o {os.path.join(dirname, "protein.gro")} ' \
               f'-water tip3p -ignh ' \
               f'-i {os.path.join(dirname, "posre.itp")} ' \
               f'-p {os.path.join(dirname, "topol.top")} -ff {pytest.ff}'

        assert run_check_subprocess(cmd, env=os.environ.copy(), ignore_error=True) is not None
        assert os.path.isfile(f'{os.path.join(dirname, "protein.gro")}')
        assert os.path.isfile(f'{os.path.join(dirname, "topol.top")}')


#@pytest.mark.skip(reason="Test directly with -k test_prepare_ligand. Takes some time")
@preparation_test
def test_prepare_ligand():
    with TemporaryDirectory(dir=os.path.curdir) as dirname:
        # list with prepared directory or empty list
        result = prepare_input_ligands(ligand_fname=pytest.ligand_mol_file,
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
        assert os.path.isfile(f'{os.path.join(result[0], "resid.txt")}')

        molid, resid = next(get_mol_resid_pair(f'{os.path.join(result[0], "resid.txt")}'))

        assert molid == '1ke7_LS3'
        assert os.path.isfile(f'{os.path.join(result[0], f"{molid}.itp")}')
        assert os.path.isfile(f'{os.path.join(result[0], f"{molid}.mol2")}')
