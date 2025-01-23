import os

import pytest

from streamd.run_md import start

md_test = pytest.mark.skipif(
    "not config.getoption('--run-md') and not config.getoption('--run-md-full')",
    reason="Only run when --run-md is given",
)
md_detailed_test = pytest.mark.skipif(
    "not config.getoption('--run-md-full')",
    reason="Only run when --run-md-full is given",
)

# a temporary full run_md test to check consistence of the full pipline and correctness of each conventional output files
# more functionally detailed tests will be added later
@md_test
@pytest.mark.filterwarnings('ignore::DeprecationWarning')
def test_run_md_full_pipline(dir_with_input_for_preparation):
    wdir = dir_with_input_for_preparation

    expected_output_files = ['md_fit.xtc', 'md_out.xtc',
                             'index.ndx', 'topol.top', 'md_out.log']

    assert not os.path.isfile(os.path.join(wdir, f"finished_complexes_test.txt"))
    assert not os.path.isdir(os.path.join(wdir, 'md_files', 'md_preparation'))
    assert not os.path.isdir(os.path.join(wdir, 'md_files', 'md_run', 'protein_HIS_1ke7_LS3'))

    start(wdir=wdir ,
          protein=os.path.join(wdir, 'protein_HIS.pdb'), lfile=os.path.join(wdir,'ligand.mol'),
          system_lfile=None,
          noignh=False,
          no_dr=False,
          forcefield_name=pytest.ff,
          npt_time_ps=10,
          nvt_time_ps=10,
          mdtime_ns=0.001,
          topol=None,
          topol_itp_list=None,
          posre_list_protein=None,
          wdir_to_continue_list=None,
          deffnm='md_out',
          tpr_prev=None, cpt_prev=None, xtc_prev=None,
          ligand_list_file_prev=None, ligand_resid='UNL',
          activate_gaussian=None, gaussian_exe=None, gaussian_basis=None, gaussian_memory=None,
          metal_resnames=None, metal_charges={}, mcpbpy_cut_off=None,
          seed=120,
          steps=None,
          hostfile=None,
          ncpu=len(os.sched_getaffinity(0)),
          mdrun_per_node=1,
          compute_device='auto',
          gpu_ids=None,
          ntmpi_per_gpu=None,
          clean_previous=False,
          not_clean_backup_files=False,
          unique_id='test',
          active_site_dist=5.0,
          save_traj_without_water=False,
          mdp_dir=None, bash_log='bash.log')


    assert os.path.isfile(os.path.join(wdir, f"finished_complexes_test.txt"))
    assert os.path.isfile(os.path.join(wdir, f"rmsd_mean_std_time-ranges_test.csv"))
    assert os.path.isfile(os.path.join(wdir, f"rmsd_mean_std_time-ranges_test.html"))

    assert os.path.isdir(os.path.join(wdir,'md_files', 'md_preparation','protein'))
    assert os.path.isdir(os.path.join(wdir,'md_files', 'md_preparation','ligands'))
    assert os.path.isdir(os.path.join(wdir, 'md_files','md_run','protein_HIS_1ke7_LS3'))
    assert os.path.isdir(os.path.join(wdir, 'md_files','md_run','protein_HIS_1ke7_LS3', 'md_analysis'))


    assert os.path.getsize(os.path.join(wdir, f"finished_complexes_test.txt")) > 0
    assert os.path.getsize(os.path.join(wdir, f"rmsd_mean_std_time-ranges_test.csv")) > 0
    assert os.path.getsize(os.path.join(wdir, f"rmsd_mean_std_time-ranges_test.html")) > 0

    for f in expected_output_files:
        assert os.path.isfile(os.path.join(wdir, 'md_files', 'md_run', 'protein_HIS_1ke7_LS3', f))

    with open(os.path.join(wdir, 'md_files','md_run','protein_HIS_1ke7_LS3', 'md_out.log')) as inp:
        data = inp.readlines()
        assert 'Finished mdrun' in data[-2]

