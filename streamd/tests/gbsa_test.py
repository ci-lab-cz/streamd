import os

import pandas as pd
import pytest

from streamd.run_gbsa import start, get_number_of_frames, get_used_number_of_frames, get_mmpbsa_start_end_interval
#,get_mmpbsa_start_end_interval, run_get_frames_from_wdir, \
#    parse_gmxMMPBSA_output, run_gbsa_task)

gbsa_full_test = pytest.mark.skipif(
    "not config.getoption('--run-gbsa') and not config.getoption('--run-gbsa-full')",
    reason="Only run when --run-gbsa is given",
)

gbsa_detailed_test = pytest.mark.skipif(
    "not config.getoption('--run-gbsa-full')",
    reason="Only run when --run-gbsa-full is given",
)


@gbsa_full_test
def test_get_frames(dir_with_streamd_output_for_gbsa):
    wdir = dir_with_streamd_output_for_gbsa

    startframe, endframe, interval = get_mmpbsa_start_end_interval(os.path.join(wdir, 'mmpbsa.in'))
    assert startframe == 1
    assert endframe == 100
    assert interval == 1

    res = get_number_of_frames(os.path.join(wdir, 'md_fit.xtc' ), env=os.environ.copy())
    assert res
    num_frames, timestep = list(map(int, res))
    assert num_frames == 6
    assert timestep == 10

    used_number_of_frames = get_used_number_of_frames(var_number_of_frames=[num_frames],
                                                      startframe=startframe,
                                                      endframe=endframe,
                                                      interval=interval)
    assert used_number_of_frames == 6


@gbsa_full_test
def test_run_gbsa_full_pipline(dir_with_streamd_output_for_gbsa):
    wdir = dir_with_streamd_output_for_gbsa

    assert not os.path.isfile( os.path.join(wdir, f"finished_gbsa_files_test.txt"))
    assert not os.path.isfile( os.path.join(wdir, f"FINAL_RESULTS_MMPBSA_test.dat"))
    assert not os.path.isfile( os.path.join(wdir, f"FINAL_RESULTS_MMPBSA_test.csv"))
    assert not os.path.isfile(os.path.join(wdir, f'GBSA_output_test.csv'))
    assert not os.path.isfile(os.path.join(wdir, f'PBSA_output_test.csv'))

    start(wdir_to_run=[wdir],
          tpr='md_out.tpr',
          xtc='md_fit.xtc',
          topol='topol.top',
          index='index.ndx',
          out_wdir=wdir,
          mmpbsa=os.path.join(wdir, 'mmpbsa.in'),
          ncpu=len(os.sched_getaffinity(0)),
          ligand_resid='UNL',
          append_protein_selection=None,
          hostfile=None,
          unique_id='test',
          bash_log='bash.log',
          gmxmmpbsa_out_files=None, clean_previous=False)

    assert os.path.isfile(os.path.join(wdir, f"finished_gbsa_files_test.txt"))
    assert os.path.isfile(os.path.join(wdir, f"FINAL_RESULTS_MMPBSA_test.dat"))
    assert os.path.isfile(os.path.join(wdir, f"FINAL_RESULTS_MMPBSA_test.csv"))
    assert os.path.isfile(os.path.join(wdir, f'GBSA_output_test.csv'))
    assert os.path.isfile(os.path.join(wdir, f'PBSA_output_test.csv'))

    assert os.path.getsize(os.path.join(wdir, f"finished_gbsa_files_test.txt")) > 0

    gbsa_out = pd.read_csv(os.path.join(wdir, f'GBSA_output_test.csv'), sep='\t')
    pbsa_out = pd.read_csv(os.path.join(wdir, f'PBSA_output_test.csv'), sep='\t')

    assert not gbsa_out.empty
    assert not pbsa_out.empty


@gbsa_detailed_test
def test_run_gbsa_full_pipline_from_files(dir_with_streamd_output_for_gbsa):
    wdir = dir_with_streamd_output_for_gbsa

    assert not os.path.isfile( os.path.join(wdir, f"finished_gbsa_files_test.txt"))
    assert not os.path.isfile( os.path.join(wdir, f"FINAL_RESULTS_MMPBSA_test.dat"))
    assert not os.path.isfile( os.path.join(wdir, f"FINAL_RESULTS_MMPBSA_test.csv"))
    assert not os.path.isfile(os.path.join(wdir, f'GBSA_output_test.csv'))
    assert not os.path.isfile(os.path.join(wdir, f'PBSA_output_test.csv'))

    start(wdir_to_run=None,
              tpr=os.path.join(wdir, 'md_out.tpr'),
              xtc=os.path.join(wdir, 'md_fit.xtc'),
              topol=os.path.join(wdir, 'topol.top'),
              index=os.path.join(wdir, 'index.ndx'),
              out_wdir=wdir,
              mmpbsa=os.path.join(wdir, 'mmpbsa.in'),
              ncpu=len(os.sched_getaffinity(0)),
              ligand_resid='UNL',
              append_protein_selection=None,
              hostfile=None,
              unique_id='test',
              bash_log='bash.log',
              gmxmmpbsa_out_files=None, clean_previous=False)

    assert os.path.isfile(os.path.join(wdir, f"finished_gbsa_files_test.txt"))
    assert os.path.isfile(os.path.join(wdir, f"FINAL_RESULTS_MMPBSA_test.dat"))
    assert os.path.isfile(os.path.join(wdir, f"FINAL_RESULTS_MMPBSA_test.csv"))
    assert os.path.isfile(os.path.join(wdir, f'GBSA_output_test.csv'))
    assert os.path.isfile(os.path.join(wdir, f'PBSA_output_test.csv'))

    assert os.path.getsize(os.path.join(wdir, f"finished_gbsa_files_test.txt")) > 0

    gbsa_out = pd.read_csv(os.path.join(wdir, f'GBSA_output_test.csv'), sep='\t')
    pbsa_out = pd.read_csv(os.path.join(wdir, f'PBSA_output_test.csv'), sep='\t')

    assert not gbsa_out.empty
    assert not pbsa_out.empty
