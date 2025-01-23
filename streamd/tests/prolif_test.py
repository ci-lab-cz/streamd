import os

import pandas as pd
import pytest

from streamd.prolif.run_prolif import start

prolif_test = pytest.mark.skipif(
    "not config.getoption('--run-prolif') and not config.getoption('--run-prolif-full')",
    reason="Only run when --run-prolif is given",
)
prolif_detailed_test = pytest.mark.skipif(
    "not config.getoption('--run-prolif-full')",
    reason="Only run when --run-prolif-full is given",
)


@pytest.mark.filterwarnings('ignore::DeprecationWarning')
@prolif_test
def test_run_prolif_full_pipline(dir_with_streamd_output_for_prolif):
    wdir = dir_with_streamd_output_for_prolif

    occupancy = 0.6

    assert not os.path.isfile(os.path.join(wdir, f"finished_prolif_files_test.txt"))
    assert not os.path.isfile(os.path.join(wdir, f"plifs.csv"))
    assert not os.path.isfile(os.path.join(wdir, f"plif.png"))
    assert not os.path.isfile(os.path.join(wdir, f"plif_framemap.png"))
    assert not os.path.isfile(os.path.join(wdir, f"plif_occupancy0.6.html"))
    assert not os.path.isfile(os.path.join(wdir, f"prolif_output_test_occupancy{occupancy}.png"))
    assert not os.path.isfile(os.path.join(wdir, f"prolif_output_test.csv"))

    start(wdir_to_run=[wdir],
          wdir_output=wdir,
          tpr='md_out.tpr',
          xtc='md_fit.xtc',
          step=1,
          append_protein_selection=None,
          protein_selection='protein',
          ligand_resid='UNL',
          hostfile=None,
          ncpu=len(os.sched_getaffinity(0)),
          n_jobs=None,
          occupancy=occupancy,
          plot_width=10, plot_height=10,
          save_viz=True, unique_id='test',
          pdb=None, verbose=True)

    assert os.path.isfile(os.path.join(wdir, f"finished_prolif_files_test.txt"))
    assert os.path.isfile(os.path.join(wdir, f"plifs.csv"))
    assert os.path.isfile(os.path.join(wdir, f"plif.png"))
    assert os.path.isfile(os.path.join(wdir, f"plif_framemap.png"))
    assert os.path.isfile(os.path.join(wdir, f"plif_occupancy0.6.html"))
    assert os.path.isfile(os.path.join(wdir, f"prolif_output_test_occupancy{occupancy}.png"))
    assert os.path.isfile(os.path.join(wdir, f"prolif_output_test.csv"))

    assert os.path.getsize(os.path.join(wdir, f"finished_prolif_files_test.txt")) > 0
    assert os.path.getsize(os.path.join(wdir, f"plif.png")) > 0
    assert os.path.getsize(os.path.join(wdir, f"plif_framemap.png")) > 0
    assert os.path.getsize(os.path.join(wdir, f"prolif_output_test_occupancy{occupancy}.png")) > 0

    plifs_out = pd.read_csv(os.path.join(wdir, f"plifs.csv"), sep='\t')
    prolif_output_test_out = pd.read_csv(os.path.join(wdir, f"prolif_output_test.csv"), sep='\t')

    assert not plifs_out.empty
    assert not prolif_output_test_out.empty



@prolif_detailed_test
@pytest.mark.filterwarnings('ignore::DeprecationWarning')
def test_run_prolif_full_pipline_from_files(dir_with_streamd_output_for_prolif):
    wdir = dir_with_streamd_output_for_prolif

    occupancy = 0.6

    assert not os.path.isfile(os.path.join(wdir, f"finished_prolif_files_test.txt"))
    assert not os.path.isfile(os.path.join(wdir, f"plifs.csv"))
    assert not os.path.isfile(os.path.join(wdir, f"plif.png"))
    assert not os.path.isfile(os.path.join(wdir, f"plif_framemap.png"))
    assert not os.path.isfile(os.path.join(wdir, f"plif_occupancy0.6.html"))
    assert not os.path.isfile(os.path.join(wdir, f"prolif_output_test_occupancy{occupancy}.png"))
    assert not os.path.isfile(os.path.join(wdir, f"prolif_output_test.csv"))

    start(wdir_to_run=None,
          wdir_output=wdir,
          tpr=os.path.join(wdir, 'md_out.tpr'),
          xtc=os.path.join(wdir, 'md_fit.xtc'),
          step=1,
          append_protein_selection=None,
          protein_selection='protein',
          ligand_resid='UNL',
          hostfile=None,
          ncpu=len(os.sched_getaffinity(0)),
          n_jobs=None,
          occupancy=occupancy,
          plot_width=10, plot_height=10,
          save_viz=True, unique_id='test',
          pdb=None, verbose=True)

    assert os.path.isfile(os.path.join(wdir, f"finished_prolif_files_test.txt"))
    assert os.path.isfile(os.path.join(wdir, f"plifs.csv"))
    assert os.path.isfile(os.path.join(wdir, f"plif.png"))
    assert os.path.isfile(os.path.join(wdir, f"plif_framemap.png"))
    assert os.path.isfile(os.path.join(wdir, f"plif_occupancy0.6.html"))
    assert os.path.isfile(os.path.join(wdir, f"prolif_output_test_occupancy{occupancy}.png"))
    assert os.path.isfile(os.path.join(wdir, f"prolif_output_test.csv"))

    assert os.path.getsize(os.path.join(wdir, f"finished_prolif_files_test.txt")) > 0
    assert os.path.getsize(os.path.join(wdir, f"plif.png")) > 0
    assert os.path.getsize(os.path.join(wdir, f"plif_framemap.png")) > 0
    assert os.path.getsize(os.path.join(wdir, f"prolif_output_test_occupancy{occupancy}.png")) > 0

    plifs_out = pd.read_csv(os.path.join(wdir, f"plifs.csv"), sep='\t')
    prolif_output_test_out = pd.read_csv(os.path.join(wdir, f"prolif_output_test.csv"), sep='\t')

    assert not plifs_out.empty
    assert not prolif_output_test_out.empty