import os

import pandas as pd
import pytest


analysis_test = pytest.mark.skipif(
    "not config.getoption('--run-analysis')",
    reason="Only run when --run-analysis is given",
)


@analysis_test
@pytest.mark.filterwarnings('ignore::DeprecationWarning')
def test_rmsd_analysis(dir_with_streamd_output_for_analysis,
                       # list_expected_system_analysis_output,
                       # list_expected_analysis_output
                       ):
    # import here to avoid bunch of not functional warnings from matplotlib
    from streamd.analysis.md_system_analysis import run_md_analysis

    list_expected_analysis_output = [
        f'gyrate_{pytest.system_name}.xvg',
        f'gyrate_{pytest.system_name}.png',
        f'rmsf_{pytest.system_name}.xvg',
        f'rmsf_{pytest.system_name}.png',
        f'rmsf_{pytest.system_name}.pdb'
    ]

    list_expected_system_analysis_output = [
        'md_fit.xtc', 'md_fit_nowater.xtc', 'md_short_forcheck.xtc',
        'frame.pdb'
    ]

    dir_with_streamd_output_files = dir_with_streamd_output_for_analysis
    md_analysis_dir = os.path.join(dir_with_streamd_output_files, 'md_analysis')
    rmsd_file = os.path.join(md_analysis_dir, f'rmsd_{pytest.system_name}.csv')

    expected_output = (
        rmsd_file,
        md_analysis_dir,
        dir_with_streamd_output_files
    )
    expected_output_system_analysis_files = [
        os.path.join(dir_with_streamd_output_files, i) for i in list_expected_system_analysis_output
    ]
    expected_output_analysis_files = [
        os.path.join(dir_with_streamd_output_files, 'md_analysis', i) for i in list_expected_analysis_output
    ]

    for i in expected_output_system_analysis_files:
        assert not os.path.isfile(i)
    for i in expected_output_analysis_files:
        assert not os.path.isfile(i)

    res = run_md_analysis((dir_with_streamd_output_files, 'md_out'),
                          mdtime_ns=0.05,
                          project_dir=pytest.streamd_directory,
                          bash_log='bash.log',
                          active_site_dist=5.0, ligand_resid='UNL',
                          save_traj_without_water=True,
                          analysis_dirname='md_analysis',
                          ligand_list_file_prev=None, env=None,
                          system_name=pytest.system_name)

    assert res == expected_output

    assert os.path.isdir(md_analysis_dir)
    assert os.path.isfile(rmsd_file)

    for i in expected_output_system_analysis_files:
        assert os.path.isfile(i)
    for i in expected_output_analysis_files:
        assert os.path.isfile(i)


@analysis_test
@pytest.mark.filterwarnings('ignore::DeprecationWarning')
# @pytest.mark.skip(reason="Ignore")
@pytest.mark.parametrize("rmsd_type_list", [
    pytest.param(['backbone', 'ligand', f'ActiveSite5.0A'], id="backbone, ligand, ActiveSite5.0A"),
    pytest.param(['backbone'], id="backbone"),
    pytest.param(['ActiveSite5.0A'], id="ActiveSite5.0A"),
])
def test_rmsd_analysis(rmsd_type_list, dir_and_rmsd_files):
    wdir, rmsd_file_list = dir_and_rmsd_files
    rmsd_expected_output = os.path.join(wdir, f'rmsd_all_systems_test.csv')
    rmsd_expected_output_ranges = os.path.join(wdir, 'rmsd_mean_std_time-ranges_test.csv')
    rmsd_expected_output_html = os.path.join(wdir, 'rmsd_mean_std_time-ranges_test.html')

    from streamd.analysis.run_analysis import run_rmsd_analysis

    assert not os.path.isfile(rmsd_expected_output)
    assert not os.path.isfile(rmsd_expected_output_ranges)
    assert not os.path.isfile(rmsd_expected_output_html)

    # test rmsd file with only protein rmsd available (protein only in water case)
    if rmsd_type_list == ['backbone']:
        for rmsd_file in rmsd_file_list:
            data = pd.read_csv(rmsd_file, sep='\t')
            data.loc[:, 'ligand_name'] = None
            data.to_csv(rmsd_file, sep='\t', index=False)
        del data

    run_rmsd_analysis(rmsd_file_list,
                      wdir=wdir,
                      unique_id='test',
                      time_ranges=None,
                      rmsd_type_list=rmsd_type_list,
                      paint_by_fname=None,
                      title=None)

    assert os.path.isfile(rmsd_expected_output)
    assert os.path.isfile(rmsd_expected_output_ranges)
    assert os.path.isfile(rmsd_expected_output_html)

    for csv_file in [rmsd_expected_output, rmsd_expected_output_ranges]:
        data = pd.read_csv(csv_file, sep='\t')
        assert not data.empty

    assert os.path.getsize(rmsd_expected_output_html) > 0

@analysis_test
@pytest.mark.filterwarnings('ignore::DeprecationWarning')
def test_rmsd_analysis_html_paintby(dir_and_rmsd_files, tmp_experimental_file_for_html_paintby):
    wdir, rmsd_file_list = dir_and_rmsd_files
    paint_by_file = tmp_experimental_file_for_html_paintby

    assert os.path.isfile(paint_by_file)

    rmsd_expected_output = os.path.join(wdir, f'rmsd_all_systems_test.csv')
    rmsd_expected_output_ranges = os.path.join(wdir, 'rmsd_mean_std_time-ranges_test.csv')
    rmsd_expected_output_html = os.path.join(wdir, 'rmsd_mean_std_time-ranges_test.html')

    assert not os.path.isfile(rmsd_expected_output)
    assert not os.path.isfile(rmsd_expected_output_ranges)
    assert not os.path.isfile(rmsd_expected_output_html)

    from streamd.analysis.run_analysis import run_rmsd_analysis

    run_rmsd_analysis(rmsd_file_list,
                      wdir=wdir,
                      unique_id='test',
                      time_ranges=None,
                      rmsd_type_list=['backbone', 'ligand'],
                      paint_by_fname=paint_by_file,
                      title=None)

    assert os.path.isfile(rmsd_expected_output)
    assert os.path.isfile(rmsd_expected_output_ranges)
    assert os.path.isfile(rmsd_expected_output_html)

    for csv_file in [rmsd_expected_output, rmsd_expected_output_ranges]:
        data = pd.read_csv(csv_file, sep='\t')
        assert not data.empty

    assert os.path.getsize(rmsd_expected_output_html) > 0

    # check experimental column in html
    line_pKi_bar = '"coloraxis":{"colorbar":{"title":{"text":"pKi"}}'
    with open(rmsd_expected_output_html) as html_file:
        html_text = html_file.read()

    assert line_pKi_bar in html_text

