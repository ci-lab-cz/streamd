"""Tests for analysis routines."""

import os
import sys
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest


analysis_test = pytest.mark.skipif(
    "not config.getoption('--run-analysis')",
    reason="Only run when --run-analysis is given",
)


def _make_test_universe(frames):
    """Create a small MDAnalysis universe with protein and ligand test atoms."""
    mda = pytest.importorskip("MDAnalysis")
    if not hasattr(mda, "Universe"):
        pytest.skip("MDAnalysis is stubbed")

    atom_resindex = [0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2]
    universe = mda.Universe.empty(
        11,
        n_residues=3,
        atom_resindex=atom_resindex,
        trajectory=True,
    )
    universe.add_TopologyAttr('name', ['N', 'CA', 'C', 'O', 'N', 'CA', 'C', 'O', 'C1', 'C2', 'C3'])
    universe.add_TopologyAttr('type', ['N', 'C', 'C', 'O', 'N', 'C', 'C', 'O', 'C', 'C', 'C'])
    universe.add_TopologyAttr('mass', [14.0, 12.0, 12.0, 16.0, 14.0, 12.0, 12.0, 16.0, 12.0, 12.0, 12.0])
    universe.add_TopologyAttr('resname', ['ALA', 'GLY', 'UNL'])
    universe.add_TopologyAttr('resid', [1, 2, 3])
    universe.add_TopologyAttr('segid', ['A'])
    universe.load_new(np.asarray(frames, dtype=np.float32), order='fac', dt=20)
    return universe


def _base_coords():
    """Return reference coordinates for the in-memory RMSD test system."""
    return np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [1.0, 1.0, 1.0],
        [10.0, 0.0, 0.0],
        [11.0, 0.0, 0.0],
        [11.0, 1.0, 0.0],
        [11.0, 1.0, 1.0],
        [0.2, 0.5, 0.2],
        [0.5, 0.7, 0.4],
        [0.7, 0.2, 0.3],
    ])


def _rmsd_with_local(frames):
    """Calculate standard and local ligand RMSD for in-memory test frames."""
    from streamd.analysis import md_system_analysis as analysis

    universe = _make_test_universe(frames)
    pocket_selection, ligand_selection = analysis._reference_pocket_selections(universe, 'UNL', 5.0)
    rmsd_df = analysis.rmsd_for_atomgroups(universe, selection1='backbone', selection2=[ligand_selection])
    if pocket_selection:
        rmsd_df = analysis._add_local_ligand_rmsd(rmsd_df, universe, pocket_selection, ligand_selection)
    return rmsd_df, ligand_selection, pocket_selection


def _rotate_translate(coords):
    """Apply a rigid-body rotation and translation to coordinates."""
    angle = np.deg2rad(90.0)
    rotation = np.array([
        [np.cos(angle), -np.sin(angle), 0.0],
        [np.sin(angle), np.cos(angle), 0.0],
        [0.0, 0.0, 1.0],
    ])
    return coords @ rotation.T + np.array([3.0, -2.0, 1.5])


def _write_rmsd_input(path, drop_columns=(), column_overrides=None):
    """Write a standard two-frame RMSD input table for convergence tests."""
    data = {
        'time(ns)': [0.0, 0.02],
        'backbone': [0.0, 0.0],
        'ligand': [1.0, 2.0],
        'ligand_local': [0.5, 1.5],
        'ligand_name': ['ligand_a', 'ligand_a'],
        'system': ['protein_a_replica1', 'protein_a_replica1'],
    }
    data.update(column_overrides or {})
    pd.DataFrame(data).drop(columns=list(drop_columns)).to_csv(
        path,
        sep='\t',
        index=False,
    )
    return path


def test_rmsd_for_atomgroups_uses_trajectory_time(monkeypatch):
    """RMSD time must come from MDAnalysis timestamps, not frame spacing."""
    from streamd.analysis import md_system_analysis as analysis

    class FakeTrajectory:
        def __getitem__(self, index):
            return index

    class FakeRMSD:
        def __init__(self, *args, **kwargs):
            self.results = SimpleNamespace(rmsd=None)

        def run(self):
            self.results.rmsd = np.array(
                [
                    [0.0, 0.0, 0.0, 0.0],
                    [1.0, 20.0, 1.0, 2.0],
                    [2.0, 40.0, 2.0, 4.0],
                ]
            )

    monkeypatch.setattr(analysis, 'rms', SimpleNamespace(RMSD=FakeRMSD))
    universe = SimpleNamespace(trajectory=FakeTrajectory())

    rmsd_df = analysis.rmsd_for_atomgroups(universe, selection1='backbone', selection2=['ligand'])

    assert rmsd_df['time(ns)'].tolist() == [0.0, 0.02, 0.04]
    assert 'frame' not in rmsd_df.columns
    assert rmsd_df['ligand'].tolist() == [0.0, 2.0, 4.0]


def test_rmsd_for_atomgroups_uses_actual_mdanalysis_timestamps():
    """In-memory trajectories with 20 ps spacing produce ns timestamps from MDAnalysis."""
    coords = _base_coords()
    rmsd_df, ligand_selection, _ = _rmsd_with_local([coords, coords, coords])

    assert rmsd_df['time(ns)'].tolist() == [0.0, 0.02, 0.04]
    assert rmsd_df[ligand_selection].tolist() == [0.0, 0.0, 0.0]


def test_global_and_local_ligand_rmsd_are_rigid_body_invariant():
    """Rigid-body motion of the full complex is removed by backbone and pocket fitting."""
    coords = _base_coords()
    moved = _rotate_translate(coords)

    rmsd_df, ligand_selection, _ = _rmsd_with_local([coords, moved])

    assert rmsd_df['backbone'].iloc[1] == pytest.approx(0.0, abs=1e-5)
    assert rmsd_df[ligand_selection].iloc[1] == pytest.approx(0.0, abs=1e-5)
    assert rmsd_df['ligand_local'].iloc[1] == pytest.approx(0.0, abs=1e-5)


def test_ligand_translation_is_detected_after_backbone_alignment():
    """A 3 A ligand-only translation remains visible after protein-backbone fitting."""
    coords = _base_coords()
    moved = coords.copy()
    moved[8:] += np.array([3.0, 0.0, 0.0])

    rmsd_df, ligand_selection, _ = _rmsd_with_local([coords, moved])

    assert rmsd_df['backbone'].iloc[1] == pytest.approx(0.0, abs=1e-5)
    assert rmsd_df[ligand_selection].iloc[1] == pytest.approx(3.0, abs=1e-5)
    assert rmsd_df['ligand_local'].iloc[1] == pytest.approx(3.0, abs=1e-5)


def test_local_pocket_alignment_removes_pocket_domain_motion():
    """When pocket and ligand move together, ligand_local remains near zero."""
    coords = _base_coords()
    moved = coords.copy()
    moved[:4] += np.array([0.0, 3.0, 0.0])
    moved[8:] += np.array([0.0, 3.0, 0.0])

    rmsd_df, ligand_selection, _ = _rmsd_with_local([coords, moved])

    assert rmsd_df[ligand_selection].iloc[1] > 0.1
    assert rmsd_df['ligand_local'].iloc[1] == pytest.approx(0.0, abs=1e-5)


def test_reference_pocket_membership_is_fixed_to_frame_zero():
    """Pocket indices are selected from frame 0 even if the ligand moves later."""
    from streamd.analysis import md_system_analysis as analysis

    coords = _base_coords()
    moved = coords.copy()
    moved[8:] += np.array([10.0, 0.0, 0.0])
    universe = _make_test_universe([coords, moved])

    universe.trajectory[1]
    pocket_selection, _ = analysis._reference_pocket_selections(universe, 'UNL', 5.0)

    assert pocket_selection == 'index 0 1 2 3'


def test_ligand_local_is_omitted_for_insufficient_pocket(caplog):
    """Invalid local pockets warn and do not prevent standard RMSD calculation."""
    from streamd.analysis import md_system_analysis as analysis

    coords = _base_coords()
    coords[:4] = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [3.0, 0.0, 0.0],
    ])
    universe = _make_test_universe([coords, coords])

    with caplog.at_level('WARNING'):
        pocket_selection, ligand_selection = analysis._reference_pocket_selections(universe, 'UNL', 5.0)

    rmsd_df = analysis.rmsd_for_atomgroups(universe, selection1='backbone', selection2=[ligand_selection])

    assert pocket_selection is None
    assert 'ligand_local' not in rmsd_df.columns
    assert 'fewer than three non-collinear pocket-backbone atoms' in caplog.text


def _make_two_ligand_universe():
    """Create a universe with two copies of the ligand resname (two binding sites)."""
    mda = pytest.importorskip("MDAnalysis")
    if not hasattr(mda, "Universe"):
        pytest.skip("MDAnalysis is stubbed")

    atom_resindex = [0, 0, 0, 1, 1, 1]
    universe = mda.Universe.empty(
        6,
        n_residues=2,
        atom_resindex=atom_resindex,
        trajectory=True,
    )
    universe.add_TopologyAttr('name', ['C1', 'C2', 'C3', 'C1', 'C2', 'C3'])
    universe.add_TopologyAttr('type', ['C', 'C', 'C', 'C', 'C', 'C'])
    universe.add_TopologyAttr('resname', ['UNL', 'UNL'])
    universe.add_TopologyAttr('resid', [1, 2])
    universe.add_TopologyAttr('segid', ['A'])
    coords = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [10.0, 0.0, 0.0],
        [11.0, 0.0, 0.0],
        [10.0, 1.0, 0.0],
    ])
    universe.load_new(np.asarray([coords], dtype=np.float32), order='fac', dt=20)
    return universe


def test_reference_pocket_skipped_for_multiple_ligand_residues(caplog):
    """Multiple ligand copies (e.g. per chain) skip the single-site pocket metrics."""
    from streamd.analysis import md_system_analysis as analysis

    universe = _make_two_ligand_universe()

    with caplog.at_level('WARNING'):
        pocket_selection, ligand_selection = analysis._reference_pocket_selections(universe, 'UNL', 5.0)

    assert pocket_selection is None
    assert ligand_selection == 'resname UNL and not name H*'
    assert 'single binding site' in caplog.text


def test_md_analysis_shell_uses_center_and_trajectory_fit_groups_without_backbone_xtc():
    """Centering uses center_group; fitted outputs use trajectory_fit_group."""
    script = os.path.join(pytest.script_directory, 'script_sh', 'md_analysis.sh')
    with open(script) as inp:
        lines = inp.readlines()

    command_lines = [line for line in lines if line.lstrip().startswith('printf')]
    center_lines = [line for line in command_lines if '-center' in line]
    fit_lines = [line for line in command_lines if '-fit rot+trans' in line]
    rmsf_lines = [line for line in command_lines if 'gmx rmsf' in line]

    assert center_lines
    assert all('"$center_group"' in line for line in center_lines)
    assert len(fit_lines) == 3
    assert all('"$index_group"' not in line for line in fit_lines)
    assert all('"$fit_group"' not in line for line in fit_lines)
    assert all('md_fit_backbone' not in line for line in command_lines)
    assert any('-o md_fit.xtc' in line and '"$trajectory_fit_group"' in line and '"System"' in line for line in fit_lines)
    assert any('-o md_fit_nowater.xtc' in line and '"$trajectory_fit_group"' in line and '"non-Water"' in line for line in fit_lines)
    assert any('-o md_out_nowater.gro' in line and '"$trajectory_fit_group"' in line and '"non-Water"' in line for line in fit_lines)
    assert rmsf_lines
    assert all('"Protein"' in line for line in rmsf_lines)
    assert all('-f md_fit.xtc' in line for line in rmsf_lines)
    assert all('md_fit_backbone.xtc' not in line for line in rmsf_lines)
    assert all('-nofit' not in line for line in rmsf_lines)


@pytest.mark.parametrize(
    'molid_resid_pairs, expected_center_group, expected_trajectory_fit_group',
    [
        pytest.param([], 1, 1, id='protein-only'),
        pytest.param([('ligand', 'UNL')], 7, 7, id='protein-ligand'),
        pytest.param([('cofactor', 'COF')], 8, 8, id='protein-cofactor'),
        pytest.param([('ligand', 'UNL'), ('cofactor', 'COF')], 9, 7, id='protein-ligand-cofactor'),
        pytest.param([('ligand', 'UNL'), ('cofactor', 'GTP'), ('metal', 'MG')], 10, 7, id='protein-ligand-multiple-cofactors'),
    ],
)
def test_resolve_analysis_groups_for_supported_systems(
    tmp_path,
    molid_resid_pairs,
    expected_center_group,
    expected_trajectory_fit_group,
):
    """Centering follows the complex; trajectory fitting follows protein plus primary ligand."""
    from streamd.analysis import md_system_analysis as analysis

    index_list = [
        'System',
        'Protein',
        'non-Water',
        'UNL',
        'COF',
        'GTP',
        'MG',
        'Protein_UNL',
        'Protein_COF',
        'Protein_UNL_COF',
        'Protein_UNL_GTP_MG',
    ]

    center_group, trajectory_fit_group, _ = analysis._resolve_analysis_groups(
        index_list=index_list,
        molid_resid_pairs=molid_resid_pairs,
        ligand_resid='UNL',
        wdir=str(tmp_path),
        bash_log='bash.log',
        env=None,
    )

    assert center_group == expected_center_group
    assert trajectory_fit_group == expected_trajectory_fit_group


def test_resolve_analysis_groups_creates_missing_combined_groups_once(tmp_path, monkeypatch):
    """Missing center groups are created from unique residues and reloaded from index.ndx."""
    from streamd.analysis import md_system_analysis as analysis

    index_states = [
        ['System', 'Protein', 'non-Water', 'UNL', 'GTP', 'Protein_UNL'],
        ['System', 'Protein', 'non-Water', 'UNL', 'GTP', 'Protein_UNL', 'Protein_UNL_GTP'],
    ]
    queries = []

    def fake_make_group_ndx(query, wdir, bash_log, env=None):
        queries.append(query)
        return True

    def fake_get_index(index_file, env=None):
        return index_states[1]

    monkeypatch.setattr(analysis, 'make_group_ndx', fake_make_group_ndx)
    monkeypatch.setattr(analysis, 'get_index', fake_get_index)

    center_group, trajectory_fit_group, index_list = analysis._resolve_analysis_groups(
        index_list=index_states[0],
        molid_resid_pairs=[
            ('ligand', 'UNL'),
            ('duplicate-ligand-entry', 'UNL'),
            ('cofactor', 'GTP'),
        ],
        ligand_resid='UNL',
        wdir=str(tmp_path),
        bash_log='bash.log',
        env=None,
    )

    assert queries == ['1|3|4']
    assert center_group == 6
    assert trajectory_fit_group == 5
    assert index_list == index_states[1]


def test_resolve_analysis_groups_aborts_on_missing_residue_group(tmp_path, monkeypatch, caplog):
    """A residue with no index group signals a broken topology and skips the system."""
    from streamd.analysis import md_system_analysis as analysis

    index_list = ['System', 'Protein', 'non-Water', 'UNL', 'Protein_UNL']

    def fail_make_group_ndx(*args, **kwargs):
        raise AssertionError('No group should be created for a residue absent from the system')

    monkeypatch.setattr(analysis, 'make_group_ndx', fail_make_group_ndx)

    with caplog.at_level('ERROR'):
        result = analysis._resolve_analysis_groups(
            index_list=list(index_list),
            molid_resid_pairs=[('cofactor', 'ZZZ'), ('ligand', 'UNL')],
            ligand_resid='UNL',
            wdir=str(tmp_path),
            bash_log='bash.log',
            env=None,
        )

    assert result is None
    assert 'ZZZ' in caplog.text


def test_run_md_analysis_propagates_center_and_trajectory_fit_groups(tmp_path, monkeypatch):
    """The shell command receives separate center_group and trajectory_fit_group variables."""
    from streamd.analysis import md_system_analysis as analysis

    (tmp_path / 'index.ndx').write_text(
        '\n'.join(
            [
                '[ System ]',
                '[ Protein ]',
                '[ non-Water ]',
                '[ UNL ]',
                '[ Protein_UNL ]',
                '[ COF ]',
                '[ Protein_UNL_COF ]',
            ]
        )
    )
    (tmp_path / 'all_ligand_resid.txt').write_text('ligand\tUNL\ncofactor\tCOF\n')

    commands = []

    def fake_run_check_subprocess(cmd, *args, **kwargs):
        commands.append(cmd)
        return True

    monkeypatch.setattr(analysis, 'run_check_subprocess', fake_run_check_subprocess)
    monkeypatch.setattr(analysis, 'md_rmsd_analysis', lambda **kwargs: str(tmp_path / 'rmsd.csv'))
    monkeypatch.setattr(analysis, 'create_last_frame_file', lambda **kwargs: None)
    monkeypatch.setattr(analysis, 'convertxvg2png', lambda *args, **kwargs: None)

    result = analysis.run_md_analysis(
        (str(tmp_path), 'md_out'),
        mdtime_ns=0.05,
        project_dir=pytest.streamd_directory,
        bash_log='bash.log',
        active_site_dist=5.0,
        ligand_resid='UNL',
        save_traj_without_water=True,
        analysis_dirname='md_analysis',
        ligand_list_file_prev=None,
        env=None,
        system_name='protein_UNL',
    )

    assert result == (str(tmp_path / 'rmsd.csv'), str(tmp_path / 'md_analysis'), str(tmp_path))
    assert len(commands) == 1
    assert 'center_group=6' in commands[0]
    assert 'trajectory_fit_group=4' in commands[0]
    assert ' fit_group=' not in commands[0]
    assert 'index_group=' not in commands[0]


def test_run_md_analysis_removes_no_water_temporary_files(tmp_path, monkeypatch):
    """Temporary no-water files are removed when they are not retained."""
    from streamd.analysis import md_system_analysis as analysis

    (tmp_path / 'index.ndx').write_text('\n'.join(['[ System ]', '[ Protein ]', '[ non-Water ]']))

    def fake_run_check_subprocess(cmd, *args, **kwargs):
        for filename in ['md_out_nowater.tpr', 'md_out_nowater.gro', 'md_fit_nowater.xtc']:
            (tmp_path / filename).write_text('temporary')
        return True

    monkeypatch.setattr(analysis, 'run_check_subprocess', fake_run_check_subprocess)
    monkeypatch.setattr(analysis, 'md_rmsd_analysis', lambda **kwargs: str(tmp_path / 'rmsd.csv'))
    monkeypatch.setattr(analysis, 'create_last_frame_file', lambda **kwargs: None)
    monkeypatch.setattr(analysis, 'convertxvg2png', lambda *args, **kwargs: None)

    result = analysis.run_md_analysis(
        (str(tmp_path), 'md_out'),
        mdtime_ns=0.05,
        project_dir=pytest.streamd_directory,
        bash_log='bash.log',
        save_traj_without_water=False,
        analysis_dirname='md_analysis',
        ligand_list_file_prev=None,
        env=None,
        system_name='protein_only',
    )

    assert result == (str(tmp_path / 'rmsd.csv'), str(tmp_path / 'md_analysis'), str(tmp_path))
    assert not (tmp_path / 'md_out_nowater.tpr').exists()
    assert not (tmp_path / 'md_out_nowater.gro').exists()
    assert not (tmp_path / 'md_fit_nowater.xtc').exists()


def test_gbsa_default_trajectory_remains_md_fit():
    """GBSA defaults must keep using the complex-fitted md_fit.xtc pathway."""
    run_gbsa = os.path.join(pytest.streamd_directory, 'run_gbsa.py')
    with open(run_gbsa) as inp:
        source = inp.read()

    assert "xtc = 'md_fit.xtc'" in source
    assert '-cg {protein_index} {ligand_index} -ct {xtc}' in source
    assert "xtc = 'md_fit_backbone.xtc'" not in source
    assert "xtc = 'md_fit_backbone_nowater.xtc'" not in source


def test_local_ligand_rmsd_uses_groupselection_column(monkeypatch):
    """ligand_local must read column 3: frame, time, fit-group RMSD, ligand RMSD."""
    from streamd.analysis import md_system_analysis as analysis

    class FakeTrajectory:
        def __getitem__(self, index):
            return index

    class FakeRMSD:
        def __init__(self, *args, **kwargs):
            self.results = SimpleNamespace(rmsd=None)

        def run(self):
            self.results.rmsd = np.array(
                [
                    [0.0, 0.0, 0.1, 1.1],
                    [1.0, 20.0, 0.2, 1.2],
                ]
            )

    monkeypatch.setattr(analysis, 'rms', SimpleNamespace(RMSD=FakeRMSD))
    rmsd_df = pd.DataFrame({'time(ns)': [0.0, 0.02]})
    universe = SimpleNamespace(trajectory=FakeTrajectory())

    result = analysis._add_local_ligand_rmsd(
        rmsd_df=rmsd_df,
        universe=universe,
        pocket_selection='index 1 2 3',
        ligand_selection='resname UNL and not name H*',
    )

    assert result['ligand_local'].tolist() == [1.1, 1.2]


def test_run_rmsd_analysis_aggregates_ligand_local(tmp_path):
    """Requested ligand and ligand_local metrics are aggregated independently."""
    from streamd.analysis.run_analysis import run_rmsd_analysis

    rmsd_file = _write_rmsd_input(tmp_path / 'rmsd.csv')

    run_rmsd_analysis(
        [str(rmsd_file)],
        wdir=str(tmp_path),
        unique_id='local',
        time_ranges=[(0.0, 0.02)],
        rmsd_type_list=['ligand', 'ligand_local'],
        paint_by_fname=None,
        title=None,
    )

    result = pd.read_csv(tmp_path / 'rmsd_mean_std_time-ranges_local.csv', sep='\t')
    assert set(result['rmsd_system']) == {'ligand', 'ligand_local'}


@pytest.mark.parametrize('rmsd_files', [[], None])
def test_run_rmsd_analysis_requires_input_files(tmp_path, rmsd_files):
    """Missing input files fail with a clear public API error."""
    from streamd.analysis.run_analysis import run_rmsd_analysis

    with pytest.raises(ValueError, match='At least one RMSD input file is required'):
        run_rmsd_analysis(
            rmsd_files,
            wdir=str(tmp_path),
            unique_id='missing-input',
            time_ranges=[(0.0, 0.02)],
            rmsd_type_list=['ligand'],
            paint_by_fname=None,
            title=None,
        )


def test_run_rmsd_analysis_cli_requires_input(monkeypatch, capsys):
    """The CLI should reject missing --input during argument parsing."""
    from streamd.analysis import run_analysis

    monkeypatch.setattr(sys, 'argv', ['run_rmsd_analysis'])

    with pytest.raises(SystemExit) as exc:
        run_analysis.main()

    assert exc.value.code == 2
    captured = capsys.readouterr()
    assert 'required' in captured.err
    assert '--input' in captured.err


@pytest.mark.parametrize('missing_base', [['time(ns)'], ['system'], ['ligand_name']])
def test_run_rmsd_analysis_requires_base_columns(tmp_path, missing_base):
    """Missing structural base columns fail with a clear per-file error."""
    from streamd.analysis.run_analysis import run_rmsd_analysis

    rmsd_file = _write_rmsd_input(tmp_path / 'rmsd.csv', drop_columns=missing_base)

    with pytest.raises(ValueError) as exc:
        run_rmsd_analysis(
            [str(rmsd_file)],
            wdir=str(tmp_path),
            unique_id='missing-base',
            time_ranges=[(0.0, 0.02)],
            rmsd_type_list=['ligand'],
            paint_by_fname=None,
            title=None,
        )

    message = str(exc.value)
    assert str(rmsd_file) in message
    for column in missing_base:
        assert column in message


def test_run_rmsd_analysis_skips_metric_missing_from_some_files(tmp_path, caplog):
    """A metric absent from only some files is kept and warned about, not fatal."""
    from streamd.analysis.run_analysis import run_rmsd_analysis

    file_with_local = _write_rmsd_input(tmp_path / 'with_local.csv')
    file_without_local = _write_rmsd_input(
        tmp_path / 'without_local.csv',
        drop_columns=['ligand_local'],
        column_overrides={
            'ligand_name': ['ligand_b', 'ligand_b'],
            'system': ['protein_b_replica1', 'protein_b_replica1'],
        },
    )

    with caplog.at_level('WARNING'):
        run_rmsd_analysis(
            [str(file_with_local), str(file_without_local)],
            wdir=str(tmp_path),
            unique_id='mixed-missing-local',
            time_ranges=[(0.0, 0.02)],
            rmsd_type_list=['ligand', 'ligand_local'],
            paint_by_fname=None,
            title=None,
        )

    result = pd.read_csv(tmp_path / 'rmsd_mean_std_time-ranges_mixed-missing-local.csv', sep='\t')
    assert set(result['rmsd_system']) == {'ligand', 'ligand_local'}
    assert 'ligand_local' in caplog.text
    assert str(file_without_local) in caplog.text


def test_run_rmsd_analysis_drops_metric_absent_from_all_files(tmp_path, caplog):
    """A metric absent from every file is dropped with a warning if others survive."""
    from streamd.analysis.run_analysis import run_rmsd_analysis

    rmsd_file = _write_rmsd_input(tmp_path / 'rmsd.csv', drop_columns=['ligand_local'])

    with caplog.at_level('WARNING'):
        run_rmsd_analysis(
            [str(rmsd_file)],
            wdir=str(tmp_path),
            unique_id='drop-all',
            time_ranges=[(0.0, 0.02)],
            rmsd_type_list=['ligand', 'ligand_local'],
            paint_by_fname=None,
            title=None,
        )

    result = pd.read_csv(tmp_path / 'rmsd_mean_std_time-ranges_drop-all.csv', sep='\t')
    assert set(result['rmsd_system']) == {'ligand'}
    assert 'ligand_local' in caplog.text


def test_run_rmsd_analysis_errors_when_no_requested_metric_present(tmp_path):
    """If none of the requested metrics exist anywhere, fail with a clear error."""
    from streamd.analysis.run_analysis import run_rmsd_analysis

    rmsd_file = _write_rmsd_input(tmp_path / 'rmsd.csv', drop_columns=['ligand', 'ligand_local'])

    with pytest.raises(ValueError, match='None of the requested'):
        run_rmsd_analysis(
            [str(rmsd_file)],
            wdir=str(tmp_path),
            unique_id='no-metric',
            time_ranges=[(0.0, 0.02)],
            rmsd_type_list=['ligand', 'ligand_local'],
            paint_by_fname=None,
            title=None,
        )


@pytest.mark.parametrize(
    'start, end, expected',
    [
        (0.0, 10.0, [(0.0, 10.0), (5, 10.0), (9.0, 10.0)]),
        (10.0, 20.0, [(10.0, 20.0), (15, 20.0), (19.0, 20.0)]),
        (0.0, 9.5, [(0.0, 9.5), (5, 9.5), (8.5, 9.5)]),
        (3.0, 10.0, [(3.0, 10.0), (7, 10.0), (9.0, 10.0)]),
        (3.2, 3.8, [(3.2, 3.8)]),
    ],
)
def test_default_time_ranges_use_whole_nanosecond_midpoint(start, end, expected):
    """Default second-half windows begin at the next whole-ns midpoint boundary."""
    from streamd.analysis.run_analysis import _default_time_ranges

    assert _default_time_ranges(start, end) == expected


@pytest.mark.parametrize(
    'start, end, expected',
    [
        (
            9.5,
            10.9,
            [
                (9.5, 10.9),
                (10, 10.9),
                (9.9, 10.9),
            ],
        ),
        (
            99.8,
            101.0,
            [
                (99.8, 101.0),
                (100.0, 101.0),
            ],
        ),
    ],
)
def test_default_time_ranges_keep_midpoint_inside_interval(start, end, expected):
    """Rounded midpoint windows must not start beyond the trajectory end."""
    from streamd.analysis.run_analysis import _default_time_ranges

    assert _default_time_ranges(start, end) == expected


@pytest.mark.parametrize(
    'start, end',
    [
        (0.0, 9.5),
        (3.0, 10.0),
        (9.5, 10.9),
        (99.8, 101.0),
    ],
)
def test_default_time_ranges_are_valid(start, end):
    """Every generated default range must stay within the available trajectory."""
    from streamd.analysis.run_analysis import _default_time_ranges

    time_ranges = _default_time_ranges(start, end)

    assert len(time_ranges) == len(set(time_ranges))

    for range_start, range_end in time_ranges:
        assert start <= range_start <= range_end <= end


def test_default_time_ranges_remove_duplicate_ranges():
    """Equivalent midpoint and final-1-ns ranges are emitted once."""
    from streamd.analysis.run_analysis import _default_time_ranges

    assert _default_time_ranges(0.0, 2.0) == [
        (0.0, 2.0),
        (1.0, 2.0),
    ]


def test_default_time_ranges_normalize_float_boundaries():
    """Floating-point artifacts are rounded in emitted default ranges."""
    from streamd.analysis.run_analysis import _default_time_ranges

    assert _default_time_ranges(0.0, 1.1) == [
        (0.0, 1.1),
        (1.0, 1.1),
        (0.1, 1.1),
    ]


def test_run_rmsd_analysis_paint_by_requires_value_column_for_ligand_system(tmp_path):
    """Protein-ligand paint-by files need a property column beyond identifiers."""
    from streamd.analysis.run_analysis import run_rmsd_analysis

    rmsd_file = _write_rmsd_input(tmp_path / 'rmsd.csv')

    paint_by_file = tmp_path / 'paint_by.csv'
    pd.DataFrame({
        'protein_name': ['protein_a'],
        'ligand_name': ['ligand_a'],
    }).to_csv(paint_by_file, sep='\t', index=False)

    with pytest.raises(ValueError) as exc:
        run_rmsd_analysis(
            [str(rmsd_file)],
            wdir=str(tmp_path),
            unique_id='paint-by-missing-value',
            time_ranges=[(0.0, 0.02)],
            rmsd_type_list=['ligand'],
            paint_by_fname=str(paint_by_file),
            title=None,
        )

    message = str(exc.value)
    assert str(paint_by_file) in message
    assert 'protein_name' in message
    assert 'ligand_name' in message


def test_run_rmsd_analysis_paint_by_requires_value_column_for_protein_only(tmp_path):
    """Protein-only paint-by files need a property column beyond protein_name."""
    from streamd.analysis.run_analysis import run_rmsd_analysis

    rmsd_file = _write_rmsd_input(
        tmp_path / 'rmsd.csv',
        column_overrides={'ligand_name': [None, None]},
    )

    paint_by_file = tmp_path / 'paint_by.csv'
    pd.DataFrame({
        'protein_name': ['protein_a'],
    }).to_csv(paint_by_file, sep='\t', index=False)

    with pytest.raises(ValueError) as exc:
        run_rmsd_analysis(
            [str(rmsd_file)],
            wdir=str(tmp_path),
            unique_id='protein-only-paint-by-missing-value',
            time_ranges=[(0.0, 0.02)],
            rmsd_type_list=['backbone'],
            paint_by_fname=str(paint_by_file),
            title=None,
        )

    message = str(exc.value)
    assert str(paint_by_file) in message
    assert 'protein_name' in message


@pytest.mark.parametrize(
    'extra_columns',
    [
        {'pKi': [7.2]},
        {'pKi': [7.2], 'activity_class': ['active']},
    ],
)
def test_run_rmsd_analysis_uses_first_paint_by_value_column(tmp_path, monkeypatch, extra_columns):
    """Valid paint-by files use the first non-identifier column for colouring."""
    from streamd.analysis import run_analysis

    rmsd_file = _write_rmsd_input(tmp_path / 'rmsd.csv')

    paint_by_file = tmp_path / 'paint_by.csv'
    paint_by_data = {
        'protein_name': ['protein_a'],
        'ligand_name': ['ligand_a'],
    }
    paint_by_data.update(extra_columns)
    pd.DataFrame(paint_by_data).to_csv(paint_by_file, sep='\t', index=False)

    plot_calls = []

    def fake_plot_rmsd_mean_std(**kwargs):
        plot_calls.append(kwargs)

    monkeypatch.setattr(run_analysis, 'plot_rmsd_mean_std', fake_plot_rmsd_mean_std)

    run_analysis.run_rmsd_analysis(
        [str(rmsd_file)],
        wdir=str(tmp_path),
        unique_id='paint-by-valid',
        time_ranges=[(0.0, 0.02)],
        rmsd_type_list=['ligand'],
        paint_by_fname=str(paint_by_file),
        title=None,
    )

    assert plot_calls[0]['paint_by_col'] == 'pKi'


@analysis_test
@pytest.mark.filterwarnings('ignore::DeprecationWarning')
def test_rmsd_analysis(dir_with_streamd_output_for_analysis,
                       # list_expected_system_analysis_output,
                       # list_expected_analysis_output
                       ):
    """Ensure RMSD analysis generates expected files."""
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

    assert not os.path.isfile(os.path.join(dir_with_streamd_output_files, 'md_fit_backbone.xtc'))
    assert not os.path.isfile(os.path.join(dir_with_streamd_output_files, 'md_fit_backbone_nowater.xtc'))

    rmsd_data = pd.read_csv(rmsd_file, sep='\t')
    expected_rmsd_columns = {
        'time(ns)',
        'backbone',
        'ligand',
        'ActiveSite5.0A',
        'ligand_local',
        'ligand_name',
        'system',
        'protein_name',
        'replica',
        'directory',
    }
    assert expected_rmsd_columns.issubset(rmsd_data.columns)


@analysis_test
@pytest.mark.filterwarnings('ignore::DeprecationWarning')
# @pytest.mark.skip(reason="Ignore")
@pytest.mark.parametrize("rmsd_type_list", [
    pytest.param(['backbone', 'ligand', f'ActiveSite5.0A'], id="backbone, ligand, ActiveSite5.0A"),
    pytest.param(['backbone'], id="backbone"),
    pytest.param(['ActiveSite5.0A'], id="ActiveSite5.0A"),
])
def test_rmsd_analysis(rmsd_type_list, dir_and_rmsd_files):
    """Run summary RMSD analysis with optional time ranges."""
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
    """Validate painting of RMSD plots by experimental data."""
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
