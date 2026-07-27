"""Tests for the ``--ligand_forcefield`` (GAFF / GAFF2) ligand selection.

These tests do not require a real AmberTools installation. External
AmberTools binaries (``antechamber``, ``parmchk2``, ``tleap``) are replaced by
tiny fake executables that record the argument list they are called with, so we
can assert exactly which flags reach each tool. The higher-level Python
propagation test mocks ``run_check_subprocess``/``prepare_tleap`` and captures
the generated commands.
"""

import logging
import os
import subprocess
from unittest import mock

import pytest

from streamd.run_md import get_parser
from streamd.preparation import ligand_preparation as lp
from streamd.preparation.ligand_preparation import (
    LIGAND_LEAPRC,
    backup_stale_ligand_files,
    ensure_no_conflicting_run,
    ensure_run_forcefield_compatible,
    has_existing_ligand_files,
    prepare_tleap,
    resolve_ligand_forcefield,
)


def _tleap_template():
    return os.path.join(pytest.script_directory, 'tleap.in')


def _script(name):
    return os.path.join(pytest.script_directory, 'script_sh', name)


def _make_fake_bin(bindir, name, body):
    """Create an executable shell stub named *name* in *bindir*."""
    path = os.path.join(bindir, name)
    with open(path, 'w') as fh:
        fh.write('#!/bin/bash\n' + body + '\n')
    os.chmod(path, 0o755)
    return path


def _run_script(script, bindir, workdir, env_extra):
    env = os.environ.copy()
    env['PATH'] = f'{bindir}{os.pathsep}' + env['PATH']
    env.update(env_extra)
    return subprocess.run(['bash', script], cwd=workdir, env=env,
                          capture_output=True, text=True)


# --------------------------------------------------------------------------- #
# 1. CLI argument tests
# --------------------------------------------------------------------------- #
class TestCli:
    def test_default_is_gaff(self):
        assert get_parser().parse_args([]).ligand_forcefield == 'gaff'

    def test_gaff_accepted(self):
        args = get_parser().parse_args(['--ligand_forcefield', 'gaff'])
        assert args.ligand_forcefield == 'gaff'

    def test_gaff2_accepted(self):
        args = get_parser().parse_args(['--ligand_forcefield', 'gaff2'])
        assert args.ligand_forcefield == 'gaff2'

    def test_invalid_value_rejected(self):
        with pytest.raises(SystemExit):
            get_parser().parse_args(['--ligand_forcefield', 'openff'])


# --------------------------------------------------------------------------- #
# 2. LEaP template rendering
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize('ff,expected', [('gaff', 'leaprc.gaff'), ('gaff2', 'leaprc.gaff2')])
def test_tleap_sources_correct_leaprc(tmp_path, ff, expected):
    out = tmp_path / 'tleap.in'
    prepare_tleap(_tleap_template(), str(out), molid='MYLIG',
                  conda_env_path='/env', ligand_forcefield=ff)
    lines = out.read_text().splitlines()
    assert lines[0] == f'source /env/dat/leap/cmd/{expected}'
    # molid substitution still works and mapping is exact (gaff must not become gaff2)
    assert 'MYLIG.frcmod' in out.read_text()
    assert 'MYLIG.mol2' in out.read_text()
    assert LIGAND_LEAPRC[ff] == expected


def test_tleap_default_is_backward_compatible(tmp_path):
    """Omitting ligand_forcefield reproduces the previous GAFF behaviour."""
    out = tmp_path / 'tleap.in'
    prepare_tleap(_tleap_template(), str(out), molid='LIG', conda_env_path='/env')
    assert out.read_text().splitlines()[0] == 'source /env/dat/leap/cmd/leaprc.gaff'


def test_tleap_invalid_value_raises(tmp_path):
    with pytest.raises(ValueError):
        prepare_tleap(_tleap_template(), str(tmp_path / 'x.in'), molid='LIG',
                      conda_env_path='/env', ligand_forcefield='openff')


# --------------------------------------------------------------------------- #
# 3. Command generation reaching the real tools (fake AmberTools binaries)
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize('ff', ['gaff', 'gaff2'])
def test_antechamber_receives_atom_type_and_bcc(tmp_path, ff):
    """ligand_mol2prep.sh must pass ``-at <ff>`` and keep ``-c bcc``."""
    bindir = tmp_path / 'bin'
    bindir.mkdir()
    work = tmp_path / 'work'
    work.mkdir()
    args_file = tmp_path / 'antechamber.args'
    _make_fake_bin(str(bindir), 'antechamber',
                   f'printf "%s\\n" "$@" > "{args_file}"\nexit 0')
    (work / 'LIG.mol').write_text('dummy')

    r = _run_script(_script('ligand_mol2prep.sh'), str(bindir), str(work), {
        'input_dirname': str(work), 'lfile': str(work / 'LIG.mol'),
        'molid': 'LIG', 'resid': 'UNL', 'charge': '0', 'dr': 'yes',
        'ligand_forcefield': ff,
    })
    assert r.returncode == 0, r.stderr
    args = args_file.read_text().splitlines()
    assert '-at' in args and args[args.index('-at') + 1] == ff
    assert '-c' in args and args[args.index('-c') + 1] == 'bcc'


@pytest.mark.parametrize('ff', ['gaff', 'gaff2'])
def test_parmchk2_receives_parameter_set(tmp_path, ff):
    """ligand_prep.sh must pass ``-s <ff>`` to parmchk2."""
    bindir = tmp_path / 'bin'
    bindir.mkdir()
    work = tmp_path / 'work'
    work.mkdir()
    args_file = tmp_path / 'parmchk2.args'
    _make_fake_bin(str(bindir), 'parmchk2',
                   f'printf "%s\\n" "$@" > "{args_file}"\nexit 0')
    # stop the chain right after parmchk2 so the fake tleap/gmx pipeline is not needed
    _make_fake_bin(str(bindir), 'tleap', 'exit 1')
    (work / 'LIG.mol2').write_text('dummy')

    r = _run_script(_script('ligand_prep.sh'), str(bindir), str(work), {
        'input_dirname': str(work), 'script_path': pytest.script_directory,
        'molid': 'LIG', 'ligand_forcefield': ff,
    })
    # the fake tleap intentionally fails, so the script stops right after parmchk2
    assert r.returncode != 0
    assert args_file.exists()
    args = args_file.read_text().splitlines()
    assert '-s' in args and args[args.index('-s') + 1] == ff
    # -frc must not be added (protein FF is not a ligand parameter set)
    assert '-frc' not in args


@pytest.mark.parametrize('ff', ['gaff', 'gaff2'])
def test_gaussian_antechamber_receives_forcefield_and_resp(tmp_path, ff):
    """The Gaussian/RESP branch must pass ``-at <ff>`` and keep ``-c resp``."""
    bindir = tmp_path / 'bin'
    bindir.mkdir()
    work = tmp_path / 'work'
    work.mkdir()
    args_file = tmp_path / 'antechamber.args'

    # record args at the first (RESP) antechamber call and stop the script there
    _make_fake_bin(str(bindir), 'antechamber', f'printf "%s\\n" "$@" >> "{args_file}"\nexit 1')
    _make_fake_bin(str(bindir), 'g09', 'exit 0')
    _make_fake_bin(str(bindir), 'python', 'exit 0')
    # prepare_Gaussian_input.py (mocked above) would create these; they must exist so the
    # `g09 < ${molid}_*_gaussian.com` input redirections succeed and the script reaches
    # the RESP antechamber step.
    (work / 'LIG_opt_gaussian.com').write_text('')
    (work / 'LIG_charg_gaussian.com').write_text('')

    r = _run_script(_script('ligand_mol2prep_by_gaussian.sh'), str(bindir), str(work), {
        'input_dirname': str(work), 'lfile': str(work / 'LIG.mol'),
        'script_path': pytest.script_directory, 'molid': 'LIG', 'resid': 'UNL',
        'charge': '0', 'gaussian_version': 'g09', 'activate_gaussian': '',
        'ligand_forcefield': ff,
    })

    assert r.returncode != 0  # fake antechamber exits 1 to stop the script here
    assert args_file.exists()
    args = args_file.read_text().splitlines()
    assert '-at' in args and args[args.index('-at') + 1] == ff
    assert '-c' in args and args[args.index('-c') + 1] == 'resp'


def test_prep_ligand_rejects_invalid_forcefield(tmp_path):
    """Unsupported force fields are rejected before any external tool is invoked.

    Invalid values are blocked by argparse ``choices`` on the CLI (see TestCli) and,
    for values reaching the API directly (e.g. via a YAML config default, which
    argparse does not re-validate), by this guard in prep_ligand.
    """
    mol = mock.MagicMock()
    with pytest.raises(ValueError):
        lp.prep_ligand(
            (mol, 'LIG', 'UNL'),
            script_path=pytest.script_directory,
            project_dir=pytest.streamd_directory,
            wdir_ligand=str(tmp_path),
            conda_env_path='/env',
            bash_log='bash.log',
            ligand_forcefield='openff',
        )


def test_shell_defaults_to_gaff_when_unset(tmp_path):
    """A missing ligand_forcefield env var falls back to gaff (backward compat)."""
    bindir = tmp_path / 'bin'
    bindir.mkdir()
    work = tmp_path / 'work'
    work.mkdir()
    args_file = tmp_path / 'antechamber.args'
    _make_fake_bin(str(bindir), 'antechamber',
                   f'printf "%s\\n" "$@" > "{args_file}"\nexit 0')
    (work / 'LIG.mol').write_text('dummy')

    r = _run_script(_script('ligand_mol2prep.sh'), str(bindir), str(work), {
        'input_dirname': str(work), 'lfile': str(work / 'LIG.mol'),
        'molid': 'LIG', 'resid': 'UNL', 'charge': '0', 'dr': 'yes',
        # ligand_forcefield intentionally omitted
    })
    assert r.returncode == 0, r.stderr
    args = args_file.read_text().splitlines()
    assert args[args.index('-at') + 1] == 'gaff'


# --------------------------------------------------------------------------- #
# 4. Python-level propagation: prep_ligand -> all three tools
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize('ff', ['gaff', 'gaff2'])
def test_prep_ligand_propagates_to_all_stages(tmp_path, monkeypatch, ff):
    """A partial implementation (missing any of the three) fails this test."""
    calls = []
    tleap_kwargs = {}

    def fake_run(cmd, *args, **kwargs):
        calls.append(cmd[0] if isinstance(cmd, (list, tuple)) else cmd)
        return True

    def fake_prepare_tleap(*args, **kwargs):
        tleap_kwargs.update(kwargs)

    monkeypatch.setattr(lp, 'run_check_subprocess', fake_run)
    monkeypatch.setattr(lp, 'prepare_tleap', fake_prepare_tleap)
    monkeypatch.setattr(lp, 'reorder_hydrogens', lambda m: m)
    monkeypatch.setattr(lp.Chem, 'AddHs', lambda m, **k: m)
    monkeypatch.setattr(lp.Chem, 'MolToMolFile', lambda *a, **k: None)
    monkeypatch.setattr(lp.Chem, 'MolFromSmarts', lambda s: s)
    monkeypatch.setattr(lp.rdmolops, 'GetFormalCharge', lambda m: 0)

    mol = mock.MagicMock()
    mol.HasSubstructMatch.return_value = False  # not boron-containing

    res = lp.prep_ligand(
        (mol, 'LIG', 'UNL'),
        script_path=pytest.script_directory,
        project_dir=pytest.streamd_directory,
        wdir_ligand=str(tmp_path),
        conda_env_path='/env',
        bash_log='bash.log',
        ligand_forcefield=ff,
    )

    assert res == os.path.join(str(tmp_path), 'LIG')

    mol2prep_cmds = [c for c in calls if 'ligand_mol2prep.sh' in c]
    ligprep_cmds = [c for c in calls if 'ligand_prep.sh' in c]

    # antechamber (atom typing) branch
    assert mol2prep_cmds, 'ligand_mol2prep.sh was not invoked'
    assert f'ligand_forcefield="{ff}"' in mol2prep_cmds[0]
    # parmchk2 branch
    assert ligprep_cmds, 'ligand_prep.sh was not invoked'
    assert f'ligand_forcefield="{ff}"' in ligprep_cmds[0]
    # LEaP branch
    assert tleap_kwargs.get('ligand_forcefield') == ff


def test_supplied_mol2_skips_antechamber_and_warns(tmp_path, monkeypatch, caplog):
    """A pre-existing MOL2 is not re-typed: Antechamber is skipped and a warning is logged.

    parmchk2 (-s) and LEaP still follow the selected force field, but no atom
    typing occurs, so StreaMD must not silently claim GAFF2 conversion.
    """
    calls = []
    tleap_kwargs = {}

    def fake_run(cmd, *args, **kwargs):
        calls.append(cmd[0] if isinstance(cmd, (list, tuple)) else cmd)
        return True

    def fake_prepare_tleap(*args, **kwargs):
        tleap_kwargs.update(kwargs)

    class _Res:
        name = 'LIG'

    class _Struct:
        residues = [_Res()]

        def save(self, *a, **k):
            pass

    class _Loaded:
        def to_structure(self):
            return _Struct()

    monkeypatch.setattr(lp, 'run_check_subprocess', fake_run)
    monkeypatch.setattr(lp, 'prepare_tleap', fake_prepare_tleap)
    monkeypatch.setattr(lp.pmd, 'load_file', lambda *a, **k: _Loaded())

    mol2_in = tmp_path / 'input.mol2'
    mol2_in.write_text('dummy')
    mol = mock.MagicMock()

    with caplog.at_level(logging.WARNING):
        res = lp.prep_ligand(
            (mol, 'LIG', 'UNL'),
            script_path=pytest.script_directory,
            project_dir=pytest.streamd_directory,
            wdir_ligand=str(tmp_path),
            conda_env_path='/env',
            bash_log='bash.log',
            mol2_file=str(mol2_in),
            ligand_forcefield='gaff2',
        )

    assert res == os.path.join(str(tmp_path), 'LIG')
    # Antechamber (ligand_mol2prep.sh) must NOT be invoked for a supplied MOL2 ...
    assert not any('ligand_mol2prep.sh' in c for c in calls)
    # ... but parmchk2/LEaP still follow the selected force field
    assert any('ligand_prep.sh' in c and 'ligand_forcefield="gaff2"' in c for c in calls)
    assert tleap_kwargs.get('ligand_forcefield') == 'gaff2'
    # a warning must be emitted so the behaviour is not silent
    warnings = [r.getMessage() for r in caplog.records if r.levelno == logging.WARNING]
    assert any('as-is' in m and 'gaff2' in m for m in warnings), warnings


# --------------------------------------------------------------------------- #
# 5. Force-field metadata: no cross-parameter-set reuse of cached topologies
# --------------------------------------------------------------------------- #
def _stub_prep_calls(monkeypatch):
    """Neutralize external tools / RDKit so prep_ligand's control flow can run."""
    calls = []

    def fake_run(cmd, *a, **k):
        calls.append(cmd[0] if isinstance(cmd, (list, tuple)) else cmd)
        return True

    monkeypatch.setattr(lp, 'run_check_subprocess', fake_run)
    monkeypatch.setattr(lp, 'prepare_tleap', lambda *a, **k: None)
    monkeypatch.setattr(lp, 'reorder_hydrogens', lambda m: m)
    monkeypatch.setattr(lp.Chem, 'AddHs', lambda m, **k: m)
    monkeypatch.setattr(lp.Chem, 'MolToMolFile', lambda *a, **k: None)
    monkeypatch.setattr(lp.Chem, 'MolFromSmarts', lambda s: s)
    monkeypatch.setattr(lp.rdmolops, 'GetFormalCharge', lambda m: 0)
    return calls


def _make_prepared_ligand_dir(base, molid='LIG', forcefield=None):
    """Create a directory that looks like an already-prepared ligand."""
    ligdir = os.path.join(str(base), molid)
    os.makedirs(ligdir, exist_ok=True)
    open(os.path.join(ligdir, f'{molid}.itp'), 'w').close()
    open(os.path.join(ligdir, f'posre_{molid}.itp'), 'w').close()
    if forcefield is not None:
        lp.write_prepared_ligand_forcefield(ligdir, forcefield)
    return ligdir


def _run_prep(tmp_path, ff):
    mol = mock.MagicMock()
    mol.HasSubstructMatch.return_value = False
    return lp.prep_ligand(
        (mol, 'LIG', 'UNL'),
        script_path=pytest.script_directory,
        project_dir=pytest.streamd_directory,
        wdir_ligand=str(tmp_path),
        conda_env_path='/env',
        bash_log='bash.log',
        ligand_forcefield=ff,
    )


def test_metadata_read_write_roundtrip(tmp_path):
    assert lp.read_prepared_ligand_forcefield(str(tmp_path)) is None
    lp.write_prepared_ligand_forcefield(str(tmp_path), 'gaff2')
    assert lp.read_prepared_ligand_forcefield(str(tmp_path)) == 'gaff2'


def test_metadata_recorded_after_successful_prep(tmp_path, monkeypatch):
    _stub_prep_calls(monkeypatch)
    res = _run_prep(tmp_path, 'gaff2')
    assert lp.read_prepared_ligand_forcefield(res) == 'gaff2'


def test_skip_only_when_forcefield_matches(tmp_path, monkeypatch):
    _make_prepared_ligand_dir(tmp_path, forcefield='gaff2')
    calls = _stub_prep_calls(monkeypatch)
    res = _run_prep(tmp_path, 'gaff2')
    assert res == os.path.join(str(tmp_path), 'LIG')
    assert calls == []  # nothing re-run: cached GAFF2 topology reused


def test_mismatch_forces_reparameterization(tmp_path, monkeypatch, caplog):
    """gaff then gaff2 in the same directory must NOT reuse the GAFF topology."""
    ligdir = _make_prepared_ligand_dir(tmp_path, forcefield='gaff')
    # add a mol2 that would otherwise be silently reused (Antechamber skipped)
    open(os.path.join(ligdir, 'LIG.mol2'), 'w').close()
    calls = _stub_prep_calls(monkeypatch)
    with caplog.at_level(logging.WARNING):
        res = _run_prep(tmp_path, 'gaff2')
    assert res == ligdir
    # Antechamber atom typing + parmchk2/LEaP were re-run ...
    assert any('ligand_mol2prep.sh' in c for c in calls)
    assert any('ligand_prep.sh' in c for c in calls)
    # ... the record was replaced with the new force field ...
    assert lp.read_prepared_ligand_forcefield(ligdir) == 'gaff2'
    # ... the previous parameterization was PRESERVED (moved, not deleted) ...
    backup_dir = os.path.join(ligdir, 'backup_gaff')
    assert os.path.isdir(backup_dir)
    assert os.path.isfile(os.path.join(backup_dir, 'LIG.mol2'))
    assert os.path.isfile(os.path.join(backup_dir, 'LIG.itp'))
    assert lp.read_prepared_ligand_forcefield(backup_dir) == 'gaff'
    # ... and the switch was not silent
    assert any('gaff' in m and 'gaff2' in m for m in
               [r.getMessage() for r in caplog.records if r.levelno == logging.WARNING])


def test_legacy_topology_without_metadata_skips_for_gaff(tmp_path, monkeypatch):
    """A pre-option topology (no record) is GAFF, so a gaff request still skips."""
    _make_prepared_ligand_dir(tmp_path, forcefield=None)
    calls = _stub_prep_calls(monkeypatch)
    res = _run_prep(tmp_path, 'gaff')
    assert res == os.path.join(str(tmp_path), 'LIG')
    assert calls == []


def test_legacy_topology_without_metadata_reprepares_for_gaff2(tmp_path, monkeypatch):
    ligdir = _make_prepared_ligand_dir(tmp_path, forcefield=None)
    calls = _stub_prep_calls(monkeypatch)
    res = _run_prep(tmp_path, 'gaff2')
    assert res == ligdir
    assert any('ligand_mol2prep.sh' in c for c in calls)
    assert lp.read_prepared_ligand_forcefield(ligdir) == 'gaff2'


def _make_interrupted_ligand_dir(base, molid='LIG'):
    """Simulate an interrupted GAFF run: intermediates but no .itp/posre/metadata."""
    ligdir = os.path.join(str(base), molid)
    os.makedirs(ligdir, exist_ok=True)
    open(os.path.join(ligdir, f'{molid}.mol2'), 'w').close()
    open(os.path.join(ligdir, f'{molid}.frcmod'), 'w').close()
    return ligdir


def test_interrupted_prep_reparameterizes_on_mismatch(tmp_path, monkeypatch):
    """An interrupted GAFF run (mol2/frcmod, no itp/posre) must NOT be reused for gaff2.

    The mismatch check must fire even though the completed-output condition
    (itp AND posre) is not satisfied - otherwise the stale GAFF .mol2 would be
    reused with parmchk2 -s gaff2 / leaprc.gaff2.
    """
    ligdir = _make_interrupted_ligand_dir(tmp_path)
    calls = _stub_prep_calls(monkeypatch)
    res = _run_prep(tmp_path, 'gaff2')
    assert res == ligdir
    # the stale GAFF .mol2 was moved aside ...
    backup_dir = os.path.join(ligdir, 'backup_gaff')
    assert os.path.isdir(backup_dir)
    assert os.path.isfile(os.path.join(backup_dir, 'LIG.mol2'))
    # ... so Antechamber re-runs with gaff2 atom typing instead of being skipped
    assert any('ligand_mol2prep.sh' in c for c in calls)
    assert lp.read_prepared_ligand_forcefield(ligdir) == 'gaff2'


def test_interrupted_prep_same_forcefield_resumes(tmp_path, monkeypatch):
    """Rerunning an interrupted GAFF run with gaff reuses the existing .mol2 (resume)."""
    ligdir = _make_interrupted_ligand_dir(tmp_path)
    calls = _stub_prep_calls(monkeypatch)
    res = _run_prep(tmp_path, 'gaff')
    assert res == ligdir
    # same force field: no backup, and the consistent .mol2 is reused (Antechamber skipped) ...
    assert not os.path.isdir(os.path.join(ligdir, 'backup_gaff'))
    assert not any('ligand_mol2prep.sh' in c for c in calls)
    # ... while parmchk2/LEaP still run and the force field is recorded on success
    assert any('ligand_prep.sh' in c for c in calls)
    assert lp.read_prepared_ligand_forcefield(ligdir) == 'gaff'


def test_interrupted_gaff2_with_recorded_metadata_resumes(tmp_path, monkeypatch):
    """An interrupted GAFF2 run (FF recorded before the heavy steps) resumes, not regenerates."""
    ligdir = _make_interrupted_ligand_dir(tmp_path)
    lp.write_prepared_ligand_forcefield(ligdir, 'gaff2')  # recorded at prep start, before the crash

    calls = _stub_prep_calls(monkeypatch)
    res = _run_prep(tmp_path, 'gaff2')
    assert res == ligdir
    # matching force field: no backup, and the gaff2 .mol2 is reused (Antechamber skipped)
    assert not os.path.isdir(os.path.join(ligdir, 'backup_gaff2'))
    assert not any('ligand_mol2prep.sh' in c for c in calls)
    assert any('ligand_prep.sh' in c for c in calls)
    assert lp.read_prepared_ligand_forcefield(ligdir) == 'gaff2'


def test_interrupted_gaff2_metadata_blocks_reuse_as_gaff(tmp_path, monkeypatch):
    """A recorded GAFF2 interruption must not be silently resumed as GAFF (correctness)."""
    ligdir = _make_interrupted_ligand_dir(tmp_path)
    lp.write_prepared_ligand_forcefield(ligdir, 'gaff2')

    calls = _stub_prep_calls(monkeypatch)
    res = _run_prep(tmp_path, 'gaff')
    assert res == ligdir
    # switching to gaff: the gaff2 intermediates are backed up, not reused ...
    assert os.path.isdir(os.path.join(ligdir, 'backup_gaff2'))
    assert os.path.isfile(os.path.join(ligdir, 'backup_gaff2', 'LIG.mol2'))
    # ... and Antechamber re-runs to assign gaff atom types
    assert any('ligand_mol2prep.sh' in c for c in calls)
    assert lp.read_prepared_ligand_forcefield(ligdir) == 'gaff'


def test_backup_failure_aborts_without_reusing_stale_files(tmp_path, monkeypatch):
    """If a stale file cannot be moved, preparation must abort, not reuse it."""
    ligdir = _make_interrupted_ligand_dir(tmp_path)  # GAFF intermediates, no metadata
    calls = _stub_prep_calls(monkeypatch)

    def _boom(*a, **k):
        raise OSError('device or resource busy')

    monkeypatch.setattr(lp.os, 'rename', _boom)

    with pytest.raises(RuntimeError):
        _run_prep(tmp_path, 'gaff2')

    # aborted before any (re)parameterization ran: the stale GAFF .mol2 was NOT reused
    assert not any('ligand_mol2prep.sh' in c or 'ligand_prep.sh' in c for c in calls)
    # the stale file is still present (move failed) and was not silently accepted
    assert os.path.isfile(os.path.join(ligdir, 'LIG.mol2'))


def test_run_forcefield_guard_allows_matching_run_dir(tmp_path, caplog):
    """A run directory built with the same force field is accepted and the reuse is logged."""
    lp.write_prepared_ligand_forcefield(str(tmp_path), 'gaff2')
    (tmp_path / 'md_fit.xtc').write_text('x')
    with caplog.at_level(logging.INFO):
        ensure_run_forcefield_compatible(str(tmp_path), 'gaff2')  # must not raise
    assert any('already exists' in r.getMessage() and 'gaff2' in r.getMessage()
               for r in caplog.records)


def test_run_forcefield_guard_blocks_mismatch_with_record(tmp_path):
    """An existing run (recorded gaff) is not reused for a gaff2 request."""
    lp.write_prepared_ligand_forcefield(str(tmp_path), 'gaff')
    (tmp_path / 'md_fit.xtc').write_text('x')
    with pytest.raises(ValueError):
        ensure_run_forcefield_compatible(str(tmp_path), 'gaff2')


def test_run_forcefield_guard_blocks_legacy_dir_conservatively(tmp_path):
    """A legacy run directory (content but no record) is assumed GAFF and blocks a gaff2 request."""
    (tmp_path / 'topol.top').write_text('x')  # some content, no ligand_forcefield.txt
    with pytest.raises(ValueError):
        ensure_run_forcefield_compatible(str(tmp_path), 'gaff2')


def test_run_forcefield_guard_blocks_mismatch_on_existing_run_dir(tmp_path):
    """A run directory with a mismatched record is blocked even before any MD output exists."""
    lp.write_prepared_ligand_forcefield(str(tmp_path), 'gaff2')  # non-empty, recorded gaff2
    with pytest.raises(ValueError):
        ensure_run_forcefield_compatible(str(tmp_path), 'gaff')


def test_run_forcefield_guard_allows_empty_dir(tmp_path):
    """An empty run directory (nothing prepared) is allowed."""
    ensure_run_forcefield_compatible(str(tmp_path), 'gaff2')  # must not raise


def test_no_conflicting_run_stops_if_any_replica_conflicts(tmp_path):
    """A pre-flight over all replicas stops if any one used a different force field."""
    md_run = tmp_path / 'md_run'
    r1 = md_run / 'complex_replica1'
    r2 = md_run / 'complex_replica2'
    r1.mkdir(parents=True)
    r2.mkdir(parents=True)
    lp.write_prepared_ligand_forcefield(str(r1), 'gaff2')  # matches
    lp.write_prepared_ligand_forcefield(str(r2), 'gaff')   # conflicts
    with pytest.raises(ValueError):
        ensure_no_conflicting_run(str(md_run), 'gaff2')


def test_no_conflicting_run_passes_when_all_replicas_consistent(tmp_path):
    md_run = tmp_path / 'md_run'
    for name in ('complex_replica1', 'complex_replica2'):
        d = md_run / name
        d.mkdir(parents=True)
        lp.write_prepared_ligand_forcefield(str(d), 'gaff2')
    ensure_no_conflicting_run(str(md_run), 'gaff2')  # must not raise


def test_no_conflicting_run_noop_when_no_run_dir(tmp_path):
    """No md_run directory yet -> nothing to check."""
    ensure_no_conflicting_run(str(tmp_path / 'md_run'), 'gaff2')  # must not raise


def test_resolve_forcefield_adopts_recorded_when_not_explicit(tmp_path):
    """Omitting --ligand_forcefield adopts the FF a previous run used."""
    run = tmp_path / 'md_files' / 'md_run' / 'complex_replica1'
    run.mkdir(parents=True)
    lp.write_prepared_ligand_forcefield(str(run), 'gaff2')
    assert resolve_ligand_forcefield([str(tmp_path / 'md_files')], 'gaff', explicit=False) == 'gaff2'


def test_resolve_forcefield_keeps_requested_when_explicit(tmp_path):
    """An explicit --ligand_forcefield is honored even if a different one is recorded."""
    run = tmp_path / 'md_files' / 'md_run' / 'complex_replica1'
    run.mkdir(parents=True)
    lp.write_prepared_ligand_forcefield(str(run), 'gaff2')
    assert resolve_ligand_forcefield([str(tmp_path / 'md_files')], 'gaff', explicit=True) == 'gaff'


def test_resolve_forcefield_default_when_nothing_recorded(tmp_path):
    (tmp_path / 'md_files').mkdir()
    assert resolve_ligand_forcefield([str(tmp_path / 'md_files')], 'gaff', explicit=False) == 'gaff'


def test_resolve_forcefield_ignores_backup_records(tmp_path):
    """A preserved backup_<ff> snapshot must not make resolution ambiguous."""
    lig = tmp_path / 'md_files' / 'md_preparation' / 'ligands' / 'LIG'
    (lig / 'backup_gaff').mkdir(parents=True)
    lp.write_prepared_ligand_forcefield(str(lig), 'gaff2')             # current
    lp.write_prepared_ligand_forcefield(str(lig / 'backup_gaff'), 'gaff')  # preserved old
    assert resolve_ligand_forcefield([str(tmp_path / 'md_files')], 'gaff', explicit=False) == 'gaff2'


def test_resolve_forcefield_logs_when_recorded_matches_implicit(tmp_path, caplog):
    """Omitting the flag and matching the recorded FF is logged (reuse confirmation)."""
    run = tmp_path / 'md_files' / 'md_run' / 'c_replica1'
    run.mkdir(parents=True)
    lp.write_prepared_ligand_forcefield(str(run), 'gaff')
    with caplog.at_level(logging.INFO):
        assert resolve_ligand_forcefield([str(tmp_path / 'md_files')], 'gaff', explicit=False) == 'gaff'
    assert any('Reusing previously recorded ligand force field' in r.getMessage() for r in caplog.records)


def test_resolve_forcefield_logs_when_recorded_matches_explicit(tmp_path, caplog):
    """An explicit choice that matches the recorded FF is logged."""
    run = tmp_path / 'md_files' / 'md_run' / 'c_replica1'
    run.mkdir(parents=True)
    lp.write_prepared_ligand_forcefield(str(run), 'gaff')
    with caplog.at_level(logging.INFO):
        assert resolve_ligand_forcefield([str(tmp_path / 'md_files')], 'gaff', explicit=True) == 'gaff'
    assert any('matches the record' in r.getMessage() for r in caplog.records)


def test_resolve_forcefield_ambiguous_raises(tmp_path):
    r1 = tmp_path / 'md_files' / 'md_run' / 'c_replica1'
    r2 = tmp_path / 'md_files' / 'md_run' / 'c_replica2'
    r1.mkdir(parents=True)
    r2.mkdir(parents=True)
    lp.write_prepared_ligand_forcefield(str(r1), 'gaff')
    lp.write_prepared_ligand_forcefield(str(r2), 'gaff2')
    with pytest.raises(ValueError):
        resolve_ligand_forcefield([str(tmp_path / 'md_files')], 'gaff', explicit=False)


def test_partial_ligand_cache_is_detected_and_backed_up(tmp_path):
    """An interrupted preparation (only .mol2/.frcmod, no metadata) is detected and backed up.

    Directly exercises the partial-cache branch the changelog promises: an interrupted
    run leaves intermediate files and no metadata record.
    """
    molid = "ligand"
    ligand_dir = tmp_path / molid
    ligand_dir.mkdir()

    mol2 = ligand_dir / f"{molid}.mol2"
    frcmod = ligand_dir / f"{molid}.frcmod"
    mol2.write_text("partial GAFF MOL2")
    frcmod.write_text("partial GAFF parameters")

    assert has_existing_ligand_files(str(ligand_dir), molid)

    backup_dir = backup_stale_ligand_files(str(ligand_dir), molid, "gaff")

    assert backup_dir is not None
    # removed from the active directory ...
    assert not mol2.exists()
    assert not frcmod.exists()
    # ... and preserved (moved, not deleted) in the backup subdirectory
    assert os.path.isfile(os.path.join(backup_dir, mol2.name))
    assert os.path.isfile(os.path.join(backup_dir, frcmod.name))
