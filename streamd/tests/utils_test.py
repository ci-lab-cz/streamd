"""Tests for last-frame snapshot generation helpers in :mod:`streamd.utils.utils`."""

import pytest

from streamd.utils import utils


class _FakeCompleted:
    """Minimal stand-in for ``subprocess.CompletedProcess`` with byte streams."""

    def __init__(self, stdout=b'', stderr=b'', returncode=0):
        self.stdout = stdout
        self.stderr = stderr
        self.returncode = returncode


def _patch_gmx_check(monkeypatch, stdout=b'', stderr=b'', returncode=0):
    """Make ``gmx check`` return canned output without invoking GROMACS."""
    def fake_run(cmd, *args, **kwargs):
        return _FakeCompleted(stdout=stdout, stderr=stderr, returncode=returncode)

    monkeypatch.setattr(utils.subprocess, 'run', fake_run)


# --------------------------------------------------------------------------- #
# get_last_frame_time parser
# --------------------------------------------------------------------------- #

def test_last_frame_time_representative(monkeypatch):
    """The representative gmx check line yields the final time in ps."""
    _patch_gmx_check(monkeypatch, stderr=b'Last frame      10000 time 100000.000\n')
    assert utils.get_last_frame_time('traj.xtc') == 100000.0


def test_last_frame_time_scientific_notation(monkeypatch):
    """Scientific-notation times are parsed."""
    _patch_gmx_check(monkeypatch, stderr=b'Last frame 10000 time 1.00000e+05\n')
    assert utils.get_last_frame_time('traj.xtc') == 100000.0


def test_last_frame_time_variable_whitespace(monkeypatch):
    """Tabs and irregular spacing between fields are tolerated."""
    _patch_gmx_check(monkeypatch, stderr=b'Last frame\t10000\ttime\t100000\n')
    assert utils.get_last_frame_time('traj.xtc') == 100000.0


def test_last_frame_time_from_stdout(monkeypatch):
    """The final time is found when gmx writes it to stdout."""
    _patch_gmx_check(monkeypatch, stdout=b'Last frame 10000 time 100000.000\n')
    assert utils.get_last_frame_time('traj.xtc') == 100000.0


def test_last_frame_time_from_stderr(monkeypatch):
    """The final time is found when gmx writes it to stderr."""
    _patch_gmx_check(monkeypatch, stderr=b'Last frame 10000 time 100000.000\n')
    assert utils.get_last_frame_time('traj.xtc') == 100000.0


def test_last_frame_time_uses_last_occurrence(monkeypatch):
    """When gmx prints several 'Last frame' lines, the final one wins."""
    _patch_gmx_check(
        monkeypatch,
        stderr=(
            b'Last frame 5000 time 50000.000\n'
            b'Last frame 8000 time 80000.000\n'
            b'Last frame 10000 time 100000.000\n'
        ),
    )
    assert utils.get_last_frame_time('traj.xtc') == 100000.0


def test_last_frame_time_real_gmx_check_format(monkeypatch):
    """Parse output shaped like real 'gmx check' (2025.2): 'Checking file' on stdout,
    frames and the carriage-return-overwritten 'Last frame' line on stderr.

    The 'Last frame' entry shares one physical line with many 'Reading frame ... time'
    progress entries; the parser must isolate 'Last frame' and not a 'Reading frame'.
    """
    stdout = b'Checking file /path/md_out.xtc\n'
    stderr = (
        b'Reading frame       0 time    0.000   \n'
        b'# Atoms  67796\n'
        b'Precision 0.001 (nm)\n'
        b'Reading frame       1 time   10.000   '
        b'Reading frame       2 time   20.000   '
        b'Reading frame      10 time  100.000   '
        b'Last frame         10 time  100.000   \n'
        b'\n'
        b'Item        #frames Timestep (ps)\n'
        b'Step            11    10\n'
    )
    _patch_gmx_check(monkeypatch, stdout=stdout, stderr=stderr, returncode=0)
    assert utils.get_last_frame_time('md_out.xtc') == 100.0


def test_last_frame_time_gromacs_2023_no_last_frame_line(monkeypatch):
    """Real 'gmx check' (GROMACS 2023.4) prints no 'Last frame' line.

    Only periodic 'Reading frame' progress lines and the summary 'Step' row are
    emitted. The last progress line (frame 90 @ 900 ps) is NOT the final frame:
    the trajectory has 92 frames (indices 0-91, 10 ps step), so the true final
    time is 910 ps. Reproduces the reported crash scenario.
    """
    stderr = (
        b'Reading frame       0 time    0.000   \n'
        b'# Atoms  67796\n'
        b'Precision 0.001 (nm)\n'
        b'Reading frame      90 time  900.000   \n'
        b'\n'
        b'Item        #frames Timestep (ps)\n'
        b'Step            92    10\n'
        b'Time            92    10\n'
        b'Coords          92    10\n'
        b'Box             92    10\n'
    )
    _patch_gmx_check(monkeypatch, stdout=b'Checking file md_out.xtc\n', stderr=stderr)
    assert utils.get_last_frame_time('md_out.xtc') == 910.0


def test_last_frame_time_2023_single_progress_line(monkeypatch):
    """Fallback works when only the first frame progress line is printed.

    Extends frame 0 @ 0 ps over 11 frames (10 ps step) -> final time 100 ps.
    """
    stderr = (
        b'Reading frame       0 time    0.000   \n'
        b'\n'
        b'Item        #frames Timestep (ps)\n'
        b'Step            11    10\n'
    )
    _patch_gmx_check(monkeypatch, stderr=stderr)
    assert utils.get_last_frame_time('md_out.xtc') == 100.0


def test_last_frame_time_missing_line_raises(monkeypatch):
    """No 'Last frame' line and no summary table raises a clear error with gmx output."""
    _patch_gmx_check(monkeypatch, stdout=b'Reading frame 0 time 0.000\n')
    with pytest.raises(RuntimeError) as excinfo:
        utils.get_last_frame_time('traj.xtc')
    assert 'Reading frame 0 time 0.000' in str(excinfo.value)


def test_last_frame_time_failed_execution_raises(monkeypatch):
    """A failed gmx check execution raises with the captured output."""
    _patch_gmx_check(
        monkeypatch,
        stderr=b'Fatal error: cannot read file traj.xtc\n',
        returncode=1,
    )
    with pytest.raises(RuntimeError) as excinfo:
        utils.get_last_frame_time('traj.xtc')
    assert 'Fatal error' in str(excinfo.value)


def test_last_frame_time_nonzero_returncode_raises_even_with_valid_line(monkeypatch):
    """A nonzero return code fails clearly even if a parseable line is present."""
    _patch_gmx_check(
        monkeypatch,
        stderr=b'Last frame 10000 time 100000.000\nError: trajectory is truncated\n',
        returncode=1,
    )
    with pytest.raises(RuntimeError) as excinfo:
        utils.get_last_frame_time('traj.xtc')
    message = str(excinfo.value)
    assert 'return code 1' in message
    assert 'trajectory is truncated' in message


# --------------------------------------------------------------------------- #
# create_last_frame_file command generation
# --------------------------------------------------------------------------- #

def _capture_cmd(monkeypatch, last_time):
    """Stub time lookup and command runner; return a list capturing the cmd."""
    commands = []
    monkeypatch.setattr(utils, 'get_last_frame_time', lambda **kwargs: last_time)
    monkeypatch.setattr(
        utils, 'run_check_subprocess',
        lambda cmd, *a, **k: commands.append(cmd) or True,
    )
    return commands


def test_create_last_frame_file_visualization_command(monkeypatch):
    """last_frame.pdb generation adds PBC/centering and the correct dump time."""
    commands = _capture_cmd(monkeypatch, last_time=100000.0)

    utils.create_last_frame_file(
        wdir='/wd', tpr='md.tpr', xtc='md.xtc',
        out_file='last_frame.pdb', bash_log='bash.log', env=None,
        center_group=6, index='index.ndx',
    )

    assert len(commands) == 1
    cmd = commands[0]
    assert '-dump 100000' in cmd
    assert '-pbc mol' in cmd
    assert '-center' in cmd
    assert '-ur compact' in cmd
    # selection input: center group first, then System
    assert cmd.index('"6"') < cmd.index('"System"')


def test_create_start_frame_file_dumps_time_zero_with_pbc(monkeypatch):
    """start_frame.pdb uses -dump 0 and the same visualization flags as last_frame."""
    commands = []
    # dump_time is passed explicitly, so get_last_frame_time must not be needed;
    # make it raise to prove it is never consulted for the start frame.
    def _boom(**kwargs):
        raise AssertionError('get_last_frame_time should not be called when dump_time is set')
    monkeypatch.setattr(utils, 'get_last_frame_time', _boom)
    monkeypatch.setattr(
        utils, 'run_check_subprocess',
        lambda cmd, *a, **k: commands.append(cmd) or True,
    )

    utils.create_last_frame_file(
        wdir='/wd', tpr='md.tpr', xtc='md.xtc',
        out_file='start_frame.pdb', bash_log='bash.log', env=None,
        center_group=6, index='index.ndx', dump_time=0,
    )

    assert len(commands) == 1
    cmd = commands[0]
    assert '-dump 0' in cmd
    assert '-pbc mol' in cmd
    assert '-center' in cmd
    assert '-ur compact' in cmd
    assert cmd.index('"6"') < cmd.index('"System"')
    assert 'start_frame.pdb' in cmd


def test_create_last_frame_file_plain_command_no_centering(monkeypatch):
    """Without a center group (continuation .gro) no PBC/centering is applied."""
    commands = _capture_cmd(monkeypatch, last_time=100000.0)

    utils.create_last_frame_file(
        wdir='/wd', tpr='md.tpr', xtc='md.xtc',
        out_file='md_out.gro', bash_log='bash.log', env=None,
    )

    assert len(commands) == 1
    cmd = commands[0]
    assert '-dump 100000' in cmd
    assert '-pbc mol' not in cmd
    assert '-center' not in cmd
    assert '-ur compact' not in cmd


def test_create_last_frame_file_continuation_gro_uses_final_time_no_pbc(monkeypatch):
    """Continuation .gro dumps the actual final time but gets no visualization flags."""
    # Real parser over gmx check output: 10001-frame trajectory, final time 100000 ps.
    _patch_gmx_check(monkeypatch, stderr=b'Last frame      10000 time 100000.000\n')
    commands = []
    monkeypatch.setattr(
        utils, 'run_check_subprocess',
        lambda cmd, *a, **k: commands.append(cmd) or True,
    )

    # No center_group / index -> continuation .gro path (run_md.py call site).
    utils.create_last_frame_file(
        wdir='/wd', tpr='md.tpr', xtc='md.xtc',
        out_file='md_out.gro', bash_log='bash.log', env=None,
    )

    assert len(commands) == 1
    cmd = commands[0]
    # actual final time, not the frame count
    assert '-dump 100000' in cmd
    assert '-dump 10001' not in cmd
    # no visualization-oriented PBC/centering flags on the production .gro
    assert '-pbc mol' not in cmd
    assert '-center' not in cmd
    assert '-ur compact' not in cmd
    # writes only the System group (single selection line)
    assert '"System"' in cmd
    assert 'md_out.gro' in cmd


def test_create_last_frame_file_regression_uses_time_not_frame_count(monkeypatch):
    """A 10001-frame / 100000 ps trajectory dumps at time 100000, not frame 10001."""
    # Real parser over representative gmx check output (last frame index 10000).
    _patch_gmx_check(monkeypatch, stderr=b'Last frame      10000 time 100000.000\n')
    commands = []
    monkeypatch.setattr(
        utils, 'run_check_subprocess',
        lambda cmd, *a, **k: commands.append(cmd) or True,
    )

    utils.create_last_frame_file(
        wdir='/wd', tpr='md.tpr', xtc='md.xtc',
        out_file='last_frame.pdb', bash_log='bash.log', env=None,
        center_group=6, index='index.ndx',
    )

    assert len(commands) == 1
    cmd = commands[0]
    assert '-dump 100000' in cmd
    assert '-dump 10001' not in cmd


# --------------------------------------------------------------------------- #
# check_to_continue_simulation_time
# --------------------------------------------------------------------------- #

def test_check_to_continue_stops_when_target_reached(monkeypatch):
    """Continuation is stopped once the trajectory reaches the desired time."""
    _patch_gmx_check(monkeypatch, stderr=b'Last frame  20000 time 200000.000\n')
    assert utils.check_to_continue_simulation_time('md.xtc', 200000, env=None) is False


def test_check_to_continue_continues_when_target_not_reached(monkeypatch):
    """Continuation proceeds while the trajectory is shorter than requested."""
    _patch_gmx_check(monkeypatch, stderr=b'Last frame  5000 time 50000.000\n')
    assert utils.check_to_continue_simulation_time('md.xtc', 200000, env=None) is True


def test_check_to_continue_uses_last_frame_time_not_frame_count(monkeypatch):
    """The decision follows the real final time, not number_of_frames * timestep.

    Emulates a merged trajectory whose true final time is 150000 ps (< target, so
    it must continue). The 'Step' table row (30001 frames x 10 ps) would make the
    old ``frames * timestep`` heuristic compute 300010 ps and wrongly stop.
    """
    stderr = (
        b'Reading frame       0 time    0.000   '
        b'Last frame      30000 time 150000.000   \n'
        b'\n'
        b'Item        #frames Timestep (ps)\n'
        b'Step         30001    10\n'
    )
    _patch_gmx_check(monkeypatch, stderr=stderr)
    # true final time 150000 ps < 200000 ps target -> must continue
    assert utils.check_to_continue_simulation_time('md.xtc', 200000, env=None) is True


def test_check_to_continue_gromacs_2023_no_last_frame_line(monkeypatch):
    """Continuation decision works on GROMACS 2023.x output lacking a 'Last frame' line.

    Reconstructed final time is 910 ps (< 200000 ps target) -> must continue.
    """
    stderr = (
        b'Reading frame       0 time    0.000   \n'
        b'Reading frame      90 time  900.000   \n'
        b'\n'
        b'Item        #frames Timestep (ps)\n'
        b'Step            92    10\n'
    )
    _patch_gmx_check(monkeypatch, stderr=stderr)
    assert utils.check_to_continue_simulation_time('md.xtc', 200000, env=None) is True


def test_check_to_continue_defaults_to_true_when_time_unknown(monkeypatch):
    """If the final time cannot be determined, continuation is not blocked."""
    _patch_gmx_check(monkeypatch, stderr=b'Fatal error: cannot read file\n', returncode=1)
    assert utils.check_to_continue_simulation_time('md.xtc', 200000, env=None) is True
