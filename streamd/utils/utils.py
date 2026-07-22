"""Generic helper functions used across the StreaMD project."""

from contextlib import contextmanager
from glob import glob
import logging
import os
import re
import shlex
import shutil
import subprocess
from tempfile import mkdtemp
from typing import Iterable, Tuple, Dict, Any

import argparse
import yaml

import MDAnalysis as mda


def parse_with_config(parser: argparse.ArgumentParser, cli_args: Iterable[str]) -> Tuple[argparse.Namespace, Dict[str, Any]]:
    """Parse ``cli_args`` using ``parser`` honouring an optional ``--config`` file.

    Parameters
    ----------
    parser:
        Configured ``argparse.ArgumentParser`` that includes a ``--config`` option.
    cli_args:
        Sequence of command-line arguments excluding the program name.

    Returns
    -------
    Tuple[argparse.Namespace, Dict[str, Any]]
        The parsed arguments and the subset of configuration values that were
        applied as defaults.
    """

    args, _ = parser.parse_known_args(cli_args)
    config_args: Dict[str, Any] = {}
    if getattr(args, "config", None):
        with open(args.config) as fh:
            config_args = yaml.safe_load(fh) or {}
        if not isinstance(config_args, dict):
            raise ValueError("Config file must contain key-value pairs")
        valid_keys = {action.dest for action in parser._actions}
        cli_dests = {
            action.dest
            for action in parser._actions
            if any(opt in cli_args for opt in action.option_strings)
        }
        config_args = {
            k: v for k, v in config_args.items() if k in valid_keys and k not in cli_dests
        }

        # Convert configuration values according to parser specifications so that
        # defaults set from YAML mimic CLI parsing behaviour.  This handles
        # ``nargs`` options as well as applying any ``type`` conversions.
        for action in parser._actions:
            dest = action.dest
            if dest not in config_args:
                continue
            value = config_args[dest]

            #argparse represents store_true/store_false actions with nargs set to 0

            if action.nargs not in (None, '?', 0):
                if isinstance(value, str):
                    value = value.split()
                elif not isinstance(value, (list, tuple)):
                    value = [value]
                if action.type:
                    value = [action.type(v) for v in value]
            else:
                if action.type:
                    value = action.type(value)

            config_args[dest] = value

        parser.set_defaults(**config_args)

    return parser.parse_args(cli_args), config_args


def filepath_type(x: str, ext=None, check_exist=True, exist_type='file', create_dir=False):
    """Validate file paths and optionally ensure existence and extension.

    :param x: Path to validate.
    :param ext: Iterable of allowed file extensions without leading dots.
    :param check_exist: If ``True``, ensure the path exists.
    :param exist_type: Expectation for existing path: ``'file'`` or ``'dir'``.
    :param create_dir: Create the directory if it does not exist and ``x`` is a directory path.
    :return: Absolute path to the validated file or directory.
    """
    # str - parse config arguments when used numbers
    value = os.path.abspath(str(x)) if x else x
    if create_dir:
        os.makedirs(value, exist_ok=True)
    if check_exist:
        if exist_type == 'file' and not os.path.isfile(value):
            raise FileExistsError(f'{value} does not exist')
        if exist_type == 'dir' and not os.path.isdir(value):
            raise NotADirectoryError(f'{value} directory does not exist')
    if ext and os.path.splitext(value)[1].replace('.', '') not in ext:
        raise FileExistsError(f'File {value} has a wrong extension. Allowed are {", ".join(ext)}')
    return value


def get_index(index_file, env=None):
    """Return group names from a GROMACS index file.

    :param index_file: Path to the ``index.ndx`` file.
    :param env: Optional environment variables for subprocess calls.
    :return: List of index group names.
    """
    with open(index_file) as input:
        data = input.read()

    if not data:
        create_ndx(index_file, env=env)
        with open(index_file) as input:
            data = input.read()

    groups = [i.strip() for i in re.findall(r'\[(.*)\]', data)]
    return groups


def make_group_ndx(query, wdir, bash_log, env=None):
    """Create a new index group using ``gmx make_ndx``.

    :param query: ``make_ndx`` selection string.
    :param wdir: Working directory containing ``index.ndx``.
    :param bash_log: Log file capturing shell output.
    :param env: Optional environment variables for subprocess calls.
    :return: ``True`` on success, ``False`` otherwise.
    """
    cmd = f'''
        cd {wdir}
        gmx make_ndx -f solv_ions.gro -n index.ndx << INPUT  >> {os.path.join(wdir,bash_log)} 2>&1
        {query}
        q
        INPUT
        '''
    if not run_check_subprocess(cmd, key=wdir, log=os.path.join(wdir,bash_log), env=env):
       return False

    return True

def create_ndx(index_file, env=None):
    """Generate an index file with ``gmx make_ndx`` if absent.

    :param index_file: Path where the index file should exist.
    :param env: Optional environment variables for subprocess calls.
    :return: ``None``.
    """
    wdir = os.path.dirname(index_file)
    cmd = f'''
        cd {wdir}
        gmx make_ndx -f {os.path.join(wdir, "solv_ions.gro")} -o {index_file} << INPUT 
        q
        INPUT
        '''
    run_check_subprocess(cmd, key=wdir, log=None, env=env)

def get_mol_resid_pair(fname):
    """Yield molecule and residue identifiers from a text file.

    :param fname: File containing tab-separated molecule and residue pairs.
    :return: Generator yielding ``(molid, resid)`` tuples.
    """
    with open(fname) as inp:
        data = inp.readlines()
    for molid_resid_pair in data:
        pair = [i.strip() for i in molid_resid_pair.split('\t')]
        if pair:
            molid, resid = pair
            yield molid, resid

def run_check_subprocess(cmd, key=None, log=None, env=None, ignore_error=False):
    """Run a shell command and log failures.

    :param cmd: Command string executed via ``subprocess.check_output``.
    :param key: Identifier used in log messages.
    :param log: Path to a log file to mention when the command fails.
    :param env: Optional environment variables for the subprocess.
    :param ignore_error: Suppress raising/logging errors if ``True``.
    :return: ``True`` if the command succeeded, ``False`` otherwise.
    """
    try:
        subprocess.check_output(cmd, shell=True, env=env, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        if log:
            logging.warning(f'Failed run for {key}. Check log {log}\n')
        if not ignore_error:
            logging.exception(f'Failed run for {key}. Error:{e}', stack_info=True)
        return False
    return True

def get_protein_resid_set(protein_fname):
    """Return set of residue names present in a protein structure.

    :param protein_fname: Path to a protein PDB or GRO file.
    :return: Set of residue names found in the structure.
    """
    protein = mda.Universe(protein_fname)
    protein_resid_set = set(protein.residues.resnames.tolist())
    return protein_resid_set

def get_number_of_frames(xtc, env):
    """Return number of frames and timestep from an XTC trajectory.

    :param xtc: Path to the trajectory ``.xtc`` file.
    :param env: Optional environment variables for subprocess calls.
    :return: Tuple of ``(frames, timestep)`` or ``None`` if parsing fails.
    """
    res = subprocess.run(f'gmx check -f {xtc}', shell=True, capture_output=True, env=env)

    res_parsed = re.findall(r'Step\s*(\d*)\s*(\d*)\n', res.stderr.decode("utf-8"))
    if res_parsed:
        frames, timestep = res_parsed[0]
        # starts with 0
        logging.info(f'{xtc} has {int(frames)} frames ({timestep} timestep)')
        return int(frames), int(timestep) if timestep else None
    else:
        logging.warning(f'Failed to read number of frames of {xtc} trajectory. {res}')
        return None

def get_last_frame_time(xtc, env=None):
    """Return the final trajectory time in ps parsed from ``gmx check``.

    Two output shapes are handled:

    * Some GROMACS versions print an explicit ``Last frame  N time  T`` line; when
      present its time is used directly (this is exact, including for merged
      trajectories with a non-uniform output interval).
    * ``gmx check`` in other versions (e.g. GROMACS 2023.x) prints **no**
      ``Last frame`` line, only periodic ``Reading frame N time T`` progress lines
      plus a summary table (``Step  #frames  timestep``). The last progress line is
      not necessarily the final frame, because gmx prints them at intervals. The
      final time is then reconstructed by anchoring on the last observed
      ``(index, time)`` and extending to the last frame (``#frames - 1``) with the
      reported timestep.

    :param xtc: Path to the trajectory ``.xtc`` file.
    :param env: Optional environment variables for subprocess calls.
    :return: Final trajectory time in picoseconds as ``float``.
    :raises RuntimeError: If ``gmx check`` fails (nonzero return code) or the final
        time cannot be parsed from its output.
    """
    # Time literal: integer, decimal or scientific notation.
    time_re = r'([-+]?\d*\.?\d+(?:[eE][+-]?\d+)?)'
    res = subprocess.run(f'gmx check -f {shlex.quote(str(xtc))}', shell=True,
                         capture_output=True, env=env)
    # gmx may write the relevant lines to stdout or stderr depending on the version.
    combined_output = res.stdout.decode('utf-8', errors='replace') + '\n' + \
                      res.stderr.decode('utf-8', errors='replace')
    if res.returncode != 0:
        raise RuntimeError(
            f'"gmx check" failed for {xtc} (return code {res.returncode}). '
            f'GROMACS output:\n{combined_output}')
    # Preferred path: an explicit "Last frame" line (variable whitespace tolerated).
    # gmx may emit several progress lines; use the last reported "Last frame" time.
    last_frame = re.findall(rf'Last frame\s+\d+\s+time\s+{time_re}', combined_output)
    if last_frame:
        return float(last_frame[-1])
    # Fallback: reconstruct from "Reading frame N time T" progress lines and the
    # "Step  #frames  timestep" summary row (GROMACS versions with no "Last frame").
    progress = re.findall(rf'(?:Reading|Last) frame\s+(\d+)\s+time\s+{time_re}',
                          combined_output)
    table = re.search(r'Step\s+(\d+)\s+(\d+)', combined_output)
    if progress and table:
        last_idx, last_time = int(progress[-1][0]), float(progress[-1][1])
        nframes, timestep = int(table.group(1)), int(table.group(2))
        # Frames are 0-indexed: extend from the last printed frame to the final one.
        return last_time + max(0, (nframes - 1) - last_idx) * timestep
    raise RuntimeError(
        f'Cannot determine the final trajectory time of {xtc} from "gmx check". '
        f'GROMACS output:\n{combined_output}')

def backup_and_replace(src_file, target_file, copy=False):
    backup_prev_files(target_file)
    if not copy:
        shutil.move(src_file, target_file)
    else:
        shutil.copy(src_file, target_file)


def backup_prev_files(file_to_backup, wdir=None, copy=False):
    """Rename or copy an existing file to avoid overwriting.

    :param file_to_backup: File path that may be overwritten.
    :param wdir: Directory where backups are stored. Defaults to the file's directory.
    :param copy: If ``True``, copy instead of moving the file.
    :return: ``None``.
    """
    if not os.path.isfile(file_to_backup):
        return None
    if wdir is None:
        wdir = os.path.dirname(file_to_backup)
    n = len(glob(os.path.join(wdir, f'#{os.path.basename(file_to_backup)}.*#'))) + 1
    new_f = os.path.join(wdir, f'#{os.path.basename(file_to_backup)}.{n}#')
    if not copy:
        shutil.move(file_to_backup, new_f)
    else:
        shutil.copy(file_to_backup, new_f)
    logging.warning(f'Backup previous file {file_to_backup} to {new_f}')

def check_to_continue_simulation_time(xtc, new_mdtime_ps, env):
    """Check whether a trajectory reached the desired simulation time.

    The current length is taken from the actual final frame time reported by
    ``gmx check`` (:func:`get_last_frame_time`), not from
    ``number_of_frames * timestep``. The latter assumes a single uniform output
    interval and can be wrong for concatenated/continued trajectories whose parts
    used different output intervals (or contain a duplicate frame at the join).

    :param xtc: Path to the trajectory ``.xtc`` file.
    :param new_mdtime_ps: Desired simulation time in picoseconds.
    :param env: Optional environment variables for subprocess calls.
    :return: ``False`` if the desired time is already reached, otherwise ``True``.
    """
    try:
        current_time_ps = get_last_frame_time(xtc=xtc, env=env)
    except RuntimeError as e:
        # If the current length cannot be determined, do not block continuation.
        logging.warning(f'Could not determine the current length of {xtc}; '
                        f'continuation will proceed. {e}')
        return True
    logging.info(f'The length of the found trajectory is {current_time_ps / 1000} ns. '
                 f'The desired length is {new_mdtime_ps / 1000} ns.')
    if current_time_ps >= new_mdtime_ps:
        logging.warning(f'The desired length of the found simulation trajectory {xtc} has been already reached. '
                        f'Calculations will be interrupted.')
        return False
    return True

def merge_parts_of_simulation(start_xtc, part_xtc, new_xtc, wdir, bash_log, env=None):
    """Concatenate trajectory parts using ``gmx trjcat``.

    :param start_xtc: Initial trajectory file.
    :param part_xtc: Additional trajectory fragment to append.
    :param new_xtc: Output trajectory file name.
    :param wdir: Working directory containing the trajectory files.
    :param bash_log: Log file capturing shell output.
    :param env: Optional environment variables for subprocess calls.
    :return: ``None``.
    """
    cmd =f'''
cd {wdir}
gmx trjcat -f {start_xtc} {part_xtc} -o {new_xtc} -tu fs >> {os.path.join(wdir,bash_log)} 2>&1
    '''
    return run_check_subprocess(cmd, key=part_xtc, log=os.path.join(wdir, bash_log), env=env)

def create_last_frame_file(wdir, tpr, xtc, out_file, bash_log, env, center_group=None, index=None, dump_time=None):
    """Extract a single trajectory frame into a structure file.

    By default the *final* frame is written: it is selected by the final
    trajectory *time* (ps) reported by :func:`get_last_frame_time`, because
    ``gmx trjconv -dump`` expects a time, not a frame index. Pass ``dump_time`` to
    extract a different frame instead (e.g. ``0`` for the first/start frame).

    When ``center_group`` and ``index`` are provided (used for ``last_frame.pdb`` and
    ``start_frame.pdb``) the snapshot is made visualization-ready: molecules are
    reconstructed across periodic boundaries (``-pbc mol``), ``center_group`` is centered
    (``-center``) and the whole ``System`` is written in a compact unit cell
    (``-ur compact``). No rotational fit is applied. Without those arguments the raw frame
    of ``System`` is written unchanged (used for the continuation ``.gro``).

    :param wdir: directory containing the trajectory files.
    :param tpr: Path to the run-input ``.tpr`` file.
    :param xtc: Path to the trajectory ``.xtc`` file.
    :param out_file: Path of the structure file to write.
    :param bash_log: Log file capturing shell output.
    :param env: Optional environment variables for subprocess calls.
    :param center_group: Index group (name or numeric id) to center on; enables PBC
        reconstruction, centering and compact unit-cell output when given together with
        ``index``.
    :param index: Path to the ``index.ndx`` file defining ``center_group``.
    :param dump_time: Trajectory time in ps to extract. Defaults to the final frame time
        determined from ``gmx check``.
    :return: ``None``.
    """
    if dump_time is None:
        last_time = get_last_frame_time(xtc=xtc, env=env)
        # keep an integer literal when the time is integral (e.g. 100000, not 100000.0)
        dump_time = int(last_time) if float(last_time).is_integer() else last_time
    log_path = os.path.join(wdir, bash_log)
    if center_group is not None and index is not None:
        cmd = f'''
    cd {shlex.quote(wdir)}
    printf '%s\\n%s\\n' "{center_group}" "System" | gmx trjconv -s {shlex.quote(tpr)} -f {shlex.quote(xtc)} -o {shlex.quote(out_file)} -dump {dump_time} -pbc mol -center -ur compact -n {shlex.quote(index)} >> {shlex.quote(log_path)} 2>&1
        '''
    else:
        cmd = f'''
    cd {shlex.quote(wdir)}
    printf '%s\\n' "System" | gmx trjconv -s {shlex.quote(tpr)} -f {shlex.quote(xtc)} -o {shlex.quote(out_file)} -dump {dump_time} >> {shlex.quote(log_path)} 2>&1
        '''
    run_check_subprocess(cmd, key=out_file, log=log_path, env=env)

@contextmanager
def temporary_directory_debug(remove=True, suffix=None, dir=None):
    """Create a temporary directory removed after use.

    :param remove: Whether to remove the directory upon exit.
    :param suffix: Optional suffix for the temporary directory name.
    :param dir: Parent directory in which to create the temporary directory.
    :return: Context manager yielding the path to the temporary directory.
    """
    if dir is None:
        dir = os.path.curdir
    path = os.path.abspath(mkdtemp(dir=dir, suffix=suffix))
    # os.makedirs(path, exist_ok=True)
    try:
        yield path
    finally:
        if remove:
            try:
                shutil.rmtree(path)
            except OSError as e:
                logging.error(f'\nCould not remove the tmp directory: {path}. Error: {e}\n')
                # subprocess.call(f'rm -r {path}', shell=True)
                shutil.rmtree(path, ignore_errors=True)
