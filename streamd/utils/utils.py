from glob import glob
import logging
import os
import re
import shutil
import subprocess

import MDAnalysis as mda


def filepath_type(x, ext=None, check_exist=True, exist_type='file', create_dir=False):
    value = os.path.abspath(x) if x else x
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
    with open(index_file) as input:
        data = input.read()

    if not data:
        create_ndx(index_file, env=env)
        with open(index_file) as input:
            data = input.read()

    groups = [i.strip() for i in re.findall(r'\[(.*)\]', data)]
    return groups


def make_group_ndx(query, wdir, bash_log, env=None):
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
    wdir = os.path.dirname(index_file)
    cmd = f'''
        cd {wdir}
        gmx make_ndx -f {os.path.join(wdir, "solv_ions.gro")} -o {index_file} << INPUT 
        q
        INPUT
        '''
    run_check_subprocess(cmd, key=wdir, log=None, env=env)

def get_mol_resid_pair(fname):
    with open(fname) as inp:
        data = inp.readlines()
    for molid_resid_pair in data:
        pair = [i.strip() for i in molid_resid_pair.split('\t')]
        if pair:
            molid, resid = pair
            yield molid, resid

def run_check_subprocess(cmd, key=None, log=None, env=None, ignore_error=False):
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
    protein = mda.Universe(protein_fname)
    protein_resid_set = set(protein.residues.resnames.tolist())
    return protein_resid_set

def get_number_of_frames(xtc, env):
    res = subprocess.run(f'gmx check -f {xtc}', shell=True, capture_output=True, env=env)
    res_parsed = re.findall('Step[ ]*([0-9]*)[ ]*([0-9]*)\n', res.stderr.decode("utf-8"))
    if res_parsed:
        frames, timestep = res_parsed[0]
        # starts with 0
        logging.info(f'{xtc} has {int(frames)} frames')
        return int(frames), int(timestep)
    else:
        logging.warning(f'Failed to read number of frames of {xtc} trajectory. {res}')
        return None

def backup_prev_files(file_to_backup, wdir=None, copy=False):
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
    current_number_of_frames, timestep = get_number_of_frames(xtc=xtc, env=env)
    if current_number_of_frames and timestep:
        time_ns = (current_number_of_frames*timestep-timestep)/1000
        logging.info(f'The length of the found trajectory is {time_ns} ns. Should be continued until {new_mdtime_ps/1000} ns.')
        if current_number_of_frames * timestep >= new_mdtime_ps:
            logging.warning(f'The desired length of the found simulation trajectory {xtc} has been already reached. '
                            f'Calculations will be interrupted.')
            return False
    return True

def merge_parts_of_simulation(start_xtc, part_xtc, new_xtc, wdir, bash_log, env=None):
    cmd =f'''
cd {wdir}
gmx trjcat -f {start_xtc} {part_xtc} -o {new_xtc} -tu fs >> {os.path.join(wdir,bash_log)} 2>&1
    '''
    run_check_subprocess(cmd, key=part_xtc, log=os.path.join(wdir, bash_log), env=env)
