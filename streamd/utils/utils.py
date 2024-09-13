import logging
import os
import re
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
    if ext and os.path.splitext(value)[1].lstrip('.') not in ext:
        raise FileExistsError(f'File {value} has a wronh extension. Allowed are {", ".join(ext)}')
    return value


def get_index(index_file):
    with open(index_file) as input:
        data = input.read()

    if not data:
        create_ndx(index_file)
        with open(index_file) as input:
            data = input.read()

    groups = [i.strip() for i in re.findall('\[(.*)\]', data)]
    return groups


def make_group_ndx(query, wdir, bash_log):
    cmd = f'''
        cd {wdir}
        gmx make_ndx -f solv_ions.gro -n index.ndx << INPUT  >> {os.path.join(wdir,bash_log)} 2>&1
        {query}
        q
        INPUT
        '''
    if not run_check_subprocess(cmd, key=wdir, log=os.path.join(wdir,bash_log)):
       return False

    return True

def create_ndx(index_file):
    wdir = os.path.dirname(index_file)
    cmd = f'''
        cd {wdir}
        gmx make_ndx -f {os.path.join(wdir, "solv_ions.gro")} -o {index_file} << INPUT 
        q
        INPUT
        '''
    run_check_subprocess(cmd, key=wdir, log=None)

def get_mol_resid_pair(fname):
    with open(fname) as inp:
        data = inp.readlines()
    for molid_resid_pair in data:
        pair = [i.strip() for i in molid_resid_pair.split('\t')]
        if pair:
            molid, resid = pair
            yield molid, resid

def run_check_subprocess(cmd, key, log, env=None):
    try:
        subprocess.check_output(cmd, shell=True, env=env, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        logging.exception(f'{key}. {f"Check log {log}" if log else ""}\nError:{e}', stack_info=True)
        return False
    return True

def get_protein_resid_set(protein_fname):
    protein = mda.Universe(protein_fname)
    protein_resid_set = set(protein.residues.resnames.tolist())
    return protein_resid_set