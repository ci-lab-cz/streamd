import logging
import os
import re
import subprocess


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
    groups = [i.strip() for i in re.findall('\[(.*)\]', data)]
    return groups


def make_group_ndx(query, wdir):
    try:
        subprocess.check_output(f'''
        cd {wdir}
        gmx make_ndx -f solv_ions.gro -n index.ndx << INPUT
        {query}
        q
        INPUT
        ''', shell=True)
    except subprocess.CalledProcessError as e:
        logging.error(f'{wdir}\t{e}\n')
        return False
    return True


def get_mol_resid_pair(fname):
    with open(fname) as inp:
        data = inp.readlines()
    for molid_resid_pair in data:
        pair = [i.strip() for i in molid_resid_pair.split('\t')]
        if pair:
            molid, resid = pair
            yield molid, resid

def run_check_subprocess(cmd, key):
    try:
        subprocess.check_output(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        logging.exception(f'{key}\nError:{e}', stack_info=True)
        return False
    return True
