import os
import re
import subprocess
import logging


def filepath_type(x):
    if x:
        return os.path.abspath(x)
    else:
        return x


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
