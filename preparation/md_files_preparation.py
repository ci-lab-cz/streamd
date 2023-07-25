import logging
import os
import shutil
import subprocess
from glob import glob

from utils.utils import get_index, make_group_ndx, get_mol_resid_pair


def add_ligands_to_topol(all_itp_list, all_posres_list, all_resids, topol):
    def get_n_line_of_protein_in_mol_section(topol_file):
        with open(topol_file) as input:
            data = input.read()
        molecules_section = data[data.find('[ molecules ]'):]
        last_protein_n_line = None
        for n, line in enumerate(molecules_section.split('\n'), 1):
            if line.startswith('Protein'):
                last_protein_n_line = n
        return last_protein_n_line

    itp_include_list, posres_include_list, resid_include_list = [], [], []
    for itp, posres, resid in zip(all_itp_list, all_posres_list, all_resids):
        itp_include_list.append(f'; Include {resid} topology\n'
                                f'#include "{itp}"\n')
        posres_include_list.append(f'; {resid} position restraints\n#ifdef POSRES_{resid}\n'
                                   f'#include "{posres}"\n#endif\n')
        resid_include_list.append(f'{resid}             1')

    edit_topology_file(topol, pattern="; Include forcefield parameters",
                       add='\n'.join(itp_include_list), how='after', n=3)
    # reverse order since add before pattern
    edit_topology_file(topol, pattern="; Include topology for ions",
                       add='\n'.join(posres_include_list[::-1]), how='before')
    # ligand molecule should be after protein (or all protein chains listed)
    n_line_after_molecules_section = get_n_line_of_protein_in_mol_section(topol_file=topol)
    edit_topology_file(topol, pattern='[ molecules ]', add='\n'.join(resid_include_list),
                       how='after', n=n_line_after_molecules_section)


def edit_mdp(md_file, pattern, replace):
    new_mdp = []
    with open(md_file) as inp:
        for line in inp.readlines():
            if line.startswith(pattern):
                new_mdp.append(f'{replace.strip()}\n')
            else:
                new_mdp.append(line)
    with open(md_file, 'w') as out:
        out.write(''.join(new_mdp))


def edit_topology_file(topol_file, pattern, add, how='before', n=0):
    with open(topol_file) as input:
        data = input.read()

    if n == 0:
        data = data.replace(pattern, f'{add}\n{pattern}' if how == 'before' else f'{pattern}\n{add}')
    else:
        data = data.split('\n')
        ind = data.index(pattern)
        data.insert(ind + n, add)
        data = '\n'.join(data)

    with open(topol_file, 'w') as output:
        output.write(data)


def prep_md_files(wdir_var_ligand, protein_name, wdir_system_ligand_list, wdir_protein, wdir_md, clean_previous=False):
    '''

    :param wdir_var_ligand:
    :param protein_name:
    :param wdir_system_ligand_list:
    :param wdir_protein:
    :param wdir_md:
    :param clean_previous:
    :return: wdir to run md, dict with the list of md files which will be add to topol.top or to complex.gro
    Important: topol.top requires the same order of all additional files, so the same order should be preserved in all lists of the output dict
    '''

    def copy_md_files_to_wdir(files, wdir_copy_to):
        for file in files:
            shutil.copy(file, os.path.join(wdir_copy_to, os.path.basename(file)))

    if wdir_var_ligand:
        var_lig_molid, var_lig_resid = next(get_mol_resid_pair(os.path.join(wdir_var_ligand, 'resid.txt')))
        wdir_md_cur = os.path.join(wdir_md, f'{protein_name}_{var_lig_molid}')
    else:
        var_lig_molid, var_lig_resid = None, None
        wdir_md_cur = os.path.join(wdir_md, protein_name)

    if clean_previous and os.path.isdir(wdir_md_cur):
        shutil.rmtree(wdir_md_cur)

    os.makedirs(wdir_md_cur, exist_ok=True)

    # topol for protein for all chains
    copy_md_files_to_wdir(glob(os.path.join(wdir_protein, '*.itp')), wdir_copy_to=wdir_md_cur)
    copy_md_files_to_wdir([os.path.join(wdir_protein, "topol.top")], wdir_copy_to=wdir_md_cur)

    # prep ligands and cofactor
    # copy ligand.itp to md_wdir_cur. Will be edited
    md_files_dict = {'itp': [], 'gro': [], 'posres': [], 'molid': [], 'resid': []}
    if wdir_var_ligand:
        md_files_dict['itp'].append(os.path.join(wdir_md_cur, f'{var_lig_molid}.itp'))
        md_files_dict['gro'].append(os.path.join(wdir_var_ligand, f'{var_lig_molid}.gro'))
        md_files_dict['posres'].append(os.path.join(wdir_var_ligand, f'posre_{var_lig_molid}.itp'))
        md_files_dict['molid'].append(var_lig_molid)
        md_files_dict['resid'].append(var_lig_resid)

        copy_md_files_to_wdir([os.path.join(wdir_var_ligand, f'{var_lig_molid}.itp')], wdir_copy_to=wdir_md_cur)

    # copy system_lig itp to ligand_md_wdir
    if wdir_system_ligand_list:
        for wdir_system_ligand in wdir_system_ligand_list:
            system_lig_molid, system_lig_resid = next(get_mol_resid_pair(os.path.join(wdir_system_ligand, 'resid.txt')))
            copy_md_files_to_wdir([os.path.join(wdir_system_ligand, f'{system_lig_molid}.itp')],
                                  wdir_copy_to=wdir_md_cur)
            md_files_dict['itp'].append(os.path.join(wdir_md_cur, f'{system_lig_molid}.itp'))
            md_files_dict['gro'].append(os.path.join(wdir_system_ligand, f'{system_lig_molid}.gro'))
            md_files_dict['posres'].append(os.path.join(wdir_system_ligand, f'posre_{system_lig_molid}.itp'))
            md_files_dict['molid'].append(system_lig_molid)
            md_files_dict['resid'].append(system_lig_resid)

    return wdir_md_cur, md_files_dict


def prepare_mdp_files(wdir_md_cur, all_resids, script_path, nvt_time_ps, npt_time_ps, mdtime_ns):
    if not os.path.isfile(os.path.join(wdir_md_cur, 'index.ndx')):
        try:
            subprocess.check_output(f'''
            cd {wdir_md_cur}
            gmx make_ndx -f solv_ions.gro << INPUT
            q
            INPUT
            ''', shell=True)

        except subprocess.CalledProcessError as e:
            logging.exception(f'{wdir_md_cur}\n{e}', stack_info=True)
            return None

    index_list = get_index(os.path.join(wdir_md_cur, 'index.ndx'))
    # make couple_index_group
    couple_group_ind = '|'.join([str(index_list.index(i)) for i in ['Protein'] + all_resids])
    couple_group = '_'.join(['Protein'] + all_resids)

    non_couple_group = f'!{couple_group}'

    for mdp_fname in ['nvt.mdp', 'npt.mdp', 'md.mdp']:
        mdp_file = os.path.join(script_path, mdp_fname)
        shutil.copy(mdp_file, wdir_md_cur)
        md_fname = os.path.basename(mdp_file)

        edit_mdp(md_file=os.path.join(wdir_md_cur, md_fname),
                 pattern='tc-grps',
                 replace=f'tc-grps                 = {couple_group} {non_couple_group}; two coupling groups')
        steps = 0
        if md_fname == 'nvt.mdp':
            steps = int(nvt_time_ps * 1000 / 2)
        if md_fname == 'npt.mdp':
            steps = int(npt_time_ps * 1000 / 2)
        if md_fname == 'md.mdp':
            # picoseconds=mdtime*1000; femtoseconds=picoseconds*1000; steps=femtoseconds/2
            steps = int(mdtime_ns * 1000 * 1000 / 2)

        edit_mdp(md_file=os.path.join(wdir_md_cur, md_fname),
                 pattern='nsteps',
                 replace=f'nsteps                  = {steps}        ;')

    if couple_group not in index_list:
        if not make_group_ndx(couple_group_ind, wdir_md_cur):
            return None
    if non_couple_group not in index_list:
        index_list = get_index(os.path.join(wdir_md_cur, 'index.ndx'))
        non_couple_group_ind = f'!{index_list.index(couple_group)}'
        if not make_group_ndx(non_couple_group_ind, wdir_md_cur):
            return None

    return wdir_md_cur
