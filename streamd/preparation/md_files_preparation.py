import logging
import os
import shutil
from glob import glob


from streamd.utils.utils import get_index, make_group_ndx, get_mol_resid_pair, create_ndx

def get_couple_groups(mdp_file):
    with open(mdp_file) as inp:
        data = inp.readlines()

    tcl_group_string = None
    for line in data:
        if line.startswith('tc-grps'):
            tcl_group_string = line.strip()
            break

    if tcl_group_string:
        # tc-grps                 = Protein_UNL Water_and_ions; two coupling groups - more accurate' to catch Protein_UNL and Water_and_ions
        groups = tcl_group_string.split('=')[-1].split(';')[0].split()
        logging.warning(f'Create coupling groups: {groups} from user-defined mdp file {mdp_file}')

        return groups


def create_couple_group_in_index_file(couple_group, index_file, wdir, env, bash_log):
    '''

    :param couple_group: [Protein, Water, UNL]
    :param index_file: will be modified. To add new group if not exists
    :param env:
    :param bash_log:
    :return: None
    '''
    index_list = get_index(index_file, env=env)

    if couple_group not in index_list:
        couple_group_gmx =  couple_group.replace('_&_','&').replace('_','|')

        # obtain indexes of all groups
        couple_group_dict = {i: str(index_list.index(i)) for i in couple_group_gmx.replace('&','_').replace('|','_').replace('!','').split('_')}

        # replace group names by indexes to avoid duplicated names
        couple_group_ind = couple_group_gmx
        for group, index in couple_group_dict.items():
            couple_group_ind = couple_group_ind.replace(group, index)

        logging.info(f'Creating temperature coupling for {couple_group} group: {couple_group_ind}')
        if not make_group_ndx(couple_group_ind, wdir, bash_log=bash_log, env=env):
            logging.ERROR(f'Could not create a coupling group: {couple_group}. Calculations will be interrupted. '
                          f'Check {os.path.join(wdir,bash_log)}')
            return None

    return couple_group



def check_if_info_already_added_to_topol(topol, string):
    with open(topol) as inp:
        data = inp.read()
    return string in data # True or False

def add_ligands_to_topol(all_itp_list, all_posres_list, all_resids, topol):
    itp_posres_include_list, resid_include_list = [], []
    for itp, posres, resid in zip(all_itp_list, all_posres_list, all_resids):
        itp_posres_include_list.append(f'; Include {resid} topology\n'
                                f'#include "{itp}"\n\n'
                                f'; {resid} position restraints\n#ifdef POSRES\n'
                                   f'#include "{posres}"\n#endif\n')

        resid_include_list.append(f'{resid}             1')

    if not check_if_info_already_added_to_topol(topol, '\n'.join(itp_posres_include_list)):
        # edit_topology_file(topol, pattern="; Include forcefield parameters",
        #                add='\n'.join(itp_include_list), how='after', n=3)
        edit_topology_file(topol, pattern="; Include water topology",
                       add='\n'.join(itp_posres_include_list), how='before')

    if not check_if_info_already_added_to_topol(topol, '\n'.join(resid_include_list)):
        # ligand molecule should be after protein (or all protein chains listed)
        edit_topology_file(topol, pattern=None, add='\n'.join(resid_include_list), n=-1)


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

def edit_topology_file(topol_file, pattern, add, how='before', n=0, count_pattern=1):
    with open(topol_file) as input:
        data = input.read()

    if pattern is not None:
        if n!=0:
            data = data.split('\n')
            ind = data.index(pattern)
            if how == 'after':
                data.insert(ind + n, add)
            else:
                data.insert(ind - n, add)
            data = '\n'.join(data)

        elif count_pattern:
            arr = data.split(pattern)
            part1 = pattern.join(arr[:count_pattern])
            part2 = pattern.join(arr[count_pattern:])
            if how == 'before':
                data = part1 + add + '\n' + pattern + part2
            else:
                data = part1 + pattern + '\n' + add + part2

    if pattern is None:
        data = data.split('\n')
        data.insert(n, add)
        data = '\n'.join(data)

    with open(topol_file, 'w') as output:
        output.write(data)


def prep_md_files(wdir_var_ligand, protein_name, wdir_system_ligand_list, wdir_protein, wdir_md, clean_previous=False):
    '''

    :param wdir_var_ligand:
    :param protein_name:
    :param wdir_system_ligand_list:
    :param wdir_protein: None (if mcpbpy procedure is applied) or path where topol and gro files are stored
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
    if wdir_protein is not None:
        copy_md_files_to_wdir(glob(os.path.join(wdir_protein, '*.itp')), wdir_copy_to=wdir_md_cur)
        # don't rewrite existed topol.top
        if not os.path.isfile(os.path.join(wdir_md_cur, "topol.top")):
            copy_md_files_to_wdir([os.path.join(wdir_protein, "topol.top")], wdir_copy_to=wdir_md_cur)

    # prep ligands and cofactor
    # copy ligand.itp to md_wdir_cur. Will be edited
    md_files_dict = {'itp_orig': [], 'itp': [], 'gro': [], 'posres': [], 'molid': [], 'resid': []}
    if wdir_var_ligand:
        #copy posres
        copy_md_files_to_wdir([os.path.join(wdir_var_ligand, f'posre_{var_lig_molid}.itp')], wdir_copy_to=wdir_md_cur)
        md_files_dict['itp_orig'].append(os.path.join(wdir_var_ligand, f'{var_lig_molid}.itp'))
        md_files_dict['itp'].append(f'{var_lig_molid}.itp')
        md_files_dict['posres'].append(f'posre_{var_lig_molid}.itp')
        md_files_dict['gro'].append(os.path.join(wdir_var_ligand, f'{var_lig_molid}.gro'))
        md_files_dict['molid'].append(var_lig_molid)
        md_files_dict['resid'].append(var_lig_resid)

    if wdir_system_ligand_list:
        for wdir_system_ligand in wdir_system_ligand_list:
            system_lig_molid, system_lig_resid = next(get_mol_resid_pair(os.path.join(wdir_system_ligand, 'resid.txt')))
            copy_md_files_to_wdir([os.path.join(wdir_system_ligand, f'posre_{system_lig_molid}.itp')],
                                  wdir_copy_to=wdir_md_cur)
            md_files_dict['itp_orig'].append(os.path.join(wdir_system_ligand, f'{system_lig_molid}.itp'))
            md_files_dict['itp'].append(f'{system_lig_molid}.itp')
            md_files_dict['posres'].append(f'posre_{system_lig_molid}.itp')
            md_files_dict['gro'].append(os.path.join(wdir_system_ligand, f'{system_lig_molid}.gro'))
            md_files_dict['molid'].append(system_lig_molid)
            md_files_dict['resid'].append(system_lig_resid)

    return wdir_md_cur, md_files_dict


def prepare_mdp_files(wdir_md_cur, all_resids, nvt_time_ps, npt_time_ps,
                      mdtime_ns, user_mdp_files, bash_log, seed, explicit_args=(), env=None):
    '''

    :param wdir_md_cur:
    :param all_resids:
    :param nvt_time_ps:
    :param npt_time_ps:
    :param mdtime_ns:
    :param user_mdp_files:
    :param bash_log:
    :param seed:
    :param explicit_args: list of mdp file names to change if time and mdp files were explicitly provided by user
    :param env:
    :return:
    '''
    if not os.path.isfile(os.path.join(wdir_md_cur, 'index.ndx')):
        create_ndx(os.path.join(wdir_md_cur, 'index.ndx'), env=env)


    # Check if default coupling groups will be used
    default_couple_group, default_non_couple_group = None, None
    if not all([i in user_mdp_files for i in ['nvt.mdp', 'npt.mdp', 'md.mdp']]):
        default_couple_group = create_couple_group_in_index_file(couple_group='_'.join(['Protein'] + all_resids),
                                      index_file=os.path.join(wdir_md_cur, 'index.ndx'),
                                      wdir=wdir_md_cur, env=env, bash_log=bash_log)
        if not default_couple_group:
            return None

        default_non_couple_group = f'!{default_couple_group}'
        index_list = get_index(os.path.join(wdir_md_cur, 'index.ndx'), env=env)

        if default_non_couple_group not in index_list:
            non_couple_group_ind = f'!{index_list.index(default_couple_group)}'
            if not make_group_ndx(non_couple_group_ind, wdir_md_cur, bash_log=bash_log, env=env):
                return None

        logging.info(f"Default temperature coupling groups '{default_couple_group}' and '{default_non_couple_group}' "
                        f" will be used for { {'nvt.mdp', 'npt.mdp', 'md.mdp'} - set(user_mdp_files)} mdp files")

    # couple_group_ind = '|'.join([str(index_list.index(i)) for i in ['Protein'] + all_resids])


    for mdp_fname in ['nvt.mdp', 'npt.mdp', 'md.mdp']:
        mdp_file = os.path.join(wdir_md_cur, mdp_fname)

        if mdp_fname in user_mdp_files:
            # create user-defined coupling groups
            user_coupling_groups = get_couple_groups(mdp_file)
            if user_coupling_groups:
                for user_group in user_coupling_groups:
                    couple_group = create_couple_group_in_index_file(couple_group=user_group,
                                                                     index_file=os.path.join(wdir_md_cur, 'index.ndx'),
                                                                     wdir=wdir_md_cur, env=env, bash_log=bash_log)
                    if user_group != couple_group:
                        return None

        else:
            # use default parameters
            edit_mdp(md_file=mdp_file,
                     pattern='tc-grps',
                     replace=f'tc-grps                 = {default_couple_group} {default_non_couple_group}; two coupling groups')

        # check if seed were preset explicitly and need to be changed in user mdp files
        if mdp_fname == 'nvt.mdp' and ('seed' in explicit_args or mdp_fname not in user_mdp_files):
            logging.info(f'Set seed: {seed} in nvt.mdp')
            edit_mdp(md_file=mdp_file,
                     pattern='gen_seed',
                     replace=f'gen_seed                = {seed}        ;')

        # check if time were preset explicitly and need to be changed in user mdp files
        if mdp_fname in explicit_args or mdp_fname not in user_mdp_files:
            steps = 0
            if mdp_fname == 'nvt.mdp':
                steps = int(nvt_time_ps * 1000 / 2)
            if mdp_fname == 'npt.mdp':
                steps = int(npt_time_ps * 1000 / 2)
            if mdp_fname == 'md.mdp':
                # picoseconds=mdtime*1000; femtoseconds=picoseconds*1000; steps=femtoseconds/2
                steps = int(mdtime_ns * 1000 * 1000 / 2)

            logging.info(f'Set {steps/1000/1000*2} ns = {steps/1000 *2} ps = {steps} steps in {mdp_fname}')
            edit_mdp(md_file=mdp_file,
                     pattern='nsteps',
                     replace=f'nsteps                  = {steps}        ;')


    return wdir_md_cur
