import argparse
import os
import sys
import shutil
from glob import glob
from multiprocessing import cpu_count
import math
import time
import random
import subprocess
from rdkit import Chem
from dask.distributed import Client
import logging
from datetime import datetime
import re

class RawTextArgumentDefaultsHelpFormatter(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass


def filepath_type(x):
    if x:
        return os.path.abspath(x)
    else:
        return x

def init_dask_cluster(n_tasks_per_node, ncpu, hostfile=None):
    '''

    :param n_tasks_per_node: number of task on a single server
    :param ncpu: number of cpu on a single server
    :param hostfile:
    :return:
    '''
    if hostfile:
        with open(hostfile) as f:
            hosts = [line.strip() for line in f]
            n_servers = sum(1 if line.strip() else 0 for line in f)
    else:
        n_servers = 1

    n_workers = n_servers * n_tasks_per_node
    n_threads = math.ceil(ncpu / n_tasks_per_node)
    if hostfile is not None:
        cmd = f'dask ssh --hostfile {hostfile} --nworkers {n_workers} --nthreads {n_threads} &'
        subprocess.check_output(cmd, shell=True)
        time.sleep(10)
        dask_client = Client(hosts[0] + ':8786', connection_limit=2048)
    else:
        dask_client = Client(n_workers=n_workers, threads_per_worker=n_threads)   # to run dask on a single server

    dask_client.forward_logging(level=logging.INFO)
    return dask_client


def make_all_itp(fileitp_list, out_file):
    atom_type_list = []
    start_columns = None
    #'[ atomtypes ]\n; name    at.num    mass    charge ptype  sigma      epsilon\n'
    for f in fileitp_list:
        with open(f) as input:
            data = input.read()
        start = data.find('[ atomtypes ]')
        end = data.find('[ moleculetype ]') - 1
        atom_type_list.extend(data[start:end].split('\n')[2:])
        if start_columns is None:
            start_columns = data[start:end].split('\n')[:2]
        new_data = data[:start] + data[end + 1:]
        with open(f, 'w') as itp_ouput:
            itp_ouput.write(new_data)

    atom_type_uniq = [i for i in set(atom_type_list) if i]
    with open(out_file, 'w') as ouput:
        ouput.write('\n'.join(start_columns)+'\n')
        ouput.write('\n'.join(atom_type_uniq)+'\n')


def complex_preparation(protein_gro, ligand_gro_list, out_file):
    atoms_list = []
    with open(protein_gro) as input:
        prot_data = input.readlines()
        atoms_list.extend(prot_data[2:-1])

    for f in ligand_gro_list:
        with open(f) as input:
            data = input.readlines()
        atoms_list.extend(data[2:-1])

    n_atoms = len(atoms_list)
    with open(out_file, 'w') as output:
        output.write(prot_data[0])
        output.write(f'{n_atoms}\n')
        output.write(''.join(atoms_list))
        output.write(prot_data[-1])


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
            logging.error(f'{wdir_md_cur}\t{e}\n')
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
        # check if works test
        index_list = get_index(os.path.join(wdir_md_cur, 'index.ndx'))
        non_couple_group_ind = f'!{index_list.index(couple_group)}'
        if not make_group_ndx(non_couple_group_ind, wdir_md_cur):
            return None

    return wdir_md_cur

def run_complex_prep(var_lig_data, wdir_protein, system_lig_data,
                     protein_name, wdir_system_ligand, wdir_md,
                     script_path, project_dir, mdtime_ns,
                     npt_time_ps, nvt_time_ps, clean_previous=False):

    wdir_md_cur, all_itp_list, all_gro_list, all_posres_list, all_lig_molids, all_resids = \
        prep_md_files(var_lig_data=var_lig_data, protein_name=protein_name, system_lig_data=system_lig_data,
                                wdir_protein=wdir_protein, wdir_md=wdir_md, clean_previous=clean_previous)

    protein_gro = os.path.join(wdir_protein, f'{protein_name}.gro')

    if var_lig_data:
        wdir_var_ligand_cur, var_lig_molid, var_lig_resid = var_lig_data
    else:
        wdir_var_ligand_cur, var_lig_molid, var_lig_resid = None, None, None


    if all_itp_list:
        if not os.path.isfile(os.path.join(wdir_md_cur, "all.itp")):
            # make all itp and edit current itps
            make_all_itp(all_itp_list, out_file=os.path.join(wdir_md_cur, 'all.itp'))
        else:
            logging.warning(f'{wdir_md_cur}. Prepared itp files exist. Skip topol preparation step\n')

        add_ligands_to_topol(all_itp_list, all_posres_list, all_resids, topol=os.path.join(wdir_md_cur, "topol.top"))
        # copy molid-resid pairs for variable ligand and all system ligands
        edit_topology_file(topol_file=os.path.join(wdir_md_cur, "topol.top"), pattern="; Include forcefield parameters",
                           add=f'; Include all topology\n#include "{os.path.join(wdir_md_cur, "all.itp")}"\n',
                           how='after', n=3)

        with open(os.path.join(wdir_md_cur, 'all_ligand_resid.txt'), 'w') as out:
            if wdir_var_ligand_cur and os.path.isfile(os.path.join(wdir_var_ligand_cur, 'resid.txt')):
                with open(os.path.join(wdir_var_ligand_cur, 'resid.txt')) as inp:
                    out.write(inp.read())
            if wdir_system_ligand and os.path.isfile(os.path.join(wdir_system_ligand, 'all_resid.txt')):
                with open(os.path.join(wdir_system_ligand, 'all_resid.txt')) as inp:
                    out.write(inp.read())

    # complex
    if not os.path.isfile(os.path.join(wdir_md_cur, 'complex.gro')):
        complex_preparation(protein_gro=protein_gro,
                            ligand_gro_list=all_gro_list,
                            out_file=os.path.join(wdir_md_cur, 'complex.gro'))
    else:
        logging.warning(f'{wdir_md_cur}. Prepared complex file exists. Skip complex preparation step\n')

    for mdp_fname in ['ions.mdp','minim.mdp']:
        mdp_file = os.path.join(script_path, mdp_fname)
        shutil.copy(mdp_file, wdir_md_cur)

    if not os.path.isfile(os.path.join(wdir_md_cur, 'solv_ions.gro')):
        try:
            subprocess.check_output(f'wdir={wdir_md_cur} bash {os.path.join(project_dir, "solv_ions.sh")}', shell=True)
        except subprocess.CalledProcessError as e:
            logging.error(f'{wdir_md_cur}\t{e}\n')
            return None
    else:
        logging.warning(f'{wdir_md_cur}. Prepared solv_ions.gro file exists. Skip solvation and ion preparation step\n')

    if not prepare_mdp_files(wdir_md_cur=wdir_md_cur, all_resids=all_resids,
                      script_path=script_path, nvt_time_ps=nvt_time_ps,
                      npt_time_ps=npt_time_ps, mdtime_ns=mdtime_ns):
        return None

    return wdir_md_cur


def edit_topology_file(topol_file, pattern, add, how='before', n=0):
    with open(topol_file) as input:
        data = input.read()

    if n == 0:
        data = data.replace(pattern, f'{add}\n{pattern}' if how == 'before' else f'{pattern}\n{add}')
    else:
        data = data.split('\n')
        ind = data.index(pattern)
        data.insert(ind+n, add)
        data = '\n'.join(data)

    with open(topol_file, 'w') as output:
        output.write(data)

def get_n_line_of_protein_in_mol_section(topol_file):
    with open(topol_file) as input:
        data = input.read()
    molecules_section = data[data.find('[ molecules ]'):]
    last_protein_n_line = None
    for n, line in enumerate(molecules_section.split('\n'), 1):
        if line.startswith('Protein'):
            last_protein_n_line = n
    return last_protein_n_line


def add_ligands_to_topol(all_itp_list, all_posres_list, all_resids, topol):
    itp_include_list, posres_include_list, resid_include_list = [], [], []
    for itp, posres, resid in zip(all_itp_list, all_posres_list, all_resids):
        itp_include_list.append(f'; Include {resid} topology\n'
                                f'#include "{itp}"\n')
        posres_include_list.append(f'; {resid} position restraints\n#ifdef POSRES_{resid}\n'
                                  f'#include "{posres}"\n#endif\n')
        resid_include_list.append(f'{resid}             1')

    edit_topology_file(topol, pattern="; Include forcefield parameters",
                add='\n'.join(itp_include_list), how='after', n=3)
    #reverse order since add before pattern
    edit_topology_file(topol, pattern="; Include topology for ions",
                add='\n'.join(posres_include_list[::-1]), how='before')
    # ligand molecule should be after protein (or all protein chains listed)
    n_line_after_molecules_section = get_n_line_of_protein_in_mol_section(topol_file=topol)
    edit_topology_file(topol, pattern='[ molecules ]', add='\n'.join(resid_include_list),
                       how='after', n=n_line_after_molecules_section)


def prep_ligand(mol, script_path, project_dir, wdir_ligand, addH=True):
    molid = mol.GetProp('_Name')
    resid = mol.GetProp('resid')

    wdir_ligand_cur = os.path.join(wdir_ligand, molid)
    os.makedirs(wdir_ligand_cur, exist_ok=True)

    mol_file = os.path.join(wdir_ligand_cur, f'{molid}.mol')

    if os.path.isfile(os.path.join(wdir_ligand_cur, f'{molid}.itp')) and os.path.isfile(os.path.join(wdir_ligand_cur, f'posre_{molid}.itp'))\
            and os.path.isfile(os.path.join(wdir_ligand_cur, 'resid.txt')):
        logging.warning(f'{molid}.itp and posre_{molid}.itp files already exist. Mol preparation step will be skipped for such molecule\n')
        with open(os.path.join(wdir_ligand_cur, 'resid.txt')) as inp:
            molid, resid = inp.read().strip().split('\t')

        return wdir_ligand_cur, molid, resid

    if addH:
        mol = Chem.AddHs(mol, addCoords=True)

    Chem.MolToMolFile(mol, mol_file)

    try:
        subprocess.check_output(f'script_path={script_path} lfile={mol_file} input_dirname={wdir_ligand_cur} name={resid} bash {os.path.join(project_dir, "lig_prep.sh")}',
                                shell=True)
    except subprocess.CalledProcessError as e:
        logging.error(f'{molid}\t{e}\n')
        return None

    # create log for molid resid corresponding
    with open(os.path.join(wdir_ligand_cur, 'resid.txt'), 'w') as out:
        out.write(f'{molid}\t{resid}\n')
    with open(os.path.join(wdir_ligand, 'all_resid.txt'), 'a') as out:
        out.write(f'{molid}\t{resid}\n')

    return wdir_ligand_cur, molid, resid

def prep_md_files(var_lig_data, protein_name, system_lig_data, wdir_protein, wdir_md, clean_previous=False):
    def copy_md_files_to_wdir(files, wdir_copy_to):
        for file in files:
            shutil.copy(file, os.path.join(wdir_copy_to, os.path.basename(file)))

    if var_lig_data:
        wdir_var_ligand_cur, var_lig_molid, var_lig_resid = var_lig_data
        wdir_md_cur = os.path.join(wdir_md, f'{protein_name}_{var_lig_molid}')
    else:
        wdir_md_cur = os.path.join(wdir_md, protein_name)

    if clean_previous and os.path.isdir(wdir_md_cur):
        shutil.rmtree(wdir_md_cur)

    os.makedirs(wdir_md_cur, exist_ok=True)

    # topol for protein for all chains
    copy_md_files_to_wdir(glob(os.path.join(wdir_protein,'*.itp')), wdir_copy_to=wdir_md_cur)
    copy_md_files_to_wdir([os.path.join(wdir_protein, "topol.top")], wdir_copy_to=wdir_md_cur)

    # prep ligands and cofactor
    # copy ligand.itp to md_wdir_cur. Will be edited
    if var_lig_data:
        wdir_var_ligand_cur, var_lig_molid, var_lig_resid = var_lig_data
        all_itp_list, all_gro_list, all_posres_list, all_lig_molids, all_resids = [os.path.join(wdir_md_cur, f'{var_lig_molid}.itp')],\
                                                     [os.path.join(wdir_var_ligand_cur, f'{var_lig_molid}.gro')],\
                                                     [os.path.join(wdir_var_ligand_cur, f'posre_{var_lig_molid}.itp')],\
                                                     [var_lig_molid], \
                                                     [var_lig_resid]
        copy_md_files_to_wdir([os.path.join(wdir_var_ligand_cur, f'{var_lig_molid}.itp')], wdir_copy_to=wdir_md_cur)
    else:
        all_itp_list, all_gro_list, all_posres_list, all_lig_molids, all_resids = [], [], [], [], []

    # copy system_lig itp to ligand_md_wdir
    for wdir_system_ligand_cur, system_lig_molid, system_lig_resid in system_lig_data:
        copy_md_files_to_wdir([os.path.join(wdir_system_ligand_cur, f'{system_lig_molid}.itp')], wdir_copy_to=wdir_md_cur)
        all_itp_list.append(os.path.join(wdir_md_cur, f'{system_lig_molid}.itp'))
        all_gro_list.append(os.path.join(wdir_system_ligand_cur, f'{system_lig_molid}.gro'))
        all_posres_list.append(os.path.join(wdir_system_ligand_cur, f'posre_{system_lig_molid}.itp'))
        all_lig_molids.append(system_lig_molid)
        all_resids.append(system_lig_resid)

    return wdir_md_cur, all_itp_list, all_gro_list, all_posres_list, all_lig_molids, all_resids


def supply_mols(fname, set_resid=None):
    def create_random_resid():
        # gro
        ascii_uppercase_digits = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789'
        return ''.join(random.choices(ascii_uppercase_digits, k=3))

    def add_resid(mol, n, input_fname, set_resid, used_resids):
        if set_resid is None:
            cur_resid = create_random_resid()
            while cur_resid == 'UNL' or cur_resid in used_resids:
                cur_resid = create_random_resid()
            mol.SetProp('resid', cur_resid)
        else:
            mol.SetProp('resid', set_resid)

        if not mol.HasProp('_Name'):
            mol.SetProp('_Name', f'{input_fname}_ID{n}')
        return mol

    used_resids = []

    if fname.endswith('.sdf'):
        for n, mol in enumerate(Chem.SDMolSupplier(fname, removeHs=False)):
            if mol:
                mol = add_resid(mol, n, input_fname=os.path.basename(fname).strip('.sdf'),
                                set_resid=set_resid,
                                used_resids=used_resids)
                used_resids.append(mol.GetProp('resid'))
                yield mol

    if fname.endswith('.mol'):
        mol = Chem.MolFromMolFile(fname, removeHs=False)
        if mol:
            mol = add_resid(mol, n=1, input_fname=os.path.basename(fname).strip('.mol'),
                            set_resid=set_resid,
                            used_resids=used_resids)
            used_resids.append(mol.GetProp('resid'))
            yield mol


def calc_dask(func, main_arg, dask_client, dask_report_fname=None, **kwargs):
    main_arg = iter(main_arg)
    Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)
    if dask_client is not None:
        from dask.distributed import as_completed, performance_report
        # https://stackoverflow.com/a/12168252/895544 - optional context manager
        from contextlib import contextmanager
        none_context = contextmanager(lambda: iter([None]))()
        with (performance_report(filename=dask_report_fname) if dask_report_fname is not None else none_context):
            nworkers = len(dask_client.scheduler_info()['workers'])
            futures = []
            for i, arg in enumerate(main_arg, 1):
                futures.append(dask_client.submit(func, arg, **kwargs))
                if i == nworkers * 2:  # you may submit more tasks then workers (this is generally not necessary if you do not use priority for individual tasks)
                    break
            seq = as_completed(futures, with_results=True)
            for i, (future, results) in enumerate(seq, 1):
                yield results
                del future
                try:
                    arg = next(main_arg)
                    new_future = dask_client.submit(func, arg, **kwargs)
                    seq.add(new_future)
                except StopIteration:
                    continue

def run_equilibration(wdir, project_dir):
    if os.path.isfile(os.path.join(wdir, 'npt.gro')) and os.path.isfile(os.path.join(wdir, 'npt.cpt')):
        logging.warning(f'{wdir}. Checkpoint files after Equilibration exist. '
                         f'Equilibration step will be skipped ')
        return wdir

    try:
        subprocess.check_output(f'wdir={wdir} bash {os.path.join(project_dir, "equlibration.sh")}', shell=True)
    except subprocess.CalledProcessError as e:
        logging.error(f'{wdir}\t{e}\n')
        return None
    return wdir

def run_simulation(wdir, project_dir):
    if os.path.isfile(os.path.join(wdir, 'md_out.tpr')) and os.path.isfile(os.path.join(wdir, 'md_out.cpt'))\
            and os.path.isfile(os.path.join(wdir, 'md_out.xtc')) and os.path.isfile(os.path.join(wdir, 'md_out.log')):
        logging.warning(f'{wdir}. md_out.xtc and md_out.tpr exist. '
                        f'MD simulation step will be skipped. '
                        f'You can rerun the script and use --wdir_to_continue {wdir} --md_time time_in_ns to extend current trajectory.')
        return wdir
    try:
        subprocess.check_output(f'wdir={wdir} bash {os.path.join(project_dir, "md.sh")}', shell=True)
    except subprocess.CalledProcessError as e:
        logging.error(f'{wdir}\t{e}\n')
        return None
    return wdir

def md_lig_rmsd_analysis(molid, resid, wdir, tu):
    index_list = get_index(os.path.join(wdir, 'index.ndx'))
    if f'{resid}_&_!H*' not in index_list:
        if not make_group_ndx(query=f'{index_list.index(resid)} & ! a H*', wdir=wdir):
            return None
        index_list = get_index(os.path.join(wdir, 'index.ndx'))
    index_ligand_noH = index_list.index(f'{resid}_&_!H*')

    try:
        subprocess.check_output(f'''
        cd {wdir}
        gmx rms -s md_out.tpr -f md_fit.xtc -o rmsd_{molid}.xvg -n index.ndx  -tu {tu} <<< "Backbone  {index_ligand_noH}"''', shell=True)
    except subprocess.CalledProcessError as e:
        logging.error(f'{wdir}\t{e}\n')


def run_md_analysis(wdir, deffnm, mdtime_ns, project_dir):
    index_list = get_index(os.path.join(wdir, 'index.ndx'))

    if os.path.isfile(os.path.join(wdir, 'all_ligand_resid.txt')):
        with open(os.path.join(wdir, 'all_ligand_resid.txt')) as input:
            ligands_data = input.readlines()
    else:
        ligands_data = []

    if ligands_data:
        resid = 'UNL'
        if f'Protein_{resid}' not in index_list:
            if not make_group_ndx(query=f'"Protein"|{index_list.index(resid)}', wdir=wdir):
                return None
            index_list = get_index(os.path.join(wdir, 'index.ndx'))

        index_group = index_list.index(f'Protein_{resid}')
    else:
        index_group = index_list.index('Protein')

    tu = 'ps' if mdtime_ns <= 10 else 'ns'
    dtstep = 50 if mdtime_ns <= 10 else 100

    try:
        subprocess.check_output(f'wdir={wdir} index_group={index_group} tu={tu} dtstep={dtstep} deffnm={deffnm} bash {os.path.join(project_dir, "md_analysis.sh")}', shell=True)
    except subprocess.CalledProcessError as e:
        logging.error(f'{wdir}\t{e}\n')
        return None

    # molid resid pairs
    # calc rmsd
    for cur_ligand in ligands_data:
        cur_ligand_molid, cur_ligand_resid = cur_ligand.split('\t')
        cur_ligand_molid, cur_ligand_resid = cur_ligand_molid.strip(), cur_ligand_resid.strip()
        md_lig_rmsd_analysis(molid=cur_ligand_molid, resid=cur_ligand_resid, wdir=wdir, tu=tu)

    return wdir

def get_prev_last_step(md_log):
    with open(md_log) as inp:
        data = inp.read()
    last_step = re.findall('Writing checkpoint, step ([0-9]*) [A-Za-z0-9: ]*\n', data)
    if last_step:
        return int(last_step[0])


def continue_md_from_dir(wdir_to_continue, deffnm_prev, deffnm_next, mdtime_ns, project_dir):
    tpr = os.path.join(wdir_to_continue, f'{deffnm_prev}.tpr')
    cpk = os.path.join(wdir_to_continue, f'{deffnm_prev}.cpk')
    xtc = os.path.join(wdir_to_continue, f'{deffnm_prev}.xtc')
    md_log = os.path.join(wdir_to_continue, f'{deffnm_prev}.log')
    return continue_md(tpr=tpr, cpk=cpk, xtc=xtc, md_log=md_log, wdir=wdir_to_continue,
                       mdtime_ns=mdtime_ns, deffnm_prev=deffnm_prev, deffnm_next=deffnm_next, project_dir=project_dir)


def continue_md(tpr, cpk, xtc, md_log, wdir, mdtime_ns, deffnm_prev, deffnm_next, project_dir):
    last_step = get_prev_last_step(md_log)
    new_steps = int(mdtime_ns * 1000 * 1000 / 2) - last_step
    if new_steps <= 0:
        logging.error('Fail to continue MD simulation. New time equals or less then previously calculated simulation')
        return None
    try:
        subprocess.check_output(f'wdir={wdir} tpr={tpr} cpk={cpk} new_mdsteps={new_steps} '
                                f'deffnm_prev={deffnm_prev} deffnm_next={deffnm_next} xtc={xtc}  bash {os.path.join(project_dir, "continue_md.sh")}', shell=True)
    except subprocess.CalledProcessError as e:
        logging.error(f'{wdir}\t{e}\n')
        return None

    return wdir



def main(protein, wdir, lfile=None, system_lfile=None,
         forcefield_num=6, addH=True,
         gromacs_version="GROMACS/2021.4-foss-2020b-PLUMED-2.7.3",
         mdtime_ns=1, npt_time_ps=100, nvt_time_ps=100,
         topol=None, posre_protein=None,
         wdir_to_continue_list=None, deffnm_prev='md_out',
         hostfile=None, ncpu=1):
    '''

    :param protein:
    :param lfile:
    :param mdtime_ns:
    :param system_lfile:
    :param wdir:
    :param md_param:
    :param gromacs_version:
    :param hostfile:
    :param ncpu:
    :param topol:
    :param posre_protein:
    :param forcefield_num:
    :return:
    '''
    # TODO: add docstring with description of args


    try:
        subprocess.check_output(f'module load {gromacs_version}', shell=True)
    except subprocess.CalledProcessError as e:
        logging.error(e)
        return None

    project_dir = os.path.dirname(os.path.abspath(__file__))
    script_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'scripts')

    if (wdir_to_continue_list is None) and (protein is not None):
        # create dirs
        wdir_protein = os.path.join(wdir, 'md_preparation', 'protein')
        wdir_ligand = os.path.join(wdir, 'md_preparation', 'var_lig')
        wdir_system_ligand = os.path.join(wdir, 'md_preparation', 'system_lig')

        wdir_md = os.path.join(wdir, 'md_preparation', 'md_files')

        os.makedirs(wdir_md, exist_ok=True)
        os.makedirs(wdir_protein, exist_ok=True)
        os.makedirs(wdir_ligand, exist_ok=True)
        os.makedirs(wdir_system_ligand, exist_ok=True)

        if not os.path.isfile(protein):
            raise FileExistsError(f'{protein} does not exist')

        pname, p_ext = os.path.splitext(os.path.basename(protein))
        if not os.path.isfile(f'{os.path.join(wdir_protein, pname)}.gro'):
            if p_ext != '.gro' or topol is None or posre_protein is None:
                try:
                    subprocess.check_output(f'gmx pdb2gmx -f {protein} -o {os.path.join(wdir_protein, pname)}.gro -water tip3p -ignh '
                          f'-i {os.path.join(wdir_protein, "posre.itp")} '
                          f'-p {os.path.join(wdir_protein, "topol.top")}'
                          f'<<< {forcefield_num}', shell=True)
                except subprocess.CalledProcessError as e:
                    logging.error(e)
                    return None
            else:
                if not os.path.isfile(os.path.join(wdir_protein, protein)):
                    shutil.copy(protein, os.path.join(wdir_protein, protein))
                if not os.path.isfile(os.path.join(wdir_protein, 'topol.top')):
                    shutil.copy(topol, os.path.join(wdir_protein, 'topol.top'))
                if not os.path.isfile(os.path.join(wdir_protein, 'posre.itp')):
                    shutil.copy(posre_protein, os.path.join(wdir_protein, 'posre.itp'))
        else:
            logging.warning(f'{os.path.join(wdir_protein, pname)}.gro protein file exists. '
                            f'Protein preparation step will be skipped.')


        # Part 1. Preparation. Run on each cpu
        dask_client = init_dask_cluster(hostfile=hostfile, n_tasks_per_node=ncpu, ncpu=ncpu)
        try:
            system_lig_data = [] # wdir_ligand_cur, molid, resid
            if system_lfile is not None:
                if not os.path.isfile(system_lfile):
                    raise FileExistsError(f'{system_lfile} does not exist')

                mols = supply_mols(system_lfile, set_resid=None)
                for res in calc_dask(prep_ligand, mols, dask_client,
                                                  script_path=script_path, project_dir=project_dir,
                                                  wdir_ligand=wdir_system_ligand,
                                                  addH=addH):
                    if not res:
                        logging.error(f'Error with system ligand (cofactor) preparation. The calculation will be interrupted\n')
                        return
                    system_lig_data.append(res)

            var_lig_data = [] # wdir_ligand_cur, molid, resid
            if lfile is not None:
                if not os.path.isfile(lfile):
                    raise FileExistsError(f'{lfile} does not exist')

                mols = supply_mols(lfile, set_resid='UNL')
                for res in calc_dask(prep_ligand, mols, dask_client,
                                      script_path=script_path, project_dir=project_dir,
                                      wdir_ligand=wdir_ligand, addH=addH):
                    if res:
                        var_lig_data.append(res)

            # make all.itp and create complex
            var_complex_prepared_dirs = []
            if var_lig_data:
                for res in calc_dask(run_complex_prep, var_lig_data, dask_client,
                                     system_lig_data=system_lig_data,
                                     protein_name=pname,
                                     wdir_protein=wdir_protein,
                                     wdir_system_ligand=wdir_system_ligand,
                                     clean_previous=clean_previous, wdir_md=wdir_md,
                                     script_path=script_path, project_dir=project_dir, mdtime_ns=mdtime_ns,
                                     npt_time_ps=npt_time_ps, nvt_time_ps=nvt_time_ps):
                    if res:
                        var_complex_prepared_dirs.append(res)
            else:
                res = run_complex_prep(var_lig_data=[], wdir_protein=wdir_protein,
                                     system_lig_data=system_lig_data,
                                     protein_name=pname,
                                     wdir_system_ligand=wdir_system_ligand,
                                     clean_previous=clean_previous, wdir_md=wdir_md,
                                     script_path=script_path, project_dir=project_dir, mdtime_ns=mdtime_ns,
                                     npt_time_ps=npt_time_ps, nvt_time_ps=nvt_time_ps)
                if res:
                    var_complex_prepared_dirs.append(res)

        finally:
            dask_client.shutdown()

        # Part 2. Equilibration and MD simulation. Run on all cpu
        dask_client = init_dask_cluster(hostfile=hostfile, n_tasks_per_node=1, ncpu=ncpu)
        try:
            var_eq_dirs = []
            #os.path.dirname(var_lig)
            for res in calc_dask(run_equilibration, var_complex_prepared_dirs, dask_client, project_dir=project_dir):
                if res:
                    var_eq_dirs.append(res)

            var_md_dirs = []
            # os.path.dirname(var_lig)
            for res in calc_dask(run_simulation, var_eq_dirs, dask_client, project_dir=project_dir):
                if res:
                    var_md_dirs.append(res)

        finally:
            dask_client.shutdown()

        deffnm = 'md_out'

    else: # continue prev md]
        dask_client = init_dask_cluster(hostfile=hostfile, n_tasks_per_node=1, ncpu=ncpu)
        try:
            var_md_dirs = []
            deffnm = f'{deffnm_prev}_{mdtime_ns}'
            # os.path.dirname(var_lig)
            for res in calc_dask(continue_md_from_dir, wdir_to_continue_list, dask_client,
                                 deffnm_prev=deffnm_prev, deffnm_next=deffnm, mdtime_ns=mdtime_ns, project_dir=project_dir):
                if res:
                    var_md_dirs.append(res)

        finally:
            dask_client.shutdown()

    logging.info(f'Continue of simulation of {var_md_dirs} were successfully finished\n')

    # Part 3. MD Analysis. Run on each cpu
    dask_client = init_dask_cluster(hostfile=hostfile, n_tasks_per_node=ncpu, ncpu=ncpu)
    try:
        var_md_analysis_dirs = []
        # os.path.dirname(var_lig)
        for res in calc_dask(run_md_analysis, var_md_dirs,
                             dask_client, deffnm=deffnm, mdtime_ns=mdtime_ns, project_dir=project_dir):
            if res:
                var_md_analysis_dirs.append(res)
    finally:
        dask_client.shutdown()

    logging.info(f'\nAnalysis of md simulation of {var_md_analysis_dirs} were successfully finished\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=''' ''')
    parser.add_argument('-p', '--protein', metavar='FILENAME', required=False, type=filepath_type,
                        help='input file with compound. Supported formats: *.pdb or gro')
    parser.add_argument('-d', '--wdir', metavar='WDIR', default=None, type=filepath_type,
                        help='')
    parser.add_argument('-l', '--ligand', metavar='FILENAME', required=False, type=filepath_type,
                        help='input file with compound. Supported formats: *.mol or sdf or gro')
    parser.add_argument('--cofactor', metavar='FILENAME', default=None, type=filepath_type,
                        help='input file with compound. Supported formats: *.mol or sdf or gro')
    parser.add_argument('--not_add_H', action='store_true', default=False,
                        help='disable to add hydrogens to molecules before simulation.')
    parser.add_argument('--clean_previous_md', action='store_true', default=False,
                        help='Remove all previous prepared for the current system MD files.\n'
                             'Prepared ligand and protein files will be used if it exists.')
    parser.add_argument('--hostfile', metavar='FILENAME', required=False, type=str, default=None,
                        help='text file with addresses of nodes of dask SSH cluster. The most typical, it can be '
                             'passed as $PBS_NODEFILE variable from inside a PBS script. The first line in this file '
                             'will be the address of the scheduler running on the standard port 8786. If omitted, '
                             'calculations will run on a single machine as usual.')
    parser.add_argument('-c', '--ncpu', metavar='INTEGER', required=False, default=cpu_count(), type=int,
                        help='number of CPU per server. Use all cpus by default.')
    parser.add_argument('--topol', metavar='topol.top', required=False, default=None, type=filepath_type,
                        help='Required if gro file of the protein is provided')
    parser.add_argument('--posre', metavar='posre.itp', required=False, default=None, type=filepath_type,
                        help='Required if gro file of the protein is provided')
    ##
    parser.add_argument('--md_time', metavar='ns', required=False, default=1, type=float,
                        help='Time of MD simulation in ns')
    parser.add_argument('--npt_time', metavar='ps', required=False, default=100, type=int,
                        help='Set up NPT time equilibration in ns')
    parser.add_argument('--nvt_time', metavar='ps', required=False, default=100, type=int,
                        help='Set up NVT time equilibration in ns')
    #continue md
    parser.add_argument('--tpr', metavar='FILENAME', required=False, default=None, type=filepath_type,
                        help='TPR file from the previous MD simulation')
    parser.add_argument('--cpt', metavar='FILENAME', required=False, default=None, type=filepath_type,
                        help='CPT file from previous simulation')
    parser.add_argument('--md_log', metavar='FILENAME', required=False, default=None, type=filepath_type,
                        help='Log file from previous simulation. Example: md_out.log')
    parser.add_argument('--wdir_to_continue', metavar='DIRNAME', required=False, default=None, nargs='+', type=filepath_type,
                        help='wdir for previous simulation. Should consist of: tpr, cpt, xtc, md_log files')
    parser.add_argument('--deffnm', metavar='preffix for md files', required=False, default='md_out',
                        help='preffix for previous md files. Used only if wdir_to_continue is used')


    args = parser.parse_args()

    if args.wdir is None:
        wdir = os.getcwd()
    else:
        wdir = args.wdir

    log_file = os.path.join(wdir, f'log_{os.path.basename(str(args.protein))[:-4]}_{os.path.basename(str(args.ligand))[:-4]}{datetime.now().strftime("%d-%m-%Y-%H-%M-%S")}.log')

    logging.basicConfig(format='%(asctime)s - %(name)s - %(levelname)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S',
                        level=logging.INFO,
                        handlers=[logging.FileHandler(log_file),
                                  logging.StreamHandler()])

    logging.getLogger('distributed').setLevel('WARNING')
    logging.getLogger('distributed.worker').setLevel('WARNING')
    logging.getLogger('distributed.core').setLevel('WARNING')
    logging.getLogger('distributed.comm').setLevel('WARNING')
    logging.getLogger('bockeh').setLevel('WARNING')


    logging.info(args)
    try:
        main(protein=args.protein, lfile=args.ligand, addH=not args.not_add_H,
         clean_previous=args.clean_previous_md, system_lfile=args.cofactor,
         topol=args.topol, posre_protein=args.posre,  mdtime_ns=args.md_time, npt_time_ps=args.npt_time, nvt_time_ps=args.nvt_time,
         wdir_to_continue_list=args.wdir_to_continue, deffnm_prev=args.deffnm,
         hostfile=args.hostfile, ncpu=args.ncpu, wdir=wdir)
    finally:
        logging.shutdown()
