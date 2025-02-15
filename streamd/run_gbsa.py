#!/usr/bin/env python3
# ==============================================================================
# author          : Aleksandra Ivanova, Olena Mokshyna, Pavel Polishchuk
# date            : 14-06-2024
# version         :
# python_version  :
# copyright       :
# license         : MIT
# ==============================================================================
import argparse
import logging
import math
import os
import pathlib
import re
import tempfile
import shutil
import subprocess
from datetime import datetime
from functools import partial
from multiprocessing import Pool

import pandas as pd

from streamd.utils.dask_init import init_dask_cluster, calc_dask
from streamd.utils.utils import (get_index, make_group_ndx, filepath_type, run_check_subprocess,
                                 get_number_of_frames)

logging.getLogger('distributed').setLevel('CRITICAL')
logging.getLogger('asyncssh').setLevel('CRITICAL')
logging.getLogger('distributed.worker').setLevel('CRITICAL')
logging.getLogger('distributed.core').setLevel('CRITICAL')
logging.getLogger('distributed.comm').setLevel('CRITICAL')
logging.getLogger('distributed.nanny').setLevel('CRITICAL')
logging.getLogger('bockeh').setLevel('CRITICAL')


def run_gbsa_task(wdir, tpr, xtc, topol, index, mmpbsa, np, ligand_resid, append_protein_selection,
                  unique_id, env, bash_log, clean_previous):

    def calc_gbsa(wdir, tpr, xtc, topol, index, mmpbsa, np, protein_index,
                  ligand_index, unique_id, env, bash_log):
        output = os.path.join(wdir, f"FINAL_RESULTS_MMPBSA_{unique_id}.dat")
        with tempfile.TemporaryDirectory(dir=wdir) as tmpdirname:
            logging.info(f'tmp intermediate dir: {tmpdirname}')
            cmd = f'cd {tmpdirname}; mpirun -np {np} gmx_MMPBSA MPI -O -i {mmpbsa} ' \
                  f' -cs {tpr} -ci {index} -cg {protein_index} {ligand_index} -ct {xtc} -cp {topol} -nogui ' \
                  f'-o {output} ' \
                  f'-eo {os.path.join(wdir, f"FINAL_RESULTS_MMPBSA_{unique_id}.csv")}' \
                  f' >> {os.path.join(wdir, bash_log)} 2>&1'

            if not run_check_subprocess(cmd, key=xtc, log=os.path.join(wdir, bash_log), env=env):
                run_check_subprocess(f'ls {tmpdirname}', key=tmpdirname, log=os.path.join(wdir, bash_log), env=env)
                return None

        return output

    if clean_previous:
        clean_temporary_gmxMMBPSA_files(wdir)

    if not os.path.isfile(tpr) or not os.path.isfile(xtc) or not os.path.isfile(topol) or not os.path.isfile(index):
        logging.warning(f'{wdir} cannot run gbsa. Check if there are missing files: {tpr} {xtc} {topol} {index}')
        return None

    index_list = get_index(index, env=env)
    if append_protein_selection is None:
        protein_index = index_list.index('Protein')
    else:
        add_group_ids = {}
        for i in append_protein_selection:
            if i in index_list:
                add_group_ids[str(index_list.index(i))] = i
            else:
                logging.warning(f'{wdir} Could not find resname {i}. It will not be used in gbsa calculation. Check your query carefully.')
        if add_group_ids:
            query = f"{index_list.index('Protein')}|{'|'.join(add_group_ids.keys())}"
            name_query = f"Protein_{'_'.join(add_group_ids.values())}"
            if name_query not in index_list:
                if not make_group_ndx(query, wdir, bash_log=bash_log, env=env):
                    return None
                index_list = get_index(index, env=env)

            protein_index = index_list.index(name_query)
            logging.info(f'{name_query} selection will be used as a protein system')
        else:
            protein_index = index_list.index('Protein')

    logging.info(f'{protein_index} number of index selection will be used as a protein system')

    ligand_index = index_list.index(ligand_resid)

    output = calc_gbsa(wdir=wdir, tpr=tpr, xtc=xtc, topol=topol,
                       index=index, mmpbsa=mmpbsa,
                       np=np, protein_index=protein_index,
                       ligand_index=ligand_index,
                       unique_id=unique_id,
                       env=env,
                       bash_log=bash_log)

    if os.path.isfile(os.path.join(wdir, 'gmx_MMPBSA.log')):
        shutil.copy(os.path.join(wdir, 'gmx_MMPBSA.log'), os.path.join(wdir, f'gmx_MMPBSA_{unique_id}.log'))

    return output


def run_gbsa_from_wdir(wdir, tpr, xtc, topol, index, mmpbsa, np, ligand_resid,
                       append_protein_selection, unique_id, env, bash_log, clean_previous):
    tpr = os.path.join(wdir, tpr)
    xtc = os.path.join(wdir, xtc)
    topol = os.path.join(wdir, topol)
    index = os.path.join(wdir, index)
    return run_gbsa_task(wdir=wdir, tpr=tpr, xtc=xtc,
                         topol=topol, index=index, mmpbsa=mmpbsa,
                         np=np, ligand_resid=ligand_resid, append_protein_selection=append_protein_selection,
                         unique_id=unique_id,
                         env=env, bash_log=bash_log, clean_previous=clean_previous)


def clean_temporary_gmxMMBPSA_files(wdir, prefix="_GMXMMPBSA_"):
    # remove intermediate files
    try:
        subprocess.check_output(f'cd {wdir}; gmx_MMPBSA --clean -prefix "{prefix}" > /dev/null 2>&1 ', shell=True)
    except subprocess.CalledProcessError as e:
        logging.error(f'{wdir}\t{e}\n')
        return None


def parse_gmxMMPBSA_output(fname):
    def get_IE_values(IE_parsed_out):
        IE_res = {}
        IE_columns = [i.strip() for i in IE_parsed_out[0][0].split('  ') if i]
        IE_values = [i.strip() for i in IE_parsed_out[0][1].split('  ') if i]
        for n, i in enumerate(IE_columns):
            IE_res[f'IE_{i}'] = IE_values[n]
        return IE_res

    def get_delta_total_values(delta_total_columns_re, delta_total_values_re):
        delta_total_res = {}
        delta_total_columns = [i.strip() for i in delta_total_columns_re[0].split('  ') if i]
        delta_total_values = [i.strip() for i in delta_total_values_re[0].split('  ') if i]
    #['Energy Component' - skip, 'Average', 'SD(Prop.)', 'SD', 'SEM(Prop.)', 'SEM']

        for n, i in enumerate(delta_total_columns[1:]):
            delta_total_res[f'ΔTOTAL_{i}'] = delta_total_values[n]
        return delta_total_res

    def get_Gbinding_values(Gbind_parsed_out):
        Gbinding_res = {}
        G_values = [i.strip() for i in Gbind_parsed_out[0].split(' ') if i and i != '+/-']
        Gbinding_res['ΔGbinding'] = G_values[0]
        Gbinding_res['ΔGbinding+/-'] = G_values[1]
        return Gbinding_res

    with open(fname) as inp:
        data = inp.read()
    IE_GB = re.findall(
        r'Energy Method[ ]*?Entropy[ ]*?(σ\(Int. Energy\)[ ]*?Average[ ]*?SD[ ]*?SEM)\n[-]*?\nGB[ ]*?IE[ ]*?([0-9-\. ]*)\n',
        data)
    IE_PB = re.findall(
        r'Energy Method[ ]*?Entropy[ ]*?(σ\(Int. Energy\)[ ]*?Average[ ]*?SD[ ]*?SEM)\n[-]*?\nPB[ ]*?IE[ ]*?([0-9-\. ]*)\n',
        data)

    delta_total_columns_GB = re.findall(
    r'GENERALIZED BORN:[A-Z0-9\w\W\n]+?Delta \(Complex - Receptor - Ligand\):\n([A-Za-z \(\).]+)\n',
    data)
    delta_total_columns_PB = re.findall(
    r'POISSON BOLTZMANN:[A-Z0-9\w\W\n]+?Delta \(Complex - Receptor - Ligand\):\n([A-Za-z \(\).]+)\n',
    data)

    # delta_total_GB = re.findall('GENERALIZED BORN:[A-Z0-9\w\W\n]+?ΔTOTAL[ ]+([-.0-9]+)[ ]+([-.0-9]+)[ ]+([-.0-9]+)[ ]+([-.0-9]+)[ ]+', data)
    delta_total_GB = re.findall(r'GENERALIZED BORN:[A-Z0-9\w\W\n]+?ΔTOTAL[ ]+([-.0-9 ]+)\n', data)
    delta_total_PB = re.findall(r'POISSON BOLTZMANN:[A-Z0-9\w\W\n]+?ΔTOTAL[ ]+([-.0-9 ]+)\n', data)

    # G_binding_GB = re.findall('GENERALIZED BORN:[A-Z0-9\w\W\n]+?Using Interaction Entropy Approximation:\nΔG binding[ =]+([0-9-.]+)[ +/\-]+?([0-9-.]+)\n', data)
    G_binding_GB = re.findall(
        r'GENERALIZED BORN:[A-Z0-9\w\W\n]+?Using Interaction Entropy Approximation:\nΔG binding[ =]+([0-9+\-\./ ]+)\n',
        data)
    G_binding_PB = re.findall(
        r'POISSON BOLTZMANN:[A-Z0-9\w\W\n]+?Using Interaction Entropy Approximation:\nΔG binding[ =]+([0-9+\-\./ ]+)\n',
        data)
    name = pathlib.PurePath(fname).parent.name
    out_res = {'GBSA': {'Name': name, 'directory': os.path.dirname(fname)},
               'PBSA': {'Name': name, 'directory': os.path.dirname(fname)}}

    if G_binding_GB:
        out_res['GBSA'].update(get_Gbinding_values(G_binding_GB))
    if G_binding_PB:
        out_res['PBSA'].update(get_Gbinding_values(G_binding_PB))

    if delta_total_columns_GB and delta_total_GB:
        out_res['GBSA'].update(get_delta_total_values(
            delta_total_columns_re=delta_total_columns_GB,
            delta_total_values_re=delta_total_GB))
    if delta_total_columns_PB and delta_total_PB:
        out_res['PBSA'].update(get_delta_total_values(
            delta_total_columns_re=delta_total_columns_PB,
            delta_total_values_re=delta_total_PB))

    if IE_GB:
        out_res['GBSA'].update(get_IE_values(IE_GB))
    if IE_PB:
        out_res['PBSA'].update(get_IE_values(IE_PB))

    return out_res

def run_get_frames_from_wdir(wdir, xtc, env):
    return get_number_of_frames(os.path.join(wdir, xtc), env=env)


def get_mmpbsa_start_end_interval(mmpbsa):
    with open(mmpbsa) as inp:
        mmpbsa_data = inp.read()
    logging.info(f'{mmpbsa}:\n{mmpbsa_data}')

    startframe, endframe, interval = None, None, None
    for line in mmpbsa_data.split('\n'):
        if line.startswith('#'):
            continue
        line = line.strip()
        if 'startframe' in line:
            startframe = re.findall('startframe[ ]*=[ ]*([0-9]*)', line)
        if 'endframe' in line:
            endframe = re.findall('endframe[ ]*=[ ]*([0-9]*)', line)
        if 'interval' in line:
            interval = re.findall('interval[ ]*=[ ]*([0-9]*)', line)

    startframe = int(startframe[0]) if startframe else 1
    endframe = int(endframe[0]) if endframe else 9999999  # default value of gmxMMBPSA
    interval = int(interval[0]) if interval else 1

    return startframe, endframe, interval

def get_used_number_of_frames(var_number_of_frames, startframe, endframe, interval):
    return math.ceil((min(min(var_number_of_frames), endframe) - (startframe - 1)) / interval)

def start(wdir_to_run, tpr, xtc, topol, index, out_wdir, mmpbsa, ncpu, ligand_resid,
          append_protein_selection, hostfile, unique_id, bash_log,
          gmxmmpbsa_out_files=None, clean_previous=False):
    '''

    :param wdir_to_run: list
    :param tpr: path to file
    :param xtc: path to file
    :param topol: path to file
    :param index: path to file
    :param out_wdir: dir path
    :param mmpbsa: path to file
    :param ncpu: iny
    :param ligand_resid: str
    :param append_protein_selection: None or str
    :param hostfile: None or path to file
    :param unique_id: str (unique id for output files)
    :param bash_log: file name
    :param gmxmmpbsa_out_files: pathes to files
    :param clean_previous: bool
    :return:
    '''
    dask_client, cluster, pool = None, None, None
    var_gbsa_out_files = []
    if gmxmmpbsa_out_files is None:
        # gmx_mmpbsa requires that the run must have at least as many frames as processors. Thus we get and use the min number of used frames as NP
        if not mmpbsa:
            mmpbsa = os.path.join(out_wdir, f'mmpbsa_{unique_id}.in')
            project_dir = os.path.dirname(os.path.abspath(__file__))
            shutil.copy(os.path.join(project_dir, 'scripts', 'gbsa', 'mmpbsa.in'), mmpbsa)
            logging.warning(f'No mmpbsa.in file was set up. Template will be used. Created file: {mmpbsa}.')

        startframe, endframe, interval = get_mmpbsa_start_end_interval(mmpbsa)

        if wdir_to_run is not None:
            var_number_of_frames = []
            with Pool(ncpu) as pool:
                for res in pool.imap_unordered(partial(run_get_frames_from_wdir,
                                    xtc=xtc, env=os.environ.copy()), wdir_to_run):
                    if res:
                        var_number_of_frames.append(res[0])

            if not var_number_of_frames:
                logging.error(f'Could not parse number of frames from all xtc files: {wdir_to_run}. '
                              f'Calculations will be interrupted')
                return None

            used_number_of_frames = get_used_number_of_frames(var_number_of_frames=var_number_of_frames,
                                                              startframe=startframe,
                                                              endframe=endframe,
                                                              interval=interval)
            n_tasks_per_node = ncpu // min(ncpu, used_number_of_frames)
            logging.info(f'{used_number_of_frames} frames will be used')
            logging.info(f'{min(ncpu, used_number_of_frames)} NP will be used')
            # run energy calculation
            try:
                dask_client, cluster = init_dask_cluster(hostfile=hostfile, n_tasks_per_node=n_tasks_per_node,
                                                         ncpu=ncpu)
                var_gbsa_out_files = []
                for res in calc_dask(run_gbsa_from_wdir, wdir_to_run, dask_client=dask_client,
                                     tpr=tpr, xtc=xtc, topol=topol, index=index,
                                     mmpbsa=mmpbsa, np=min(ncpu, used_number_of_frames),
                                     ligand_resid=ligand_resid,
                                     append_protein_selection=append_protein_selection,
                                     unique_id=unique_id, env=os.environ.copy(),
                                     bash_log=bash_log, clean_previous=clean_previous):
                    if res:
                        var_gbsa_out_files.append(res)
            finally:
                if dask_client:
                    dask_client.retire_workers(dask_client.scheduler_info()['workers'],
                                               close_workers=True, remove=True)
                    dask_client.shutdown()
                if cluster:
                    cluster.close()

        elif tpr is not None and xtc is not None and topol is not None and index is not None:
            number_of_frames, _ = get_number_of_frames(xtc, env=os.environ.copy())
            used_number_of_frames = math.ceil((min(number_of_frames, endframe) - (startframe - 1)) / interval)
            logging.info(f'{min(ncpu, used_number_of_frames)} NP will be used')
            if used_number_of_frames <= 0:
                logging.error('Used number of frames are less or equal than 0. Run will be interrupted')
                raise ValueError
            res = run_gbsa_task(wdir=os.path.dirname(xtc), tpr=tpr, xtc=xtc, topol=topol, index=index, mmpbsa=mmpbsa,
                          np=min(ncpu, used_number_of_frames), ligand_resid=ligand_resid, append_protein_selection=append_protein_selection,
                          unique_id=unique_id, env=os.environ.copy(),
                          bash_log=bash_log, clean_previous=clean_previous)
            var_gbsa_out_files.append(res)

    else:
        var_gbsa_out_files = gmxmmpbsa_out_files

    # collect energies
    if var_gbsa_out_files:
        GBSA_output_res, PBSA_output_res = [], []
        with Pool(ncpu) as pool:
            for res in pool.imap_unordered(parse_gmxMMPBSA_output, var_gbsa_out_files):
                if res:
                    GBSA_output_res.append(res['GBSA'])
                    PBSA_output_res.append(res['PBSA'])

        pd_gbsa = pd.DataFrame(GBSA_output_res).sort_values('Name')
        pd_pbsa = pd.DataFrame(PBSA_output_res).sort_values('Name')

        if list(pd_gbsa.columns) != ['Name']:
            pd_gbsa.to_csv(os.path.join(out_wdir, f'GBSA_output_{unique_id}.csv'), sep='\t', index=False)
        if list(pd_pbsa.columns) != ['Name']:
            pd_pbsa.to_csv(os.path.join(out_wdir, f'PBSA_output_{unique_id}.csv'), sep='\t', index=False)

        finished_complexes_file = os.path.join(out_wdir, f"finished_gbsa_files_{unique_id}.txt")
        with open(finished_complexes_file, 'w') as output:
            output.write("\n".join(var_gbsa_out_files))

        logging.info(
            f'gmxMMPBSA energy calculation of {len(var_gbsa_out_files)} were successfully finished.\n'
            f'Successfully finished complexes have been saved in {finished_complexes_file} file')


def main():
    parser = argparse.ArgumentParser(description='''Run MM-GBSA/MM-PBSA calculation using gmx_MMPBSA tool''')
    parser.add_argument('-i', '--wdir_to_run', metavar='DIRNAME', required=False, default=None, nargs='+',
                        type=partial(filepath_type, exist_type='dir'),
                        help='''single or multiple directories for simulations.
                             Should consist of: tpr, xtc, ndx files''')
    parser.add_argument('--topol', metavar='topol.top', required=False, default=None, type=filepath_type,
                        help='topol file from the the MD simulation. Will be ignored if --wdir_to_run is used')
    parser.add_argument('--tpr', metavar='md_out.tpr', required=False, default=None, type=filepath_type,
                        help='tpr file from the the MD simulation. Will be ignored if --wdir_to_run is used')
    parser.add_argument('--xtc', metavar='md_fit.xtc', required=False, default=None, type=filepath_type,
                        help='xtc file of the simulation. Trajectory should have no PBC and be fitted on the Protein_Ligand group. '
                             'Will be ignored if --wdir_to_run is used')
    parser.add_argument('--index', metavar='index.ndx', required=False, default=None, type=filepath_type,
                        help='Gromacs index file from the simulation. Will be ignored if --wdir_to_run is used')
    parser.add_argument('-m', '--mmpbsa', metavar='mmpbsa.in', required=False, type=filepath_type,
                        help='MMPBSA input file. If not set up default template will be used.')
    parser.add_argument('-d', '--wdir', metavar='WDIR', default=None,
                        type=partial(filepath_type, check_exist=False, create_dir=True),
                        help='Working directory for program output. If not set the current directory will be used.')
    parser.add_argument('--out_files', nargs='+', default=None, type=filepath_type,
                        help='gmxMMPBSA out files (FINAL*.dat) to parse. If set will be used over other variables.')
    parser.add_argument('--hostfile', metavar='FILENAME', required=False, type=str, default=None,
                        help='text file with addresses of nodes of dask SSH cluster. The most typical, it can be '
                             'passed as $PBS_NODEFILE variable from inside a PBS script. The first line in this file '
                             'will be the address of the scheduler running on the standard port 8786. If omitted, '
                             'calculations will run on a single machine as usual.')
    parser.add_argument('-c', '--ncpu', metavar='INTEGER', required=False, default=len(os.sched_getaffinity(0)), type=int,
                        help='number of CPU per server. Use all available cpus by default.')
    parser.add_argument('--ligand_id', metavar='UNL', default='UNL', help='Ligand residue ID')
    parser.add_argument('-a', '--append_protein_selection', metavar='STRING', required=False, default=None,
                        nargs = '*', help='residue IDs whuch will be included in the protein system (cofactors).'
                             'Example: ZN MG')
    parser.add_argument('--clean_previous', action='store_true', default=False,
                        help=' Clean previous temporary gmxMMPBSA files')
    parser.add_argument('-o','--out_suffix', default=None,
                        help='Unique suffix for output files. By default, start-time_unique-id.'
                             'Unique suffix is used to separate outputs from different runs.')

    args = parser.parse_args()

    if args.wdir is None:
        wdir = os.getcwd()
    else:
        wdir = args.wdir

    out_time = f'{datetime.now().strftime("%d-%m-%Y-%H-%M-%S")}'
    if args.out_suffix:
        unique_id = args.out_suffix
    else:
        import secrets
        out_suffix = secrets.token_hex(3)
        unique_id = f'{out_time}_unique-id-{out_suffix}'

    log_file = os.path.join(wdir, f'log_mmpbsa_{unique_id}.log')
    bash_log = os.path.join(wdir, f'log_mmpbsa_bash_{unique_id}.log')

    if args.wdir_to_run is not None:
        tpr = 'md_out.tpr'
        xtc = 'md_fit.xtc'
        topol = 'topol.top'
        index = 'index.ndx'
    else:
        tpr = args.tpr
        xtc = args.xtc
        topol = args.topol
        index = args.index

    logging.basicConfig(format='%(asctime)s - %(levelname)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S',
                        level=logging.INFO,
                        handlers=[logging.FileHandler(log_file),
                                  logging.StreamHandler()])


    logging.info(args)
    try:
        start(tpr=tpr, xtc=xtc, topol=topol,
              index=index, out_wdir=wdir, wdir_to_run=args.wdir_to_run,
              mmpbsa=args.mmpbsa, ncpu=args.ncpu, unique_id=unique_id,
              gmxmmpbsa_out_files=args.out_files, ligand_resid=args.ligand_id,
              append_protein_selection=args.append_protein_selection,
              hostfile=args.hostfile, bash_log=bash_log, clean_previous=args.clean_previous)
    finally:
        logging.shutdown()
