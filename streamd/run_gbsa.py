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
import re
import shutil
import subprocess
from datetime import datetime
from functools import partial
from multiprocessing import cpu_count

import pandas as pd

from streamd.utils.dask_init import init_dask_cluster, calc_dask
from streamd.utils.utils import get_index, make_group_ndx, filepath_type, run_check_subprocess


def run_gbsa_task(wdir, tpr, xtc, topol, index, mmpbsa, np, ligand_resid, append_protein_selection, out_time,
                  env, bash_log, clean_previous):
    def calc_gbsa(wdir, tpr, xtc, topol, index, mmpbsa, np, protein_index, ligand_index, out_time, env, bash_log):
        output = os.path.join(wdir, f"FINAL_RESULTS_MMPBSA_{out_time}.dat")
        cmd = f'cd {wdir}; mpirun -np {np} gmx_MMPBSA MPI -O -i {mmpbsa} ' \
              f' -cs {tpr} -ci {index} -cg {protein_index} {ligand_index} -ct {xtc} -cp {topol} -nogui ' \
              f'-o {output} ' \
              f'-eo {os.path.join(wdir, f"FINAL_RESULTS_MMPBSA_{out_time}.csv")}' \
              f' >> {os.path.join(wdir, bash_log)} 2>&1'
        if not run_check_subprocess(cmd, key=xtc, log=os.path.join(wdir, bash_log), env=env):
            return None
        return output

    if clean_previous:
        clean_temporary_gmxMMBPSA_files(wdir)

    if not os.path.isfile(tpr) or not os.path.isfile(xtc) or not os.path.isfile(topol) or not os.path.isfile(index):
        logging.warning(f'{wdir} cannot run gbsa. Check if there are missing files: {tpr} {xtc} {topol} {index}')
        return None

    index_list = get_index(index)
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
                if not make_group_ndx(query, wdir, bash_log=bash_log):
                    return None
                index_list = get_index(index)

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
                       out_time=out_time,
                       env=env,
                       bash_log=bash_log)

    if os.path.isfile(os.path.join(wdir, 'gmx_MMPBSA.log')):
        shutil.copy(os.path.join(wdir, 'gmx_MMPBSA.log'), os.path.join(wdir, f'gmx_MMPBSA_{out_time}.log'))

    clean_temporary_gmxMMBPSA_files(wdir)

    return output


def run_gbsa_from_wdir(wdir, tpr, xtc, topol, index, mmpbsa, np, ligand_resid, append_protein_selection, out_time,
                       env, bash_log, clean_previous):
    tpr = os.path.join(wdir, tpr)
    xtc = os.path.join(wdir, xtc)
    topol = os.path.join(wdir, topol)
    index = os.path.join(wdir, index)
    return run_gbsa_task(wdir, tpr, xtc, topol, index, mmpbsa, np, ligand_resid, append_protein_selection, out_time,
                         env, bash_log, clean_previous)


def clean_temporary_gmxMMBPSA_files(wdir):
    # remove intermediate files
    try:
        subprocess.check_output(f'cd {wdir}; gmx_MMPBSA --clean', shell=True)
    except subprocess.CalledProcessError as e:
        logging.error(f'{wdir}\t{e}\n')
        return None


def parse_gmxMMPBSA_output(fname):
    def get_IE_values(IE_parsed_out):
        IE_res = {}
        IE_columns = [i.strip() for i in IE_parsed_out[0][0].split('  ') if i]
        IE_values = [i.strip() for i in IE_parsed_out[0][1].split('  ') if i]
        for n, i in enumerate(IE_columns):
            IE_res[f'IE{i}'] = IE_values[n]
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
        'Energy Method[ ]*?Entropy[ ]*?(σ\(Int. Energy\)[ ]*?Average[ ]*?SD[ ]*?SEM)\n[-]*?\nGB[ ]*?IE[ ]*?([0-9\. ]*)\n',
        data)
    IE_PB = re.findall(
        'Energy Method[ ]*?Entropy[ ]*?(σ\(Int. Energy\)[ ]*?Average[ ]*?SD[ ]*?SEM)\n[-]*?\nPB[ ]*?IE[ ]*?([0-9\. ]*)\n',
        data)

    delta_total_columns_GB = re.findall(
    'GENERALIZED BORN:[A-Z0-9\w\W\n]+?Delta \(Complex - Receptor - Ligand\):\n([A-Za-z \(\).]+)\n',
    data)
    delta_total_columns_PB = re.findall(
    'POISSON BOLTZMANN:[A-Z0-9\w\W\n]+?Delta \(Complex - Receptor - Ligand\):\n([A-Za-z \(\).]+)\n',
    data)

    # delta_total_GB = re.findall('GENERALIZED BORN:[A-Z0-9\w\W\n]+?ΔTOTAL[ ]+([-.0-9]+)[ ]+([-.0-9]+)[ ]+([-.0-9]+)[ ]+([-.0-9]+)[ ]+', data)
    delta_total_GB = re.findall('GENERALIZED BORN:[A-Z0-9\w\W\n]+?ΔTOTAL[ ]+([-.0-9 ]+)\n', data)
    delta_total_PB = re.findall('POISSON BOLTZMANN:[A-Z0-9\w\W\n]+?ΔTOTAL[ ]+([-.0-9 ]+)\n', data)

    # G_binding_GB = re.findall('GENERALIZED BORN:[A-Z0-9\w\W\n]+?Using Interaction Entropy Approximation:\nΔG binding[ =]+([0-9-.]+)[ +/\-]+?([0-9-.]+)\n', data)
    G_binding_GB = re.findall(
        'GENERALIZED BORN:[A-Z0-9\w\W\n]+?Using Interaction Entropy Approximation:\nΔG binding[ =]+([0-9+\-\./ ]+)\n',
        data)
    G_binding_PB = re.findall(
        'POISSON BOLTZMANN:[A-Z0-9\w\W\n]+?Using Interaction Entropy Approximation:\nΔG binding[ =]+([0-9+\-\./ ]+)\n',
        data)

    out_res = {'GBSA': {'Name': fname}, 'PBSA': {'Name': fname}}
    if IE_GB:
        out_res['GBSA'].update(get_IE_values(IE_GB))
    if IE_PB:
        out_res['PBSA'].update(get_IE_values(IE_PB))
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

    return out_res


def get_number_of_frames(xtc, env):
    res = subprocess.run(f'gmx check -f {xtc}', shell=True, capture_output=True, env=env)
    frames = re.findall('Step[ ]*([0-9]*)[ ]*[0-9]*\n', res.stderr.decode("utf-8"))
    if frames:
        logging.info(f'{xtc} has {frames[0]} frames')
        return int(frames[0])


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


def start(wdir_to_run, tpr, xtc, topol, index, out_wdir, mmpbsa, ncpu, ligand_resid, append_protein_selection,
          hostfile, out_time, bash_log,
          gmxmmpbsa_out_files=None, clean_previous=False):
    dask_client, cluster = None, None
    var_gbsa_out_files = []
    if gmxmmpbsa_out_files is None:
        # gmx_mmpbsa requires that the run must have at least as many frames as processors. Thus we get and use the min number of used frames as NP
        if not mmpbsa:
            mmpbsa = os.path.join(out_wdir, f'mmpbsa_{out_time}.in')
            project_dir = os.path.dirname(os.path.abspath(__file__))
            shutil.copy(os.path.join(project_dir, 'scripts', 'gbsa', 'mmpbsa.in'), mmpbsa)
            logging.warning(f'No mmpbsa.in file was set up. Template will be used. Created file: {mmpbsa}.')

        startframe, endframe, interval = get_mmpbsa_start_end_interval(mmpbsa)
        if wdir_to_run is not None:
            try:
                dask_client, cluster = init_dask_cluster(hostfile=hostfile,
                                                         n_tasks_per_node=min(ncpu, len(wdir_to_run)),
                                                         use_multi_servers=True if len(wdir_to_run) > ncpu else False,
                                                         ncpu=ncpu)
                var_number_of_frames = []
                for res in calc_dask(run_get_frames_from_wdir, wdir_to_run, dask_client=dask_client,
                                     xtc=xtc, env=os.environ.copy()):
                    if res:
                        var_number_of_frames.append(res)
            finally:
                if dask_client:
                    dask_client.retire_workers(dask_client.scheduler_info()['workers'],
                                               close_workers=True, remove=True)
                    dask_client.shutdown()
                if cluster:
                    cluster.close()

            used_number_of_frames = math.ceil((min(min(var_number_of_frames), endframe) - (startframe - 1)) / interval)
            n_tasks_per_node = ncpu // min(ncpu, used_number_of_frames)
            #todo 64 2 mol 32 booked -> 34 use

            logging.info(f'{min(ncpu, used_number_of_frames)} NP will be used')
            # run energy calculation
            try:
                dask_client, cluster = init_dask_cluster(hostfile=hostfile, n_tasks_per_node=n_tasks_per_node,
                                                         ncpu=ncpu)
                var_gbsa_out_files = []
                for res in calc_dask(run_gbsa_from_wdir, wdir_to_run, dask_client=dask_client,
                                     tpr=tpr, xtc=xtc, topol=topol, index=index,
                                     mmpbsa=mmpbsa, np=min(ncpu, used_number_of_frames), ligand_resid=ligand_resid,
                                     append_protein_selection=append_protein_selection,
                                     out_time=out_time, env=os.environ.copy(),
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
            number_of_frames = get_number_of_frames(xtc, env=None)
            used_number_of_frames = math.ceil((min(number_of_frames, endframe) - (startframe - 1)) / interval)
            logging.info(f'{min(ncpu, used_number_of_frames)} NP will be used')
            if used_number_of_frames <= 0:
                logging.error('Used number of frames are less or equal than 0. Run will be interrupted')
                raise ValueError
            run_gbsa_task(wdir=os.path.dirname(xtc), tpr=tpr, xtc=xtc, topol=topol, index=index, mmpbsa=mmpbsa,
                          np=min(ncpu, used_number_of_frames), ligand_resid=ligand_resid, append_protein_selection=append_protein_selection,
                          out_time=out_time, env=os.environ.copy(),
                          bash_log=bash_log, clean_previous=clean_previous)

    else:
        var_gbsa_out_files = gmxmmpbsa_out_files

    # collect energies
    if var_gbsa_out_files:
        GBSA_output_res, PBSA_output_res = [], []
        try:
            dask_client, cluster = init_dask_cluster(hostfile=hostfile, n_tasks_per_node=len(var_gbsa_out_files),
                                                     ncpu=ncpu)
            for res in calc_dask(parse_gmxMMPBSA_output, var_gbsa_out_files, dask_client=dask_client):
                if res:
                    GBSA_output_res.append(res['GBSA'])
                    PBSA_output_res.append(res['PBSA'])
        finally:
            if dask_client:
                dask_client.retire_workers(dask_client.scheduler_info()['workers'],
                                           close_workers=True, remove=True)
                dask_client.shutdown()
            if cluster:
                cluster.close()

        pd_gbsa = pd.DataFrame(GBSA_output_res).sort_values('Name')
        pd_pbsa = pd.DataFrame(PBSA_output_res).sort_values('Name')

        if list(pd_gbsa.columns) != ['Name']:
            pd_gbsa.to_csv(os.path.join(out_wdir, f'GBSA_output_{out_time}.csv'), sep='\t', index=False)
        if list(pd_pbsa.columns) != ['Name']:
            pd_pbsa.to_csv(os.path.join(out_wdir, f'PBSA_output_{out_time}.csv'), sep='\t', index=False)

        logging.info(
            f'gmxMMPBSA energy calculation of {len(var_gbsa_out_files)} were successfully finished.\nFinished: {var_gbsa_out_files}\n')


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
    parser.add_argument('-c', '--ncpu', metavar='INTEGER', required=False, default=cpu_count(), type=int,
                        help='number of CPU per server. Use all cpus by default.')
    parser.add_argument('--ligand_id', metavar='UNL', default='UNL', help='Ligand residue ID')
    parser.add_argument('-a', '--append_protein_selection', metavar='STRING', required=False, default=None,
                        nargs = '*', help='residue IDs whuch will be included in the protein system (cofactors).'
                             'Example: ZN MG')
    parser.add_argument('--clean_previous', action='store_true', default=False,
                        help=' Clean previous temporary gmxMMPBSA files')

    args = parser.parse_args()

    if args.wdir is None:
        wdir = os.getcwd()
    else:
        wdir = args.wdir

    out_time = f'{datetime.now().strftime("%d-%m-%Y-%H-%M-%S")}'
    log_file = os.path.join(wdir, f'log_mmpbsa_{out_time}.log')
    bash_log = os.path.join(wdir, f'log_mmpbsa_bash_{out_time}.log')

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

    logging.basicConfig(format='%(asctime)s - %(name)s - %(levelname)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S',
                        level=logging.INFO,
                        handlers=[logging.FileHandler(log_file),
                                  logging.StreamHandler()])

    logging.getLogger('distributed').setLevel('CRITICAL')
    logging.getLogger('asyncssh').setLevel('CRITICAL')
    logging.getLogger('distributed.worker').setLevel('CRITICAL')
    logging.getLogger('distributed.core').setLevel('CRITICAL')
    logging.getLogger('distributed.comm').setLevel('CRITICAL')
    logging.getLogger('distributed.nanny').setLevel('CRITICAL')
    logging.getLogger('bockeh').setLevel('CRITICAL')

    logging.info(args)
    try:
        start(tpr=tpr, xtc=xtc, topol=topol,
              index=index, out_wdir=wdir, wdir_to_run=args.wdir_to_run,
              mmpbsa=args.mmpbsa, ncpu=args.ncpu, out_time=out_time,
              gmxmmpbsa_out_files=args.out_files, ligand_resid=args.ligand_id, append_protein_selection=args.append_protein_selection,
              hostfile=args.hostfile, bash_log=bash_log, clean_previous=args.clean_previous)
    finally:
        logging.shutdown()
