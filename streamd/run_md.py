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
import os
import shutil
from datetime import datetime
from functools import partial
from glob import glob
import json

from streamd.md_analysis import run_md_analysis
from streamd.preparation.complex_preparation import run_complex_preparation
from streamd.preparation.ligand_preparation import prepare_input_ligands, check_mols
from streamd.utils.dask_init import init_dask_cluster, calc_dask
from streamd.utils.utils import filepath_type, run_check_subprocess, get_protein_resid_set
from streamd.mcpbpy_md import mcbpy_md


class RawTextArgumentDefaultsHelpFormatter(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass


def run_equilibration(wdir, project_dir, bash_log, ncpu, compute_device, device_param, gpu_args, env=None):
    if os.path.isfile(os.path.join(wdir, 'npt.gro')) and os.path.isfile(os.path.join(wdir, 'npt.cpt')):
        logging.warning(f'{wdir}. Checkpoint files after Equilibration step exist. '
                        f'Equilibration step will be skipped ')
        return wdir
    cmd = (f'wdir={wdir} ncpu={ncpu} compute_device={compute_device} device_param={device_param} gpu_args={gpu_args} '
           f'bash {os.path.join(project_dir, "scripts/script_sh/equlibration.sh")} '
           f'>> {os.path.join(wdir, bash_log)} 2>&1'),
    if not run_check_subprocess(cmd, wdir, log=os.path.join(wdir, bash_log), env=env):
        return None
    return wdir


def run_simulation(wdir, project_dir, bash_log, mdtime_ns,
                   tpr, cpt, xtc, deffnm, deffnm_next, ncpu,
                   compute_device, device_param, gpu_args, env=None):
    # continue/extend simulation if checkpoint files exist
    if (tpr is not None and os.path.isfile(tpr) and cpt is not None and os.path.isfile(cpt) and xtc is not None and os.path.isfile(str(xtc))) or \
        (os.path.isfile(os.path.join(wdir, f'{deffnm}.tpr')) and os.path.isfile(os.path.join(wdir, f'{deffnm}.cpt'))
         and os.path.isfile(os.path.join(wdir, f'{deffnm}.xtc'))) :
        logging.warning(f'{wdir}. {deffnm}.xtc and {deffnm}.tpr and  {deffnm}.cpt exist. '
                        f'MD simulation will be continued until the setup simulation steps are reached.')
        if continue_md_from_dir(wdir_to_continue=wdir, tpr=tpr, cpt=cpt, xtc=xtc,
                                deffnm=deffnm, deffnm_next=deffnm_next,
                                mdtime_ns=mdtime_ns, project_dir=project_dir, bash_log=bash_log,
                                ncpu=ncpu, compute_device=compute_device, device_param=device_param,
                                gpu_args=gpu_args, env=env) is None:
            return None

        return (wdir, deffnm)
    cmd = (f'wdir={wdir} ncpu={ncpu} compute_device={compute_device} gpu_args={gpu_args} device_param={device_param} deffnm={deffnm} '
           f'bash {os.path.join(project_dir, "scripts/script_sh/md.sh")} >> {os.path.join(wdir, bash_log)} 2>&1')
    if not run_check_subprocess(cmd, wdir, log=os.path.join(wdir, bash_log), env=env):
        return None
    return (wdir, deffnm)


def continue_md_from_dir(wdir_to_continue, tpr, cpt, xtc, deffnm, deffnm_next,
                         mdtime_ns, project_dir, bash_log, ncpu, compute_device, device_param, gpu_args, env=None):
    def continue_md(tpr, cpt, xtc, wdir, new_mdtime_ps, deffnm_next, project_dir, bash_log, compute_device, env):
        cmd = f'wdir={wdir} tpr={tpr} cpt={cpt} xtc={xtc} new_mdtime_ps={new_mdtime_ps} ' \
              f'deffnm_next={deffnm_next} ncpu={ncpu} compute_device={compute_device} device_param={device_param} gpu_args={gpu_args} bash {os.path.join(project_dir, "scripts/script_sh/continue_md.sh")}' \
              f'>> {os.path.join(wdir, bash_log)} 2>&1'
        if run_check_subprocess(cmd, wdir, log=os.path.join(wdir, bash_log), env=env):
            return wdir
        return None

    def backup_prev_files(file_to_backup, wdir):
        n = len(glob(os.path.join(wdir, f'#{os.path.basename(file_to_backup)}.*#'))) + 1
        new_f = os.path.join(wdir, f'#{os.path.basename(file_to_backup)}.{n}#')
        shutil.move(file_to_backup, new_f)
        logging.warning(f'Backup previous file {file_to_backup} to {new_f}')

    if tpr is None:
        tpr = os.path.join(wdir_to_continue, f'{deffnm}.tpr')
    if cpt is None:
        cpt = os.path.join(wdir_to_continue, f'{deffnm}.cpt')
    if xtc is None:
        xtc = os.path.join(wdir_to_continue, f'{deffnm}.xtc')

    for i in [tpr, cpt, xtc]:
        if not os.path.isfile(i):
            logging.exception(
                f'No {i} file was found. Cannot continue the simulation. Calculations will be interrupted ',
                stack_info=True)
            return None

    new_mdtime_ps = int(mdtime_ns * 1000)

    if continue_md(tpr=tpr, cpt=cpt, xtc=xtc, wdir=wdir_to_continue,
                   new_mdtime_ps=new_mdtime_ps, deffnm_next=deffnm_next, project_dir=project_dir,
                   compute_device=compute_device, env=env, bash_log=bash_log):
        for f in glob(os.path.join(wdir_to_continue, f'{deffnm_next}.*')):
            # check previous existing files with the same name
            backup_prev_files(file_to_backup=os.path.join(wdir_to_continue, os.path.basename(f).replace(deffnm_next, deffnm)),
                              wdir=wdir_to_continue)
            shutil.move(f, os.path.join(wdir_to_continue, os.path.basename(f).replace(deffnm_next, deffnm)))

        return wdir_to_continue


def start(protein, wdir, lfile, system_lfile, noignh,
          forcefield_name, npt_time_ps, nvt_time_ps, mdtime_ns,
          topol, topol_itp_list, posre_list_protein,
          wdir_to_continue_list, deffnm,
          tpr_prev, cpt_prev, xtc_prev, ligand_list_file_prev, ligand_resid,
          activate_gaussian, gaussian_exe, gaussian_basis, gaussian_memory,
          metal_resnames, metal_charges, mcpbpy_cut_off,
          seed, steps, hostfile, ncpu, mdrun_per_node, compute_device, gpu_ids, ntmpi_per_gpu, clean_previous,
          not_clean_backup_files, out_time, bash_log=None):
    '''
    :param protein: protein file - pdb or gro format
    :param wdir: None or path
    :param lfile: None or file
    :param system_lfile: None or file. Mol or sdf format
    :param noignh: don't use -ignh argument (gmx pdb2gmx)
    :param forcefield_name: str
    :param clean_previous: boolean. Remove all previous md files
    :param mdtime_ns: float. Time in ns
    :param npt_time_ps: int. Time in ps
    :param nvt_time_ps: int. Time in ps
    :param topol: None or file
    :param topol_itp_list: None or list of files
    :param posre_list_protein: None or list of files
    :param wdir_to_continue_list: list of paths
    :param tpr_prev: None or file
    :param cpt_prev: None or file
    :param xtc_prev: None or file
    :param ligand_resid: UNL. Used for md analysis only if continue simulation
    :param ligand_list_file_prev: None or file
    :param deffnm_prev: md_out
    :param hostfile: None or file
    :param ncpu:
    :param compute_device:
    not_clean_log_files: boolean. Remove backup md files (starts with #)
    :return:
    '''

    project_dir = os.path.dirname(os.path.abspath(__file__))
    script_path = os.path.join(project_dir, 'scripts')
    script_mdp_path = os.path.join(script_path, 'mdp')

    wdir_md = os.path.join(wdir, 'md_files', 'md_run')

    dask_client, cluster = None, None

    # GPU calculations settings
    gpu_args = ''
    # To set where to execute (cpu or gpu) the interactions and update steps during gmx mdrun
    device_param = f"-update {compute_device} -pme {compute_device} -bonded {compute_device} -pmefft {compute_device}"

    if compute_device == 'gpu' or gpu_ids:
        ngpus = len(gpu_ids) if gpu_ids else 1
        # https://gromacs.bioexcel.eu/t/using-multiple-gpus-on-one-machine/5974
        k = ngpus * ntmpi_per_gpu
        if k > 1:
            device_param = f"{device_param} -npme 1"
        gpu_args = f"-ntmpi {k} -ntomp {(ncpu // k) // mdrun_per_node}"
        if gpu_ids:
            gpu_args = gpu_args + f" -gpu_id {','.join(gpu_ids)}"
        gpu_args = f"'{gpu_args}'"

    device_param = f"'{device_param}'"

    # Start
    if tpr_prev is None or cpt_prev is None or xtc_prev is None:
        # preparation
        if (steps is None or 1 in steps) and wdir_to_continue_list is None:
            # create dirs
            ligand_resid = 'UNL'
            pname, p_ext = os.path.splitext(os.path.basename(protein))

            wdir_protein = os.path.join(wdir, 'md_files', 'md_preparation', 'protein', pname)
            wdir_ligand = os.path.join(wdir, 'md_files', 'md_preparation', 'ligands')
            wdir_system_ligand = os.path.join(wdir, 'md_files', 'md_preparation', 'cofactors')
            wdir_metal = os.path.join(wdir, 'md_files', 'md_preparation', 'metals', pname)

            os.makedirs(wdir_md, exist_ok=True)
            os.makedirs(wdir_protein, exist_ok=True)
            os.makedirs(wdir_ligand, exist_ok=True)
            os.makedirs(wdir_system_ligand, exist_ok=True)
            os.makedirs(wdir_metal, exist_ok=True)

            # check if already exist in the working directory
            if not metal_resnames or (not gaussian_exe or not activate_gaussian):
                if not os.path.isfile(f'{os.path.join(wdir_protein, pname)}.gro') or not os.path.isfile(
                        os.path.join(wdir_protein, "topol.top")):
                    if p_ext != '.gro' or topol is None or posre_list_protein is None:
                        logging.info('Start protein preparation')
                        cmd = f'gmx pdb2gmx -f {protein} -o {os.path.join(wdir_protein, pname)}.gro -water tip3p {"-ignh" if not noignh else "-noignh"} ' \
                              f'-i {os.path.join(wdir_protein, "posre.itp")} ' \
                              f'-p {os.path.join(wdir_protein, "topol.top")} -ff {forcefield_name} >> {os.path.join(wdir_protein, bash_log)} 2>&1'
                        if not run_check_subprocess(cmd, protein, log=os.path.join(wdir_protein, bash_log)):
                            return None
                        logging.info(f'Successfully finished protein preparation\n')
                    else:
                        target_path = os.path.join(wdir_protein, os.path.basename(protein))
                        if not os.path.isfile(target_path):
                            shutil.copy(protein, target_path)
                        target_path = os.path.join(wdir_protein, 'topol.top')
                        if not os.path.isfile(target_path):
                            shutil.copy(topol, target_path)
                        # multiple chains
                        for posre_protein in posre_list_protein:
                            target_path = os.path.join(wdir_protein, os.path.basename(posre_protein))
                            if not os.path.isfile(target_path):
                                shutil.copy(posre_protein, target_path)
                        if topol_itp_list is not None:
                            if len(posre_list_protein) != len(topol_itp_list):
                                logging.exception(
                                    'The number of protein_chainX.itp files should be equal the number of posre_protein_chainX.itp files.'
                                    ' Check --topol_itp and --posre arguments')
                                return None
                            for topol_itp in topol_itp_list:
                                target_path = os.path.join(wdir_protein, os.path.basename(topol_itp))
                                if not os.path.isfile(target_path):
                                    shutil.copy(topol_itp, target_path)

                else:
                    logging.warning(f'{os.path.join(wdir_protein, pname)}.gro and topol.top files exist. '
                                    f'Protein preparation step will be skipped.')

            # Part 1. Ligand Preparation
            protein_resid_set = get_protein_resid_set(protein)
            if system_lfile is not None:
                logging.info('Start cofactor preparation')
                number_of_mols, problem_mols = check_mols(system_lfile)
                if problem_mols:
                    logging.exception(f'Cofactor molecules: {problem_mols} from {system_lfile} cannot be processed. Script will be interrupted.')
                    return None

                system_lig_wdirs = prepare_input_ligands(system_lfile, preset_resid=None, protein_resid_set=protein_resid_set, script_path=script_path,
                                                         project_dir=project_dir, wdir_ligand=wdir_system_ligand,
                                                         gaussian_exe=gaussian_exe, activate_gaussian=activate_gaussian,
                                                         gaussian_basis=gaussian_basis, gaussian_memory=gaussian_memory,
                                                         hostfile=hostfile, ncpu=ncpu, bash_log=bash_log)
                if number_of_mols != len(system_lig_wdirs):
                    logging.exception(f'Error with cofactor preparation. Only {len(system_lig_wdirs)} from {number_of_mols} preparation were finished.'
                                      f' The calculation will be interrupted')
                    return None

                logging.info(f'Successfully finished {len(system_lig_wdirs)} cofactor preparation\n')
            else:
                system_lig_wdirs = []

            if lfile is not None:
                logging.info('Start ligand preparation')
                number_of_mols, problem_mols = check_mols(lfile)
                if problem_mols:
                    logging.warning(f'Ligand molecules: {problem_mols} from {lfile} cannot be processed.'
                                    f' Such molecules will be skipped.')

                var_lig_wdirs = prepare_input_ligands(lfile, preset_resid=ligand_resid, protein_resid_set=protein_resid_set, script_path=script_path,
                                                      project_dir=project_dir, wdir_ligand=wdir_ligand,
                                                      gaussian_exe=gaussian_exe, activate_gaussian=activate_gaussian,
                                                      gaussian_basis=gaussian_basis, gaussian_memory=gaussian_memory,
                                                      hostfile=hostfile, ncpu=ncpu, bash_log=bash_log)
                if number_of_mols != len(var_lig_wdirs):
                    logging.warning(f'Problem with ligand preparation. Only {len(var_lig_wdirs)} from {number_of_mols} preparation were finished.'
                                    f' Such molecules will be skipped.')

                logging.info(f'Successfully finished {len(var_lig_wdirs)} ligand preparation\n')
            else:
                var_lig_wdirs = [[]]  # run protein in water only simulation

            if not var_lig_wdirs:
                return None
            # Part 2 Complex preparation
            try:
                # make all.itp and create complex
                logging.info('Start complex preparation')
                var_complex_prepared_dirs = []
    # check conflict
                # Part 2.1 MCPBPY Metal-Complex preparation
                if metal_resnames and gaussian_exe and activate_gaussian:
                    logging.info('Start MCPBPY procedure')
                    dask_client, cluster = init_dask_cluster(hostfile=hostfile,
                                                             n_tasks_per_node=1,
                                                             use_multi_servers=True if len(var_lig_wdirs) > ncpu else False,
                                                             ncpu=ncpu)
                    for res in calc_dask(mcbpy_md.main, var_lig_wdirs, dask_client,
                                  protein_name=pname, protein_file=protein,
                                  metal_resnames=metal_resnames, metal_charges=metal_charges,
                                  wdir_metal=wdir_metal, system_lig_wdirs=system_lig_wdirs,
                                  wdir_md=wdir_md, script_path=script_path, ncpu=ncpu,
                                  activate_gaussian=activate_gaussian, gaussian_version=gaussian_exe,
                                  gaussian_basis=gaussian_basis, gaussian_memory=gaussian_memory,
                                  bash_log=bash_log, seed=seed, nvt_time_ps=nvt_time_ps, npt_time_ps=npt_time_ps,
                                  mdtime_ns=mdtime_ns, cut_off=mcpbpy_cut_off, env=os.environ.copy()):
                        if res:
                            var_complex_prepared_dirs.append(res)

                    logging.info('MCPBPY procedure: Finish MCPBPY preparation')
                else:
                    dask_client, cluster = init_dask_cluster(hostfile=hostfile,
                                                             n_tasks_per_node=min(ncpu, len(var_lig_wdirs)),
                                                             use_multi_servers=True if len(var_lig_wdirs) > ncpu else False,
                                                             ncpu=ncpu)
                    for res in calc_dask(run_complex_preparation, var_lig_wdirs, dask_client,
                                         wdir_system_ligand_list=system_lig_wdirs,
                                         protein_name=pname, wdir_protein=wdir_protein,
                                         clean_previous=clean_previous, wdir_md=wdir_md,
                                         script_path=script_mdp_path, project_dir=project_dir, mdtime_ns=mdtime_ns,
                                         npt_time_ps=npt_time_ps, nvt_time_ps=nvt_time_ps, bash_log=bash_log, seed=seed, env=os.environ.copy()):
                        if res:
                            var_complex_prepared_dirs.append(res)

                logging.info(f'Successfully finished {len(var_complex_prepared_dirs)} complexes preparation\n')
            finally:
                if dask_client:
                    dask_client.retire_workers(dask_client.scheduler_info()['workers'],
                                               close_workers=True, remove=True)
                    dask_client.shutdown()
                if cluster:
                    cluster.close()

        else:
            var_complex_prepared_dirs = wdir_to_continue_list

        # Part 3. Equilibration and MD simulation. Run on all cpu
        var_eq_dirs = []
        var_md_dirs_deffnm = []
        if (steps is None or 2 in steps or 3 in steps) and var_complex_prepared_dirs:
            try:
                dask_client, cluster = init_dask_cluster(hostfile=hostfile,
                                                         n_tasks_per_node=mdrun_per_node,
                                                         use_multi_servers=True if len(var_complex_prepared_dirs) > 1 else False,
                                                         ncpu=ncpu)
                if steps is None or 2 in steps:
                    logging.info('Start Equilibration steps')
                    for res in calc_dask(run_equilibration, var_complex_prepared_dirs, dask_client,
                                         project_dir=project_dir, bash_log=bash_log,
                                         ncpu=ncpu//mdrun_per_node, compute_device=compute_device,
                                         device_param=device_param, gpu_args=gpu_args,
                                         env=os.environ.copy()):
                        if res:
                            var_eq_dirs.append(res)
                    logging.info(f'Successfully finished {len(var_eq_dirs)} Equilibration step\n')
                else:
                    var_eq_dirs = wdir_to_continue_list

                if steps is None or 3 in steps:
                    logging.info('Start Simulation step')
                    for res in calc_dask(run_simulation, var_eq_dirs, dask_client,
                                         project_dir=project_dir, bash_log=bash_log,
                                         mdtime_ns=mdtime_ns,
                                         tpr=tpr_prev, cpt=cpt_prev, xtc=xtc_prev,
                                         deffnm=deffnm, deffnm_next=f'{deffnm}_{out_time}',
                                         ncpu=ncpu//mdrun_per_node, compute_device=compute_device,
                                         device_param=device_param, gpu_args=gpu_args,
                                         env=os.environ.copy()):
                        if res:
                            var_md_dirs_deffnm.append(res)
                    logging.info(
                        f'Simulation of {len(var_md_dirs_deffnm)} complexes were successfully finished\nFinished: {var_md_dirs_deffnm}\n')

            finally:
                if dask_client:
                    dask_client.retire_workers(dask_client.scheduler_info()['workers'],
                                               close_workers=True, remove=True)
                    dask_client.shutdown()
                if cluster:
                    cluster.close()


        elif 4 in steps:
            var_md_dirs_deffnm = [(i, deffnm) for i in wdir_to_continue_list]

    # Part 3. MD Analysis. Run on each cpu
    if (steps is None or 4 in steps) and var_md_dirs_deffnm:
        try:
            dask_client, cluster = init_dask_cluster(hostfile=hostfile,
                                                     n_tasks_per_node=min(ncpu, len(var_md_dirs_deffnm)),
                                                     use_multi_servers=True if len(var_md_dirs_deffnm) > ncpu else False,
                                                     ncpu=ncpu)
            logging.info('Start Analysis of the simulations')
            var_md_analysis_dirs = []
            # os.path.dirname(var_lig)
            for res in calc_dask(run_md_analysis, var_md_dirs_deffnm,
                                 dask_client, mdtime_ns=mdtime_ns, project_dir=project_dir,
                                 bash_log=bash_log, ligand_resid=ligand_resid, ligand_list_file_prev=ligand_list_file_prev,
                                 env=os.environ.copy()):
                if res:
                    var_md_analysis_dirs.append(res)
        finally:
            if dask_client:
                dask_client.retire_workers(dask_client.scheduler_info()['workers'],
                                           close_workers=True, remove=True)
                dask_client.shutdown()
            if cluster:
                cluster.close()

        logging.info(
            f'Analysis of MD simulations of {len(var_md_analysis_dirs)} complexes were successfully finished\nFinished: {var_md_analysis_dirs}')

    if not not_clean_backup_files:
        if wdir_to_continue_list is None:
            for f in glob(os.path.join(wdir_md, '*', '#*#')):
                # if '.tpr.' not in f and '.xtc.' not in f:
                #     os.remove(f)
                os.remove(f)
            for f in glob(os.path.join(wdir_md, '*', '*.trr')):
                os.remove(f)
        else:
            for wdir_md in wdir_to_continue_list:
                for f in glob(os.path.join(wdir_md, '#*#')):
                    # if '.tpr.' not in f and '.xtc.' not in f:
                    #     os.remove(f)
                    os.remove(f)
                for f in glob(os.path.join(wdir_md, '*', '*.trr')):
                    os.remove(f)


def main():
    parser = argparse.ArgumentParser(description='''Run or continue MD simulation.\n
    Allowed systems: Protein, Protein-Ligand, Protein-Cofactors(multiple), Protein-Ligand-Cofactors(multiple) ''')
    parser1 = parser.add_argument_group('Standard Molecular Dynamics Simulation Run')
    parser1.add_argument('-p', '--protein', metavar='FILENAME', required=False,
                        type=partial(filepath_type, ext=('pdb', 'gro'), check_exist=True),
                        help='Input file of protein. Supported formats: *.pdb or gro')
    parser1.add_argument('-d', '--wdir', metavar='WDIR', default=None,
                        type=partial(filepath_type, check_exist=False, create_dir=True),
                        help='Working directory. If not set the current directory will be used.')
    parser1.add_argument('-l', '--ligand', metavar='FILENAME', required=False,
                        type=partial(filepath_type, ext=('mol', 'sdf', 'mol2')),
                        help='Input file with compound(s). Supported formats: *.mol or sdf or mol2')
    parser1.add_argument('--cofactor', metavar='FILENAME', default=None,
                        type=partial(filepath_type, ext=('mol', 'sdf', 'mol2')),
                        help='Input file with compound(s). Supported formats: *.mol or sdf or mol2')
    parser1.add_argument('--clean_previous_md', action='store_true', default=False,
                        help='Remove a production MD simulation directory if it exists to re-initialize production MD setup')
    parser1.add_argument('--hostfile', metavar='FILENAME', required=False, type=str, default=None,
                        help='Text file with addresses of nodes of dask SSH cluster. The most typical, it can be '
                             'passed as $PBS_NODEFILE variable from inside a PBS script. The first line in this file '
                             'will be the address of the scheduler running on the standard port 8786. If omitted, '
                             'calculations will run on a single machine as usual.')
    parser1.add_argument('-c', '--ncpu', metavar='INTEGER', required=False,
                         default=len(os.sched_getaffinity(0)), #returns set of CPUs available
                         type=int, help='Number of CPU per server. Use all available cpus by default.')
    parser1.add_argument('--mdrun_per_node', metavar='INTEGER', required=False,
                         default=1, type=int,
                         help='Number of simulations to run per 1 server.'
                              'In case of multiple simulations per node, the available CPU cores will be split evenly across these simulations.'
                              'By default, run only 1 simulation per node and use all available cpus')
    parser1.add_argument('--device', metavar='cpu', required=False, default='auto',
                         choices=['cpu','gpu','auto'], type=lambda x: str(x).lower(),
                         help='Calculate bonded and non-bonded interactions on: auto, cpu, gpu')
    parser1.add_argument('--gpu_ids', metavar='GPU ID', required=False, default=None,
                         nargs='+', type=str, help='List of unique GPU device IDs available to use. '
                                                   'Use in case of multiple GPUs usage or using exact GPU devices.')
    parser1.add_argument('--ntmpi_per_gpu', metavar='int', required=False, default=1,
                         type=int, help='The number of thread-MPI ranks to start per GPU. The default, 1, will start one rank per GPU')
    parser1.add_argument('--topol', metavar='topol.top', required=False, default=None, type=filepath_type,
                        help='Topology file (required if a gro-file is provided for the protein).'
                             'All output files obtained from gmx2pdb should preserve the original names')
    parser1.add_argument('--topol_itp', metavar='topol_chainA.itp topol_chainB.itp', required=False, nargs='+',
                        default=None, type=filepath_type,
                        help='itp files for individual protein chains (required if a gro-file is provided for the protein).'
                             'All output files obtained from gmx2pdb should preserve the original names')
    parser1.add_argument('--posre', metavar='posre.itp', required=False, nargs='+', default=None, type=filepath_type,
                        help='posre file(s) (required if a gro-file is provided for the protein).'
                             'All output files obtained from gmx2pdb should preserve the original names')
    parser1.add_argument('--protein_forcefield', metavar='amber99sb-ildn', required=False, default='amber99sb-ildn', type=str,
                        help='Force Field for protein preparation.'
                             'Available FF can be found at Miniconda3/envs/md/share/gromacs/top')
    parser1.add_argument('--noignh', required=False,  action='store_true', default=False,
                         help="By default, Streamd uses gmx pdb2gmx -ignh, which re-adds hydrogens using residue names "
                              "(the correct protonation states must be provided by user) and ignores the original hydrogens."
                              " If the --noignh argument is used, the original hydrogen atoms will be preserved during the preparation,"
                              " although there may be problems with recognition of atom names by GROMACS.")
    parser1.add_argument('--md_time', metavar='ns', required=False, default=1, type=float,
                        help='Time of MD simulation in ns')
    parser1.add_argument('--npt_time', metavar='ps', required=False, default=100, type=int,
                        help='Time of NPT equilibration in ps')
    parser1.add_argument('--nvt_time', metavar='ps', required=False, default=100, type=int,
                        help='Time of NVT equilibration in ps')
    parser1.add_argument('--seed', metavar='int', required=False, default=-1, type=int,
                        help='seed')
    parser1.add_argument('--not_clean_backup_files', action='store_true', default=False,
                        help='Not to remove all backups of md files')
    parser1.add_argument('--steps', default=None, nargs='*', type=int,
                        help='Run a particular step(s) of the StreaMD run. '
                             'Options:'
                             '1 - run preparation step (protein, ligand, cofactor preparation)\n'
                             '2 - run MD equilibration step (minimization, NVT, NPT)\n'
                             '3 - run MD simulation\n'
                             '4 - run MD analysis.\n'
                             'Ex: 3 4\n'
                             'If 2/3/4 step(s) are used --wdir_to_continue argument should be used to provide\n'
                             'directories with files obtained during the step 1')
    parser1.add_argument('--wdir_to_continue', metavar='DIRNAME', required=False, default=None, nargs='+',
                         type=partial(filepath_type, exist_type='dir'),
                         help='''Single or multiple directories contain simulations created by the tool.
                            Use with steps 2,3,4 to continue run.\n'
                                    Should consist of: tpr, cpt, xtc and all_ligand_resid.txt files. 
                                    File all_ligand_resid.txt is optional and used to run md analysis for the ligands.\n
                                    If you want to continue your own simulation not created by the tool use --tpr, --cpt, --xtc and --wdir or arguments 
                                    (--ligand_list_file is optional and required to run md analysis after simulation )''')
    # continue md
    parser2 = parser.add_argument_group('Continue or Extend Molecular Dynamics Simulation')
    parser2.add_argument('--deffnm', metavar='preffix for md files', required=False, default='md_out',
                        help='''Preffix for the md files. Used to run, extend or continue the simulation.
                            If --wdir_to_continue is used files as deffnm.tpr, deffnm.cpt, deffnm.xtc will be searched from --wdir_to_continue directories''')
    parser2.add_argument('--tpr', metavar='FILENAME', required=False, default=None, type=filepath_type,
                        help='Use explicit tpr arguments to continue a non-StreaMD simulation')
    parser2.add_argument('--cpt', metavar='FILENAME', required=False, default=None, type=filepath_type,
                        help='Use explicit cpt arguments to continue a non-StreaMD simulation')
    parser2.add_argument('--xtc', metavar='FILENAME', required=False, default=None, type=filepath_type,
                        help='Use explicit xtc arguments to continue a non-StreaMD simulation')
    parser2.add_argument('--ligand_list_file', metavar='all_ligand_resid.txt', default=None, type=filepath_type,
                        help='''If you want automatic md analysis for ligands was run after continue of non-StreaMD simulation you should set ligand_list file. 
                                 Format of the file (no headers): user_ligand_id\tgromacs_ligand_id. Example: my_ligand\tUNL.
                                 Can be set up or placed into --wdir_to_continue directory(ies)''')
    parser2.add_argument('--ligand_id', metavar='UNL', default='UNL', type=str,
                        help='''If you want to run an automatic md analysis for the ligand after continue of simulation you can set ligand_id if it is not UNL default value''')
    # boron-containing molecules
    parser3 = parser.add_argument_group('Boron-containing molecules or MCPBPY usage (use together with Standard Molecular Dynamics Simulation Run arguments group)')
    # boron-containing molecules and mcpbpy
    parser3.add_argument('--activate_gaussian', metavar='module load Gaussian/09-d01', required=False, default=None,
                        help='String to load gaussian module if necessary')
    parser3.add_argument('--gaussian_exe', metavar='g09 or /apps/all/Gaussian/09-d01/g09/g09', required=False,
                        default=None,
                        help='Path to gaussian executable or alias. Required to run preparation of boron-containing compounds.')
    parser3.add_argument('--gaussian_basis', metavar='B3LYP/6-31G*', required=False,
                        default='B3LYP/6-31G*', help='Gaussian Basis')
    parser3.add_argument('--gaussian_memory', metavar='120GB', required=False,
                        default='120GB', help='Gaussian Memory Usage')
    # mcpbpy
    parser4 = parser.add_argument_group(
        'MCPBPY usage (Use together with Standard Molecular Dynamics Simulation Run and Boron-containing molecules arguments group)')
    parser4.add_argument('--metal_resnames', metavar='MN', required=False, default=None, nargs='*',
                        help='Metal residue names to run MCPB.py procedure. '
                             'Start MCPBPY procedure only if gaussian_exe and activate_gaussian arguments are set up,'
                             'Otherwise standard gmx2pdb procedure will be run.')
    parser4.add_argument('--metal_cutoff', metavar='2.8', required=False, default=2.8,
                        help='Metal residue cutoff to run MCPB.py procedure')
    parser4.add_argument('--metal_charges', metavar='{MN:2, ZN:2, CA:2}', type=json.loads, required=False,
                        default={'MN':2, 'ZN':2, 'CA':2},
                        help='Metal residue charges in dictionary format'
                             'Start MCPBPY procedure only if metal_resnames and gaussian_exe and activate_gaussian arguments are set up, '
                             'otherwise standard gmx2pdb procedure will be run')

    args = parser.parse_args()

    if args.wdir is None:
        wdir = os.getcwd()
    else:
        wdir = args.wdir

    if args.steps is not None:
        if not all([i in [1, 2, 3, 4] for i in args.steps]):
            raise ValueError(f'--steps {args.steps} argument is not valid. Please choose the combination from: 1, 2, 3, 4')
        if any([i in [2, 3, 4] for i in args.steps]) and args.wdir_to_continue is None:
            raise ValueError(f'--wdir_to_continue argument is not valid. '
                         f'If you set up --step {args.steps} you need to provide directories containing md files'
                         f' created by previous steps')

    out_time = f'{datetime.now().strftime("%d-%m-%Y-%H-%M-%S")}'
    log_file = os.path.join(wdir,
                            f'log_{os.path.basename(str(args.protein))[:-4]}_{os.path.basename(str(args.ligand))[:-4]}_{os.path.basename(str(args.cofactor))[:-4]}_'
                            f'{out_time}.log')
    bash_log = f'streamd_bash_{os.path.basename(str(args.protein))[:-4]}_{os.path.basename(str(args.ligand))[:-4]}_{os.path.basename(str(args.cofactor))[:-4]}_{out_time}.log'

    logging.basicConfig(format='%(asctime)s - %(name)s - %(levelname)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S',
                        level=logging.INFO,
                        handlers=[logging.FileHandler(log_file),
                                  logging.StreamHandler()])

    logging.getLogger('distributed').setLevel('WARNING')
    logging.getLogger('asyncssh').setLevel('WARNING')
    logging.getLogger('distributed.worker').setLevel('WARNING')
    logging.getLogger('distributed.core').setLevel('WARNING')
    logging.getLogger('distributed.comm').setLevel('WARNING')
    logging.getLogger('distributed.nanny').setLevel('CRITICAL')
    logging.getLogger('bockeh').setLevel('WARNING')

    logging.info(args)

    ncpu = min(max(0, args.ncpu), len(os.sched_getaffinity(0)))
    if ncpu != args.ncpu:
        logging.warning('The number of available CPUs are less than specified value. '
                        f'The tool will use only {ncpu} CPUs.')

    try:
        start(protein=args.protein,
              lfile=args.ligand, system_lfile=args.cofactor, noignh=args.noignh,
              topol=args.topol, topol_itp_list=args.topol_itp, posre_list_protein=args.posre,
              forcefield_name=args.protein_forcefield, npt_time_ps=args.npt_time,
              nvt_time_ps=args.nvt_time, mdtime_ns=args.md_time,
              wdir_to_continue_list=args.wdir_to_continue, deffnm=args.deffnm,
              tpr_prev=args.tpr, cpt_prev=args.cpt, xtc_prev=args.xtc,
              ligand_list_file_prev=args.ligand_list_file, ligand_resid=args.ligand_id,
              activate_gaussian=args.activate_gaussian, gaussian_exe=args.gaussian_exe,
              gaussian_basis=args.gaussian_basis, gaussian_memory=args.gaussian_memory,
              hostfile=args.hostfile, ncpu=ncpu, mdrun_per_node=args.mdrun_per_node, compute_device=args.device,
              gpu_ids=args.gpu_ids, ntmpi_per_gpu=args.ntmpi_per_gpu, wdir=wdir, seed=args.seed, steps=args.steps,
              clean_previous=args.clean_previous_md, not_clean_backup_files=args.not_clean_backup_files,
              metal_resnames=args.metal_resnames, metal_charges=args.metal_charges,
              mcpbpy_cut_off=args.metal_cutoff, out_time=out_time, bash_log=bash_log)
    finally:
        logging.shutdown()