#!/usr/bin/env python3
"""Utilities for running the StreaMD molecular dynamics pipeline.

This module orchestrates preparation, equilibration, simulation (and
continuation) and analysis steps for MD workflow.
"""
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
import sys
from datetime import datetime
from functools import partial
from glob import glob
import json
import re

from streamd.analysis.md_system_analysis import run_md_analysis
from streamd.analysis.run_analysis import run_rmsd_analysis
from streamd.preparation.complex_preparation import run_complex_preparation
from streamd.preparation.ligand_preparation import prepare_input_ligands, check_mols
from streamd.utils.dask_init import init_dask_cluster, calc_dask
from streamd.utils.utils import (
    filepath_type,
    run_check_subprocess,
    get_protein_resid_set,
    backup_prev_files,
    check_to_continue_simulation_time,
    merge_parts_of_simulation,
    parse_with_config,
)
from streamd.mcpbpy_md import mcbpy_md

logging.getLogger('distributed').setLevel('WARNING')
logging.getLogger('asyncssh').setLevel('WARNING')
logging.getLogger('MDAnalysis').setLevel('CRITICAL')
logging.getLogger('bockeh').setLevel('WARNING')


class RawTextArgumentDefaultsHelpFormatter(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass


def run_equilibration(wdir, project_dir, bash_log, ncpu, compute_device,
                      device_param, gpu_args, analysis_dirname='md_analysis', env=None):
    """Run the equilibration stage of the StreaMD pipeline.

    The function executes an external shell script to perform the
    equilibration step. If checkpoint files are already present, the step is
    skipped to avoid recomputation.

    :param wdir: Directory containing simulation files.
    :param project_dir: Base project directory with bundled scripts.
    :param bash_log: Name of the log file to capture shell output.
    :param ncpu: Number of CPU cores to use.
    :param compute_device: Target device ("cpu" or "gpu") for the run.
    :param device_param: Command line parameters for device selection.
    :param gpu_args: Additional GPU-related command line arguments.
    :param analysis_dirname: Name of the output analysis directory.
    :param env: Optional environment variables for the subprocess call.
    :return: Path to the working directory or ``None`` on failure.
    """
    if os.path.isfile(os.path.join(wdir, 'npt.gro')) and os.path.isfile(os.path.join(wdir, 'npt.cpt')):
        logging.warning(f'{wdir}. Checkpoint files after Equilibration step exist. '
                        f'Equilibration step will be skipped ')
        return wdir

    wdir_out_analysis = os.path.join(wdir, analysis_dirname)
    os.makedirs(wdir_out_analysis, exist_ok=True)
    system_name = os.path.split(wdir)[-1]

    cmd = (f'wdir={wdir} ncpu={ncpu} compute_device={compute_device} device_param={device_param} '
           f'gpu_args={gpu_args} wdir_out_analysis={wdir_out_analysis} system_name={system_name} '
           f'bash {os.path.join(project_dir, "scripts/script_sh/equlibration.sh")} '
           f'>> {os.path.join(wdir, bash_log)} 2>&1'),
    if not run_check_subprocess(cmd, wdir, log=os.path.join(wdir, bash_log), env=env):
        return None
    return wdir


def run_simulation(wdir, project_dir, bash_log, mdtime_ns,
                   tpr, cpt, xtc, deffnm, deffnm_next, ncpu,
                   compute_device, device_param, gpu_args, env=None):
    """Execute the main molecular dynamics simulation.

    The simulation is either started from scratch or continued from existing
    checkpoint files. All heavy lifting is delegated to an external shell
    script which runs the GROMACS ``mdrun`` command.

    :param wdir: Working directory of the simulation.
    :param project_dir: Base directory containing helper scripts.
    :param bash_log: Name of the log file to collect shell output.
    :param mdtime_ns: Desired simulation time in nanoseconds.
    :param tpr: Optional path to a ``.tpr`` file for continuation.
    :param cpt: Optional path to a ``.cpt`` checkpoint file.
    :param xtc: Optional trajectory file used when continuing a run.
    :param deffnm: Base name for simulation output files.
    :param deffnm_next: Base name for the next continuation step.
    :param ncpu: Number of CPU cores to use.
    :param compute_device: Target device ("cpu" or "gpu").
    :param device_param: Command line parameters describing device usage.
    :param gpu_args: Additional GPU specific arguments.
    :param env: Optional environment variables for the subprocess call.
    :return: Tuple with working directory and ``deffnm`` or ``None`` on failure.
    """
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
                         mdtime_ns, project_dir, bash_log, ncpu, compute_device,
                         device_param, gpu_args, env=None):
    """Continue a molecular dynamics run from an existing directory.

    This helper checks for required checkpoint files, validates simulation
    length, and merges any partial trajectories produced by previous
    continuation attempts. The actual continuation is executed through an
    external shell script.

    :param wdir_to_continue: Directory containing previous run data.
    :param tpr: Path to a ``.tpr`` file or ``None`` to use the default.
    :param cpt: Path to a checkpoint ``.cpt`` file or ``None``.
    :param xtc: Path to a trajectory ``.xtc`` file or ``None``.
    :param deffnm: Base name of current simulation files.
    :param deffnm_next: Base name for the next continuation step.
    :param mdtime_ns: Total desired simulation time in nanoseconds.
    :param project_dir: Base directory with helper scripts.
    :param bash_log: Log file capturing shell output.
    :param ncpu: Number of CPU cores available.
    :param compute_device: Target compute device.
    :param device_param: Device-specific command line parameters.
    :param gpu_args: GPU-specific command line parameters.
    :param env: Optional environment variables for the subprocess call.
    :return: Path to the working directory or ``None`` if continuation fails.
    """
    def continue_md(tpr, cpt, xtc, wdir, new_mdtime_ps, deffnm_next, project_dir, bash_log, compute_device, env):
        """Execute continuation step for a single run.

        :param tpr: Path to the ``.tpr`` file describing the system.
        :param cpt: Checkpoint ``.cpt`` file to restart from.
        :param xtc: Input trajectory fragment to append.
        :param wdir: Working directory where the command is executed.
        :param new_mdtime_ps: Target simulation time in picoseconds.
        :param deffnm_next: Base name for output of the next step.
        :param project_dir: Base directory with helper scripts.
        :param bash_log: Log file capturing shell output.
        :param compute_device: Selected compute device.
        :param env: Optional environment variables for the subprocess call.
        :return: Working directory if successful, otherwise ``None``.
        """
        cmd = f'wdir={wdir} tpr={tpr} cpt={cpt} xtc={xtc} new_mdtime_ps={new_mdtime_ps} ' \
              f'deffnm_next={deffnm_next} ncpu={ncpu} compute_device={compute_device} device_param={device_param} gpu_args={gpu_args} bash {os.path.join(project_dir, "scripts/script_sh/continue_md.sh")}' \
              f'>> {os.path.join(wdir, bash_log)} 2>&1'
        if run_check_subprocess(cmd, wdir, log=os.path.join(wdir, bash_log), env=env):
            return wdir
        return None

    if tpr is None:
        tpr = os.path.join(wdir_to_continue, f'{deffnm}.tpr')
    if cpt is None:
        cpt = os.path.join(wdir_to_continue, f'{deffnm}.cpt')
    if xtc is None:
        xtc = os.path.join(wdir_to_continue, f'{deffnm}.xtc')

    for i in [tpr, cpt, xtc]:
        if not os.path.isfile(i):
            logging.exception(
                f'No {i} file was found. Cannot continue the simulation. '
                f'Calculations will be interrupted.',
                stack_info=True)
            return None

    new_mdtime_ps = int(mdtime_ns * 1000)

    # check calculated time
    if not check_to_continue_simulation_time(xtc=xtc, new_mdtime_ps=new_mdtime_ps, env=env):
        return wdir_to_continue

    # check if can find unfinished or not merged continued trajectories {deffnm}_cont_
    found_already_continued_parts_simulations = glob(os.path.join(wdir_to_continue, f'{deffnm}_cont_*.xtc'))
    if found_already_continued_parts_simulations:
        if len(found_already_continued_parts_simulations) > 1:
            found_already_continued_parts_simulations = sorted(found_already_continued_parts_simulations,
                                                key=lambda f: datetime.strptime(os.path.basename(f).replace(f'{deffnm}_cont_', '').split('.')[0], "%d-%m-%Y-%H-%M-%S"))
            logging.warning(
                f'Found more than 1 continued part of the simulations: {found_already_continued_parts_simulations}.'
                f'They will be merged by the date sorted order. Trying to merge')

        else:
            logging.warning(
                f'Found a continued part of the simulation: {found_already_continued_parts_simulations[0]}.'
                f' Try to merge')

        for part_xtc in found_already_continued_parts_simulations:
            # md_out_cont.part0002.xtc and md_out_cont.tpr
            deffnm_part = os.path.basename(part_xtc).split('.part')[0]
            logging.info(deffnm_part)
            if re.findall(f'\.part[0-9]*\.xtc', part_xtc):
                # not merged part
                new_xtc = os.path.join(wdir_to_continue, f'{deffnm_part}.xtc')
                merge_parts_of_simulation(start_xtc=xtc,
                                           part_xtc=part_xtc,
                                           new_xtc=new_xtc,
                                           wdir=wdir_to_continue,
                                           bash_log=bash_log,
                                           env=env)
                backup_prev_files(file_to_backup=part_xtc, wdir=wdir_to_continue)
            # backup old part files - xtc, log, tpr, cpt, edr
            for cont_sim_file in glob(os.path.join(wdir_to_continue, f'{deffnm_part}*')):
                ext = os.path.splitext(cont_sim_file)[1]
                # simulation start files
                if ext == '.xtc':
                    user_file_to_replace = xtc
                elif ext == '.cpt':
                    user_file_to_replace = cpt
                elif ext == '.tpr':
                    user_file_to_replace = tpr
                else:
                    user_file_to_replace = os.path.join(wdir_to_continue, f'{deffnm}{ext}')

                if os.path.isfile(user_file_to_replace):
                    backup_prev_files(file_to_backup=user_file_to_replace, wdir=wdir_to_continue)
                    # replace start simulation file with new merged one
                    backup_prev_files(file_to_backup=cont_sim_file, wdir=wdir_to_continue, copy=True)
                    shutil.move(cont_sim_file, user_file_to_replace)

        # check new merged trajectory time
        logging.warning(f'{xtc} and {found_already_continued_parts_simulations} were merged successfully.')
        if not check_to_continue_simulation_time(xtc=xtc, new_mdtime_ps=new_mdtime_ps, env=env):
            return wdir_to_continue

    if continue_md(tpr=tpr, cpt=cpt, xtc=xtc, wdir=wdir_to_continue,
                   new_mdtime_ps=new_mdtime_ps, deffnm_next=deffnm_next, project_dir=project_dir,
                   compute_device=compute_device, env=env, bash_log=bash_log):
        logging.warning('Simulation extension completed successfully')
        # backup cont part files
        for f in glob(os.path.join(wdir_to_continue, f'{deffnm_next}.part*.*')):
            backup_prev_files(file_to_backup=f, wdir=wdir_to_continue)
        #
        for f in glob(os.path.join(wdir_to_continue, f'{deffnm_next}.*')):
            # check previous existing files with the same name
            if os.path.isfile(os.path.join(wdir_to_continue, os.path.basename(f).replace(deffnm_next, deffnm))):
                backup_prev_files(file_to_backup=os.path.join(wdir_to_continue, os.path.basename(f).replace(deffnm_next, deffnm)),
                              wdir=wdir_to_continue)
                backup_prev_files(file_to_backup=f, wdir=wdir_to_continue, copy=True)
            # replace old file with continued one
                shutil.move(f, os.path.join(wdir_to_continue, os.path.basename(f).replace(deffnm_next, deffnm)))

        return wdir_to_continue


def start(protein, wdir, lfile, system_lfile, noignh, no_dr,
          forcefield_name, npt_time_ps, nvt_time_ps, mdtime_ns,
          topol, topol_itp_list, posre_list_protein,
          wdir_to_continue_list, deffnm,
          tpr_prev, cpt_prev, xtc_prev, ligand_list_file_prev, ligand_resid,
          activate_gaussian, gaussian_exe, gaussian_basis, gaussian_memory,
          metal_resnames, metal_charges, mcpbpy_cut_off,
          seed, steps, hostfile, ncpu, mdrun_per_node, compute_device, gpu_ids, ntmpi_per_gpu, clean_previous,
          not_clean_backup_files, unique_id,
          active_site_dist=5.0, save_traj_without_water=False,
          explicit_args=(), mdp_dir=None, bash_log=None):
    """Run StreaMD pipeline.

    :param protein: Protein file in PDB or GRO format.
    :param wdir: Working directory or ``None`` for current path.
    :param lfile: Optional ligand file.
    :param system_lfile: Optional cofactor file in MOL or SDF format.
    :param noignh: Do not use ``-ignh`` argument in ``gmx pdb2gmx``.
    :param no_dr: Disable acdoctor mode when processing ligand file.
    :param forcefield_name: Force field to use for preparation.
    :param clean_previous: Whether to remove previous MD files.
    :param mdtime_ns: Simulation time in nanoseconds.
    :param npt_time_ps: Equilibration NPT time in picoseconds.
    :param nvt_time_ps: Equilibration NVT time in picoseconds.
    :param topol: Optional topology file.
    :param topol_itp_list: Optional list of ITP files for protein chains.
    :param posre_list_protein: Optional list of position restraint files.
    :param wdir_to_continue_list: List of directories with runs to continue.
    :param deffnm: Base name for simulation output files.
    :param tpr_prev: Optional ``.tpr`` file for continuation.
    :param cpt_prev: Optional checkpoint ``.cpt`` file.
    :param xtc_prev: Optional trajectory ``.xtc`` file.
    :param ligand_resid: Ligand residue name used for analysis on continuation.
    :param ligand_list_file_prev: Optional previous ligand list file.
    :param activate_gaussian: Command used to activate the Gaussian environment.
    :param gaussian_exe: Path to the Gaussian executable.
    :param gaussian_basis: Gaussian basis set name.
    :param gaussian_memory: Amount of memory to allocate for Gaussian jobs.
    :param metal_resnames: List of metal residue names to process.
    :param metal_charges: Mapping of metal residue names to their charges.
    :param mcpbpy_cut_off: Cut-off distance for MCPB.py metal center processing.
    :param seed: Random seed for stochastic steps.
    :param steps: Tuple selecting specific pipeline stages to run.
    :param hostfile: Optional file with Dask host addresses.
    :param ncpu: Number of CPU cores.
    :param mdrun_per_node: Number of ``mdrun`` instances to launch per node.
    :param compute_device: Compute device to use.
    :param gpu_ids: Optional list of GPU IDs.
    :param ntmpi_per_gpu: Number of thread-MPI ranks per GPU.
    :param unique_id: Unique identifier for run artifacts.
    :param bash_log: Name of log file capturing shell output.
    :param mdp_dir: Directory containing MDP files.
    :param explicit_args: Tuple of additional command-line arguments.
    :param not_clean_backup_files: Keep backup files instead of removing.
    :param active_site_dist: Distance threshold in Å for active-site analysis.
    :param save_traj_without_water: Whether to output trajectories without water.
    :return: ``None``.
    """
    project_dir = os.path.dirname(os.path.abspath(__file__))
    script_path = os.path.join(project_dir, 'scripts')
    script_mdp_path = os.path.join(script_path, 'mdp')

    wdir_md = os.path.join(wdir, 'md_files', 'md_run')
    analysis_dirname = 'md_analysis'

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
                        if not run_check_subprocess(cmd, protein, log=os.path.join(wdir_protein, bash_log), env=os.environ.copy()):
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
                                                         project_dir=project_dir, wdir_ligand=wdir_system_ligand, no_dr=no_dr,
                                                         gaussian_exe=gaussian_exe, activate_gaussian=activate_gaussian,
                                                         gaussian_basis=gaussian_basis, gaussian_memory=gaussian_memory,
                                                         hostfile=hostfile, ncpu=ncpu, bash_log=bash_log)
                if number_of_mols != len(system_lig_wdirs):
                    logging.exception(f'Error with the cofactor preparation. Only {len(system_lig_wdirs)} from {number_of_mols} preparation were finished.'
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
                                                      project_dir=project_dir, wdir_ligand=wdir_ligand, no_dr=no_dr,
                                                      gaussian_exe=gaussian_exe, activate_gaussian=activate_gaussian,
                                                      gaussian_basis=gaussian_basis, gaussian_memory=gaussian_memory,
                                                      hostfile=hostfile, ncpu=ncpu, bash_log=bash_log)
                if number_of_mols != len(var_lig_wdirs):
                    logging.warning(f'Problem with the ligand preparation. Only {len(var_lig_wdirs)} from {number_of_mols} preparation were finished.'
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
                                         npt_time_ps=npt_time_ps, nvt_time_ps=nvt_time_ps,
                                         mdp_dir=mdp_dir, explicit_args=explicit_args,
                                         bash_log=bash_log, seed=seed, env=os.environ.copy()):
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
                                         analysis_dirname=analysis_dirname,
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
                                         deffnm=deffnm, deffnm_next=f'{deffnm}_cont_{unique_id}',
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
        logging.info('Start Analysis of the simulations')
        var_md_analysis_res = []
        try:
            dask_client, cluster = init_dask_cluster(hostfile=hostfile,
                                                     n_tasks_per_node=min(ncpu, len(var_md_dirs_deffnm)),
                                                     use_multi_servers=True if len(var_md_dirs_deffnm) > ncpu else False,
                                                     ncpu=ncpu)
            for res in calc_dask(run_md_analysis, var_md_dirs_deffnm,
                                 dask_client, mdtime_ns=mdtime_ns, project_dir=project_dir,
                                 bash_log=bash_log, ligand_resid=ligand_resid,
                                 ligand_list_file_prev=ligand_list_file_prev,
                                 save_traj_without_water=save_traj_without_water,
                                 analysis_dirname=analysis_dirname,
                                 env=os.environ.copy()):
                if res:
                    # (rmsd_out_file, md_analysis_dir, md_cur_wdir)
                    var_md_analysis_res.append(res)
        finally:
            if dask_client:
                dask_client.retire_workers(dask_client.scheduler_info()['workers'],
                                           close_workers=True, remove=True)
                dask_client.shutdown()
            if cluster:
                cluster.close()

        rmsd_files = [i[0] for i in var_md_analysis_res]
        md_dirs_analyzed = [i[2] for i in var_md_analysis_res]

        rmsd_type_list = ['backbone', 'ligand', f'ActiveSite{active_site_dist}A'] if lfile else ['backbone']
        run_rmsd_analysis(rmsd_files=rmsd_files, wdir=wdir, unique_id=unique_id,
                          time_ranges=None,
                          rmsd_type_list=rmsd_type_list)

        finished_complexes_file = os.path.join(wdir, f"finished_complexes_{unique_id}.txt")
        with open(finished_complexes_file, 'w') as output:
            output.write("\n".join([i for i in md_dirs_analyzed]))

        logging.info(
            f'Analysis of MD simulations of {len(md_dirs_analyzed)} complexes were successfully finished'
            f'\nSuccessfully finished complexes have been saved in {finished_complexes_file} file')

    if not not_clean_backup_files:
        if wdir_to_continue_list is None:
            for f in glob(os.path.join(wdir_md, '*', '#*#')):
                # if '.tpr.' not in f and '.xtc.' not in f:
                #     os.remove(f)
                os.remove(f)
            for f in glob(os.path.join(wdir_md, '*', analysis_dirname, '#*#')):
                os.remove(f)
            for f in glob(os.path.join(wdir_md, '*', '*.trr')):
                os.remove(f)
        else:
            for wdir_md in wdir_to_continue_list:
                for f in glob(os.path.join(wdir_md, '#*#')):
                    os.remove(f)
                for f in glob(os.path.join(wdir_md, analysis_dirname, '#*#')):
                    os.remove(f)

def main():
    parser = argparse.ArgumentParser(description='''Run or continue MD simulation.\n
    Allowed systems: Protein, Protein-Ligand, Protein-Cofactors(multiple), Protein-Ligand-Cofactors(multiple) ''')
    parser.add_argument('--config', metavar='FILENAME', required=False,
                        type=partial(filepath_type, ext=("yml", "yaml")),
                        help='Path to YAML configuration file with default arguments')
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
                        help='Time of MD simulation in ns. Default: 1 ns.')
    parser1.add_argument('--npt_time', metavar='ps', required=False, default=1000, type=int,
                        help='Time of NPT equilibration in ps. Default: 1000 ps.')
    parser1.add_argument('--nvt_time', metavar='ps', required=False, default=1000, type=int,
                        help='Time of NVT equilibration in ps. Default: 1000 ps.')
    parser1.add_argument('--seed', metavar='int', required=False, default=-1, type=int,
                        help='seed')
    parser1.add_argument('--no_dr', action='store_true', default=False,
                         help='Turn off the acdoctor mode and do not check/diagnose problems in the input ligand file '
                              'in the next attempt if the regular antechamber run for ligand preparation fails (ligand_mol2prep.sh script related issues). '
                              'Use this argument carefully and ensure that you provide valid structures')
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
    parser1.add_argument('--mdp_dir', metavar='Path to a directory with specific mdp files', required=False,
                         default=None, type=partial(filepath_type, exist_type='dir'),
                         help='To use specific MD settings, the user can provide a path to a directory '
                              'that contains any of the following .mdp files: '
                              'ions.mdp, minim.mdp, nvt.mdp, npt.mdp, md.mdp. '
                              'Missing .mdp files will be replaced by default StreaMD files. '
                              'Provided .mdp files will be used as templates, '
                              'although the system StreaMD parameters '
                              '(seed, nvt_time, npt_time, md_time, and tc-grps (can not be changed by user)) will override the ones provided. '
                              'Warning: The names of the files must be strictly preserved.')
    parser1.add_argument('--save_traj_without_water', action='store_true', default=False,
                         help='Save additional md_out_nowater.tpr and md_fit_nowater.xtc files '
                              'for more memory efficient analysis.')
    parser1.add_argument('--wdir_to_continue', metavar='DIRNAME', required=False, default=None, nargs='+',
                         type=partial(filepath_type, exist_type='dir'),
                         help='''Single or multiple directories contain simulations created by the tool.
                            Use with steps 2,3,4 to continue run.\n'
                                    Should consist of: tpr, cpt, xtc and all_ligand_resid.txt files. 
                                    File all_ligand_resid.txt is optional and used to run md analysis for the ligands.\n
                                    If you want to continue your own simulation not created by the tool use --tpr, --cpt, --xtc and --wdir or arguments 
                                    (--ligand_list_file is optional and required to run md analysis after simulation )''')
    parser.add_argument('-o','--out_suffix', default=None,
                        help='User unique suffix for output files')
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

    args, config_args = parse_with_config(parser, sys.argv[1:])

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

    # check if user explicitly specified md, nvt or npt time and if so change user-defined mdp files
    explicit_args = []
    if '--nvt_time' in sys.argv or 'nvt_time' in config_args:
        explicit_args.append('nvt.mdp')
    if '--npt_time' in sys.argv or 'npt_time' in config_args:
        explicit_args.append('npt.mdp')
    if '--md_time' in sys.argv or 'md_time' in config_args:
        explicit_args.append('md.mdp')
    if '--seed' in sys.argv or 'seed' in config_args:
        explicit_args.append('seed')

    out_time = f'{datetime.now().strftime("%d-%m-%Y-%H-%M-%S")}'

    if args.out_suffix:
        unique_id = args.out_suffix
    else:
        unique_id = out_time

    log_file = os.path.join(wdir,
                            f'log_{os.path.basename(str(args.protein))[:-4]}_{os.path.basename(str(args.ligand))[:-4]}_{os.path.basename(str(args.cofactor))[:-4]}_'
                            f'{unique_id}.log')
    bash_log = f'streamd_bash_{os.path.basename(str(args.protein))[:-4]}_{os.path.basename(str(args.ligand))[:-4]}_{os.path.basename(str(args.cofactor))[:-4]}_{unique_id}.log'

    logging.basicConfig(format='%(asctime)s  - %(name)s - %(levelname)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S',
                        level=logging.INFO,
                        handlers=[logging.FileHandler(log_file),
                                  logging.StreamHandler()])

    logging.info(args)

    ncpu = min(max(0, args.ncpu), len(os.sched_getaffinity(0)))
    if ncpu != args.ncpu:
        logging.warning('The number of available CPUs are less than specified value. '
                        f'The tool will use only {ncpu} CPUs.')

    try:
        start(protein=args.protein,
              lfile=args.ligand, system_lfile=args.cofactor, noignh=args.noignh, no_dr=args.no_dr,
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
              mcpbpy_cut_off=args.metal_cutoff, unique_id=unique_id,
              save_traj_without_water=args.save_traj_without_water,
              mdp_dir=args.mdp_dir,
              explicit_args=explicit_args,
              bash_log=bash_log)
    finally:
        logging.shutdown()