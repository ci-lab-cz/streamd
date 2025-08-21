#!/usr/bin/env python3
# ==============================================================================
# author          : Aleksandra Ivanova, Olena Mokshyna, Pavel Polishchuk
# date            : 14-06-2024
# version         :
# python_version  :
# copyright       :
# license         : MIT
# ==============================================================================
"""Run MM-GBSA/MM-PBSA calculations using the gmx_MMPBSA tool."""

import argparse
import logging
import math
import os
import pathlib
import re
import shutil
import subprocess
from datetime import datetime
from functools import partial
from multiprocessing import Pool
import sys

import pandas as pd

from streamd.utils.dask_init import init_dask_cluster, calc_dask
from streamd.utils.utils import (
    get_index,
    make_group_ndx,
    filepath_type,
    run_check_subprocess,
    get_number_of_frames,
    temporary_directory_debug,
    parse_with_config,
)

logging.getLogger('distributed').setLevel('CRITICAL')
logging.getLogger('asyncssh').setLevel('CRITICAL')
logging.getLogger('distributed.worker').setLevel('CRITICAL')
logging.getLogger('distributed.core').setLevel('CRITICAL')
logging.getLogger('distributed.comm').setLevel('CRITICAL')
logging.getLogger('distributed.nanny').setLevel('CRITICAL')
logging.getLogger('bockeh').setLevel('CRITICAL')


def run_gbsa_task(wdir, tpr, xtc, topol, index, mmpbsa, np, ligand_resid, append_protein_selection,
                  unique_id, env, bash_log, clean_previous, debug, decomp=False):
    """Run a GBSA calculation for a prepared simulation directory.

    :param wdir: Working directory containing the simulation files.
    :param tpr: Path to the ``.tpr`` file.
    :param xtc: Path to the trajectory ``.xtc`` file.
    :param topol: Path to the topology ``.top`` file.
    :param index: Path to the GROMACS index file.
    :param mmpbsa: gmx_MMPBSA input file.
    :param np: Number of MPI processes to use.
    :param ligand_resid: Residue name of the ligand in the index file.
    :param append_protein_selection: Additional selections to merge with the protein group.
    :param unique_id: Identifier appended to output files.
    :param env: Optional environment variables for subprocess calls.
    :param bash_log: Name of the log file capturing shell output.
    :param clean_previous: Whether to remove previous GBSA outputs.
    :param debug: If ``True``, keep intermediate files for debugging.
    :param decomp: If ``True``, perform per-residue decomposition and save additional files.
    :return: Path to the generated GBSA results file or ``None`` on failure.
    """

    def calc_gbsa(wdir, tpr, xtc, topol, index, mmpbsa, np, protein_index,
                  ligand_index, unique_id, env, bash_log, debug, decomp=False):
        """Launch gmx_MMPBSA for a single trajectory.

        :param wdir: Working directory for temporary files and logs.
        :param tpr: Path to the ``.tpr`` file.
        :param xtc: Path to the trajectory ``.xtc`` file.
        :param topol: Topology file path.
        :param index: Index file path.
        :param mmpbsa: gmx_MMPBSA input file.
        :param np: Number of MPI processes.
        :param protein_index: Index of the protein selection.
        :param ligand_index: Index of the ligand selection.
        :param unique_id: Identifier appended to output files.
        :param env: Optional environment variables for subprocess calls.
        :param bash_log: Log file capturing shell output.
        :param debug: Keep intermediate files if ``True``.
        :param decomp: If ``True``, perform per-residue decomposition and save additional files.
        :return: Path to the generated results file or ``None`` on failure.
        """
        output_file = os.path.join(wdir, f"FINAL_RESULTS_MMPBSA_{unique_id}.dat")
        output_decomp = None
        with temporary_directory_debug(dir=wdir, remove=not debug, suffix=f'_gbsa_{unique_id}') as tmpdirname:
            logging.info(f'tmp intermediate dir: {tmpdirname}')
            cmd = (
                f'cd {tmpdirname}; mpirun -np {np} gmx_MMPBSA MPI -O -i {mmpbsa} '
                f' -cs {tpr} -ci {index} -cg {protein_index} {ligand_index} -ct {xtc} -cp {topol} -nogui '
                f'-o {output_file} '
                f'-eo {os.path.join(wdir, f"FINAL_RESULTS_MMPBSA_{unique_id}.csv")} '
            )
            if decomp:
                output_decomp = os.path.join(wdir, f"FINAL_DECOMP_MMPBSA_{unique_id}.dat")
                logging.info('Run Decomposition Analysis')
                cmd += (
                    f'-do {output_decomp} '
                    f'-deo {os.path.join(wdir, f"FINAL_DECOMP_MMPBSA_{unique_id}.csv")} '
                )
            cmd += f'>> {os.path.join(wdir, bash_log)} 2>&1'

            if not run_check_subprocess(cmd, key=xtc, log=os.path.join(wdir, bash_log), env=env):
                run_check_subprocess(f'ls {tmpdirname}', key=tmpdirname, log=os.path.join(wdir, bash_log), env=env)
                return None

        output = {'out': output_file,
                  'out_decomp': output_decomp}
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
                logging.warning(
                    f'{wdir} Could not find resname {i}. It will not be used in gbsa calculation. Check your query carefully.')
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

    output = calc_gbsa(
        wdir=wdir,
        tpr=tpr,
        xtc=xtc,
        topol=topol,
        index=index,
        mmpbsa=mmpbsa,
        np=np,
        protein_index=protein_index,
        ligand_index=ligand_index,
        unique_id=unique_id,
        env=env,
        bash_log=bash_log,
        debug=debug,
        decomp=decomp,
    )

    if os.path.isfile(os.path.join(wdir, 'gmx_MMPBSA.log')):
        shutil.copy(os.path.join(wdir, 'gmx_MMPBSA.log'), os.path.join(wdir, f'gmx_MMPBSA_{unique_id}.log'))

    return output


def run_gbsa_from_wdir(wdir, tpr, xtc, topol, index, mmpbsa, np, ligand_resid,
                       append_protein_selection, unique_id, env, bash_log, clean_previous, debug, decomp=False):
    """Execute GBSA using file paths relative to the working directory.

    :param wdir: Base working directory containing simulation files.
    :param tpr: Relative path to the ``.tpr`` file inside ``wdir``.
    :param xtc: Relative path to the ``.xtc`` trajectory file inside ``wdir``.
    :param topol: Relative path to the topology file inside ``wdir``.
    :param index: Relative path to the index file inside ``wdir``.
    :param mmpbsa: gmx_MMPBSA input file path.
    :param np: Number of MPI processes to use.
    :param ligand_resid: Residue name of the ligand in the index file.
    :param append_protein_selection: Additional selections to merge with the protein group.
    :param unique_id: Identifier appended to output files.
    :param env: Optional environment variables for subprocess calls.
    :param bash_log: Name of the log file capturing shell output.
    :param clean_previous: Whether to remove previous GBSA outputs.
    :param debug: If ``True``, keep intermediate files for debugging.
    :param decomp: If ``True``, perform per-residue decomposition and save additional files.
    :return: Path to the generated GBSA results file or ``None`` on failure.
    """
    tpr = os.path.join(wdir, tpr)
    xtc = os.path.join(wdir, xtc)
    topol = os.path.join(wdir, topol)
    index = os.path.join(wdir, index)
    return run_gbsa_task(
        wdir=wdir,
        tpr=tpr,
        xtc=xtc,
        topol=topol,
        index=index,
        mmpbsa=mmpbsa,
        np=np,
        ligand_resid=ligand_resid,
        append_protein_selection=append_protein_selection,
        unique_id=unique_id,
        env=env,
        bash_log=bash_log,
        clean_previous=clean_previous,
        debug=debug,
        decomp=decomp,
    )


def clean_temporary_gmxMMBPSA_files(wdir, prefix="_GMXMMPBSA_"):
    """Remove intermediate files produced by gmx_MMPBSA.

    :param wdir: Working directory in which to run the clean command.
    :param prefix: Prefix used by gmx_MMPBSA to tag temporary files.
    :return: ``None``.
    """
    try:
        subprocess.check_output(f'cd {wdir}; gmx_MMPBSA --clean -prefix "{prefix}" > /dev/null 2>&1 ', shell=True)
    except subprocess.CalledProcessError as e:
        logging.error(f'{wdir}\t{e}\n')
        return None


def parse_gmxMMPBSA_output(fname):
    """Parse energies from a gmx_MMPBSA result file.

    :param fname: Path to the ``FINAL_RESULTS_MMPBSA`` output file.
    :return: Nested dictionaries with GBSA and PBSA energy terms.
    """

    def get_IE_values(IE_parsed_out):
        """Map interaction energy fields to values.

        :param IE_parsed_out: Parsed regex groups for interaction energies.
        :return: Dictionary with ``IE_`` prefixed energy components.
        """
        IE_res = {}
        IE_columns = [i.strip() for i in IE_parsed_out[0][0].split('  ') if i]
        IE_values = [i.strip() for i in IE_parsed_out[0][1].split('  ') if i]
        for n, i in enumerate(IE_columns):
            IE_res[f'IE_{i}'] = IE_values[n]
        return IE_res

    def get_delta_total_values(delta_total_columns_re, delta_total_values_re):
        """Collect total energy terms from regex matches.

        :param delta_total_columns_re: Regex match with column names.
        :param delta_total_values_re: Regex match with column values.
        :return: Dictionary with ``ΔTOTAL_`` prefixed energy components.
        """
        delta_total_res = {}
        delta_total_columns = [i.strip() for i in delta_total_columns_re[0].split('  ') if i]
        delta_total_values = [i.strip() for i in delta_total_values_re[0].split('  ') if i]
        # ['Energy Component' - skip, 'Average', 'SD(Prop.)', 'SD', 'SEM(Prop.)', 'SEM']

        for n, i in enumerate(delta_total_columns[1:]):
            delta_total_res[f'ΔTOTAL_{i}'] = delta_total_values[n]
        return delta_total_res

    def get_Gbinding_values(Gbind_parsed_out):
        """Parse binding free energy and its error.

        :param Gbind_parsed_out: Regex match with binding energy values.
        :return: Dictionary with ``ΔGbinding`` and its uncertainty.
        """
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


def parse_gmxMMPBSA_decomp_dat(fname):
    """Parse averaged residue decomposition from gmx_MMPBSA ``.dat`` files.

    ``gmx_MMPBSA`` produces ``FINAL_DECOMP_MMPBSA`` text files containing
    averaged per-residue energy contributions for the complex, receptor, ligand
    and their deltas.  Each section also distinguishes between total, sidechain
    and backbone contributions and may report values for both Generalized Born
    and Poisson Boltzmann models.  This helper converts such a file into a tidy
    :class:`pandas.DataFrame` with columns combining the energy term and
    statistical measure, e.g. ``Internal Avg.`` or ``Electrostatic Std. Dev.``.

    Parameters
    ----------
    fname : str | os.PathLike
        Path to ``FINAL_DECOMP_MMPBSA`` ``.dat`` output.

    Returns
    -------
    pandas.DataFrame
        Parsed residue contributions.  If the file cannot be parsed an empty
        DataFrame is returned.
    """

    sections = {"Complex", "Receptor", "Ligand", "DELTAS"}
    contributions = {
        "Total Energy Decomposition": "Total",
        "Sidechain Energy Decomposition": "Sidechain",
        "Backbone Energy Decomposition": "Backbone",
    }

    data = []
    region = None
    contrib = None
    header_top = None
    header_sub = None
    method = None

    with open(fname) as fh:
        for line in fh:
            stripped = line.strip()
            if not stripped:
                continue

            if stripped.startswith("Energy Decomposition Analysis"):
                if "Generalized Born" in stripped:
                    method = "GB"
                elif "Poisson Boltzmann" in stripped:
                    method = "PB"
                region = None
                contrib = None
                header_top = None
                header_sub = None
                continue

            if stripped.endswith(":") and stripped[:-1] in sections:
                region = stripped[:-1]
                contrib = None
                header_top = None
                header_sub = None
                continue

            key = stripped.rstrip(":")
            if key in contributions:
                contrib = contributions[key]
                header_top = None
                header_sub = None
                continue

            if contrib and header_top is None:
                header_top = [h.strip() for h in stripped.split(",")]
                continue

            if contrib and header_sub is None:
                header_sub = [h.strip() for h in stripped.split(",")]
                header = []
                current = None
                for top, sub in zip(header_top, header_sub):
                    if top:
                        current = top
                    if current == "Residue":
                        header.append("Residue")
                        continue
                    name = current
                    if sub:
                        name = f"{current} {sub}"
                    header.append(name)
                continue

            if header_sub and region and contrib and method:
                parts = [p.strip() for p in stripped.split(",")]
                if len(parts) == len(header):
                    row = dict(zip(header, parts))
                    row["Region"] = region
                    row["Contribution"] = contrib
                    row["Method"] = method
                    data.append(row)

    df = pd.DataFrame(data)
    if df.empty:
        return df

    numeric_cols = [
        c
        for c in df.columns
        if c not in {"Residue", "Region", "Contribution", "Method"}
    ]
    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")


    return df


def run_get_frames_from_wdir(wdir, xtc, env):
    """Return number of trajectory frames for a working directory.

    :param wdir: Directory containing the trajectory file.
    :param xtc: Name of the trajectory ``.xtc`` file relative to ``wdir``.
    :param env: Optional environment variables for subprocess calls.
    :return: Tuple of frame count and timestep or ``None``.
    """
    return get_number_of_frames(os.path.join(wdir, xtc), env=env)


def get_mmpbsa_params(mmpbsa):
    """Extract start, end and interval parameters from an MMPBSA input file.

    :param mmpbsa: Path to the gmx_MMPBSA input file.
    :return: Tuple of ``(startframe, endframe, interval, decomp: True/False)``.
    """
    decomp = False

    with open(mmpbsa) as inp:
        mmpbsa_data = inp.read()
    logging.info(f'{mmpbsa}:\n{mmpbsa_data}')

    startframe, endframe, interval = None, None, None
    for line in mmpbsa_data.split('\n'):
        if line.startswith('#'):
            continue
        line = line.strip().lower()
        if 'startframe' in line:
            startframe = re.findall('startframe[ ]*=[ ]*([0-9]*)', line)
        if 'endframe' in line:
            endframe = re.findall('endframe[ ]*=[ ]*([0-9]*)', line)
        if 'interval' in line:
            interval = re.findall('interval[ ]*=[ ]*([0-9]*)', line)
        if line.startswith('&decomp'):
            decomp = True
            logging.info(
                "Detected '&decomp' section in %s; enabling per-residue decomposition.",
                mmpbsa)

    startframe = int(startframe[0]) if startframe else 1
    endframe = int(endframe[0]) if endframe else 9999999  # default value of gmxMMBPSA
    interval = int(interval[0]) if interval else 1

    return startframe, endframe, interval, decomp


def get_used_number_of_frames(var_number_of_frames, startframe, endframe, interval):
    """Calculate number of frames used based on limits and stride.

    :param var_number_of_frames: Number of frames present in the trajectory.
    :param startframe: Starting frame for analysis (1-indexed).
    :param endframe: Last frame to include in analysis.
    :param interval: Stride between frames.
    :return: Number of frames that will be used.
    """
    return math.ceil((min(min(var_number_of_frames), endframe) - (startframe - 1)) / interval)


def start(wdir_to_run, tpr, xtc, topol, index, out_wdir, mmpbsa_file, ncpu, ligand_resid,
          append_protein_selection, hostfile, unique_id, bash_log,
          gmxmmpbsa_dat_files=None, gmxmmpbsa_decomp_dat_files=None,
          clean_previous=False, debug=False):
    """Start GBSA calculations and aggregate resulting energies.

    Parameters
    ----------
    wdir_to_run : list[str] | None
        Working directories containing trajectories; overrides individual file arguments.
    tpr, xtc, topol, index : str | None
        Paths to required simulation files when ``wdir_to_run`` is not used.
    out_wdir : str
        Directory to write output tables.
    mmpbsa : str
        Input file for gmx_MMPBSA.
    ncpu : int
        Number of CPU cores to use.
    ligand_resid : str
        Residue identifier for the ligand.
    append_protein_selection : str | None
        Extra selection appended to protein group.
    hostfile : str | None
        Hostfile for distributed runs.
    unique_id : str
        Unique identifier used in output filenames.
    bash_log : str
        Name of log file capturing gmx_MMPBSA output.
    gmxmmpbsa_dat_files : list[str] | None
        Precomputed gmx_MMPBSA result files to parse instead of running calculations.
    clean_previous : bool
        Remove previous outputs before running.
    debug : bool
        Enable verbose debugging output.
    decomp : bool
        Perform per-residue decomposition analysis and save additional files.

    Returns
    -------
    list[str] | None
        Paths to generated GBSA output files or ``None`` on failure.
    """
    dask_client, cluster, pool = None, None, None
    var_gbsa_dat_files = []
    decomp_dat_files = []

    # Automatically enable decomposition if the mmpbsa input defines a ``&decomp`` block.

    if gmxmmpbsa_dat_files is None and gmxmmpbsa_decomp_dat_files is None:
        # gmx_mmpbsa requires that the run must have at least as many frames as processors. Thus we get and use the min number of used frames as NP
        if not mmpbsa_file:
            mmpbsa_file = os.path.join(out_wdir, f'mmpbsa_{unique_id}.in')
            project_dir = os.path.dirname(os.path.abspath(__file__))
            shutil.copy(os.path.join(project_dir, 'scripts', 'gbsa', 'mmpbsa.in'), mmpbsa_file)
            logging.warning(f'No mmpbsa.in file was provided. Template will be used. Created file: {mmpbsa_file}.')

        startframe, endframe, interval, decomp = get_mmpbsa_params(mmpbsa_file)

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
                var_gbsa_dat_files = []
                for res in calc_dask(run_gbsa_from_wdir, wdir_to_run, dask_client=dask_client,
                                     tpr=tpr, xtc=xtc, topol=topol, index=index,
                                     mmpbsa=mmpbsa_file, np=min(ncpu, used_number_of_frames),
                                     ligand_resid=ligand_resid,
                                     append_protein_selection=append_protein_selection,
                                     unique_id=unique_id, env=os.environ.copy(),
                                     bash_log=bash_log, clean_previous=clean_previous,
                                     debug=debug, decomp=decomp):
                    if res:
                        if res.get('out', None) is not None:
                            var_gbsa_dat_files.append(res['out'])
                        if res.get('out_decomp', None) is not None:
                            decomp_dat_files.append(res['out_decomp'])
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
            res = run_gbsa_task(wdir=os.path.dirname(xtc), tpr=tpr, xtc=xtc, topol=topol, index=index,
                                mmpbsa=mmpbsa_file, np=min(ncpu, used_number_of_frames),
                                ligand_resid=ligand_resid,
                                append_protein_selection=append_protein_selection,
                                unique_id=unique_id, env=os.environ.copy(), bash_log=bash_log,
                                clean_previous=clean_previous, debug=debug, decomp=decomp)
            if res:
                if res.get('out', None):
                    var_gbsa_dat_files.append(res['out'])
                if res.get('out_decomp', None):
                    decomp_dat_files.append(res['out_decomp'])

    else:
        print('work')
        if gmxmmpbsa_dat_files is not None:
            var_gbsa_dat_files = gmxmmpbsa_dat_files
        if gmxmmpbsa_decomp_dat_files is not None:
            decomp_dat_files = gmxmmpbsa_decomp_dat_files

    # collect energies
    if var_gbsa_dat_files:
        GBSA_output_res, PBSA_output_res = [], []
        with Pool(ncpu) as pool:
            for res in pool.imap_unordered(parse_gmxMMPBSA_output, var_gbsa_dat_files):
                if res:
                    if res['GBSA']:
                        GBSA_output_res.append(res['GBSA'])
                    if res['PBSA']:
                        PBSA_output_res.append(res['PBSA'])

        if GBSA_output_res:
            pd_gbsa = pd.DataFrame(GBSA_output_res).sort_values('Name')
            if list(pd_gbsa.columns) != ['Name']:
                pd_gbsa.to_csv(os.path.join(out_wdir, f'GBSA_output_{unique_id}.csv'), sep='\t', index=False)

        if PBSA_output_res:
            pd_pbsa = pd.DataFrame(PBSA_output_res).sort_values('Name')
            if list(pd_pbsa.columns) != ['Name']:
                pd_pbsa.to_csv(os.path.join(out_wdir, f'PBSA_output_{unique_id}.csv'), sep='\t', index=False)

        finished_complexes_file = os.path.join(out_wdir, f"finished_gbsa_files_{unique_id}.txt")
        with open(finished_complexes_file, "w") as output:
            output.write("\n".join(var_gbsa_dat_files))

        logging.info(
            f"gmxMMPBSA energy calculation of {len(var_gbsa_dat_files)} were successfully finished.\n"
            f"Successfully finished complexes have been saved in {finished_complexes_file} file"
        )

    if decomp_dat_files:
        decomp_dat_dfs = []
        for f in decomp_dat_files:
            if os.path.isfile(f):
                df = parse_gmxMMPBSA_decomp_dat(f)
                if not df.empty:
                    df.insert(0, "Name", pathlib.PurePath(f).parent.name)
                    decomp_dat_dfs.append(df)
        if decomp_dat_dfs:
            all_decomp_dat = pd.concat(decomp_dat_dfs, ignore_index=True)
            gb_decomp_dat = all_decomp_dat[all_decomp_dat["Method"] == "GB"]
            if not gb_decomp_dat.empty:
                gb_decomp_dat.to_csv(
                    os.path.join(out_wdir, f"GBSA_decomp_avg_{unique_id}.csv"),
                    sep="\t",
                    index=False,
                )
            pb_decomp_dat = all_decomp_dat[all_decomp_dat["Method"] == "PB"]
            if not pb_decomp_dat.empty:
                pb_decomp_dat.to_csv(
                    os.path.join(out_wdir, f"PBSA_decomp_avg_{unique_id}.csv"),
                    sep="\t",
                    index=False,
                )


    logging.info(f"gmxMMPBSA energy calculation of {len(var_gbsa_dat_files)} complexes were successfully finished.\n")

    finished_gbsa_dat_file = os.path.join(out_wdir, f"finished_gbsa_files_{unique_id}.txt")
    if var_gbsa_dat_files:
        with open(finished_gbsa_dat_file, "w") as output:
            output.write("\n".join(var_gbsa_dat_files))
        logging.info(f"Successfully finished complexes have been saved in {finished_gbsa_dat_file} file")

    if decomp_dat_files:
        logging.info(f"gmxMMPBSA decomposition analysis of {len(decomp_dat_files)} complexes were successfully finished.\n")
        finished_decomp_dat_file = os.path.join(out_wdir, f"finished_decomp_files_{unique_id}.txt")
        with open(finished_decomp_dat_file, "w") as output:
            output.write("\n".join(decomp_dat_files))
        logging.info(f"Successfully finished decomposition outputs have been saved in {finished_decomp_dat_file} file")


def main():
    """CLI entry point for GBSA calculations."""
    parser = argparse.ArgumentParser(description='''Run MM-GBSA/MM-PBSA calculation using gmx_MMPBSA tool''')
    parser.add_argument('--config', metavar='FILENAME', required=False,
                        type=partial(filepath_type, ext=("yml", "yaml")),
                        help='Path to YAML configuration file with default arguments')
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
                        help='gmxMMPBSA out files (FINAL_RESULTS*.dat) to parse. If set will be used over other variables.')
    parser.add_argument('--out_decomp_files', nargs='+', default=None, type=filepath_type,
                        help='gmxMMPBSA decomposition out dat files (FINAL_DECOMP*.dat) to parse. If set will be used over other variables.')
    parser.add_argument('--hostfile', metavar='FILENAME', required=False, type=str, default=None,
                        help='text file with addresses of nodes of dask SSH cluster. The most typical, it can be '
                             'passed as $PBS_NODEFILE variable from inside a PBS script. The first line in this file '
                             'will be the address of the scheduler running on the standard port 8786. If omitted, '
                             'calculations will run on a single machine as usual.')
    parser.add_argument('-c', '--ncpu', metavar='INTEGER', required=False, default=len(os.sched_getaffinity(0)),
                        type=int, help='number of CPU per server. Use all available cpus by default.')
    parser.add_argument('--ligand_id', metavar='UNL', default='UNL', help='Ligand residue ID')
    parser.add_argument('-a', '--append_protein_selection', metavar='STRING', required=False, default=None,
                        nargs='*', help='residue IDs whuch will be included in the protein system (cofactors).'
                                        'Example: ZN MG')
    parser.add_argument('--clean_previous', action='store_true', default=False,
                        help=' Clean previous temporary gmxMMPBSA files')
    parser.add_argument('--debug', action='store_true', default=False,
                        help=' Save all temporary gmxMMPBSA files')
    parser.add_argument('-o', '--out_suffix', default=None,
                        help='Unique suffix for output files. By default, start-time_unique-id.'
                             'Unique suffix is used to separate outputs from different runs.')
    args, _ = parse_with_config(parser, sys.argv[1:])

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
              mmpbsa_file=args.mmpbsa, ncpu=args.ncpu, unique_id=unique_id,
              gmxmmpbsa_dat_files=args.out_files,
              gmxmmpbsa_decomp_dat_files=args.out_decomp_files,
              ligand_resid=args.ligand_id,
              append_protein_selection=args.append_protein_selection,
              hostfile=args.hostfile, bash_log=bash_log,
              clean_previous=args.clean_previous, debug=args.debug)
    finally:
        logging.shutdown()
