#!/usr/bin/env python3

import argparse
from datetime import datetime
import os
import shutil
from functools import partial
from glob import glob
import logging
import pathlib

import MDAnalysis as mda
import pandas as pd
import prolif as plf
from prolif.plotting.barcode import Barcode
from prolif.plotting.network import LigNetwork
import matplotlib.pyplot as plt

from streamd.utils.dask_init import init_dask_cluster, calc_dask
from streamd.utils.utils import filepath_type
from streamd.prolif.prolif2png import convertprolif2png
from streamd.prolif.prolif_frame_map import convertplifbyframe2png
plt.ioff()

class RawTextArgumentDefaultsHelpFormatter(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass


def backup_output(output):
    if os.path.isfile(output):
        all_outputs = glob(os.path.join(os.path.dirname(output), f'#{os.path.basename(output)}*#'))
        n = len(all_outputs) + 1
        shutil.move(output, os.path.join(os.path.dirname(output), f'#{os.path.basename(output)}.{n}#'))


def run_prolif_task(tpr, xtc, protein_selection, ligand_selection, step, verbose, output, n_jobs,
                    occupancy = 0.6, save_viz=True, dpi=300, plot_width=15, plot_height=8, pdb=None):
    '''

    :param tpr:
    :param xtc:
    :param protein_selection:
    :param ligand_selection:
    :param step:
    :param verbose:
    :param output:
    :param n_jobs:
    :param save_pics: save barcode in png and network in html
    :param dpi:
    :param plot_width:  in inches
    :param plot_height: in inches
    :return: pandas dataframe
    '''
    u = mda.Universe(tpr, xtc, in_memory=False, in_memory_step=1)

    protein = u.atoms.select_atoms(protein_selection)
    ligand = u.atoms.select_atoms(ligand_selection)

    if pdb:
        u1 = mda.Universe(pdb)
        protein_pdb = u1.atoms.select_atoms(protein_selection)
        if len(protein.residues.resids) == len(protein_pdb.residues.resids):
            protein.residues.resids = protein_pdb.residues.resids
        if len(protein.segments.segids) == len(protein_pdb.segments.segids):
            protein.segments.segids = protein_pdb.segments.segids

    fp = plf.Fingerprint(['Hydrophobic', 'HBDonor', 'HBAcceptor', 'Anionic', 'Cationic', 'CationPi', 'PiCation',
                          'PiStacking', 'MetalAcceptor'])
    fp.run(u.trajectory[::step] if step > 1 else u.trajectory, ligand, protein, progress=verbose, n_jobs=n_jobs)

    df = fp.to_dataframe()
    df.columns = ['.'.join(item.strip().lower() for item in items[1:]) for items in df.columns]
    df = df.reindex(sorted(df.columns), axis=1)
    df.to_csv(output, sep='\t')

    if save_viz:
        # barcode
        Barcode.from_fingerprint(fp).display(figsize=(plot_width, plot_height)).figure.savefig(f'{output.rstrip(".csv")}.png', dpi=dpi)
        # Net
        LigNetwork.from_fingerprint(fp, ligand_mol=ligand.convert_to('rdkit'), threshold=occupancy).save(f'{output.rstrip(".csv")}_occupancy{occupancy}.html')
        convertplifbyframe2png(plif_out_file=output, plot_width=plot_width, plot_height=plot_height)

    return df


def run_prolif_from_wdir(wdir, tpr, xtc, protein_selection, ligand_selection, step, verbose, output,
                         plot_width, plot_height, save_viz, pdb, n_jobs, occupancy):
    tpr = os.path.join(wdir, tpr)
    xtc = os.path.join(wdir, xtc)
    if pdb:
        pdb = os.path.join(wdir, pdb)
    output = os.path.join(wdir, output)
    backup_output(output)

    if not os.path.isfile(tpr) or not os.path.isfile(xtc):
        print(f'{wdir}: cannot run prolif. Check if there are missing files: {tpr} {xtc}. Skip such directory')
        return None

    run_prolif_task(tpr=tpr, xtc=xtc, protein_selection=protein_selection,
                    ligand_selection=ligand_selection, step=step, verbose=verbose, output=output,
                    plot_width=plot_width, plot_height=plot_height, save_viz=save_viz, occupancy=occupancy,
                    pdb=pdb, n_jobs=n_jobs)
    return output


def collect_outputs(output_list, output):
    df_list = []
    for i in output_list:
        df = pd.read_csv(i, sep='\t')
        # save dirname - protein_ligand pair
        df['Name'] = pathlib.PurePath(i).parent.name
        df_list.append(df)

    df_aggregated = pd.concat(df_list)
    df_aggregated = df_aggregated.fillna(False).sort_values('Frame')
    amino_acids = df_aggregated.columns.drop(['Name', 'Frame']).to_list()
    # sort by number and type of interaction
    amino_acids.sort(key=lambda x: (int(x.split('.')[0][3:]), x.split('.')[1]))
    sorted_columns = ['Name', 'Frame'] + amino_acids
    df_aggregated.loc[:, sorted_columns].to_csv(output, sep='\t', index=False)


def start(wdir_to_run, wdir_output, tpr, xtc, step, append_protein_selection,
          protein_selection, ligand_resid, hostfile, ncpu, n_jobs,
          occupancy, plot_width, plot_height, save_viz, unique_id, pdb, verbose):
    '''

    :param wdir_to_run: list
    :param wdir_output: path to dirn
    :param tpr: path to file
    :param xtc: path to file
    :param step: int
    :param append_protein_selection: str
    :param protein_selection: str
    :param ligand_resid: str
    :param hostfile: Nonen or path to file
    :param ncpu: int
    :param n_jobs: int
    :param occupancy: float
    :param plot_width: float
    :param plot_height: float
    :param save_viz: bool
    :param unique_id: str
    :param pdb: None or path to file (protein.pdb for renumbering)
    :param verbose: bool
    :return:
    '''
    output = 'plifs.csv'
    output_aggregated = os.path.join(wdir_output, f'prolif_output_{unique_id}.csv')

    if append_protein_selection is not None:
        protein_selection = f'({protein_selection}) or ({append_protein_selection})'

    ligand_selection = f'resname {ligand_resid}'

    if wdir_to_run is not None:
        dask_client, cluster = None, None
        #n_jobs_per_task = n_jobs if n_jobs <= ncpu else ncpu
        if n_jobs is None:
            # limit to 12 https://github.com/chemosim-lab/ProLIF/issues/110
            n_jobs_per_task = min(12, ncpu // len(wdir_to_run))
        else:
            n_jobs_per_task = n_jobs

        logging.info(f'Allocating {n_jobs_per_task} n_jobs per each task.')

        try:
            dask_client, cluster = init_dask_cluster(hostfile=hostfile,
                                                     n_tasks_per_node=min(len(wdir_to_run), ncpu//n_jobs_per_task),
                                                     use_multi_servers=True if len(wdir_to_run) > ncpu else False,
                                                     ncpu=ncpu)
            var_prolif_out_files = []
            for res in calc_dask(run_prolif_from_wdir, wdir_to_run, dask_client=dask_client,
                                 tpr=tpr, xtc=xtc, protein_selection=protein_selection,
                                 ligand_selection=ligand_selection, step=step, verbose=verbose, output=output,
                                 plot_width=plot_width, plot_height=plot_height, save_viz=save_viz, pdb=pdb,
                                 n_jobs=n_jobs_per_task, occupancy=occupancy):
                if res:
                    var_prolif_out_files.append(res)
        finally:
            if dask_client:
                dask_client.retire_workers(dask_client.scheduler_info()['workers'],
                                           close_workers=True, remove=True)
                dask_client.shutdown()
            if cluster:
                cluster.close()
    else:
        output = os.path.join(os.path.dirname(xtc), output)
        backup_output(output)

        run_prolif_task(tpr=tpr, xtc=xtc, protein_selection=protein_selection,
                        ligand_selection=ligand_selection,
                        step=step, verbose=verbose, output=output,
                        plot_width=plot_width, plot_height=plot_height, save_viz=save_viz,
                        pdb=pdb, n_jobs= min(12, ncpu), occupancy=occupancy)
        var_prolif_out_files = [output]

    backup_output(output_aggregated)
    collect_outputs(var_prolif_out_files, output=output_aggregated)

    convertprolif2png(output_aggregated, occupancy=occupancy, plot_width=plot_width, plot_height=plot_height)
    finished_complexes_file = os.path.join(wdir_output, f"finished_prolif_files_{unique_id}.txt")
    with open(finished_complexes_file, 'w') as output:
        output.write("\n".join(var_prolif_out_files))

    logging.info(
        f'ProLIF calculation of {len(var_prolif_out_files)} were successfully finished.\n'
         f'Successfully finished complexes have been saved in {finished_complexes_file} file')



def main():
    parser = argparse.ArgumentParser(description='Get protein-ligand interactions from MD trajectories using '
                                                 'ProLIF module.',
                                     formatter_class=RawTextArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--wdir_to_run', metavar='DIRNAME', required=False, default=None, nargs='+',
                        type=partial(filepath_type, exist_type='dir'),
                        help='''single or multiple directories for simulations.
                             Should consist of: md_out.tpr and md_fit.xtc files''')
    parser.add_argument('--xtc', metavar='FILENAME', required=False,
                        help='input trajectory file (XTC). Will be ignored if --wdir_to_run is used')
    parser.add_argument('--tpr', metavar='FILENAME', required=False,
                        help='input topology file (TPR). Will be ignored if --wdir_to_run is used')
    parser.add_argument('-l', '--ligand', metavar='STRING', required=False, default='UNL',
                        help='residue name of a ligand in the input trajectory.')
    parser.add_argument('-s', '--step', metavar='INTEGER', required=False, default=1, type=int,
                        help='step to take every n-th frame. ps')
    parser.add_argument('--protein_selection', metavar='STRING',
                        required=False, default='protein',
                        help='The protein selection atoms. '
                             'Example: "protein" or "protein and byres around 20.0 resname UNL"')
    parser.add_argument('-a', '--append_protein_selection', metavar='STRING', required=False, default=None,
                        help='the string which will be concatenated to the protein selection atoms. '
                             'Example: "resname ZN or resname MG".')
    parser.add_argument('-d', '--wdir', metavar='WDIR', default=None,
                        type=partial(filepath_type, check_exist=False, create_dir=True),
                        help='Working directory for program output. If not set the current directory will be used.')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='print progress.')
    parser.add_argument('--hostfile', metavar='FILENAME', required=False, type=str, default=None,
                        help='text file with addresses of nodes of dask SSH cluster. The most typical, it can be '
                             'passed as $PBS_NODEFILE variable from inside a PBS script. The first line in this file '
                             'will be the address of the scheduler running on the standard port 8786. If omitted, '
                             'calculations will run on a single machine as usual.')
    parser.add_argument('-c', '--ncpu', metavar='INTEGER', required=False, default=len(os.sched_getaffinity(0)), type=int,
                        help='number of CPU per server. By default, StreaMD utilizes all available cpus.')
    parser.add_argument('--n_jobs', metavar='INTEGER', required=False,
                         default=None, type=int,
                         help='Number of processes to run per each interaction analysis tasks. '
                              'By default, StreaMD distributes the specified number of cores (--ncpu) evenly '
                              'between the available CPUs and the number of tasks to execute (e.g., multiple directories provided via --wdir_to_run). '
                              'However, by default, the --n_jobs value is limited to 12 to avoid the bottleneck issue (https://github.com/chemosim-lab/ProLIF/issues/110) described by the ProLIF authors .'
                              'Users can override this limitation by explicitly specifying the --n_jobs argument value.')
    parser.add_argument('--width', metavar='FILENAME', default=15, type=int,
                        help='width of the output pictures')
    parser.add_argument('--height', metavar='FILENAME', default=10, type=int,
                        help='height of the output pictures')
    parser.add_argument('--occupancy', metavar='float', default=0.6, type=float,
                        help='occupancy of the unique contacts to show. '
                             'Applied for plifs_occupancyX.html (for each complex) and'
                             ' prolif_output_occupancyX.png (all systems aggregated plot)')
    parser.add_argument('--not_save_pics', default=False, action='store_true',
                        help='not create html and png files (by frames) for each unique trajectory.'
                             ' Only overall prolif png file will be created.')
    parser.add_argument('-o','--out_suffix',
                        metavar='string', default=None,
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

    log_file = os.path.join(wdir, f'log_prolif_{unique_id}.log')

    logging.basicConfig(format='%(asctime)s - %(name)s - %(levelname)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S',
                        level=logging.INFO,
                        handlers=[logging.FileHandler(log_file),
                                  logging.StreamHandler()])


    if args.wdir_to_run is not None:
        tpr = 'md_out.tpr'
        xtc = 'md_fit.xtc'
        pdb = 'frame.pdb'
    else:
        tpr = args.tpr
        xtc = args.xtc
        pdb = None

    logging.getLogger('distributed').setLevel('CRITICAL')
    logging.getLogger('distributed.core').setLevel('CRITICAL')
    logging.getLogger('asyncssh').setLevel('CRITICAL')
    logging.getLogger('distributed.worker').setLevel('CRITICAL')
    logging.getLogger('distributed.comm').setLevel('CRITICAL')
    logging.getLogger('distributed.nanny').setLevel('CRITICAL')
    logging.getLogger('bockeh').setLevel('CRITICAL')
    logging.getLogger('matplotlib.font_manager').setLevel('CRITICAL')

    logging.info(args)
    try:
        start(wdir_to_run=args.wdir_to_run, wdir_output=wdir, tpr=tpr,
          xtc=xtc, step=args.step, append_protein_selection=args.append_protein_selection,
          protein_selection=args.protein_selection, ligand_resid=args.ligand, hostfile=args.hostfile, ncpu=args.ncpu,
          n_jobs=args.n_jobs, occupancy=args.occupancy, plot_width=args.width, plot_height=args.height,
          save_viz=not args.not_save_pics, unique_id=unique_id, pdb=pdb,
          verbose=args.verbose)
    finally:
        logging.shutdown()


if __name__ == '__main__':
    main()
