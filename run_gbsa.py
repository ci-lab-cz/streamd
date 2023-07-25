import argparse
import logging
import os
import subprocess
from datetime import datetime
from multiprocessing import cpu_count

from utils.dask_init import init_dask_cluster, calc_dask
from utils.utils import get_index, filepath_type


# mpirun -np $NP gmx_MMPBSA MPI -O -i mmpbsa.in  -cs $tpr -ci $index -cg $index_protein $index_ligand -ct $xtc -cp $topol -nogui


def run_gbsa_from_wdir(wdir, tpr, xtc, topol, index, mmpbsa, deffnm, np, ligand_resid):
    if tpr is None:
        tpr = os.path.join(wdir, f'{deffnm}.tpr')
    if xtc is None:
        xtc = os.path.join(wdir, 'md_fit.xtc')
    if topol is None:
        topol = os.path.join(wdir, 'topol.top')
    if index is None:
        index = os.path.join(wdir, 'index.ndx')
    if mmpbsa is None:
        mmpbsa = os.path.join(wdir, 'mmpbsa.in')

    if ligand_resid is None:
        ligand_resid = 'UNL'

    index_list = get_index(os.path.join(wdir, 'index.ndx'))
    protein_index = index_list.index('Protein')
    ligand_index = index_list.index(ligand_resid)

    return calc_gbsa(wdir=wdir, tpr=tpr, xtc=xtc, topol=topol,
                     index=index, mmpbsa=mmpbsa,
                     np=np, protein_index=protein_index,
                     ligand_index=ligand_index,
                     wdir_out=wdir)


def calc_gbsa(wdir, tpr, xtc, topol, index, mmpbsa, np, protein_index, ligand_index, wdir_out):
    try:
        subprocess.check_output(f'cd {wdir}; mpirun -np {np} gmx_MMPBSA MPI -O -i {mmpbsa} '
                                f' -cs {tpr} -ci {index} -cg {protein_index} {ligand_index} -ct {xtc} -cp {topol} -nogui '
                                f'-o {os.path.join(wdir_out, "FINAL_RESULTS_MMPBSA.dat")} '
                                f'-eo {os.path.join(wdir_out, "FINAL_RESULTS_MMPBSA.csv")}', shell=True)
    except subprocess.CalledProcessError as e:
        logging.error(f'{xtc}\t{e}\n')
        return None
    return wdir


def main(wdir_to_run, tpr, xtc, topol, index, wdir, mmpbsa, deffnm, ncpu, ligand_resid, hostfile):
    # TODO calc number of frames - split ncpu
    dask_client = init_dask_cluster(hostfile=hostfile, n_tasks_per_node=1, ncpu=ncpu)
    try:
        var_gbsa_wdirs = []
        for res in calc_dask(run_gbsa_from_wdir, wdir_to_run, dask_client=dask_client,
                             tpr=tpr, xtc=xtc, topol=topol, index=index, deffnm=deffnm,
                             mmpbsa=mmpbsa, np=ncpu, ligand_resid=ligand_resid):
            if res:
                var_gbsa_wdirs.append(res)
    finally:
        dask_client.shutdown()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''Run GBSA/PBSA calculation using gmx_gbsa tool''')
    parser.add_argument('--wdir_to_run', metavar='DIRNAME', required=False, default=None, nargs='+', type=filepath_type,
                        help='''directories for the previous simulations. Use to extend or continue the simulation. '
                             Should consist of: tpr, cpt, xtc files''')
    parser.add_argument('--topol', metavar='topol.top', required=False, default=None, type=filepath_type,
                        help='topol file from the the MD simulation')
    parser.add_argument('--tpr', metavar='md_out.tpr', required=False, default=None, type=filepath_type,
                        help='tpr file from the the MD simulation')
    parser.add_argument('--xtc', metavar='md_fit.xtc', required=False, default=None, type=filepath_type,
                        help='xtc file of the simulation. Trajectory should have no PBC and be fitted on the Protein_Ligand group')
    parser.add_argument('--index', metavar='index.ndx', required=False, default=None, type=filepath_type,
                        help='index file from the simulation')
    parser.add_argument('-m', '--mmpbsa', metavar='mmpbsa.in', default=None, type=filepath_type,
                        help='')
    parser.add_argument('-d', '--wdir', metavar='WDIR', default=None, type=filepath_type,
                        help='Working directory. If not set the current directory will be used.')
    parser.add_argument('--hostfile', metavar='FILENAME', required=False, type=str, default=None,
                        help='text file with addresses of nodes of dask SSH cluster. The most typical, it can be '
                             'passed as $PBS_NODEFILE variable from inside a PBS script. The first line in this file '
                             'will be the address of the scheduler running on the standard port 8786. If omitted, '
                             'calculations will run on a single machine as usual.')
    parser.add_argument('-c', '--ncpu', metavar='INTEGER', required=False, default=cpu_count(), type=int,
                        help='number of CPU per server. Use all cpus by default.')
    parser.add_argument('--deffnm', metavar='preffix for md files', required=False, default='md_out',
                        help='''preffix for the previous md files. Use to extend or continue the simulation.
                        Only if wdir_to_continue is used. Use if each --tpr, --cpt, --xtc arguments are not set up. 
                        Files deffnm.tpr, deffnm.cpt, deffnm.xtc will be used from wdir_to_continue''')
    parser.add_argument('--ligand_id', metavar='UNL', required=False, help='')

    args = parser.parse_args()

    if args.wdir is None:
        wdir = os.getcwd()
    else:
        wdir = args.wdir
    # TODO log file?
    log_file = os.path.join(wdir, f'log_gbsa_{datetime.now().strftime("%d-%m-%Y-%H-%M-%S")}.log')

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
        main(tpr=args.tpr, xtc=args.xtc, topol=args.topol,
             index=args.index, wdir=wdir, wdir_to_run=args.wdir_to_run,
             mmpbsa=args.mmpbsa,
             deffnm=args.deffnm, ncpu=args.ncpu,
             ligand_resid=args.ligand_id, hostfile=args.hostfile)
    finally:
        logging.shutdown()
