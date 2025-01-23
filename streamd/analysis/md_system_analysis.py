from glob import glob
import os
import shutil
import pandas as pd
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import rms
from streamd.analysis.xvg2png import convertxvg2png
from streamd.analysis.plot_build import plot_rmsd
from streamd.utils.utils import get_index, make_group_ndx, get_mol_resid_pair, run_check_subprocess, backup_prev_files



def rmsd_for_atomgroups(universe, selection1, selection2=None):
    """Calulate the RMSD for selected atom groups.

    Parameters
    ----------
    universe: MDAnalysis.core.universe.Universe
        MDAnalysis universe.
    selection1: str
        Selection string for main atom group, also used during alignment.
    selection2: list of str, optional
        Selection strings for additional atom groups.

    Returns
    -------
    rmsd_df: pandas.core.frame.DataFrame
        DataFrame containing RMSD of the selected atom groups over time.
    """

    universe.trajectory[0]
    ref = universe
    rmsd_analysis = rms.RMSD(universe, ref, select=selection1, groupselections=selection2, in_memory=False)
    rmsd_analysis.run()
    columns = [selection1, *selection2] if selection2 else [selection1]
    rmsd_df = pd.DataFrame(np.round(rmsd_analysis.results.rmsd[:, 2:], 2), columns=columns)
    rmsd_df.index.name = "frame"
    rmsd_df = rmsd_df.reset_index()
    # transform to ns
    rmsd_df['time(ns)'] = rmsd_df['frame'] / 100
    rmsd_df = rmsd_df.drop('frame', axis='columns')
    return rmsd_df


def md_rmsd_analysis(tpr, xtc, wdir_out_analysis, system_name,
                     molid_resid_pairs,
                     ligand_resid="UNL", active_site_dist=5.0):
    #groupselections = ['protein']
    rmsd_out_file = os.path.join(wdir_out_analysis, f'rmsd_{system_name}.csv')
    universe = mda.Universe(tpr, xtc, in_memory=False, in_memory_step=1)
    groupselections = []
    molid_resid_pairs = dict(molid_resid_pairs)
    ligand_name = None
    if molid_resid_pairs:
        if ligand_resid in molid_resid_pairs.values():
            for molid, resid in molid_resid_pairs.items():
                if resid == ligand_resid:
                    ligand_name = molid
                    break
            groupselections.append(f'backbone and (around {active_site_dist} resname {ligand_resid})')

        groupselections.extend([f"resname {i} and not name H*" for i in molid_resid_pairs.values()])

    rmsd_df = rmsd_for_atomgroups(universe, selection1="backbone",
                                  selection2 = groupselections)
    del universe
    rmsd_df = rmsd_df.rename(
        {f'backbone and (around {active_site_dist} resname {ligand_resid})': f'ActiveSite{active_site_dist}A',
         f'resname {ligand_resid} and not name H*': 'ligand'}, axis='columns')

    rmsd_df = rmsd_df.rename({f"resname {i[1]} and not name H*": f"{i[0]}" for i in molid_resid_pairs.items()}, axis='columns')

    plot_rmsd(rmsd_df=rmsd_df, system_name=system_name, out=os.path.join(wdir_out_analysis, f'rmsd_{system_name}.png'))

    rmsd_df.loc[:, 'ligand_name'] = ligand_name
    rmsd_df.loc[:, 'system'] = system_name.replace(f'_{ligand_name}', '') if ligand_name else system_name
    rmsd_df.loc[:, 'directory'] = wdir_out_analysis

    rmsd_df.to_csv(rmsd_out_file, sep='\t', index=False)
    return rmsd_out_file


def run_md_analysis(var_md_dirs_deffnm, mdtime_ns, project_dir, bash_log,
                    active_site_dist=5.0, ligand_resid='UNL',
                    save_traj_without_water = False,
                    analysis_dirname = 'md_analysis',
                    ligand_list_file_prev=None, env=None, system_name=None):
    wdir, deffnm = var_md_dirs_deffnm

    # create subdir for analysis files only
    wdir_out_analysis = os.path.join(wdir, analysis_dirname)
    os.makedirs(wdir_out_analysis, exist_ok=True)

    if ligand_list_file_prev is None:
        molid_resid_pairs_fname = os.path.join(wdir, 'all_ligand_resid.txt')
    else:
        molid_resid_pairs_fname = ligand_list_file_prev

    # if md was continued
    if not os.path.isfile(os.path.join(wdir, 'index.ndx')):
        cmd = f'''
        cd {wdir}
        gmx make_ndx -f {deffnm}.gro << INPUT
        q
        INPUT
        '''
        if not run_check_subprocess(cmd, key=wdir, log=os.path.join(wdir, bash_log), env=env):
            return None

    index_list = get_index(os.path.join(wdir, 'index.ndx'), env=env)

    # choose group to fit the trajectory
    if os.path.isfile(molid_resid_pairs_fname) and os.path.getsize(molid_resid_pairs_fname) > 0:
        molid_resid_pairs = list(get_mol_resid_pair(molid_resid_pairs_fname))
        if f'Protein_{ligand_resid}' not in index_list:
            if not make_group_ndx(query=f'"Protein"|{index_list.index(ligand_resid)}', wdir=wdir,  bash_log=bash_log, env=env):
                return None
            index_list = get_index(os.path.join(wdir, 'index.ndx'), env=env)

        index_group = index_list.index(f'Protein_{ligand_resid}')
    else:
        molid_resid_pairs = []
        index_group = index_list.index('Protein')

    dtstep = 50 if mdtime_ns <= 10 else 100

    tpr = os.path.join(wdir, f'{deffnm}.tpr')
    xtc = os.path.join(wdir, f'{deffnm}.xtc')

    if not system_name:
        system_name = os.path.split(wdir)[-1]

    cmd = f'wdir={wdir} index_group={index_group} dtstep={dtstep} deffnm={deffnm} tpr={tpr} xtc={xtc} wdir_out_analysis={wdir_out_analysis} system_name={system_name} ' \
           f'bash {os.path.join(project_dir, "scripts/script_sh/md_analysis.sh")} >> {os.path.join(wdir, bash_log)} 2>&1'

    if not run_check_subprocess(cmd, key=wdir, log=os.path.join(wdir, bash_log), env=env):
        return None

    # molid resid pairs for all ligands in the MD system
    # calc rmsd
    # universe = mda.Universe(tpr, os.path.join(wdir, f'md_fit.xtc'))

    rmsd_out_file = md_rmsd_analysis(
        tpr=os.path.join(wdir, 'md_out_nowater.tpr'),
        xtc=os.path.join(wdir, f'md_fit_nowater.xtc'),
        # tpr=os.path.join(wdir, 'md_out.tpr'), xtc=os.path.join(wdir, f'md_fit.xtc'),
                     wdir_out_analysis=wdir_out_analysis,
                     system_name=system_name,
                     ligand_resid=ligand_resid,
                     molid_resid_pairs=molid_resid_pairs,
                     active_site_dist=active_site_dist)
    if not save_traj_without_water:
        os.remove(os.path.join(wdir, 'md_out_nowater.tpr'))
        os.remove(os.path.join(wdir, f'md_fit_nowater.xtc'))

    for xvg_file in glob(os.path.join(wdir_out_analysis, '*.xvg')):
        convertxvg2png(xvg_file, system_name=system_name, transform_nm_to_A=True)
    return rmsd_out_file, wdir_out_analysis, wdir
