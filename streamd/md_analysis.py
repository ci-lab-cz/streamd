from glob import glob
import os
import pandas as pd
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import rms
import seaborn as sns
import matplotlib.pyplot as plt
from streamd.scripts.xvg2png import convertxvg2png
from streamd.utils.utils import get_index, make_group_ndx, get_mol_resid_pair, run_check_subprocess

sns.set_context("paper", rc={"font.size":15,"axes.titlesize":15,"axes.labelsize":15},
                font_scale=1.5)
plt.figure(figsize=(15, 12))
plt.ioff()


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
    rmsd_analysis = rms.RMSD(universe, ref, select=selection1, groupselections=selection2)
    rmsd_analysis.run()
    columns = [selection1, *selection2] if selection2 else [selection1]
    rmsd_df = pd.DataFrame(np.round(rmsd_analysis.results.rmsd[:, 2:], 2), columns=columns)
    rmsd_df.index.name = "frame"
    rmsd_df = rmsd_df.reset_index()
    return rmsd_df


def plot_rmsd(rmsd_df, system_name, out):
    rmsd_df['frame'] = rmsd_df['frame'] /100
    plot = rmsd_df.set_index('frame').plot(title=f"RMSD of {system_name}")
    plt.ylabel("RMSD (Å)")
    plt.xlabel("Time (ns)")
    plt.legend(loc='lower right', ncol=len(rmsd_df.columns) // 2, frameon=False)
    plot.figure.savefig(out,  bbox_inches="tight")

def md_rmsd_analysis(universe, wdir, system_name,
                     molid_resid_pairs, ligand_resid="UNL", active_site_dist=5.0):
    #groupselections = ['protein']
    groupselections = []
    molid_resid_pairs = dict(molid_resid_pairs)
    ligand_name = None
    if molid_resid_pairs:
        if 'UNL' in molid_resid_pairs.values():
            groupselections.append(f'backbone and (around {active_site_dist} resname UNL)')
            for molid, resid in molid_resid_pairs.items():
                if resid == ligand_resid:
                    ligand_name = molid
                    break
        groupselections.extend([f"resname {i} and not name H*" for i in molid_resid_pairs.values()])

    rmsd_df = rmsd_for_atomgroups(universe, selection1="backbone",
                                  selection2 = groupselections)
    rmsd_df = rmsd_df.rename(
        {f'backbone and (around {active_site_dist} resname {ligand_resid})': f'ActiveSite{active_site_dist}A',
         f'resname {ligand_resid} and not name H*': 'ligand'}, axis='columns')
    rmsd_df = rmsd_df.rename({f"resname {i[1]} and not name H*": f"{i[0]}" for i in molid_resid_pairs.items()}, axis='columns')

    plot_rmsd(rmsd_df.drop(f'ActiveSite{active_site_dist}A', axis='columns'),
              system_name=system_name, out=os.path.join(wdir, f'rmsd_{system_name}.png'))

    rmsd_df.loc[:, 'ligand_name'] = ligand_name
    rmsd_df.loc[:, 'system'] = system_name.replace(f'_{ligand_name}', '')

    rmsd_df.to_csv(os.path.join(wdir, f'rmsd_{system_name}.csv'), sep='\t', index=False)
    return os.path.join(wdir, f'rmsd_{system_name}.csv')


def run_md_analysis(var_md_dirs_deffnm, mdtime_ns, project_dir, bash_log,
                    ligand_resid='UNL', ligand_list_file_prev=None, env=None):
    wdir, deffnm = var_md_dirs_deffnm
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

    index_list = get_index(os.path.join(wdir, 'index.ndx'))

    # choose group to fit the trajectory
    if os.path.isfile(molid_resid_pairs_fname) and os.path.getsize(molid_resid_pairs_fname) > 0:
        molid_resid_pairs = list(get_mol_resid_pair(molid_resid_pairs_fname))
        if f'Protein_{ligand_resid}' not in index_list:
            if not make_group_ndx(query=f'"Protein"|{index_list.index(ligand_resid)}', wdir=wdir,  bash_log=bash_log):
                return None
            index_list = get_index(os.path.join(wdir, 'index.ndx'))

        index_group = index_list.index(f'Protein_{ligand_resid}')
    else:
        molid_resid_pairs = []
        index_group = index_list.index('Protein')

    dtstep = 50 if mdtime_ns <= 10 else 100

    tpr = os.path.join(wdir, f'{deffnm}.tpr')
    xtc = os.path.join(wdir, f'{deffnm}.xtc')

    cmd = f'wdir={wdir} index_group={index_group} dtstep={dtstep} deffnm={deffnm} tpr={tpr} xtc={xtc} ' \
           f'bash {os.path.join(project_dir, "scripts/script_sh/md_analysis.sh")} >> {os.path.join(wdir, bash_log)} 2>&1'

    if not run_check_subprocess(cmd, key=wdir, log=os.path.join(wdir, bash_log), env=env):
        return None

    # molid resid pairs for all ligands in the MD system
    # calc rmsd
    universe = mda.Universe(tpr, os.path.join(wdir, f'md_fit.xtc'))
    md_rmsd_analysis(universe, wdir,
                     system_name=os.path.split(wdir)[-1],
                     ligand_resid=ligand_resid,
                     molid_resid_pairs=molid_resid_pairs,
                     active_site_dist=5.0)

    for xvg_file in glob(os.path.join(wdir, '*.xvg')):
        convertxvg2png(xvg_file, transform_nm_to_A=True)
    return wdir
