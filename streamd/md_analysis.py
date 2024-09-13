from glob import glob
import os

from streamd.scripts.xvg2png import convertxvg2png
from streamd.utils.utils import get_index, make_group_ndx, get_mol_resid_pair, run_check_subprocess


def md_lig_rmsd_analysis(molid, resid, tpr, xtc, wdir, tu, bash_log, project_dir, env=None):
    index_list = get_index(os.path.join(wdir, 'index.ndx'))
    if f'{resid}_&_!H*' not in index_list:
        if not make_group_ndx(query=f'{index_list.index(resid)} & ! a H*', wdir=wdir, bash_log=bash_log):
            return None
        index_list = get_index(os.path.join(wdir, 'index.ndx'))
    index_ligand_noH = index_list.index(f'{resid}_&_!H*')
    cmd = f'wdir={wdir} tpr={tpr} xtc={xtc}  molid={molid} tu={tu} index_ligand_noH={index_ligand_noH} ' \
          f'bash {os.path.join(project_dir, "scripts/script_sh/md_ligand_analysis.sh")} >> {os.path.join(wdir, bash_log)} 2>&1'

    run_check_subprocess(cmd, key=wdir, log=os.path.join(wdir, bash_log), env=env)


def run_md_analysis(var_md_dirs_deffnm, mdtime_ns, project_dir, bash_log, ligand_resid='UNL', ligand_list_file_prev=None, env=None):
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
        molid_resid_pairs = get_mol_resid_pair(molid_resid_pairs_fname)
        if f'Protein_{ligand_resid}' not in index_list:
            if not make_group_ndx(query=f'"Protein"|{index_list.index(ligand_resid)}', wdir=wdir,  bash_log=bash_log):
                return None
            index_list = get_index(os.path.join(wdir, 'index.ndx'))

        index_group = index_list.index(f'Protein_{ligand_resid}')
    else:
        molid_resid_pairs = []
        index_group = index_list.index('Protein')

    tu = 'ps' if mdtime_ns <= 10 else 'ns'
    dtstep = 50 if mdtime_ns <= 10 else 100

    tpr = os.path.join(wdir, f'{deffnm}.tpr')
    xtc = os.path.join(wdir, f'{deffnm}.xtc')

    cmd = f'wdir={wdir} index_group={index_group} tu={tu} dtstep={dtstep} deffnm={deffnm} tpr={tpr} xtc={xtc} ' \
           f'bash {os.path.join(project_dir, "scripts/script_sh/md_analysis.sh")} >> {os.path.join(wdir, bash_log)} 2>&1'

    if not run_check_subprocess(cmd, key=wdir, log=os.path.join(wdir, bash_log), env=env):
        return None

    # molid resid pairs for all ligands in the MD system
    # calc rmsd
    for molid, resid in molid_resid_pairs:
        md_lig_rmsd_analysis(molid=molid, resid=resid, xtc=os.path.join(wdir, f'md_fit.xtc'), tpr=tpr,
                             wdir=wdir, tu=tu, bash_log=bash_log, project_dir=project_dir, env=env)

    for xvg_file in glob(os.path.join(wdir, '*.xvg')):
        convertxvg2png(xvg_file, transform_nm_to_A=True)
    return wdir
