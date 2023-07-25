import logging
import os
import subprocess

from utils.utils import get_index, make_group_ndx, get_mol_resid_pair


def md_lig_rmsd_analysis(molid, resid, wdir, tu):
    index_list = get_index(os.path.join(wdir, 'index.ndx'))
    if f'{resid}_&_!H*' not in index_list:
        if not make_group_ndx(query=f'{index_list.index(resid)} & ! a H*', wdir=wdir):
            return None
        index_list = get_index(os.path.join(wdir, 'index.ndx'))
    index_ligand_noH = index_list.index(f'{resid}_&_!H*')

    try:
        subprocess.check_output(f'''
        cd {wdir}
        gmx rms -s md_out.tpr -f md_fit.xtc -o rmsd_{molid}.xvg -n index.ndx  -tu {tu} <<< "Backbone  {index_ligand_noH}"''',
                                shell=True)
    except subprocess.CalledProcessError as e:
        logging.exception(f'{wdir}\n{e}', stack_info=True)


def run_md_analysis(wdir, deffnm, mdtime_ns, project_dir):
    index_list = get_index(os.path.join(wdir, 'index.ndx'))
    molid_resid_pairs_fname = os.path.join(wdir, 'all_ligand_resid.txt')

    # choose group to fit the trajectory
    if os.path.isfile(molid_resid_pairs_fname) and os.path.getsize(molid_resid_pairs_fname) > 0:
        molid_resid_pairs = get_mol_resid_pair(molid_resid_pairs_fname)
        ligand_resid = 'UNL'
        if f'Protein_{ligand_resid}' not in index_list:
            if not make_group_ndx(query=f'"Protein"|{index_list.index(ligand_resid)}', wdir=wdir):
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

    try:
        subprocess.check_output(
            f'wdir={wdir} index_group={index_group} tu={tu} dtstep={dtstep} deffnm={deffnm} tpr={tpr} xtc={xtc} bash {os.path.join(project_dir, "scripts/script_sh/md_analysis.sh")}',
            shell=True)
    except subprocess.CalledProcessError as e:
        logging.exception(f'{wdir}\n{e}', stack_info=True)
        return None

    # molid resid pairs for all ligands in the MD system
    # calc rmsd
    for molid, resid in molid_resid_pairs:
        md_lig_rmsd_analysis(molid=molid, resid=resid, wdir=wdir, tu=tu)

    return wdir
