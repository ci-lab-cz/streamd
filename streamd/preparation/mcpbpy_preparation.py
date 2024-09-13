from glob import glob
import logging
import os
import shutil

import MDAnalysis as mda
import parmed as pmd

from streamd.utils.utils import run_check_subprocess, get_mol_resid_pair, get_index, make_group_ndx
from streamd.preparation.ligand_preparation import prepare_gaussian_files
from streamd.preparation.md_files_preparation import check_if_info_already_added_to_topol, edit_topology_file

def convert_pdb2mol2(metal_pdb, charge_dict, bash_log_curr, env):
    '''

    :param metal_pdb: metal pdb file path
    :param charge_dict: {MN:2, ZN:2, CA:2}
    :return: metal mol2 file path
    '''
    atom = os.path.basename(metal_pdb).split('_')[0]
    charge = charge_dict.get(atom, None)
    if not charge:
        logging.error(f'MCPBPY procedure: Preparation Error: Cannot find charge for {atom} metal atom. Please set up charge for the metal atom'
                      f'Calculation will be stopped')
        return None

    logging.warning(f'INFO: MCPBPY procedure: Use {charge} as a charge for {atom} metal atom')

    metal_mol2 = metal_pdb.replace(".pdb",".mol2")

    cmd = f'metalpdb2mol2.py -i {metal_pdb} -o {metal_mol2} -c {charge} >> {bash_log_curr} 2>&1'
    if not run_check_subprocess(cmd, metal_pdb, log=bash_log_curr, env=env):
        return None
    return metal_mol2


def split_metal(protein_fname, metal_resnames, wdir):
    '''
    :param protein_fname:
    :param metal_resnames:
    :param wdir:
    :return: clean_protein pdb file path, list of each metal pdb file path
    '''
    protein = mda.Universe(protein_fname)

    me_selection = ' or '.join([f'resname {i}' for i in metal_resnames])

    metal_atoms = list(protein.select_atoms(me_selection))
    protein_only = protein.select_atoms(f'not ({me_selection})')

    protein_clean_pdb = os.path.join(wdir, f'{os.path.basename(protein_fname).strip(".pdb")}_clean.pdb')
    protein_only.write(protein_clean_pdb)

    metal_pdb_list = []
    for atom in metal_atoms:
        metal_pdb_file = os.path.join(wdir, f'{atom.resname}_{atom.id}.pdb')
        protein.select_atoms(f'resid {atom.resid}').write(metal_pdb_file)

        metal_pdb_list.append(metal_pdb_file)

    return protein_clean_pdb, metal_pdb_list

def get_new_metal_ids(protein_fname, metal_resnames):
    protein = mda.Universe(protein_fname)
    me_selection = ' or '.join([f'resname {i}' for i in metal_resnames])
    metal_atoms = list(protein.select_atoms(me_selection))
    atom_ids = {}
    for atom in metal_atoms:
        atom_ids[atom.id] = atom.resname

    return atom_ids

def merge_complex(protein_pdb, ligand_mol2_list, metal_mol2_list, wdir):
    '''

    :param protein_pdb: pdb file path
    :param ligand_mol2_list: list of mol2 file paths
    :param metal_mol2_list: list of mol2 file paths
    :return: merged complex pdb file path
    '''
    def add_resids(protein, resid_list_fname):
        complex = protein
        for resid_fname in resid_list_fname:
            resid_instance = pmd.load_file(resid_fname).to_structure()
            complex = complex + resid_instance

        return complex

    protein = pmd.load_file(protein_pdb)

    complex = add_resids(protein, metal_mol2_list)
    complex = add_resids(complex, ligand_mol2_list)

    complex_file = os.path.join(wdir, 'complex.pdb')
    complex.write_pdb(complex_file)

    return complex_file

def remove_allHs_from_pdb(complex_file):
    complex_mda = mda.Universe(complex_file)
    complex_mda_noH = complex_mda.atoms.select_atoms('not (type H and protein)')

    shutil.copy(complex_file, complex_file.replace('.pdb', '_originalwithHs.pdb'))
    complex_mda_noH.write(complex_file)


def copy_rename_ligand_files(lig_wdir_list, wdir, ext_list=('mol2','frcmod')):
    '''
    mol2 and frcmod file names should correspond residue id
    :return: dict
    '''
    lig_new_file_ext_dict = {i: [] for i in ext_list}
    molids_pairs_dict = {}
    for lig_wdir in lig_wdir_list:
        if lig_wdir:
            molid, resid = next(get_mol_resid_pair(os.path.join(lig_wdir, 'resid.txt')))
            molids_pairs_dict[molid] = resid
            cur_file_name = os.path.join(lig_wdir, molid)
            new_file_name = os.path.join(wdir, resid)

            for ext in ext_list:
                if not os.path.isfile(f'{new_file_name}.{ext}'):
                    shutil.copy(f'{cur_file_name}.{ext}', f'{new_file_name}.{ext}')

                lig_new_file_ext_dict[ext].append(f'{new_file_name}.{ext}')

    return lig_new_file_ext_dict, molids_pairs_dict


def prepare_protein_in(file_in_template, file_out, complex_file,variable_ion_ids, variable_ion_mol2_files,
                       variable_ligand_mol2_files, variable_ligand_frcmod_files, gaussian_version, force_field, cut_off=2.8, large_opt=1):
    with open(file_in_template) as inp:
        data = inp.read()

    new_data = data.replace('complex.pdb', complex_file)\
        .replace('variable_ion_ids', variable_ion_ids)\
        .replace('variable_ion_mol2_files', variable_ion_mol2_files)\
        .replace('gaussian_version', gaussian_version)\
        .replace('cut_off cut_off', f'cut_off {cut_off}')\
        .replace('large_opt large_opt', f'large_opt {large_opt}')\
        .replace('force_field force_field', f'force_field {force_field}')

    if not variable_ligand_mol2_files and not variable_ligand_frcmod_files:
        new_data = new_data.replace('naa_mol2files variable_ligand_mol2_files\n', '') \
            .replace('frcmod_files variable_ligand_frcmod_files\n', '')
    else:
        new_data = new_data.replace('variable_ligand_mol2_files', variable_ligand_mol2_files)\
            .replace('variable_ligand_frcmod_files', variable_ligand_frcmod_files)\

    with open(file_out, 'w') as output:
        output.write(new_data)


def set_up_gaussian_files(wdir, ncpu, gaussian_basis, gaussian_memory):
    gaussian_files = glob(os.path.join(wdir, 'protein*.com'))
    for gau_file in gaussian_files:
        chk = gau_file.replace('.com','.chk')
        if os.path.isfile(chk):
            logging.warning(f'INFO: MCPBPY procedure: the gaussian checkpoint file exists: {chk}. Gaussian Calculation will'
                         f'begin from the last completed point in the previous calculation')
            prepare_gaussian_files(file_template=gau_file, file_out=gau_file, ncpu=ncpu,
                               gaussian_basis=gaussian_basis, gaussian_memory=gaussian_memory, opt_restart=True)
        else:
            prepare_gaussian_files(file_template=gau_file, file_out=gau_file, ncpu=ncpu,
                                   gaussian_basis=gaussian_basis, gaussian_memory=gaussian_memory, opt_restart=False)


def run_MCPBPY(protein_in_file, wdir, s, bash_log, env):
    cmd = f'cd {wdir}; MCPB.py -i {protein_in_file} -s {s} >> {bash_log} 2>&1'
    if not run_check_subprocess(cmd, wdir, log=bash_log, env=env):
        return None
    return wdir


def run_gaussian_calculation(wdir, gaussian_version, activate_gaussian, bash_log, env):
    def run_task(gau_cmd, activate_gaussian, log, wdir, env, check_only_if_exist=False):
        def check_gau_log_file(log):
            if os.path.isfile(log):
                with open(log) as checkpoint:
                    gau_log = checkpoint.read()
                if 'Normal termination of' in gau_log:
                    return True
            return False
        if (check_only_if_exist and not os.path.isfile(log)) or not check_gau_log_file(log):
            cmd = (f'cd {wdir}; {activate_gaussian};'
                   f'{gau_cmd}')
            if not run_check_subprocess(cmd, wdir, log=bash_log, env=env):
                return None
        else:
            logging.warning(f'INFO MCPBPY procedure: the gaussian calculation was finished: {wdir} {log}. Skip this step')
        return log

    small_opt_log = os.path.join(wdir, 'protein_small_opt.log')
    small_fc_log = os.path.join(wdir, 'protein_small_fc.log')
    small_opt_fchk = os.path.join(wdir, 'protein_small_opt.fchk')
    large_mk_log = os.path.join(wdir, 'protein_large_mk.log')

    logging.warning(f'INFO: MCPBPY procedure: start Gaussian geometry optimization for the small model{wdir}')
    if not run_task(gau_cmd=f'{gaussian_version} < protein_small_opt.com > {small_opt_log}', log=small_opt_log,
                    activate_gaussian=activate_gaussian, wdir=wdir, env=env):
        return None
    logging.warning(f'INFO: MCPBPY procedure: start Gaussian force constant calculation for the small model {wdir}')
    if not run_task(gau_cmd=f'{gaussian_version} < protein_small_fc.com > {small_fc_log}', log=small_fc_log,
                    activate_gaussian=activate_gaussian, wdir=wdir, env=env):
        return None
    logging.warning(f'INFO: MCPBPY procedure: generation fchk file for the small model {wdir}')
    if not run_task(gau_cmd=f'formchk protein_small_opt.chk {small_opt_fchk}', log=small_opt_fchk,
                    activate_gaussian=activate_gaussian, wdir=wdir, check_only_if_exist=True, env=env):
        return None
    logging.warning(f'INFO: MCPBPY procedure: start Gaussian geometry optimization for the large model {wdir}')
    if not run_task(gau_cmd=f'{gaussian_version} < protein_large_mk.com > {large_mk_log}', log=large_mk_log,
                    activate_gaussian=activate_gaussian, wdir=wdir, env=env):
        return None

    return wdir


def run_tleap(wdir, bash_log, env):
    #TODO check
    cmd = f'cd {wdir}; tleap -s -f protein_tleap.in > protein_tleap.out >> {bash_log} 2>&1'
    # cmd = f'cd {wdir}; tleap -s -f protein_tleap.in > protein_tleap.out'
    if not run_check_subprocess(cmd, wdir, log=bash_log, env=env):
        return None
    return wdir


def get_renamed_mcpbpy_residues(complex_original, complex_mcpbpy):
    '''

    :param complex_original:
    :param complex_mcpbpy:
    :return: dict {mcpbpy_resname:orig_resname}
    '''
    pdb_orig = mda.Universe(complex_original)
    pdb_mcpbpy = mda.Universe(complex_mcpbpy)

    residues_orig_list = pdb_orig.residues.resnames
    residues_mcpbpy_list = pdb_mcpbpy.residues.resnames

    diff_dict = {}
    for res_orig, res_mcpbpy in zip(residues_orig_list, residues_mcpbpy_list):
        if res_orig != res_mcpbpy:
            diff_dict[res_mcpbpy] = res_orig

    return diff_dict

def amber2gmx(complex_original, complex_mcpbpy, prmtop, inpcrd, wdir):
    '''
    MCPBPY renames refined amino acids. Here we rename them back and transform files from amber to gromacs format
    :param complex_original:
    :param complex_mcpbpy:
    :param prmtop:
    :param inpcrd:
    :param wdir:
    :return:
    '''

    topol_top, solv_ions_gro = os.path.join(wdir, 'topol.top'), os.path.join(wdir, 'solv_ions.gro')

    diff_residues_dict = get_renamed_mcpbpy_residues(complex_original=complex_original,
                                                     complex_mcpbpy=complex_mcpbpy)
    parm = pmd.load_file(prmtop, inpcrd)

    for res in parm.residues:
        if res.name in diff_residues_dict:
            res.name = diff_residues_dict[res.name]

    parm.save(topol_top, format='gromacs')
    parm.save(solv_ions_gro)

def create_posre(all_resids, wdir, bash_log, env):
    index_list = get_index(os.path.join(wdir, 'index.ndx'))
    couple_group_ind = '|'.join([str(index_list.index(i)) for i in ['Protein-H'] + all_resids])
    couple_group = '_'.join(['Protein-H'] + all_resids)

    if couple_group not in index_list:
        if not make_group_ndx(couple_group_ind, wdir, bash_log=bash_log):
            return None
        index_list = get_index(os.path.join(wdir, 'index.ndx'))

    cmd = (f'cd {wdir}; gmx genrestr -f solv_ions.gro -o posre.itp -fc 1000 1000 1000 -n index.ndx <<< {index_list.index(couple_group)}')

    if not run_check_subprocess(cmd, wdir, log=bash_log, env=env):
        return None
    return wdir

def add_restraints_to_topol(topol):
    string = '''; system1 position restraints
#ifdef POSRES
#include "posre.itp"
#endif'''
    if not check_if_info_already_added_to_topol(topol, string):
        edit_topology_file(topol_file=topol, pattern='\n[ moleculetype ]', add=string, count_pattern=2)



