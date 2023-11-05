from glob import glob
import logging
import os
import shutil

import MDAnalysis as mda
import parmed as pmd

from streamd.utils.utils import run_check_subprocess
from streamd.utils.utils import get_mol_resid_pair
from streamd.preparation.ligand_preparation import prepare_gaussian_files

def convert_pdb2mol2(metal_pdb, charge_dict):
    '''

    :param metal_pdb: metal pdb file path
    :param charge_dict: {MN:2, ZN:2, CA:2}
    :return: metal mol2 file path
    '''
    atom = os.path.basename(metal_pdb).split('_')[0]
    charge = charge_dict.get(atom, None)
    if not charge:
        logging.error(f'MCPBPY Preparation Error: Cannot find charge for {atom} Metal atom. Please set up charge for the Metal atom. '
                      f'Calculation will be stopped')
        return None

    metal_mol2 = metal_pdb.replace(".pdb",".mol2")

    cmd = f'metalpdb2mol2.py -i {metal_pdb} -o {metal_mol2} -c {charge}'
    if not run_check_subprocess(cmd, metal_pdb):
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


def copy_rename_ligand_files(all_lig_wdirs, dir_to_copy):
    '''
    mol2 and frcmod file names should correspond residue id
    :return: dict
    '''
    def copy_molid_to_resid_file(lig_wdir_list, wdir, ext_list=('mol2','frcmod')):
        new_file_ext_dict = {i: [] for i in ext_list}
        for lig_wdir in lig_wdir_list:
            molid, resid = list(get_mol_resid_pair(os.path.join(lig_wdir, 'resid.txt')))[0]
            cur_file_name = os.path.join(lig_wdir, molid)
            new_file_name = os.path.join(wdir, resid)

            for ext in ext_list:
                shutil.copy(f'{cur_file_name}.{ext}', f'{new_file_name}.{ext}')
                new_file_ext_dict[ext].append(f'{new_file_name}.{ext}')

        return new_file_ext_dict

    lig_new_file_ext_dict = copy_molid_to_resid_file(all_lig_wdirs, dir_to_copy)

    return lig_new_file_ext_dict


def prepare_protein_in(file_in_template, file_out, complex_file,variable_ion_ids, variable_ion_mol2_files,
                       variable_ligand_mol2_files, variable_ligand_frcmod_files, gaussian_version, force_field, cut_off=2.8, large_opt=1):
    with open(file_in_template) as inp:
        data = inp.read()

    new_data = data.replace('complex.pdb', complex_file)\
        .replace('variable_ion_ids', variable_ion_ids)\
        .replace('variable_ion_mol2_files', variable_ion_mol2_files)\
        .replace('variable_ligand_mol2_files', variable_ligand_mol2_files)\
        .replace('variable_ligand_frcmod_files', variable_ligand_frcmod_files)\
        .replace('gaussian_version', gaussian_version)\
        .replace('cut_off cut_off', f'cut_off {cut_off}')\
        .replace('large_opt large_opt', f'large_opt {large_opt}')\
        .replace('force_field force_field', f'force_field {force_field}')

    with open(file_out, 'w') as output:
        output.write(new_data)

    return file_out


def set_up_gaussian_files(wdir, ncpu, gaussian_basis, gaussian_memory):
    gaussian_files = glob(os.path.join(wdir, 'protein*.com'))
    for gau_file in gaussian_files:
        prepare_gaussian_files(file_template=gau_file, file_out=gau_file, ncpu=ncpu,
                               gaussian_basis=gaussian_basis, gaussian_memory=gaussian_memory)


def run_MCPBPY(protein_in_file, wdir, s, bash_log):
    cmd = f'cd {wdir}; MCPB.py -i {protein_in_file} -s {s} >> {bash_log} 2>&1'
    if not run_check_subprocess(cmd, wdir):
        return None
    return wdir


def run_gaussian_calculation(wdir, gaussian_version, activate_gaussian, bash_log):
    def run_task(gau_cmd, activate_gaussian, log, wdir, bash_log, check_only_if_exist=False):
        def check_gau_log_file(log):
            if os.path.isfile(log):
                with open(log) as checkpoint:
                    gau_log = checkpoint.read()
                if 'Normal termination of' in gau_log:
                    return True
            return False

        if (check_only_if_exist and not os.path.isfile(log)) or not check_gau_log_file(small_opt_log):
            cmd = (f'cd {wdir}; {activate_gaussian};'
                   f'{gau_cmd}'
                   f' >> {bash_log} 2>&1')
            if not run_check_subprocess(cmd, wdir):
                return None
        else:
            logging.info(f'MCPB.py gaussian calculation was finished: {log}. Skip this step.')
        return log

    small_opt_log = os.path.join(wdir, 'protein_small_opt.log')
    small_fc_log = os.path.join(wdir, 'protein_small_fc.log')
    small_opt_fchk = os.path.join(wdir, 'protein_small_opt.fchk')
    large_mk_log = os.path.join(wdir, 'protein_large_mk.log')

    if not run_task(gau_cmd=f'{gaussian_version} < protein_small_opt.com > protein_small_opt.log', log=small_opt_log,
                    activate_gaussian=activate_gaussian, bash_log=bash_log, wdir=wdir):
        return None
    if not run_task(gau_cmd=f'{gaussian_version} < protein_small_fc.com > protein_small_fc.log', log=small_fc_log,
                    activate_gaussian=activate_gaussian, bash_log=bash_log, wdir=wdir):
        return None
    if not run_task(gau_cmd=f'formchk protein_small_opt.chk protein_small_opt.fchk', log=small_opt_fchk,
                    activate_gaussian=activate_gaussian, bash_log=bash_log, wdir=wdir, check_only_if_exist=True):
        return None
    if not run_task(gau_cmd=f'{gaussian_version} < protein_large_mk.com > protein_large_mk.log', log=large_mk_log,
                    activate_gaussian=activate_gaussian, bash_log=bash_log, wdir=wdir):
        return None

    return wdir


def run_tleap(wdir, bash_log):
    cmd = f'cd {wdir}; tleap -s -f protein_tleap.in > protein_tleap.out >> {bash_log} 2>&1'
    if not run_check_subprocess(cmd, wdir):
        return None
    return wdir


def amber2gmx(wdir):
    parm = pmd.load_file(os.path.join(wdir, 'protein_solv.prmtop'), os.path.join(wdir, 'protein_solv.inpcrd'))
    parm.save(os.path.join(wdir, 'topol.top'), format='gromacs')
    parm.save(os.path.join(wdir, 'solv_ions.gro'))

    return wdir