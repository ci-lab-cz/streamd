import itertools
import logging
import os
import shutil
import re
from glob import glob
import re

import parmed as pmd

from rdkit import Chem
from rdkit.Chem import rdmolops

from streamd.utils.dask_init import init_dask_cluster, calc_dask
from streamd.utils.utils import run_check_subprocess

def reorder_hydrogens(mol):
    """
    Reorders the atoms in the molecule so that all hydrogens bonded to a heavy atom
    are listed immediately after that heavy atom.
    """
    new_order = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 1:  # Not a hydrogen
            new_order.append(atom.GetIdx())
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 1:  # Is a hydrogen
                    new_order.append(neighbor.GetIdx())
    mol = rdmolops.RenumberAtoms(mol, new_order)
    return mol

def supply_mols_tuple(fname, preset_resid=None, protein_resid_set=None):
    def generate_resid(protein_resid_list):
        ascii_uppercase_digits = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
        for i in itertools.product(ascii_uppercase_digits, repeat=3):
            resid = ''.join(i)
            if resid != 'UNL' and resid not in protein_resid_list:
                yield resid

    def add_ids(mol, n, input_fname, resid):
        mol.SetProp('resid', resid)
        if not mol.HasProp('_Name') or not mol.GetProp('_Name'):
            mol.SetProp('_Name', f'{input_fname}_{n}')
        return mol

    resid_generator = generate_resid(protein_resid_list=protein_resid_set)

    if fname.endswith('.sdf'):
        for n, mol in enumerate(Chem.SDMolSupplier(fname, removeHs=False), 1):
            if mol:
                if preset_resid is None:
                    resid = next(resid_generator)
                else:
                    resid = preset_resid

                mol = add_ids(mol, n, input_fname=os.path.basename(fname).rstrip('.sdf'), resid=resid)
                yield (mol, mol.GetProp('_Name'), resid)

    if fname.endswith('.mol') or fname.endswith('.mol2'):
        if fname.endswith('.mol2'):
            mol = pmd.load_file(fname).to_structure()
            mol = mol.rdkit_mol
        else:
            mol = Chem.MolFromMolFile(fname, removeHs=False)

        if mol:
            if preset_resid is None:
                resid = next(resid_generator)
            else:
                resid = preset_resid
            mol = add_ids(mol, n=1, input_fname=os.path.basename(fname).rstrip('.mol2').rstrip('.mol'), resid=resid)
            yield (mol, mol.GetProp('_Name'), resid)


def check_mols(fname):
    '''

    :param fname:
    :return: number of mols, list of problem_mols
    '''
    def check_if_problem(mol, n):
        molid = mol.GetProp('_Name') if mol.HasProp('_Name') else n
        try:
            Chem.SanitizeMol(mol)
        except:
            return molid
        return None
        
    problem_mols = []
    number_of_mols = 0
    if fname.endswith('.sdf'):
        n = 0
        for n, mol in enumerate(Chem.SDMolSupplier(fname, removeHs=False, sanitize=False), 1):
            problem_molid = check_if_problem(mol, n)
            if problem_molid:
                problem_mols.append(problem_molid)
        number_of_mols = n

    if fname.endswith('.mol') or fname.endswith('.mol2'):
        if fname.endswith('.mol2'):
            mol = pmd.load_file(fname).to_structure()
            mol = mol.rdkit_mol
        else:
            mol = Chem.MolFromMolFile(fname, removeHs=False, sanitize=False)
        problem_molid = check_if_problem(mol, 1)
        if problem_molid:
            problem_mols.append(problem_molid)
        number_of_mols = 1

    return number_of_mols, problem_mols


def make_all_itp(fileitp_input_list, fileitp_output_list, out_file):
    atom_type_list = []
    start_columns = None
    # '[ atomtypes ]\n; name    at.num    mass    charge ptype  sigma      epsilon\n'
    for itp_input, itp_output in zip(fileitp_input_list, fileitp_output_list):
        with open(itp_input) as input:
            data = input.read()
        start = data.find('[ atomtypes ]')
        end = data.find('[ moleculetype ]') - 1
        atom_type_list.extend(data[start:end].split('\n')[2:])
        if start_columns is None:
            start_columns = data[start:end].split('\n')[:2]
        new_data = data[:start] + data[end + 1:]
        with open(itp_output, 'w') as output:
            output.write(new_data)

    atom_type_uniq = [i for i in set(atom_type_list) if i]
    with open(out_file, 'w') as output:
        output.write('\n'.join(start_columns) + '\n')
        output.write('\n'.join(atom_type_uniq) + '\n')

def prepare_tleap(tleap_template, tleap, molid, conda_env_path):
    with open(tleap_template) as inp:
        data = inp.read()
    new_data = data.replace('env_path', conda_env_path).replace('ligand', molid)
    with open(tleap, 'w') as output:
        output.write(new_data)


def prepare_gaussian_files(file_template, file_out, ncpu, opt_restart=False, gaussian_basis='B3LYP/6-31G*',
                           gaussian_memory='60GB'):
    with open(file_template) as inp:
        data = inp.read()
    standard_basis = 'B3LYP/6-31G\*'
    new_data = re.sub('%NProcShared=[0-9]*', f'%NProcShared={ncpu}', data)
    new_data = re.sub('%Mem=[0-9a-zA-Z]*', f'%Mem={gaussian_memory}', new_data)
    if 'SCF=' not in new_data:
        new_data = re.sub(standard_basis, f'{gaussian_basis} SCF=XQC', new_data)
    else:
        new_data = re.sub(standard_basis, gaussian_basis, new_data)
    if opt_restart and 'Opt' in new_data and 'Opt=Restart' not in new_data:
        new_data = re.sub('Opt', 'Opt=Restart', new_data)

    with open(file_out, 'w') as output:
        output.write(new_data)


def prep_ligand(mol_tuple, script_path, project_dir, wdir_ligand, conda_env_path, bash_log, gaussian_exe=None,
                activate_gaussian=None, gaussian_basis='B3LYP/6-31G*', gaussian_memory='60GB', ncpu=1,
                mol2_file=None, env=None):
    mol, molid, resid = mol_tuple

    wdir_ligand_cur = os.path.join(wdir_ligand, molid)
    os.makedirs(wdir_ligand_cur, exist_ok=True)

    if os.path.isfile(os.path.join(wdir_ligand_cur, f'{molid}.itp')) and os.path.isfile(
            os.path.join(wdir_ligand_cur, f'posre_{molid}.itp')):
        logging.warning(
            f'{molid}.itp and posre_{molid}.itp files already exist. '
            f'Mol preparation step will be skipped for such molecule')
        if not os.path.isfile(os.path.join(wdir_ligand_cur, 'resid.txt')):
            with open(os.path.join(wdir_ligand_cur, 'resid.txt'), 'w') as out:
                out.write(f'{molid}\t{resid}\n')

        return wdir_ligand_cur

    if not mol2_file or not os.path.isfile(mol2_file):
        mol2_file = os.path.join(wdir_ligand_cur, f'{molid}.mol2')
        mol_file = os.path.join(wdir_ligand_cur, f'{molid}.mol')
        # if addH:
        mol = Chem.AddHs(mol, addCoords=True)
        # reorder hydrogens for Gromacs GPU update functionality
        mol = reorder_hydrogens(mol)

        Chem.MolToMolFile(mol, mol_file)

        charge = rdmolops.GetFormalCharge(mol)

        # generate mol2
        if not os.path.isfile(mol2_file):
            # boron atom
            if mol.HasSubstructMatch(Chem.MolFromSmarts("[#5]")):
                if gaussian_exe:
                    for file in glob(os.path.join(script_path, 'com', '*.com')):
                        prepare_gaussian_files(file_template=file,
                                               file_out=os.path.join(wdir_ligand_cur, os.path.basename(file)),
                                               ncpu=ncpu,
                                               gaussian_basis=gaussian_basis, gaussian_memory=gaussian_memory)
                    cmd = f'script_path={script_path} lfile={mol_file} input_dirname={wdir_ligand_cur} ' \
                          f'resid={resid} molid={molid} charge={charge} gaussian_version={gaussian_exe} ' \
                          f'activate_gaussian="{activate_gaussian if activate_gaussian else ""}" ' \
                          f'bash {os.path.join(project_dir, "scripts/script_sh/ligand_mol2prep_by_gaussian.sh")} ' \
                          f' >> {os.path.join(wdir_ligand_cur, bash_log)} 2>&1'
                    if not run_check_subprocess(cmd, molid, log=os.path.join(wdir_ligand_cur, bash_log), env=env):
                        return None
                else:
                    return None
            else:
                cmd = f'script_path={script_path} lfile={mol_file} input_dirname={wdir_ligand_cur} ' \
                      f'resid={resid} molid={molid} charge={charge} bash {os.path.join(project_dir, "scripts/script_sh/ligand_mol2prep.sh")} ' \
                      f' >> {os.path.join(wdir_ligand_cur, bash_log)} 2>&1',
                if not run_check_subprocess(cmd, molid, log=os.path.join(wdir_ligand_cur, bash_log), env=env):
                    return None
    else:
        mol2 = pmd.load_file(mol2_file).to_structure()
        mol2.residues[0].name = 'UNL'
        mol2.save(os.path.join(wdir_ligand_cur, f'{molid}.mol2'))
        logging.info(f'No mol2 file will be generated. {mol2_file} will be used instead')

    prepare_tleap(os.path.join(script_path, 'tleap.in'), tleap=os.path.join(wdir_ligand_cur, 'tleap.in'),
                  molid=molid, conda_env_path=conda_env_path)
    cmd = f'script_path={script_path} input_dirname={wdir_ligand_cur} ' \
          f'molid={molid} bash {os.path.join(project_dir, "scripts/script_sh/ligand_prep.sh")} ' \
          f' >> {os.path.join(wdir_ligand_cur, bash_log)} 2>&1'
    if not run_check_subprocess(cmd, molid, log=os.path.join(wdir_ligand_cur, bash_log), env=env):
        return None

    # create log for molid resid corresponding
    with open(os.path.join(wdir_ligand_cur, 'resid.txt'), 'w') as out:
        out.write(f'{molid}\t{resid}\n')
    return wdir_ligand_cur


def prepare_input_ligands(ligand_fname, preset_resid, protein_resid_set, script_path, project_dir, wdir_ligand,
                          gaussian_exe, activate_gaussian, gaussian_basis, gaussian_memory,
                          hostfile, ncpu, bash_log):
    '''

    :param ligand_fname:
    :param preset_resid:
    :param script_path:
    :param project_dir:
    :param wdir_system_ligand:
    :param gaussian_exe: str or None
    :param activate_gaussian: str or None
    :param hostfile:
    :param ncpu:
    :param bash_log:
    :return:
    '''
    lig_wdirs = []

    if ligand_fname.endswith('.mol2'):
        mol_tuple = next(supply_mols_tuple(ligand_fname, preset_resid=preset_resid, protein_resid_set=protein_resid_set))
        res = prep_ligand(mol_tuple=mol_tuple,
                          script_path=script_path, project_dir=project_dir,
                          wdir_ligand=wdir_ligand, conda_env_path=os.environ["CONDA_PREFIX"],
                          gaussian_exe=gaussian_exe, activate_gaussian=activate_gaussian,
                          gaussian_basis=gaussian_basis, gaussian_memory=gaussian_memory,
                          ncpu=ncpu, bash_log=bash_log, mol2_file=ligand_fname)
        if res:
            lig_wdirs.append(res)

    else:
        standard_mols, boron_containing_mols = [], []
        for mol_tuple in supply_mols_tuple(ligand_fname, preset_resid=preset_resid, protein_resid_set=protein_resid_set):
            mol = mol_tuple[0]
            if mol.HasSubstructMatch(Chem.MolFromSmarts("[#5]")):
                boron_containing_mols.append(mol_tuple)
            else:
                standard_mols.append(mol_tuple)

        dask_client, cluster = None, None
        # prepare boron-containig mols
        if boron_containing_mols:
            if gaussian_exe:
                try:
                    dask_client, cluster = init_dask_cluster(hostfile=hostfile,
                                                             n_tasks_per_node=1,
                                                             use_multi_servers=True if len(boron_containing_mols) > 1 else False,
                                                             ncpu=ncpu)
                    for res in calc_dask(prep_ligand, boron_containing_mols, dask_client,
                                         script_path=script_path, project_dir=project_dir,
                                         wdir_ligand=wdir_ligand, conda_env_path=os.environ["CONDA_PREFIX"],
                                         gaussian_exe=gaussian_exe, activate_gaussian=activate_gaussian,
                                         gaussian_basis=gaussian_basis, gaussian_memory=gaussian_memory,
                                         ncpu=ncpu, bash_log=bash_log,
                                         env=os.environ.copy()):
                        if res:
                            lig_wdirs.append(res)
                finally:
                    if dask_client:
                        dask_client.shutdown()
                    if cluster:
                        cluster.close()
            else:
                logging.warning(
                    f'There are molecules from {ligand_fname} which have Boron atom and to prepare such molecules you need to set up Gaussian.'
                    f' Please restart the run again and use --gaussian_exe arguments')

        if standard_mols:
            try:
                dask_client, cluster = init_dask_cluster(hostfile=hostfile,
                                                         n_tasks_per_node=min(ncpu, len(standard_mols)),
                                                         use_multi_servers= True if len(standard_mols) > ncpu else False,
                                                         ncpu=ncpu)
                for res in calc_dask(prep_ligand, standard_mols, dask_client,
                                     script_path=script_path, project_dir=project_dir,
                                     wdir_ligand=wdir_ligand, conda_env_path=os.environ["CONDA_PREFIX"],
                                     ncpu=ncpu, bash_log=bash_log,
                                     env=os.environ.copy()):
                    if res:
                        lig_wdirs.append(res)
            finally:
                if dask_client:
                    dask_client.shutdown()
                if cluster:
                    cluster.close()

    return lig_wdirs
