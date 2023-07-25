import itertools
import logging
import os
import subprocess

from rdkit import Chem


def supply_mols(fname, preset_resid=None):
    def generate_resid():
        ascii_uppercase_digits = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789'
        for i in itertools.product(ascii_uppercase_digits, repeat=3):
            resid = ''.join(i)
            if resid != 'UNL':
                yield resid

    def add_ids(mol, n, input_fname, resid):
        mol.SetProp('resid', resid)
        if not mol.HasProp('_Name'):
            mol.SetProp('_Name', f'{input_fname}_{n}')
        return mol

    resid_generator = generate_resid()

    if fname.endswith('.sdf'):
        for n, mol in enumerate(Chem.SDMolSupplier(fname, removeHs=False)):
            if mol:
                if preset_resid is None:
                    resid = next(resid_generator)
                else:
                    resid = preset_resid

                mol = add_ids(mol, n, input_fname=os.path.basename(fname).strip('.sdf'),
                              resid=resid)

                yield mol

    if fname.endswith('.mol'):
        mol = Chem.MolFromMolFile(fname, removeHs=False)
        if mol:
            if preset_resid is None:
                resid = next(resid_generator)
            else:
                resid = preset_resid
            mol = add_ids(mol, n=1, input_fname=os.path.basename(fname).strip('.mol'),
                          resid=resid)
            yield mol


def make_all_itp(fileitp_list, out_file):
    atom_type_list = []
    start_columns = None
    # '[ atomtypes ]\n; name    at.num    mass    charge ptype  sigma      epsilon\n'
    for f in fileitp_list:
        with open(f) as input:
            data = input.read()
        start = data.find('[ atomtypes ]')
        end = data.find('[ moleculetype ]') - 1
        atom_type_list.extend(data[start:end].split('\n')[2:])
        if start_columns is None:
            start_columns = data[start:end].split('\n')[:2]
        new_data = data[:start] + data[end + 1:]
        with open(f, 'w') as itp_ouput:
            itp_ouput.write(new_data)

    atom_type_uniq = [i for i in set(atom_type_list) if i]
    with open(out_file, 'w') as ouput:
        ouput.write('\n'.join(start_columns) + '\n')
        ouput.write('\n'.join(atom_type_uniq) + '\n')


def prep_ligand(mol, script_path, project_dir, wdir_ligand):
    molid = mol.GetProp('_Name')
    resid = mol.GetProp('resid')

    wdir_ligand_cur = os.path.join(wdir_ligand, molid)
    os.makedirs(wdir_ligand_cur, exist_ok=True)

    if os.path.isfile(os.path.join(wdir_ligand_cur, f'{molid}.itp')) and os.path.isfile(
            os.path.join(wdir_ligand_cur, f'posre_{molid}.itp')) \
            and os.path.isfile(os.path.join(wdir_ligand_cur, 'resid.txt')):
        logging.warning(
            f'{molid}.itp and posre_{molid}.itp files already exist.'
            f'Mol preparation step will be skipped for such molecule\n')
        return wdir_ligand_cur

    mol_file = os.path.join(wdir_ligand_cur, f'{molid}.mol')
    # if addH:
    mol = Chem.AddHs(mol, addCoords=True)
    Chem.MolToMolFile(mol, mol_file)

    try:
        subprocess.check_output(
            f'script_path={script_path} lfile={mol_file} input_dirname={wdir_ligand_cur} '
            f'name={resid} bash {os.path.join(project_dir, "scripts/script_sh/ligand_prep.sh")}',
            shell=True)
    except subprocess.CalledProcessError as e:
        logging.exception(f'{molid}\n{e}', stack_info=True)
        return None

    # create log for molid resid corresponding
    with open(os.path.join(wdir_ligand_cur, 'resid.txt'), 'w') as out:
        out.write(f'{molid}\t{resid}\n')
    # fname = os.path.join(wdir_ligand, 'all_resid.txt')
    # with daskLock(fname):
    #     with open(fname, 'a') as f:
    #         f.write(f'{molid}\t{resid}\n')

    return wdir_ligand_cur
