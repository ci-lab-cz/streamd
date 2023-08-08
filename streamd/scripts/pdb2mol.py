import argparse
import sys
from rdkit import Chem
from rdkit.Chem import AllChem
import os

def read_smi(fname, sep="\t"):
    with open(fname) as f:
        for line in f:
            items = line.strip().split(sep)
            if len(items) == 1:
                yield items[0], None
            else:
                yield items[0], items[1]


def main(fname, smi, mol_id=None, preserveH=True):
    mol_block = None
    if fname.endswith('.pdbqt'):
        with open(fname) as inp:
            pdbqt_block = inp.read()
        mol = Chem.MolFromPDBBlock('\n'.join([i[:66] for i in pdbqt_block.split('MODEL')[1].split('\n')]), removeHs=False,
                               sanitize=False)
    else:
        mol = Chem.MolFromPDBFile(fname, removeHs=False, sanitize=False)

    if mol:
        try:
            # if direction of polar Hs is not important
            if not preserveH:
                try:
                    mol = Chem.RemoveHs(mol)
                except Chem.rdchem.AtomValenceException:
                    sys.stderr.write('Warning. Could not remove Hs from original molecule')
                template_mol = Chem.RemoveHs(Chem.MolFromSmiles(smi))
            else:
                template_mol = Chem.AddHs(Chem.MolFromSmiles(smi))
            mol = AllChem.AssignBondOrdersFromTemplate(template_mol, mol)

            Chem.SanitizeMol(mol)
            # Chem.AssignStereochemistry(mol, cleanIt=True, force=True, flagPossibleStereoCenters=True)
            if mol_id is not None:
                mol.SetProp('_Name', mol_id)
            if not preserveH:
                mol = Chem.AddHs(mol, addCoords=True)
            mol_block = Chem.MolToMolBlock(mol)
        except Exception as e:
            sys.stderr.write(f'Could not assign bond orders while parsing PDB. {e} \n')
    return mol_block


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''Returns mol file from pdb. Use smi which was obtained from Marvin at pH=7.4''')
    parser.add_argument('-i', '--input', metavar='FILENAME', required=True,
                        help='input file with compound. Supported formats: *.pdb.')
    parser.add_argument('--smiles', metavar='fname', required=True,
                        help='input smile file as the template. Important: Use smi which was obtained from Marvin at pH=7.4')
    parser.add_argument('-o', '--output', metavar='FILENAME', required=False, default=None,
                        help='Output text file. If omitted output will be in stdout.')
    parser.add_argument('--preserveH', action='store_true', default=False,
                        help='Output text file. If omitted output will be in stdout.')

    args = parser.parse_args()

    molblock_fixed = None
    for smi, mol_id in read_smi(args.smiles,sep='\t'):
        if mol_id is not None:
            if os.path.basename(args.input)[:-4].lower() != mol_id.lower():
                continue
            sys.stdout.write(mol_id+': start preparation\n')
        molblock_fixed = main(args.input, smi=smi, mol_id=mol_id, preserveH=args.preserveH)
        break

    if molblock_fixed is None:
        sys.stderr.write('Error or Could not find the smi with the same mol_id as name of the pdb file. Termination\n')
        exit()

    if args.output is None:
        sys.stdout.write(str(molblock_fixed))
    else:
        with open(args.output, 'w') as output:
            output.write(molblock_fixed)


