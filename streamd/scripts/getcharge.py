import argparse
import sys
from rdkit import Chem
from rdkit.Chem import rdmolops


def main(fname):
    mol = Chem.MolFromMolFile(fname)
    if mol:
        charge = rdmolops.GetFormalCharge(mol)
    else:
        sys.stderr.write('Molecule from file {} cannot be parsed'.format(fname))
        return None
    return charge


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''Returns the formal charge for the molecule using RDKiT''')
    parser.add_argument('-i', '--input', metavar='FILENAME', required=True,
                        help='input file with compound. Supported formats: *.mol.')
    parser.add_argument('-o', '--output', metavar='FILENAME', required=False, default=None,
                        help='Output text file. If omitted output will be in stdout.')
    args = parser.parse_args()

    formal_charge = main(args.input)
    if args.output is None:
        sys.stdout.write(str(formal_charge))
    else:
        with open(args.output, 'a') as out:
            out.write('\t'.join([args.input, formal_charge]))
