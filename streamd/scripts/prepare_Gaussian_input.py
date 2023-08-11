import argparse
from rdkit import Chem
from rdkit.Chem import rdmolops, Descriptors

def main(fname, opt_param, charges_param, freq_param, out):
    mol = Chem.MolFromMolFile(fname, removeHs=False)
    title_smi = Chem.MolToSmiles(mol)

    charge = rdmolops.GetFormalCharge(mol)
    spin = Descriptors.NumRadicalElectrons(mol) + 1

    coords_as_string_list = []
    # setting from template
    with open(opt_param) as inp_opt, open(charges_param) as inp_charg, open(freq_param) as inp_freq:
        opt_job_specification = inp_opt.read().strip()
        charg_job_specification = inp_charg.read().strip()
        freq_job_specification = inp_freq.read().strip()

    for i, atom in enumerate(mol.GetAtoms()):
        positions = mol.GetConformer().GetAtomPosition(i)
        coords_as_string_list.append(','.join([atom.GetSymbol(), str(positions.x), str(positions.y), str(positions.z)+'\n']))

    out_opt, out_charg, out_freq, gout = f'{out}_opt_gaussian.com', f'{out}_charg_gaussian.com', f'{out}_freq_gaussian.com',  f'{out}.gout'
    molecule_specification = f'  {charge}  {spin}\n{"".join(coords_as_string_list)}\n'
    with open(out_opt, 'w') as out1, open (out_charg, 'w') as out2, open (out_freq, 'w') as out3:
        out1.write(f'{opt_job_specification}\n\n{title_smi}\n\n{molecule_specification}\n')
        out2.write(f'{charg_job_specification}\n\n{title_smi}\n\n{molecule_specification}\n{gout}\n\n')
        out3.write(f'{freq_job_specification}\n\n{title_smi}\n\n{molecule_specification}\n')



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''Returns the Gaussian calculation files''')
    parser.add_argument('-i', '--input', metavar='FILENAME', required=True,
                        help='input file with compound. Supported formats: *.mol.')
    parser.add_argument('--opt_param', metavar='FILENAME', required=True,
                        help='input Gaussian parameters for optimization.')
    parser.add_argument('--charges_param', metavar='FILENAME', required=True,
                        help='input Gaussian parameters for charges calculation.')
    parser.add_argument('--freq_param', metavar='FILENAME', required=True,
                        help='input Gaussian parameters for frequency calculation. Futher will be used for force field calculation by CartHess2FC.py.')
    parser.add_argument('-o', '--output', metavar='FILENAME', required=True,
                        help='Output text file')
    args = parser.parse_args()

    main(fname=args.input, opt_param=args.opt_param, charges_param=args.charges_param, freq_param=args.freq_param, out=args.output)