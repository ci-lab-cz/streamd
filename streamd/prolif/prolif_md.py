#!/usr/bin/env python3

import argparse
import MDAnalysis as mda
import prolif as plf
from rdkit import Chem
import pandas as pd
from multiprocessing import Pool
from functools import  partial
class RawTextArgumentDefaultsHelpFormatter(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass


def calc_plif(lig, n_frame, protein_group):
    #https://userguide.mdanalysis.org/stable/examples/analysis/custom_parallel_analysis.html
    protein_group.universe.trajectory[n_frame]
    # lig = Chem.AddHs(lig, addCoords=True)
    fp = plf.Fingerprint(
        ['Hydrophobic', 'HBDonor', 'HBAcceptor', 'Anionic', 'Cationic', 'CationPi', 'PiCation',
         'PiStacking', 'MetalAcceptor'])
    frame, time_ps = protein_group.universe.trajectory.frame, protein_group.universe.trajectory.time
    print("Frame: {0:5d}, Time: {1:8.3f} ps".format(frame, time_ps))
    fp.run_from_iterable([plf.Molecule.from_rdkit(lig)], plf.Molecule.from_mda(protein_group), n_jobs=1)
    df_tec = fp.to_dataframe()
    df_tec.columns = ['.'.join(item.strip().lower() for item in items[1:]) for items in df_tec.columns]
    df_tec = df_tec.reset_index()
    df_tec['Frame'] = frame
   # df_tec['Time_ps'] = time_ps
    if n_frame != protein_group.universe.trajectory.frame:
        print('!!!!PROBLEM', n_frame, protein_group.universe.trajectory.frame)
        raise ValueError
    return df_tec



def main():
    parser = argparse.ArgumentParser(description='Get protein-ligand interactions from MD trajectories using '
                                                 'ProLIF module.',
                                     formatter_class=RawTextArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--trajectory', metavar='FILENAME', required=True,
                        help='input trajectory file (XTC).')
    parser.add_argument('-t', '--topology', metavar='FILENAME', required=False,
                        help='input topology file (TPR).')
    parser.add_argument('-l', '--ligand', metavar='STRING', required=False, default='UNL',
                        help='residue name of a ligand in the input trajectory.')
    parser.add_argument('--ligand_traj', metavar='STRING', required=False, default=None,
                        help='ligand the input trajectory. SDF')
    parser.add_argument('--pdb', metavar='FILENAME', required=False, default=None,
                        help='input correct topology for the right numeration')
    parser.add_argument('-s', '--step', metavar='INTEGER', required=False, default=1, type=int,
                        help='step to take every n-th frame. ps')
    parser.add_argument('-a', '--append_protein_selection', metavar='STRING', required=False, default=None,
                        help='the string which will be concatenated to the protein selection atoms. '
                             'Example: "resname ZN or resname MG".')
    parser.add_argument('-o', '--output', metavar='FILENAME', required=True,
                        help='output text file name')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='print progress.')
    parser.add_argument('-c', default=1, type=int,
                        help='ncpu')
#   parser.add_argument('--renum_chain', action='store_true', default=False)


    args = parser.parse_args()
    u = mda.Universe(args.topology, args.trajectory)
    u.add_TopologyAttr('chainID')
 
    if args.append_protein_selection is None:
        protein_selection = 'protein'
    else:
        protein_selection = f'protein or {args.append_protein_selection}'
    prot = u.atoms.select_atoms(protein_selection)
    if args.pdb:
        u1 = mda.Universe(args.pdb)
        print('Use topology from pdb')

    # segements_starts = {}
    # for atom in prot.atoms:
    #     if '_chain' in atom.segid.lower():
    #         if atom.segid not in segements_starts:
    #             segements_starts[atom.segid] = atom.resnum - 1
    #         chain = atom.segid.split('_')[-1]
    #         atom.chainID = chain
    #         atom.residue.resid = atom.residue.resnum - segements_starts[atom.segid]
    #     else:
    #         atom.chainID = 'X'

    if args.ligand_traj is None:
        lig = u.atoms.select_atoms(f'resname {args.ligand}')
        fp = plf.Fingerprint(['Hydrophobic', 'HBDonor', 'HBAcceptor', 'Anionic', 'Cationic', 'CationPi', 'PiCation',
                              'PiStacking', 'MetalAcceptor'])
        fp.run(u.trajectory[::args.step], lig, prot, progress=args.verbose)
        df = fp.to_dataframe()
        df.columns = ['.'.join(item.strip().lower() for item in items[1:]) for items in df.columns]
        df = df.reindex(sorted(df.columns), axis=1)
        df.to_csv(args.output, sep='\t')
    else:
        lig_suppl = [i for i in Chem.SDMolSupplier(args.ligand_traj, removeHs=False)]
        if len(lig_suppl) != len(u.trajectory):
            print('Error. Lenght of ligand traj doesnot equal the lenght of full trajectory')
            raise ValueError
        with Pool(args.c) as worker_pool:
            result = worker_pool.starmap(partial(calc_plif, protein_group=prot), [(lig, ts.frame) for lig, ts in zip(lig_suppl[::args.step], u.trajectory[::args.step])])

        df = pd.concat(result).fillna(False)
        df = df.reindex(sorted(df.columns), axis=1)
        df.to_csv(args.output, sep='\t', index=False)


if __name__ == '__main__':
    main()
