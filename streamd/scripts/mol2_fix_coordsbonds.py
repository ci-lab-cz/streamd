import argparse


def main(filemol, filemol2, output):
    # work with text instead of mol because mol2 is incorrect
    with open(filemol) as inp:
        mol_data = inp.readlines()
    with open(filemol2) as inp2:
        mol2_data = inp2.readlines()
    # mol file field spec
    mol_atom_field_len = int(mol_data[3][:3])
    mol_bond_field_len = int(mol_data[3][3:6])

    mol2_atom_field_len = int(mol2_data[2][:5])
    mol2_bond_field_len = int(mol2_data[2][5:11])

    new_mol2_field_spec = str(mol_atom_field_len).rjust(5)+str(mol_bond_field_len).rjust(6)+mol2_data[2][11:]

    # indexes
    mol_atom_startindex = 4
    mol_atom_endindex = 4 + mol_atom_field_len
    mol_bond_startindex = mol_atom_endindex
    mol_bond_endindex = mol_bond_startindex+mol_bond_field_len
    # mol2 indexes
    mol2_atom_startindex = 8
    mol2_atom_endindex = 8+mol2_atom_field_len
    mol2_bond_startindex = mol2_atom_endindex+1
    mol2_bond_endindex = mol2_bond_startindex+mol2_bond_field_len
    # mol2_atom_startindex = mol2_data.index('@<TRIPOS>ATOM\n')+1
    # mol2_atom_endindex = mol2_data.index('@<TRIPOS>BOND\n')
    # mol2_bond_startindex = mol2_data.index('@<TRIPOS>BOND\n')+1
    # mol2_bond_endindex = mol2_data.index('@<TRIPOS>SUBSTRUCTURE\n')

    # fix coords
    fix_coords = []
    for line1, line2 in zip(mol_data[mol_atom_startindex:mol_atom_endindex], mol2_data[mol2_atom_startindex:mol2_atom_endindex]):
        x_mol, y_mol, z_mol = line1[:10].strip(), line1[10:20].strip(), line1[20:30].strip()
        # 11 symbols #x_mol2, y_mol2, z_mol2 = line2[16:27].strip(), line2[27:38].strip(), line2[38:49].strip()
        x_mol2, y_mol2, z_mol2 = x_mol.rjust(11), y_mol.rjust(11), z_mol.rjust(11)
        new_line2 = line2[:16]+ x_mol2+y_mol2+z_mol2+line2[49:]
        fix_coords.append(new_line2)

    # Bond fix
    fix_bonds = []
    for n, line1 in enumerate(mol_data[mol_bond_startindex:mol_bond_endindex], 1):
        atom1_mol, atom2_mol, bond_type_mol, bond_stereo = line1[:3].strip(), line1[3:6].strip(),\
                                                       line1[6:9].strip(), line1[9:].strip()
        new_line2 = str(n).rjust(6)+ atom1_mol.rjust(6)+atom2_mol.rjust(6)+bond_type_mol.rjust(2).ljust(5)+'\n'
        fix_bonds.append(new_line2)

    new_mol2 = ''.join([*mol2_data[:2], new_mol2_field_spec, *mol2_data[3:8],  *fix_coords,
                        *mol2_data[mol2_atom_endindex:mol2_bond_startindex], *fix_bonds, *mol2_data[mol2_bond_endindex:]])
    if output is None:
        output = filemol2[:-5]+'_correct.mol2'
    with open(output, 'w') as out:
        out.write(new_mol2)


    # use dict because of ar mol2 type
    # mol_bond, mol2_bond = {}, {}
    # for line1 in mol_data[mol_bond_startindex:mol_bond_endindex]:
    #     atom1_mol, atom2_mol, bond_type_mol, bond_stereo = line1[:3].strip(), line1[3:6].strip(),\
    #                                                        line1[6:9].strip(), line1[9:].strip()
    #     mol_bond[tuple(sorted([atom1_mol, atom2_mol]))] = bond_type_mol
    #
    # for line2 in mol2_data[mol2_bond_startindex:mol2_bond_endindex]:
    #     bond_id_mol2, atom1_mol2, atom2_mol2, bond_type_mol2 = line2[:6].strip(), line2[6:12].strip(), \
    #                                                        line2[12:18].strip(), line2[18:24].strip()
    #     mol2_bond[tuple(sorted([atom1_mol2, atom2_mol2]))] = bond_type_mol2
    #
    # new_bonds = []
    # # ar type in mol2
    # n = 1
    # for bond, bond_type in mol_bond.items():
    #     if bond in mol2_bond:
    #         if mol2_bond[bond] == 'ar':# and mol_bond[bond] == '2':
    #             bond_type = ' ar'
    #     atom1, atom2 = bond
    #     new_bond_line = str(n).rjust(6)+ atom1.rjust(6)+atom2.rjust(6)+bond_type.rjust(2).ljust(5)+'\n'
    #     new_bonds.append(new_bond_line)
    #     n+=1


    # new_mol2 = ''.join([*mol2_data[:mol2_startindex],*new_coords,*mol2_data[mol2_endindex:]])
    # with open(filemol2.split('.mol2')[0]+'_correctcoords.mol2', 'w') as output:
    #     output.write(new_mol2)
    # print('Finish')

# def change_bonds(filemol, filemol2):
#     with open(filemol) as inp:
#         mol_data = inp.readlines()
#     with open(filemol2) as out:
#         mol2_data = out.readlines()
#     # new_mol2 = ''.join([*mol2_data[:mol2_startindex], *new_bonds, *mol2_data[mol2_endindex:]])
    # with open(filemol2.split('.mol2')[0]+'_correctbond.mol2', 'w') as output:
    #     output.write(new_mol2)
    # print('Finish')



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''Returns the Gaussian calculation files''')
    parser.add_argument('--mol', metavar='FILENAME', required=True,
                        help='input file with compound. Supported formats: *.mol.')
    parser.add_argument('--mol2', metavar='FILENAME', required=True,
                        help='input file with compound. Supported formats: *.mol2.')
    parser.add_argument('-o', '--output', metavar='FILENAME', required=False, default=None,
                        help='Output text file')
    args = parser.parse_args()

    main(filemol=args.mol,filemol2=args.mol2,output=args.output)


