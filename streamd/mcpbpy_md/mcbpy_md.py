import os
from streamd.preparation import mcpbpy_preparation


def main(protein, metal_resnames, metal_charges, wdir_metal,
         system_lig_wdirs, lig_wdir, wdir_md, script_path, ncpu, activate_gaussian,
         gaussian_version, gaussian_basis, gaussian_memory, bash_log, cut_off=2.8):

    # split
    protein_clean_pdb, metal_pdb_list = mcpbpy_preparation.split_metal(protein_fname=protein,
                                                                       metal_resnames=metal_resnames, wdir=wdir_metal)
    # pdb to mol2
    metal_mol2_list = []
    for metal_pdb in metal_pdb_list:
        metal_mol2 = mcpbpy_preparation.convert_pdb2mol2(metal_pdb=metal_pdb, charge_dict=metal_charges)
        if metal_mol2 is not None:
            metal_mol2_list.append(metal_mol2)
        else:
            return None

    # sort mol2 by atom id - order important
    metal_mol2_list.sort(key=lambda x: int(os.path.basename(x).strip('.mol2').split('_')[1]))
    metal_atomid_list = [os.path.basename(i).strip('.mol2').split('_')[1] for i in metal_mol2_list]

    # merge
    # copy molid.mol2, molid.frcmod as resid for MCPBPY usage
    all_lig_new_file_ext_dict = mcpbpy_preparation.copy_rename_ligand_files(system_lig_wdirs+[lig_wdir], wdir_md)

    complex_file = mcpbpy_preparation.merge_complex(protein_clean_pdb,
                                                    ligand_mol2_list=all_lig_new_file_ext_dict['mol2'],
                                                    metal_mol2_list=metal_mol2_list, wdir=wdir_md)
    # prepapre protein.in
    protein_in_file = mcpbpy_preparation.prepare_protein_in(file_in_template=os.path.join(script_path, 'mcpbpy_scripts', 'protein.in'),
                       file_out=os.path.join(wdir_md, 'protein.in'),
                       complex_file=complex_file,
                       variable_ion_ids=' '.join(metal_atomid_list),
                       variable_ion_mol2_files=' '.join(metal_mol2_list),
                        # path/UNL.mol2 is not reconised
                       variable_ligand_mol2_files=' '.join([os.path.basename(i) for i in all_lig_new_file_ext_dict['mol2']]),
                       variable_ligand_frcmod_files=' '.join([os.path.basename(i) for i in all_lig_new_file_ext_dict['frcmod']]),
                       gaussian_version=gaussian_version,
                       force_field='ff99SB',
                       cut_off=cut_off)

    if not mcpbpy_preparation.run_MCPBPY(protein_in_file=protein_in_file, wdir=wdir_md, s=1, bash_log=bash_log):
        return None

    mcpbpy_preparation.set_up_gaussian_files(wdir=wdir_md, ncpu=ncpu, gaussian_basis=gaussian_basis, gaussian_memory=gaussian_memory)
    if not mcpbpy_preparation.run_gaussian_calculation(wdir=wdir_md, gaussian_version=gaussian_version, activate_gaussian=activate_gaussian, bash_log=bash_log):
        return None

    if not mcpbpy_preparation.run_MCPBPY(protein_in_file=protein_in_file, wdir=wdir_md, s=2, bash_log=bash_log):
        return None
    if not mcpbpy_preparation.run_MCPBPY(protein_in_file=protein_in_file, wdir=wdir_md, s=3, bash_log=bash_log):
        return None
    if not mcpbpy_preparation.run_MCPBPY(protein_in_file=protein_in_file, wdir=wdir_md, s=4, bash_log=bash_log):
        return None

    if not mcpbpy_preparation.run_tleap(wdir=wdir_md, bash_log=bash_log):
        return None

    #parmed -i mcpbpy_parmed.in -p 1OKL_solv.prmtop

    mcpbpy_preparation.amber2gmx(wdir=wdir_md)






