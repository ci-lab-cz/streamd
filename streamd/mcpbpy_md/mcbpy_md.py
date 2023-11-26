import logging
import os
import shutil

from streamd.preparation.md_files_preparation import prepare_mdp_files
from streamd.preparation import mcpbpy_preparation
from streamd.utils.utils import get_mol_resid_pair


def main(wdir_var_ligand, protein_name, protein_file, metal_resnames, metal_charges,
         system_lig_wdirs, wdir_metal, wdir_md, script_path, ncpu, activate_gaussian,
         gaussian_version, gaussian_basis, gaussian_memory, bash_log,
         nvt_time_ps, npt_time_ps, mdtime_ns, cut_off=2.8):

    if wdir_var_ligand:
        var_lig_molid, var_lig_resid = next(get_mol_resid_pair(os.path.join(wdir_var_ligand, 'resid.txt')))
        wdir_md_cur = os.path.join(wdir_md, f'{protein_name}_{var_lig_molid}')
    else:
        wdir_md_cur = os.path.join(wdir_md, protein_name)

    os.makedirs(wdir_md_cur, exist_ok=True)
    bash_log_curr = os.path.join(wdir_md_cur, bash_log)

    # copy molid.mol2, molid.frcmod as resid for MCPBPY usage and get molid_resid pairs
    if system_lig_wdirs or wdir_var_ligand:
        all_lig_new_file_ext_dict, molids_pairs_dict = mcpbpy_preparation.copy_rename_ligand_files(
            system_lig_wdirs + [wdir_var_ligand], wdir_md_cur)
    else:
        all_lig_new_file_ext_dict, molids_pairs_dict = {'mol2': [], 'frcmod': []}, {}

    if not ( os.path.isfile(os.path.join(wdir_md_cur, 'protein_solv.prmtop')) and os.path.isfile(os.path.join(wdir_md_cur, 'protein_solv.inpcrd')) and
        os.path.isfile(os.path.join(wdir_md_cur, 'complex.pdb')) and os.path.isfile(os.path.join(wdir_md_cur, 'protein_mcpbpy.pdb'))):

        # split metal atoms and the protein
        protein_clean_pdb, metal_pdb_list = mcpbpy_preparation.split_metal(protein_fname=protein_file,
                                                                           metal_resnames=metal_resnames, wdir=wdir_metal)
        # metal pdb to mol2
        metal_mol2_list = []
        for metal_pdb in metal_pdb_list:
            metal_mol2 = mcpbpy_preparation.convert_pdb2mol2(metal_pdb=metal_pdb, charge_dict=metal_charges)
            if metal_mol2 is not None:
                metal_mol2_list.append(metal_mol2)
            else:
                return None

        # sort mol2 by atom id - order important for complex merge
        metal_mol2_list.sort(key=lambda x: int(os.path.basename(x).strip('.mol2').split('_')[1]))

        # merge protein, metal atoms, cofactors, ligand
        complex_file = mcpbpy_preparation.merge_complex(protein_clean_pdb,
                                                        ligand_mol2_list=all_lig_new_file_ext_dict['mol2'],
                                                        metal_mol2_list=metal_mol2_list, wdir=wdir_md_cur)
        # prepapre protein.in
        protein_in_file = os.path.join(wdir_md_cur, 'protein.in')
        if not os.path.isfile(protein_in_file):
            mcpbpy_preparation.prepare_protein_in(file_in_template=os.path.join(script_path, 'mcpbpy_scripts', 'protein.in'),
                                                                file_out=protein_in_file,
                                                                complex_file=complex_file,
                                                                variable_ion_ids=' '.join(metal_atomid_list),
                                                                variable_ion_mol2_files=' '.join(metal_mol2_list),
                                                                # path/UNL.mol2 is not recognised
                                                                variable_ligand_mol2_files=' '.join([os.path.basename(i) for i in all_lig_new_file_ext_dict['mol2']]),
                                                                variable_ligand_frcmod_files=' '.join([os.path.basename(i) for i in all_lig_new_file_ext_dict['frcmod']]),
                                                                gaussian_version=gaussian_version,
                                                                force_field='ff99SB',
                                                                cut_off=cut_off)
        else:
            logging.warning(
                f'INFO: MCPBPY procedure: {protein_in_file} already exist. Use already prepared {protein_in_file} file')
        if not (os.path.isfile(os.path.join(wdir_md_cur, 'protein_small_opt.com')) and
                os.path.isfile(os.path.join(wdir_md_cur, 'protein_small_fc.com')) and
                os.path.isfile(os.path.join(wdir_md_cur, 'protein_large_mk.com'))):
            if not mcpbpy_preparation.run_MCPBPY(protein_in_file=protein_in_file, wdir=wdir_md_cur, s=1, bash_log=bash_log_curr):
                return None
        else:
            logging.warning(
                f'INFO: MCPBPY procedure: {os.path.join(wdir_md_cur, "protein_small_opt.com")}'
                f'{os.path.join(wdir_md_cur, "protein_small_fc.com")}'
                f'{os.path.join(wdir_md_cur, "protein_large_mk.com")}'
                f' already exist. Use already prepared Gaussian files')


        mcpbpy_preparation.set_up_gaussian_files(wdir=wdir_md_cur, ncpu=ncpu, gaussian_basis=gaussian_basis, gaussian_memory=gaussian_memory)

        logging.warning(f'INFO: MCPBPY procedure: Start gaussian calculation: {wdir_md_cur}')
        if not mcpbpy_preparation.run_gaussian_calculation(wdir=wdir_md_cur, gaussian_version=gaussian_version, activate_gaussian=activate_gaussian, bash_log=bash_log_curr):
            return None
        logging.warning(f'INFO: MCPBPY procedure: Finish gaussian calculation successfully: {wdir_md_cur}')

        if not mcpbpy_preparation.run_MCPBPY(protein_in_file=protein_in_file, wdir=wdir_md_cur, s=2, bash_log=bash_log_curr):
            return None
        if not mcpbpy_preparation.run_MCPBPY(protein_in_file=protein_in_file, wdir=wdir_md_cur, s=3, bash_log=bash_log_curr):
            return None
        if not mcpbpy_preparation.run_MCPBPY(protein_in_file=protein_in_file, wdir=wdir_md_cur, s=4, bash_log=bash_log_curr):
            return None

        if not mcpbpy_preparation.run_tleap(wdir=wdir_md_cur, bash_log=bash_log_curr):
            logging.warning(f'MCPBPY procedure: tleap failed: {wdir_md_cur}. '
                            f'Try to remove protein connected hydrogens in {os.path.join(wdir_md_cur, "protein_mcpbpy.pdb")} and readd it by tleap')
            mcpbpy_preparation.remove_allHs_from_pdb(os.path.join(wdir_md_cur, 'protein_mcpbpy.pdb'))
            if not mcpbpy_preparation.run_tleap(wdir=wdir_md_cur, bash_log=bash_log_curr):
                return None
    else:
        logging.warning(f'INFO: MCPBPY procedure: {wdir_md_cur}. {os.path.isfile(os.path.join(wdir_md_cur, "protein_solv.prmtop"))}'
                        f'{os.path.isfile(os.path.join(wdir_md_cur, "protein_solv.inpcrd"))} already exist. '
                        f'Skip MCPBPY calculation steps.')

    if not(os.path.isfile(os.path.join(wdir_md_cur, 'solv_ions.gro')) and os.path.isfile(
            os.path.join(wdir_md_cur, 'topol.top'))):
        mcpbpy_preparation.amber2gmx(complex_original=os.path.join(wdir_md_cur, 'complex.pdb'),
                                                complex_mcpbpy=os.path.join(wdir_md_cur,'protein_mcpbpy.pdb'),
                                                prmtop=os.path.join(wdir_md_cur, 'protein_solv.prmtop'),
                                                inpcrd=os.path.join(wdir_md_cur, 'protein_solv.inpcrd'), wdir=wdir_md_cur)

    for mdp_fname in ['ions.mdp', 'minim.mdp']:
        mdp_file = os.path.join(script_path, 'mdp', mdp_fname)
        shutil.copy(mdp_file, wdir_md_cur)

    if not prepare_mdp_files(wdir_md_cur=wdir_md_cur, all_resids=list(molids_pairs_dict.values()),
                             script_path=os.path.join(script_path, 'mdp'), nvt_time_ps=nvt_time_ps,
                             npt_time_ps=npt_time_ps, mdtime_ns=mdtime_ns):
        return None

    return wdir_md_cur




