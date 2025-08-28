"""Utilities for constructing protein–ligand complexes for MD simulations."""

import os
import shutil
import logging

from streamd.preparation.ligand_preparation import make_all_itp
from streamd.preparation.md_files_preparation import prep_md_files, add_ligands_to_topol, \
    edit_topology_file, prepare_mdp_files, check_if_info_already_added_to_topol
from streamd.utils.utils import run_check_subprocess


def complex_preparation(protein_gro, ligand_gro_list, out_file):
    """Merge protein and ligand GRO files into a single complex structure."""
    atoms_list = []
    with open(protein_gro) as input:
        prot_data = input.readlines()
        atoms_list.extend(prot_data[2:-1])

    for f in ligand_gro_list:
        with open(f) as input:
            data = input.readlines()
        atoms_list.extend(data[2:-1])

    n_atoms = len(atoms_list)
    with open(out_file, 'w') as output:
        output.write(prot_data[0])
        output.write(f'{n_atoms}\n')
        output.write(''.join(atoms_list))
        output.write(prot_data[-1])


def run_complex_preparation(wdir_var_ligand, wdir_system_ligand_list,
                            protein_name, wdir_protein, wdir_md, script_path, project_dir,
                            mdtime_ns, npt_time_ps, nvt_time_ps, clean_previous, seed, bash_log,
                            box_type='cubic', box_distance=1.0,
                            mdp_dir=None, explicit_args=(), env=None):
    """Prepare all input files for MD by combining protein and ligand data.

    :param box_type: Shape of the simulation box used by ``gmx editconf``.
    :param box_distance: Minimum distance in nm between solute and box edge.
    """
    wdir_md_cur, md_files_dict = prep_md_files(wdir_var_ligand=wdir_var_ligand, protein_name=protein_name,
                                               wdir_system_ligand_list=wdir_system_ligand_list,
                                               wdir_protein=wdir_protein,
                                               wdir_md=wdir_md, clean_previous=clean_previous)

    mdp_files_default = ['ions.mdp', 'minim.mdp', 'nvt.mdp', 'npt.mdp', 'md.mdp']

    protein_gro = os.path.join(wdir_protein, f'{protein_name}.gro')
    # ligands and cofactors
    if md_files_dict['itp']:
        # make all itp and edit itps
        if not os.path.isfile(os.path.join(wdir_md_cur, "all.itp")) or not all([os.path.isfile(os.path.join(wdir_md_cur, i)) for i in md_files_dict['itp']]):
            make_all_itp(fileitp_input_list=md_files_dict['itp_orig'], fileitp_output_list=[os.path.join(wdir_md_cur, j) for j in md_files_dict['itp']], out_file=os.path.join(wdir_md_cur, 'all.itp'))
        else:
            logging.warning(f'{wdir_md_cur}. Prepared itp files exist. Skip ligand all itp preparation step')

        # add ligands info to topology if there is no necessary row
        add_ligands_to_topol(md_files_dict['itp'], md_files_dict['posres'], md_files_dict['resid'],
                             topol=os.path.join(wdir_md_cur, "topol.top"))
        if not check_if_info_already_added_to_topol(os.path.join(wdir_md_cur, "topol.top"),
                                                    '; Include all topology\n#include "all.itp"\n'):
            edit_topology_file(topol_file=os.path.join(wdir_md_cur, "topol.top"),
                               pattern="; Include forcefield parameters",
                               add='; Include all topology\n#include "all.itp"\n',
                               how='after', n=3)

        # create file with molid resid for each ligand in the current system
        with open(os.path.join(wdir_md_cur, 'all_ligand_resid.txt'), 'w') as out:
            molid_resid_pairs = ['\t'.join(map(str, i)) for i in zip(md_files_dict['molid'], md_files_dict['resid'])]
            out.write('\n'.join(molid_resid_pairs))

    # complex
    if not os.path.isfile(os.path.join(wdir_md_cur, 'complex.gro')):
        complex_preparation(protein_gro=protein_gro,
                            ligand_gro_list=md_files_dict['gro'],
                            out_file=os.path.join(wdir_md_cur, 'complex.gro'))
    else:
        logging.warning(f'{wdir_md_cur}. Prepared complex file exists. Skip complex preparation step')


    user_mdp_files = []
    for mdp_fname in mdp_files_default:
        if mdp_dir and os.path.isfile(os.path.join(mdp_dir, mdp_fname)):
            mdp_file = os.path.join(mdp_dir, mdp_fname)
            user_mdp_files.append(mdp_fname)
            logging.warning(f'Use user provided mdp file: {mdp_file}')
        else:
            mdp_file = os.path.join(script_path, mdp_fname)

        shutil.copy(mdp_file, wdir_md_cur)

    # if not mdp_dir:
    #     for mdp_fname in ['ions.mdp', 'minim.mdp']:
    #         mdp_file = os.path.join(script_path, mdp_fname)
    #         shutil.copy(mdp_file, wdir_md_cur)

    if not os.path.isfile(os.path.join(wdir_md_cur, 'solv_ions.gro')):
        cmd = (f'wdir={wdir_md_cur} box_type={box_type} box_distance={box_distance} '
               f'bash {os.path.join(project_dir, "scripts/script_sh/solv_ions.sh")} '
               f'>> {os.path.join(wdir_md_cur, bash_log)} 2>&1')
        if not run_check_subprocess(cmd=cmd, key=wdir_md_cur, log=os.path.join(wdir_md_cur, bash_log), env=env):
            return None
    else:
        logging.warning(f'{wdir_md_cur}. Prepared solv_ions.gro file exists. Skip solvation and ion preparation step')

    if not prepare_mdp_files(wdir_md_cur=wdir_md_cur, all_resids=md_files_dict['resid'],
                             nvt_time_ps=nvt_time_ps,
                             npt_time_ps=npt_time_ps, mdtime_ns=mdtime_ns,
                             user_mdp_files=user_mdp_files,
                             explicit_args=explicit_args,
                             bash_log=bash_log, seed=seed, env=env):
        return None

    return wdir_md_cur
