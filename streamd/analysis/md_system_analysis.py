"""MD trajectory analysis helpers for StreaMD."""

import logging
from glob import glob
import os
import re
import pandas as pd
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import rms
from streamd.analysis.xvg2png import convertxvg2png
from streamd.analysis.plot_build import plot_rmsd
from streamd.utils.utils import get_index, make_group_ndx, get_mol_resid_pair, run_check_subprocess, create_last_frame_file


def rmsd_for_atomgroups(universe, selection1='backbone', selection2=None):
    """Calculate globally protein-backbone-aligned, by default, RMSD for atom groups.

    All reported group RMSDs are calculated after alignment using
    ``selection1``.  In StreaMD this is normally the protein backbone.

    Parameters
    ----------
    universe: MDAnalysis.core.universe.Universe
        MDAnalysis universe.
    selection1: str
        Selection string for main atom group, also used during alignment.
    selection2: list of str, optional
        Selection strings for additional atom groups.

    Returns
    -------
    pandas.DataFrame
        DataFrame containing RMSD of the selected atom groups over time.
    """
    universe.trajectory[0]
    ref = universe
    rmsd_analysis = rms.RMSD(universe, ref, select=selection1, groupselections=selection2, in_memory=False)
    rmsd_analysis.run()
    columns = [selection1, *(selection2 or [])]
    results = rmsd_analysis.results.rmsd
    rmsd_df = pd.DataFrame(np.round(results[:, 2:], 2), columns=columns)
    # MDAnalysis reports GROMACS trajectory time in ps; StreaMD CSVs use ns.
    times_ps = results[:, 1]
    if len(times_ps) > 1 and not np.all(np.diff(times_ps) > 0):
        logging.warning(
            'Trajectory time is not strictly increasing (%.3f ... %.3f ps); '
            'time(ns) values may be unreliable. Check that the trajectory carries '
            'per-frame time information.', times_ps[0], times_ps[-1])
    rmsd_df.insert(0, 'time(ns)', times_ps / 1000.0)
    return rmsd_df


def _unique_residue_names(molid_resid_pairs):
    """Return residue names from molecule/residue pairs, preserving first-seen order."""
    resids = []
    for _, resid in molid_resid_pairs:
        if resid not in resids:
            resids.append(resid)
    return resids


def _ensure_index_group(group_name, members, index_list, wdir, bash_log, env=None):
    """Return an index group id, creating the group from member groups when needed.

    The group is looked up by its exact (canonical) name and created when absent, so
    the group StreaMD references downstream is always the deterministically-named one
    built from the current residue order. In normal operation the same
    all_ligand_resid.txt yields the same name across runs, so re-runs reuse the group
    rather than duplicating it.
    """
    if group_name not in index_list:
        missing = [member for member in members if member not in index_list]
        if missing:
            logging.error(
                'Cannot build the %s index group: group(s) %s are missing from index.ndx. '
                'This usually means a residue listed in all_ligand_resid.txt is absent from '
                'the simulated topology/coordinates (e.g. a ligand or cofactor that was not '
                'correctly built into the system). Analysis of this system will be skipped.',
                group_name, missing)
            return None, index_list
        query = '|'.join(str(index_list.index(member)) for member in members)
        if not make_group_ndx(query=query, wdir=wdir, bash_log=bash_log, env=env):
            return None, index_list
        index_list = get_index(os.path.join(wdir, 'index.ndx'), env=env)

    if group_name not in index_list:
        logging.error('Required index group %s failed to be found/created in index.ndx', group_name)
        return None, index_list
    return index_list.index(group_name), index_list


def _resolve_analysis_groups(index_list, molid_resid_pairs, ligand_resid, wdir, bash_log, env=None):
    """Resolve separate centering and trajectory-fitting index groups for analysis."""
    if molid_resid_pairs:
        resids = _unique_residue_names(molid_resid_pairs)
        # Every residue listed in all_ligand_resid.txt must have an auto-generated
        # index group; a missing one means the residue is absent from the simulated
        # topology/coordinates, so _ensure_index_group fails and this system is
        # skipped rather than analysed as a silently incomplete complex.
        center_members = ['Protein'] + resids
        center_group_name = '_'.join(center_members)
        center_group, index_list = _ensure_index_group(
            group_name=center_group_name,
            members=center_members,
            index_list=index_list,
            wdir=wdir,
            bash_log=bash_log,
            env=env,
        )
        if ligand_resid in resids:
            trajectory_fit_members = ['Protein', ligand_resid]
            trajectory_fit_group_name = '_'.join(trajectory_fit_members)
        else:
            trajectory_fit_members = center_members
            trajectory_fit_group_name = center_group_name
        if trajectory_fit_group_name == center_group_name:
            trajectory_fit_group = center_group
        else:
            trajectory_fit_group, index_list = _ensure_index_group(
                group_name=trajectory_fit_group_name,
                members=trajectory_fit_members,
                index_list=index_list,
                wdir=wdir,
                bash_log=bash_log,
                env=env,
            )
    else:
        if 'Protein' not in index_list:
            logging.error('Required center group Protein was not found in index.ndx. '
                          f'Check if {os.path.join(wdir, "index.ndx")} file was corrupted.')
            return None
        center_group = index_list.index('Protein')
        trajectory_fit_group = center_group

    if center_group is None or trajectory_fit_group is None:
        return None
    return center_group, trajectory_fit_group, index_list


def _index_selection(atomgroup):
    """Build an MDAnalysis index selection string for a fixed atom group."""
    return 'index ' + ' '.join(str(i) for i in atomgroup.indices)


def _has_minimum_fit_atoms(atomgroup):
    """Return whether an atom group can define a stable rotational fit."""
    if len(atomgroup) < 3:
        return False
    positions = getattr(atomgroup, 'positions', None)
    if positions is None:
        return True
    return np.linalg.matrix_rank(positions - positions[0]) >= 2


def _reference_pocket_selections(universe, ligand_resid, active_site_dist):
    """Define fixed reference-frame pocket and ligand selections for local RMSD."""
    ligand_selection = f'resname {ligand_resid} and not name H*'
    universe.trajectory[0]
    ligand = universe.select_atoms(ligand_selection)
    if len(ligand) == 0:
        logging.warning('Cannot define active site: no ligand atoms found for resname %s', ligand_resid)
        return None, ligand_selection

    if len(ligand.residues) > 1:
        # Multiple copies of the ligand (e.g. one per chain in a multimer) would be
        # merged into a single cross-site pocket, and the local fit would span the
        # copies' separate binding sites, so ActiveSite/ligand_local become
        # meaningless. Skip the pocket metrics; the global 'ligand' RMSD still runs.
        # TODO (option B): support multiple binding sites by defining a pocket and
        # computing ActiveSite/ligand_local per ligand copy (per chain/segid) and
        # emitting per-copy columns instead of skipping. This needs schema/selection
        # changes in run_analysis (rmsd_type, _resolve_rmsd_columns) and the plots.
        logging.warning(
            'Found %d residues for resname %s (resids %s), e.g. copies in different '
            'chains; active-site and ligand_local RMSD are only defined for a single '
            'binding site and are skipped for this system. The global "ligand" RMSD is '
            'still reported.',
            len(ligand.residues), ligand_resid, list(ligand.residues.resids),
        )
        return None, ligand_selection

    pocket_atoms = universe.select_atoms(
        f'protein and around {active_site_dist} group ligand',
        ligand=ligand,
    )
    pocket_backbone = pocket_atoms.residues.atoms.select_atoms('backbone')
    if len(pocket_backbone) == 0:
        logging.warning(
            'Cannot define active-site RMSD: no pocket-backbone atoms within %.2f A of %s in the reference frame',
            active_site_dist,
            ligand_resid,
        )
        return None, ligand_selection

    if not _has_minimum_fit_atoms(pocket_backbone):
        logging.warning(
            'Cannot define active-site RMSD: fewer than three non-collinear pocket-backbone atoms for %s',
            ligand_resid,
        )
        return None, ligand_selection

    return _index_selection(pocket_backbone), ligand_selection


def _add_local_ligand_rmsd(rmsd_df, universe, pocket_selection, ligand_selection):
    """Append ligand RMSD after local fitting on the reference-defined pocket."""
    try:
        universe.trajectory[0]
        local_rmsd = rms.RMSD(
            universe,
            universe,
            select=pocket_selection,
            groupselections=[ligand_selection],
            in_memory=False,
        )
        local_rmsd.run()
        rmsd_df.loc[:, 'ligand_local'] = np.round(local_rmsd.results.rmsd[:, 3], 2)
    except (ValueError, RuntimeError, np.linalg.LinAlgError) as exc:
        logging.warning('Skipping ligand_local RMSD after local pocket alignment: %s', exc)
    return rmsd_df


def _parse_system_name(system: str) -> tuple[str, int]:
    """Split a system name into its base name and numeric replica id."""
    match = re.match(r"(.*)_replica(\d+)$", system)
    if match:
        return match.group(1), int(match.group(2))
    return system, 1


def _system_metadata(system_name, ligand_name):
    """Return output system, protein, and replica metadata for an RMSD table."""
    protein_name, replica = _parse_system_name(system_name)
    replica_suffix = system_name[len(protein_name):]

    if ligand_name:
        ligand_suffix = f'_{ligand_name}'
        if protein_name.endswith(ligand_suffix):
            protein_name = protein_name[:-len(ligand_suffix)]

    return f'{protein_name}{replica_suffix}', protein_name, replica


def _load_universe_with_gro_fallback(topology, trajectory):
    """Load an MDAnalysis universe, falling back from unsupported TPR to matching GRO."""
    try:
        return mda.Universe(topology, trajectory, in_memory=False, in_memory_step=1)
    except ValueError as exc:
        error_message = str(exc)
        gro = f'{os.path.splitext(topology)[0]}.gro'
        if (
            os.path.splitext(topology)[1] != '.tpr'
            or 'tpx version' not in error_message
            or 'does not support' not in error_message
            or not os.path.isfile(gro)
        ):
            raise
        logging.warning('Cannot parse %s with MDAnalysis; using %s as topology instead', topology, gro)
        return mda.Universe(gro, trajectory, in_memory=False, in_memory_step=1)


def md_rmsd_analysis(tpr, xtc, wdir_out_analysis, system_name,
                     molid_resid_pairs,
                     ligand_resid="UNL", active_site_dist=5.0):
    """Calculate RMSD profiles for a single MD system.

    ``backbone``, ``CA`` and the group RMSDs (ligand, ActiveSite, cofactors) are all
    calculated after global protein-backbone alignment; ``backbone``
    is the fit group and ``CA`` its Cα subset.  ``ActiveSite{active_site_dist}A`` is a
    reference-defined pocket-backbone RMSD after that global alignment.
    ``ligand_local``, when present, is ligand heavy-atom RMSD after local alignment
    on the same reference-defined pocket backbone.

    Parameters
    ----------
    tpr : str | os.PathLike
        Path to the MDAnalysis topology file of the simulation.
    xtc : str | os.PathLike
        Path to the trajectory ``.xtc`` file.
    wdir_out_analysis : str | os.PathLike
        Directory where analysis outputs will be stored.
    system_name : str
        Identifier of the analysed system used in file names.
    molid_resid_pairs : Iterable[tuple[str, str]]
        Pairs mapping molecule identifiers to residue names.
    ligand_resid : str, optional
        Residue name of the ligand to analyse.
    active_site_dist : float, optional
        Distance in Ångström defining the active-site region around the ligand.

    Returns
    -------
    str
        Path to the generated ``rmsd_*.csv`` file.
    """
    # groupselections = ['protein']
    rmsd_out_file = os.path.join(wdir_out_analysis, f'rmsd_{system_name}.csv')
    universe = _load_universe_with_gro_fallback(tpr, xtc)
    groupselections = []
    molid_resid_pairs = dict(molid_resid_pairs)
    mol_residue_names = _unique_residue_names(molid_resid_pairs.items())
    ligand_name = None
    active_site_selection = None
    ligand_selection = f'resname {ligand_resid} and not name H*'
    if molid_resid_pairs:
        if ligand_resid in molid_resid_pairs.values():
            for molid, resid in molid_resid_pairs.items():
                if resid == ligand_resid:
                    ligand_name = molid
                    break
            active_site_selection, ligand_selection = _reference_pocket_selections(
                universe=universe,
                ligand_resid=ligand_resid,
                active_site_dist=active_site_dist,
            )
            if active_site_selection:
                groupselections.append(active_site_selection)

        for resid in mol_residue_names:
            selection = f"resname {resid} and not name H*"
            if len(universe.select_atoms(selection)) == 0:
                logging.warning('Skipping RMSD for resname %s: no heavy atoms found', resid)
                continue
            groupselections.append(selection)

    # Superpose on the protein backbone and report every RMSD in that
    # one frame: 'backbone' is the fit group, 'CA' is the Cα subset, and the ligand/
    # ActiveSite/cofactor RMSDs are all measured relative to the same backbone alignment.
    protein_ca_selection = "protein and name CA"
    rmsd_df = rmsd_for_atomgroups(universe, selection1="backbone",
                                  selection2=[protein_ca_selection] + groupselections)
    if active_site_selection and ligand_selection in groupselections:
        rmsd_df = _add_local_ligand_rmsd(
            rmsd_df=rmsd_df,
            universe=universe,
            pocket_selection=active_site_selection,
            ligand_selection=ligand_selection,
        )
    del universe
    rmsd_df = rmsd_df.rename(
        {protein_ca_selection: 'CA',
         active_site_selection: f'ActiveSite{active_site_dist}A',
         f'resname {ligand_resid} and not name H*': 'ligand'}, axis='columns')

    rmsd_df = rmsd_df.rename({f"resname {i[1]} and not name H*": f"{i[0]}" for i in molid_resid_pairs.items()}, axis='columns')

    # Plot the backbone protein-RMSD line (not CA); the two are near-identical, and CA
    # is still written to the CSV. Ligand/ActiveSite/cofactor lines are still plotted.
    plot_rmsd(rmsd_df=rmsd_df.drop(columns=['CA'], errors='ignore'),
              system_name=system_name, out=os.path.join(wdir_out_analysis, f'rmsd_{system_name}.png'))

    rmsd_df.loc[:, 'ligand_name'] = ligand_name
    base_system, protein_name, replica = _system_metadata(system_name, ligand_name)
    rmsd_df.loc[:, 'system'] = base_system
    rmsd_df.loc[:, 'replica'] = replica
    rmsd_df.loc[:, 'protein_name'] = protein_name
    rmsd_df.loc[:, 'directory'] = wdir_out_analysis
    rmsd_df = rmsd_df[['time(ns)', *[column for column in rmsd_df.columns if column != 'time(ns)']]]

    rmsd_df.to_csv(rmsd_out_file, sep='\t', index=False)
    return rmsd_out_file


def run_md_analysis(var_md_dirs_deffnm, mdtime_ns, project_dir, bash_log,
                    active_site_dist=5.0, ligand_resid='UNL',
                    save_traj_without_water = False,
                    analysis_dirname = 'md_analysis',
                    ligand_list_file_prev=None, env=None, system_name=None):
    """Prepare trajectories and execute RMSD analysis.

    Parameters
    ----------
    var_md_dirs_deffnm : tuple[str, str]
        Tuple of working directory and deffnm base name.
    mdtime_ns : float
        Target simulation length in nanoseconds.
    project_dir : str | os.PathLike
        Base directory containing helper scripts.
    bash_log : str
        Name of log file capturing shell output.
    active_site_dist : float, optional
        Distance in Ångström to define the active site around the ligand.
    ligand_resid : str, optional
        Residue name identifying the ligand.
    save_traj_without_water : bool, optional
        Keep intermediate trajectories without water if ``True``.
    analysis_dirname : str, optional
        Name of the directory to store analysis files.
    ligand_list_file_prev : str | None, optional
        Optional file listing ligand/residue pairs from a previous run.
    env : dict | None, optional
        Environment variables for subprocess calls.
    system_name : str | None, optional
        Custom name for the MD system.

    Returns
    -------
    tuple[str, str, str] | None
        Tuple of RMSD output file path, analysis directory and working
        directory, or ``None`` on failure.
    """
    wdir, deffnm = var_md_dirs_deffnm

    # create subdir for analysis files only
    wdir_out_analysis = os.path.join(wdir, analysis_dirname)
    os.makedirs(wdir_out_analysis, exist_ok=True)

    if ligand_list_file_prev is None:
        molid_resid_pairs_fname = os.path.join(wdir, 'all_ligand_resid.txt')
    else:
        molid_resid_pairs_fname = ligand_list_file_prev

    # if md was continued
    if not os.path.isfile(os.path.join(wdir, 'index.ndx')):
        cmd = f'''
        cd {wdir}
        gmx make_ndx -f {deffnm}.gro << INPUT
        q
        INPUT
        '''
        if not run_check_subprocess(cmd, key=wdir, log=os.path.join(wdir, bash_log), env=env):
            return None

    index_list = get_index(os.path.join(wdir, 'index.ndx'), env=env)

    # choose groups to center and fit the trajectory
    if os.path.isfile(molid_resid_pairs_fname) and os.path.getsize(molid_resid_pairs_fname) > 0:
        molid_resid_pairs = list(get_mol_resid_pair(molid_resid_pairs_fname))
    else:
        molid_resid_pairs = []

    analysis_groups = _resolve_analysis_groups(
        index_list=index_list,
        molid_resid_pairs=molid_resid_pairs,
        ligand_resid=ligand_resid,
        wdir=wdir,
        bash_log=bash_log,
        env=env,
    )
    if analysis_groups is None:
        return None
    center_group, trajectory_fit_group, index_list = analysis_groups

    dtstep = 50 if mdtime_ns <= 10 else 100

    tpr = os.path.join(wdir, f'{deffnm}.tpr')
    xtc = os.path.join(wdir, f'{deffnm}.xtc')

    if not system_name:
        system_name = os.path.split(wdir)[-1]

    cmd = f'wdir={wdir} center_group={center_group} trajectory_fit_group={trajectory_fit_group} dtstep={dtstep} deffnm={deffnm} tpr={tpr} xtc={xtc} wdir_out_analysis={wdir_out_analysis} system_name={system_name} ' \
           f'bash {os.path.join(project_dir, "scripts/script_sh/md_analysis.sh")} >> {os.path.join(wdir, bash_log)} 2>&1'

    if not run_check_subprocess(cmd, key=wdir, log=os.path.join(wdir, bash_log), env=env):
        return None

    # molid resid pairs for all ligands in the MD system
    # calc rmsd
    # universe = mda.Universe(tpr, os.path.join(wdir, f'md_fit.xtc'))

    rmsd_out_file = md_rmsd_analysis(
        tpr=os.path.join(wdir, 'md_out_nowater.tpr'),
        xtc=os.path.join(wdir, 'md_fit_nowater.xtc'),
        # tpr=os.path.join(wdir, 'md_out.tpr'), xtc=os.path.join(wdir, 'md_fit.xtc'),
                     wdir_out_analysis=wdir_out_analysis,
                     system_name=system_name,
                     ligand_resid=ligand_resid,
                     molid_resid_pairs=molid_resid_pairs,
                     active_site_dist=active_site_dist)
    if not save_traj_without_water:
        for path in [
            os.path.join(wdir, 'md_out_nowater.tpr'),
            os.path.join(wdir, 'md_out_nowater.gro'),
            os.path.join(wdir, 'md_fit_nowater.xtc'),
        ]:
            if os.path.isfile(path):
                os.remove(path)

    for xvg_file in glob(os.path.join(wdir_out_analysis, '*.xvg')):
        convertxvg2png(xvg_file, system_name=system_name, transform_nm_to_A=True)

    create_last_frame_file(wdir=wdir,
                           tpr=tpr, xtc=xtc,
                           out_file=os.path.join(wdir, 'last_frame.pdb'),
                           bash_log=bash_log, env=env)

    return rmsd_out_file, wdir_out_analysis, wdir
