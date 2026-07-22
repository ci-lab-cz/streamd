#!/usr/bin/env python3

"""Run ProLIF interaction analysis over MD trajectories."""

import argparse
from datetime import datetime
import os
import shutil
from functools import partial
from glob import glob
import logging
import pathlib
import sys

import MDAnalysis as mda
import pandas as pd
import prolif as plf
from prolif.plotting.barcode import Barcode
from prolif.plotting.network import LigNetwork
import matplotlib.pyplot as plt

from streamd.utils.dask_init import init_dask_cluster, calc_dask
from streamd.utils.utils import filepath_type, parse_with_config
from streamd.prolif.prolif2png import convertprolif2png
from streamd.prolif.prolif_frame_map import convertplifbyframe2png
plt.ioff()

class RawTextArgumentDefaultsHelpFormatter(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass

def backup_output(output):
    """Rename existing output to avoid overwriting previous results."""
    if os.path.isfile(output):
        all_outputs = glob(os.path.join(os.path.dirname(output), f'#{os.path.basename(output)}*#'))
        n = len(all_outputs) + 1
        shutil.move(output, os.path.join(os.path.dirname(output), f'#{os.path.basename(output)}.{n}#'))


# ProLIF only detects interactions within this distance of the ligand
# (Fingerprint.vicinity_cutoff). Any binding-site cutoff must exceed it to keep
# the fingerprints identical to a full-protein analysis.
PROLIF_VICINITY_CUTOFF = 6.0


def restrict_protein_to_binding_site(u, protein_selection, ligand_selection, trajectory, cutoff):
    """Restrict the protein to residues near the ligand over the whole trajectory.

    ProLIF converts the *entire* protein AtomGroup to an RDKit molecule on every
    frame (bond perception + sanitization), which dominates the runtime. Because
    it only detects interactions within ``vicinity_cutoff`` (6 A) of the ligand,
    handing it the union of protein residues seen within a larger ``cutoff`` in
    any analysed frame yields identical fingerprints while converting far fewer
    atoms per frame.

    :param u: MDAnalysis Universe (already renumbered if a pdb was provided).
    :param protein_selection: MDAnalysis selection string for the protein.
    :param ligand_selection: MDAnalysis selection string for the ligand.
    :param trajectory: the (possibly stepped) trajectory iterator to scan.
    :param cutoff: distance in Angstrom; must be > 6 A to preserve results.
    :return: an AtomGroup restricted to the binding-site residues, or None if the
        selection turned out empty (caller should fall back to the full protein).
    """
    pocket = u.select_atoms(
        f'byres (({protein_selection}) and around {cutoff} ({ligand_selection}))',
        updating=True)
    resindices = set()
    current_frame = u.trajectory.frame
    for _ in trajectory:
        resindices.update(int(i) for i in pocket.residues.resindices)
    u.trajectory[current_frame]  # restore the frame pointer for downstream steps

    if not resindices:
        return None
    return u.residues[sorted(resindices)].atoms


def load_ligand_template(path):
    """Load a reference ligand structure (SDF/MOL) that carries correct bond orders.

    Used as a template for MDAnalysis' TemplateInferrer so that ProLIF assigns the ligand
    bond orders from a trusted structure instead of inferring them from the topology. This
    is usually unnecessary - see the ``ligand_sdf`` note in :func:`run_prolif_task`.

    :param path: path to a .sdf or .mol file.
    :return: an RDKit Mol with heavy atoms only (hydrogens removed, as required by
        MDAnalysis' TemplateInferrer / AssignBondOrdersFromTemplate), or None if it
        could not be read.
    """
    from rdkit import Chem
    ext = os.path.splitext(path)[1].lower()
    try:
        # removeHs=True: the template must be a heavy-atom graph - TemplateInferrer strips
        # hydrogens from the target before matching, so an H-bearing template fails to match.
        if ext == '.sdf':
            supplier = Chem.SDMolSupplier(path, removeHs=True, sanitize=True)
            return next((m for m in supplier if m is not None), None)
        if ext == '.mol':
            return Chem.MolFromMolFile(path, removeHs=True, sanitize=True)
    except Exception:
        return None
    return None


def resolve_ligand_converter_kwargs(ligand, ligand_selection, ligand_sdf, xtc):
    """Determine the RDKit converter kwargs used to interpret the ligand chemistry.

    For all-atom topologies ProLIF infers bond orders correctly from the TPR alone, which is why
    the resulting SMILES is logged as a sanity check. A ``ligand_sdf`` template is USUALLY
    UNNECESSARY: supply it only for the rare cases where the inference is wrong (united-atom /
    coarse-grained ligands, or an odd SMILES below).

    :param ligand: the ligand MDAnalysis AtomGroup.
    :param ligand_selection: MDAnalysis selection string for the ligand (used for logging).
    :param ligand_sdf: optional path to a reference ligand structure (.sdf/.mol) with correct bond
        orders, used as a template. See the ``ligand_sdf`` note in :func:`run_prolif_task`.
    :param xtc: trajectory path (used for logging only).
    :return: a tuple ``(lig_converter_kwargs, converter_kwargs)`` where ``lig_converter_kwargs`` is
        the dict of kwargs for converting the ligand AtomGroup to RDKit (for the network plot) and
        ``converter_kwargs`` is the ``(ligand_kwargs, protein_kwargs)`` tuple ProLIF's
        ``Fingerprint.run`` expects, or None to let ProLIF infer bond orders from the topology.
    """
    from rdkit import Chem
    lig_converter_kwargs = {}
    converter_kwargs = None
    try:
        topo_heavy = ligand.convert_to.rdkit().GetNumHeavyAtoms()
        topo_smiles = Chem.MolToSmiles(Chem.RemoveHs(ligand.convert_to.rdkit()))
    except Exception as e:
        topo_heavy, topo_smiles = None, None
        logging.warning(f'{xtc}: could not convert the ligand from the topology ({e}).')

    if ligand_sdf:
        if not os.path.isfile(ligand_sdf):
            logging.warning(f'{xtc}: ligand template "{ligand_sdf}" not found; using bond orders '
                            f'inferred from the topology instead (usually fine).')
        else:
            template = load_ligand_template(ligand_sdf)
            if template is None:
                logging.warning(f'{xtc}: could not parse ligand template "{ligand_sdf}" (expected .sdf/.mol); '
                                f'using topology inference instead (usually fine).')
            elif topo_heavy is not None and template.GetNumAtoms() != topo_heavy:
                logging.warning(f'{xtc}: ligand template "{os.path.basename(ligand_sdf)}" has '
                                f'{template.GetNumAtoms()} heavy atoms but the ligand has {topo_heavy}; '
                                f'template ignored, using topology inference instead.')
            else:
                from MDAnalysis.converters.RDKitInferring import TemplateInferrer
                kw = {'inferrer': TemplateInferrer(template=template)}
                try:
                    tmpl_smiles = Chem.MolToSmiles(Chem.RemoveHs(ligand.convert_to.rdkit(**kw)))
                    lig_converter_kwargs, converter_kwargs = kw, (kw, {})
                    logging.info(f'{xtc}: ligand "{ligand_selection}" bond orders assigned from template '
                                 f'{os.path.basename(ligand_sdf)}; SMILES: {tmpl_smiles}.')
                except Exception as e:
                    logging.warning(f'{xtc}: ligand template "{os.path.basename(ligand_sdf)}" did not match '
                                    f'the ligand ({e}); using topology inference instead.')

    if converter_kwargs is None and topo_smiles is not None:
        logging.info(f'{xtc}: ligand "{ligand_selection}" interpreted as SMILES: {topo_smiles} '
                     f'(inferred from topology). Usually correct without a template; pass --ligand_sdf '
                     f'only if it looks wrong.')

    return lig_converter_kwargs, converter_kwargs


def run_prolif_task(tpr, xtc, protein_selection, ligand_selection, step, verbose, output, n_jobs,
                    occupancy = 0.6, save_viz=True, dpi=300, plot_width=15, plot_height=8, pdb=None,
                    binding_site_cutoff=12.0, parallel_strategy='chunk', ligand_sdf=None,
                    water_bridge=False, water_selection='resname SOL', water_bridge_order=1,
                    water_cutoff=8.0):
    """Compute protein–ligand interaction fingerprints for a single trajectory.

    :param tpr:
    :param xtc:
    :param protein_selection:
    :param ligand_selection:
    :param step:
    :param verbose:
    :param output:
    :param n_jobs:
    :param save_pics: save barcode in png and network in html
    :param dpi:
    :param plot_width:  in inches
    :param plot_height: in inches
    :param binding_site_cutoff: restrict the protein to residues within this distance (A) of the
        ligand in any analysed frame. ProLIF re-converts the whole protein to RDKit every frame,
        so this is the main speed lever. Results stay identical while the cutoff stays above the
        ProLIF vicinity cutoff (6 A). Set to <= 0 to analyse the full protein selection.
    :param parallel_strategy: ProLIF multiprocessing strategy ('chunk', 'queue' or None for ProLIF
        auto), only used when n_jobs > 1. 'chunk' distributes trajectory chunks across workers, so the
        per-frame RDKit conversion is parallelised too; 'queue' converts frames in a producer thread on
        the main process and streams them to the workers, avoiding repeated pickling of large
        trajectories but potentially serialising the conversion.
    :param ligand_sdf: optional path to a reference ligand structure (.sdf/.mol) with correct bond
        orders, used as a template for the ligand conversion. USUALLY UNNECESSARY: for standard
        all-atom topologies ProLIF infers the ligand chemistry correctly from the TPR alone (the
        interpreted SMILES is logged so it can be checked). Provide it only for united-atom /
        coarse-grained ligands, or if the logged SMILES is wrong.
    :param water_bridge: also compute water-mediated (water-bridge) interactions
    :param water_selection: MDAnalysis selection string for water molecules (used only if water_bridge)
    :param water_bridge_order: maximum number of water molecules bridging the ligand and protein
    :param water_cutoff: only waters within this distance (A) of the ligand each frame are considered
        for water-bridge analysis. This is the main speed lever for water-bridge (ProLIF otherwise
        converts every water to RDKit every frame). Results stay identical while the cutoff exceeds
        the ProLIF vicinity cutoff (6 A); it is automatically widened for higher bridge orders. Set
        to <= 0 to consider all waters.
    :return: pandas dataframe
    """
    u = mda.Universe(tpr, xtc, in_memory=False, in_memory_step=1)

    protein = u.atoms.select_atoms(protein_selection)
    ligand = u.atoms.select_atoms(ligand_selection)

    if pdb:
        u1 = mda.Universe(pdb)
        protein_pdb = u1.atoms.select_atoms(protein_selection)
        if len(protein.residues.resids) == len(protein_pdb.residues.resids):
            protein.residues.resids = protein_pdb.residues.resids
        if len(protein.segments.segids) == len(protein_pdb.segments.segids):
            protein.segments.segids = protein_pdb.segments.segids

    trajectory = u.trajectory[::step] if step > 1 else u.trajectory

    if binding_site_cutoff and binding_site_cutoff > 0:
        # A water bridge links the ligand to a protein residue through one or more waters,
        # so the participating residue can sit further from the ligand than a direct contact.
        # Add margin per bridging water to keep water-bridge results identical too.
        effective_cutoff = binding_site_cutoff
        if water_bridge:
            effective_cutoff += 2 * PROLIF_VICINITY_CUTOFF * water_bridge_order
        if effective_cutoff <= PROLIF_VICINITY_CUTOFF:
            logging.warning(f'{xtc}: binding_site_cutoff ({binding_site_cutoff}) is not above the ProLIF '
                            f'vicinity cutoff ({PROLIF_VICINITY_CUTOFF} A); analysing the full protein instead '
                            f'to avoid dropping contacts.')
        else:
            binding_site = restrict_protein_to_binding_site(u, protein_selection, ligand_selection,
                                                            trajectory, effective_cutoff)
            if binding_site is None or binding_site.n_atoms == 0:
                logging.warning(f'{xtc}: no protein residues found within {effective_cutoff} A of the ligand; '
                                f'analysing the full protein selection instead.')
            else:
                logging.info(f'{xtc}: restricting ProLIF protein to {binding_site.n_residues} binding-site '
                             f'residues ({binding_site.n_atoms} atoms) within {effective_cutoff} A of the ligand '
                             f'(full protein: {protein.n_residues} residues).')
                protein = binding_site

    interactions = ['Hydrophobic', 'HBDonor', 'HBAcceptor', 'Anionic', 'Cationic', 'CationPi', 'PiCation',
                    'PiStacking', 'MetalAcceptor', 'XBDonor', 'XBAcceptor']
    parameters = None
    if water_bridge:
        water_all = u.atoms.select_atoms(water_selection)
        if len(water_all) == 0:
            logging.warning(f'{xtc}: water-bridge analysis was requested but no atoms matched the water selection '
                            f'"{water_selection}". Water-mediated contacts will be skipped for this trajectory.')
        else:
            # Water-bridge is enormously expensive because ProLIF converts every water molecule to
            # RDKit on every frame (in the ligand-water, water-protein and, for order>=2, water-water
            # passes). A bridging water must hydrogen-bond the ligand, so it is always close to it;
            # restricting the water to an updating (per-frame) selection near the ligand keeps only the
            # waters that can actually bridge and leaves results identical. The cutoff must exceed the
            # ProLIF vicinity cutoff (6 A) and grows with the bridge order (each extra water is one more
            # ~vicinity-length hop away from the ligand).
            if water_cutoff and water_cutoff > 0:
                effective_water_cutoff = water_cutoff + PROLIF_VICINITY_CUTOFF * (water_bridge_order - 1)
                water = u.select_atoms(
                    f'byres (({water_selection}) and around {effective_water_cutoff} ({ligand_selection}))',
                    updating=True)
                logging.info(f'{xtc}: water-bridge restricted to waters within {effective_water_cutoff} A of the '
                             f'ligand each frame (out of {water_all.n_residues} waters); set --water_cutoff 0 '
                             f'to consider all waters.')
            else:
                water = water_all
            interactions.append('WaterBridge')
            parameters = {'WaterBridge': {'water': water, 'order': water_bridge_order}}

    lig_converter_kwargs, converter_kwargs = resolve_ligand_converter_kwargs(
        ligand, ligand_selection, ligand_sdf, xtc)

    fp = plf.Fingerprint(interactions, parameters=parameters)
    fp.run(trajectory, ligand, protein, progress=verbose, n_jobs=n_jobs,
           parallel_strategy=parallel_strategy, converter_kwargs=converter_kwargs)

    df = fp.to_dataframe()
    df.columns = ['.'.join(item.strip().lower() for item in items[1:]) for items in df.columns]
    df = df.reindex(sorted(df.columns), axis=1)
    df.to_csv(output, sep='\t')

    if save_viz:
        # barcode
        Barcode.from_fingerprint(fp).display(figsize=(plot_width, plot_height)).figure.savefig(f'{output.rstrip(".csv")}.png', dpi=dpi)
        # Net
        net = LigNetwork.from_fingerprint(fp, ligand_mol=ligand.convert_to.rdkit(**lig_converter_kwargs),
                                          threshold=occupancy)
        net_output = f'{output.rstrip(".csv")}_occupancy{occupancy}.html'
        try:
            net.save(net_output, show_interaction_data=True)
        except TypeError:
            # show_interaction_data is only supported by newer ProLIF versions
            logging.warning(f'{xtc}: the installed ProLIF version does not support "show_interaction_data"; '
                            f'saving the interaction network without per-interaction data.')
            net.save(net_output)
        convertplifbyframe2png(plif_out_file=output, plot_width=plot_width, plot_height=plot_height, point_size=5)

    return df


def run_prolif_from_wdir(wdir, tpr, xtc, protein_selection, ligand_selection, step, verbose, output,
                         plot_width, plot_height, save_viz, pdb, n_jobs, occupancy,
                         binding_site_cutoff=12.0, parallel_strategy='chunk', ligand_sdf=None,
                         water_bridge=False, water_selection='resname SOL', water_bridge_order=1,
                         water_cutoff=8.0):
    """Execute ProLIF analysis using paths relative to a directory."""
    tpr = os.path.join(wdir, tpr)
    xtc = os.path.join(wdir, xtc)
    if pdb:
        pdb = os.path.join(wdir, pdb)
    output = os.path.join(wdir, output)
    backup_output(output)

    if not os.path.isfile(tpr) or not os.path.isfile(xtc):
        print(f'{wdir}: cannot run prolif. Check if there are missing files: {tpr} {xtc}. Skip such directory')
        return None

    run_prolif_task(tpr=tpr, xtc=xtc, protein_selection=protein_selection,
                    ligand_selection=ligand_selection, step=step, verbose=verbose, output=output,
                    plot_width=plot_width, plot_height=plot_height, save_viz=save_viz, occupancy=occupancy,
                    pdb=pdb, n_jobs=n_jobs, binding_site_cutoff=binding_site_cutoff,
                    parallel_strategy=parallel_strategy, ligand_sdf=ligand_sdf,
                    water_bridge=water_bridge, water_selection=water_selection,
                    water_bridge_order=water_bridge_order, water_cutoff=water_cutoff)
    return output


def collect_outputs(output_list, output):
    """Concatenate individual ProLIF CSV outputs into one file."""
    df_list = []
    for i in output_list:
        df = pd.read_csv(i, sep='\t')
        # save dirname - protein_ligand pair
        df['Name'] = pathlib.PurePath(i).parent.name
        df['directory'] = os.path.dirname(i)
        df_list.append(df)

    df_aggregated = pd.concat(df_list)
    df_aggregated = df_aggregated.fillna(False).sort_values('Frame')
    amino_acids = df_aggregated.columns.drop(['Name', 'directory', 'Frame']).to_list()
    # sort by number and type of interaction
    amino_acids.sort(key=lambda x: (int(x.split('.')[0][3:]), x.split('.')[1]))
    sorted_columns = ['Name', 'directory', 'Frame'] + amino_acids
    df_aggregated.loc[:, sorted_columns].to_csv(output, sep='\t', index=False)


def start(wdir_to_run, wdir_output, tpr, xtc, step, append_protein_selection,
          protein_selection, ligand_resid, hostfile, ncpu, n_jobs,
          occupancy, plot_width, plot_height, save_viz, unique_id, pdb, verbose,
          binding_site_cutoff=12.0, parallel_strategy='chunk', ligand_sdf=None,
          water_bridge=False, water_selection='resname SOL', water_bridge_order=1,
          water_cutoff=8.0, show_percentage=True):
    """Run ProLIF across multiple directories and aggregate results.
    :param wdir_to_run: list
    :param wdir_output: path to dirn
    :param tpr: path to file
    :param xtc: path to file
    :param step: int
    :param append_protein_selection: str
    :param protein_selection: str
    :param ligand_resid: str
    :param hostfile: Nonen or path to file
    :param ncpu: int
    :param n_jobs: int
    :param occupancy: float
    :param plot_width: float
    :param plot_height: float
    :param save_viz: bool
    :param unique_id: str
    :param pdb: None or path to file (protein.pdb for renumbering)
    :param verbose: bool
    :param binding_site_cutoff: float, restrict protein to residues within this distance (A) of the ligand
    :param parallel_strategy: str or None, ProLIF multiprocessing strategy ('chunk'/'queue'/None-for-auto)
    :param ligand_sdf: None or path to a reference ligand .sdf/.mol (bond-order template; usually unnecessary)
    :param water_bridge: bool, compute water-mediated (water-bridge) interactions
    :param water_selection: str, MDAnalysis selection for water molecules
    :param water_bridge_order: int, maximum number of bridging water molecules
    :param water_cutoff: float, only waters within this distance (A) of the ligand each frame are used
    :param show_percentage: bool, draw the occupancy percentage as a label above each dot in the aggregated plot
    :return:
    """
    output = 'plifs.csv'
    output_aggregated = os.path.join(wdir_output, f'prolif_output_{unique_id}.csv')

    if append_protein_selection is not None:
        protein_selection = f'({protein_selection}) or ({append_protein_selection})'

    ligand_selection = f'resname {ligand_resid}'

    if wdir_to_run is not None:
        dask_client, cluster = None, None
        #n_jobs_per_task = n_jobs if n_jobs <= ncpu else ncpu
        if n_jobs is None:
            # limit to 12 https://github.com/chemosim-lab/ProLIF/issues/110
            n_jobs_per_task = min(12, ncpu // len(wdir_to_run))
        else:
            n_jobs_per_task = n_jobs

        logging.info(f'Allocating {n_jobs_per_task} n_jobs per each task.')

        try:
            dask_client, cluster = init_dask_cluster(hostfile=hostfile,
                                                     n_tasks_per_node=min(len(wdir_to_run), ncpu//n_jobs_per_task),
                                                     use_multi_servers=True if len(wdir_to_run) > ncpu else False,
                                                     ncpu=ncpu)
            var_prolif_out_files = []
            for res in calc_dask(run_prolif_from_wdir, wdir_to_run, dask_client=dask_client,
                                 tpr=tpr, xtc=xtc, protein_selection=protein_selection,
                                 ligand_selection=ligand_selection, step=step, verbose=verbose, output=output,
                                 plot_width=plot_width, plot_height=plot_height, save_viz=save_viz, pdb=pdb,
                                 n_jobs=n_jobs_per_task, occupancy=occupancy,
                                 binding_site_cutoff=binding_site_cutoff,
                                 parallel_strategy=parallel_strategy, ligand_sdf=ligand_sdf,
                                 water_bridge=water_bridge, water_selection=water_selection,
                                 water_bridge_order=water_bridge_order, water_cutoff=water_cutoff):
                if res:
                    var_prolif_out_files.append(res)
        finally:
            if dask_client:
                dask_client.retire_workers(dask_client.scheduler_info()['workers'],
                                           close_workers=True, remove=True)
                dask_client.shutdown()
            if cluster:
                cluster.close()
    else:
        output = os.path.join(os.path.dirname(xtc), output)
        backup_output(output)

        run_prolif_task(tpr=tpr, xtc=xtc, protein_selection=protein_selection,
                        ligand_selection=ligand_selection,
                        step=step, verbose=verbose, output=output,
                        plot_width=plot_width, plot_height=plot_height, save_viz=save_viz,
                        pdb=pdb, n_jobs= min(12, ncpu), occupancy=occupancy,
                        binding_site_cutoff=binding_site_cutoff,
                        parallel_strategy=parallel_strategy, ligand_sdf=ligand_sdf,
                        water_bridge=water_bridge, water_selection=water_selection,
                        water_bridge_order=water_bridge_order, water_cutoff=water_cutoff)
        var_prolif_out_files = [output]

    backup_output(output_aggregated)
    collect_outputs(var_prolif_out_files, output=output_aggregated)

    convertprolif2png(output_aggregated, occupancy=occupancy,
                      plot_width=plot_width, plot_height=plot_height,
                      base_size=12, point_size=3, show_percentage=show_percentage)
    finished_complexes_file = os.path.join(wdir_output, f"finished_prolif_files_{unique_id}.txt")
    with open(finished_complexes_file, 'w') as output:
        output.write("\n".join(var_prolif_out_files))

    logging.info(
        f'ProLIF calculation of {len(var_prolif_out_files)} were successfully finished.\n'
         f'Successfully finished complexes have been saved in {finished_complexes_file} file')



def main():
    """CLI entry point for ProLIF analysis."""
    parser = argparse.ArgumentParser(description='Get protein-ligand interactions from MD trajectories using '
                                                 'ProLIF module.',
                                     formatter_class=RawTextArgumentDefaultsHelpFormatter)
    parser.add_argument('--config', metavar='FILENAME', required=False,
                        type=partial(filepath_type, ext=("yml", "yaml")),
                        help='Path to YAML configuration file with default arguments')
    parser.add_argument('-i', '--wdir_to_run', metavar='DIRNAME', required=False, default=None, nargs='+',
                        type=partial(filepath_type, exist_type='dir'),
                        help='''single or multiple directories for simulations.
                             Should consist of: md_out.tpr and md_fit.xtc files''')
    parser.add_argument('--xtc', metavar='FILENAME', required=False,
                        help='input trajectory file (XTC). Will be ignored if --wdir_to_run is used')
    parser.add_argument('--tpr', metavar='FILENAME', required=False,
                        help='input topology file (TPR). Will be ignored if --wdir_to_run is used')
    parser.add_argument('-l', '--ligand', metavar='STRING', required=False, default='UNL',
                        help='residue name of a ligand in the input trajectory.')
    parser.add_argument('-s', '--step', metavar='INTEGER', required=False, default=1, type=int,
                        help='step to take every n-th frame. ps')
    parser.add_argument('--protein_selection', metavar='STRING',
                        required=False, default='protein',
                        help='The protein selection atoms. Normally keep the default "protein": cropping '
                             'the protein to the binding site is done automatically and safely by '
                             '--binding_site_cutoff (a trajectory-wide union). Use this option to change which '
                             'atoms count as the protein (e.g. limit to a chain, add cofactor, etc.).')
    parser.add_argument('-a', '--append_protein_selection', metavar='STRING', required=False, default=None,
                        help='the string which will be concatenated to the protein selection atoms. '
                             'Example: "resname ZN or resname MG".')
    parser.add_argument('-d', '--wdir', metavar='WDIR', default=None,
                        type=partial(filepath_type, check_exist=False, create_dir=True),
                        help='Working directory for program output. If not set the current directory will be used.')
    parser.add_argument('-v', '--verbose', action='store_true', default=True,
                        help='print progress.')
    parser.add_argument('--hostfile', metavar='FILENAME', required=False, type=str, default=None,
                        help='text file with addresses of nodes of dask SSH cluster. The most typical, it can be '
                             'passed as $PBS_NODEFILE variable from inside a PBS script. The first line in this file '
                             'will be the address of the scheduler running on the standard port 8786. If omitted, '
                             'calculations will run on a single machine as usual.')
    parser.add_argument('-c', '--ncpu', metavar='INTEGER', required=False, default=len(os.sched_getaffinity(0)), type=int,
                        help='number of CPU per server. By default, StreaMD utilizes all available cpus.')
    parser.add_argument('--n_jobs', metavar='INTEGER', required=False,
                         default=None, type=int,
                         help='Number of processes to run per each interaction analysis tasks. '
                              'By default, StreaMD distributes the specified number of cores (--ncpu) evenly '
                              'between the available CPUs and the number of tasks to execute (e.g., multiple directories provided via --wdir_to_run). '
                              'However, by default, the --n_jobs value is limited to 12 to avoid the bottleneck issue (https://github.com/chemosim-lab/ProLIF/issues/110) described by the ProLIF authors .'
                              'Users can override this limitation by explicitly specifying the --n_jobs argument value.')
    parser.add_argument('--width', metavar='FILENAME', default=15, type=int,
                        help='width of the output pictures')
    parser.add_argument('--height', metavar='FILENAME', default=10, type=int,
                        help='height of the output pictures')
    parser.add_argument('--binding_site_cutoff', metavar='float', default=12.0, type=float,
                        help='Restrict the ProLIF protein selection to residues that come within this distance '
                             '(in Angstrom) of the ligand in at least one analysed frame. ProLIF re-converts the whole '
                             'protein to an RDKit molecule on every frame, so shrinking the protein is the main '
                             'speed lever and typically gives a large speedup for big proteins. Results are '
                             'identical to a full-protein analysis as long as the cutoff stays above the ProLIF '
                             'vicinity cutoff (6 A). Set to 0 to disable and analyse the full protein '
                             'selection.')
    parser.add_argument('--parallel_strategy', required=False, default='chunk',
                        choices=['chunk', 'queue', 'auto'],
                        help='ProLIF multiprocessing strategy (only used when --n_jobs > 1). "chunk" distributes '
                             'trajectory chunks across workers, so the per-frame RDKit conversion is parallelised '
                             'too. "queue" converts frames in a producer thread on the main process and streams '
                             'them to the workers, avoiding repeated pickling of large trajectories but potentially '
                             'serialising the conversion. "auto" lets ProLIF choose.')
    parser.add_argument('--ligand_sdf', metavar='FILENAME', required=False, default=None,
                        type=partial(filepath_type, ext=('sdf', 'mol'), check_exist=False),
                        help='Optional reference ligand structure (.sdf/.mol) WITH correct bond orders, used as '
                             'a template to assign the ligand bond orders for ProLIF. USUALLY UNNECESSARY: for '
                             'standard all-atom topologies the ligand chemistry (bond orders, aromaticity, '
                             'charges) is inferred correctly from the TPR alone, and the interpreted SMILES is '
                             'written to the log so you can verify it. Provide this only for united-atom / '
                             'coarse-grained ligands, or if the logged SMILES is wrong.')
    parser.add_argument('--water_bridge', default=False, action='store_true',
                        help='additionally compute water-mediated (water-bridge) interactions between the ligand '
                             'and the protein. Requires water molecules to be present in the input trajectory '
                             '(the default md_fit.xtc keeps them). This adds extra computation.')
    parser.add_argument('--water_selection', metavar='STRING', required=False, default='resname SOL',
                        help='MDAnalysis selection string for water molecules used by --water_bridge. '
                             'The GROMACS default water residue name is SOL.')
    parser.add_argument('--water_bridge_order', metavar='INTEGER', required=False, default=1, type=int,
                        help='maximum number of water molecules that can bridge the ligand and the protein '
                             '(only used with --water_bridge). Order 1 considers a single bridging water.')
    parser.add_argument('--water_cutoff', metavar='float', default=8.0, type=float,
                        help='Only waters within this distance (in Angstrom) of the ligand in a given frame are '
                             'considered for --water_bridge analysis. This is the main speed lever for '
                             'water-bridge: ProLIF otherwise converts EVERY water molecule to RDKit on every '
                             'frame, which is extremely slow for solvated systems. A bridging water must '
                             'hydrogen-bond the ligand, so results stay identical as long as the cutoff exceeds '
                             'the ProLIF vicinity cutoff (6 A); the value is widened automatically for higher '
                             '--water_bridge_order. Set to 0 to consider all waters (slow).')
    parser.add_argument('--occupancy', metavar='float', default=0.6, type=float,
                        help='occupancy of the unique contacts to show. '
                             'Applied for plifs_occupancyX.html (for each complex) and'
                             ' prolif_output_occupancyX.png (all systems aggregated plot)')
    parser.add_argument('--not_save_pics', default=False, action='store_true',
                        help='dont create html and png files (by frames) for each unique trajectory.'
                             ' Only overall prolif png plot file will be created.')
    parser.add_argument('--no-show_percentage', default=False, action='store_true',
                        help='do not show the occupancy percentage label above each dot in the aggregated '
                             'prolif_output_occupancyX.png plot (percentages are shown by default).')
    parser.add_argument('-o','--out_suffix',
                        metavar='string', default=None,
                        help='Unique suffix for output files. By default, start-time_unique-id.'
                             'Unique suffix is used to separate outputs from different runs.')
    args, _ = parse_with_config(parser, sys.argv[1:])

    if args.wdir is None:
        wdir = os.getcwd()
    else:
        wdir = args.wdir

    out_time = f'{datetime.now().strftime("%d-%m-%Y-%H-%M-%S")}'
    if args.out_suffix:
        unique_id = args.out_suffix
    else:
        import secrets
        out_suffix = secrets.token_hex(3)
        unique_id = f'{out_time}_unique-id-{out_suffix}'

    log_file = os.path.join(wdir, f'log_prolif_{unique_id}.log')

    logging.basicConfig(format='%(asctime)s - %(name)s - %(levelname)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S',
                        level=logging.INFO,
                        handlers=[logging.FileHandler(log_file),
                                  logging.StreamHandler()])


    if args.wdir_to_run is not None:
        tpr = 'md_out.tpr'
        xtc = 'md_fit.xtc'
        pdb = 'frame.pdb'
    else:
        tpr = args.tpr
        xtc = args.xtc
        pdb = None

    logging.getLogger('distributed').setLevel('CRITICAL')
    logging.getLogger('distributed.core').setLevel('CRITICAL')
    logging.getLogger('asyncssh').setLevel('CRITICAL')
    logging.getLogger('distributed.worker').setLevel('CRITICAL')
    logging.getLogger('distributed.comm').setLevel('CRITICAL')
    logging.getLogger('distributed.nanny').setLevel('CRITICAL')
    logging.getLogger('bockeh').setLevel('CRITICAL')
    logging.getLogger('matplotlib.font_manager').setLevel('CRITICAL')

    logging.info(args)
    try:
        start(wdir_to_run=args.wdir_to_run, wdir_output=wdir, tpr=tpr,
          xtc=xtc, step=args.step, append_protein_selection=args.append_protein_selection,
          protein_selection=args.protein_selection, ligand_resid=args.ligand, hostfile=args.hostfile, ncpu=args.ncpu,
          n_jobs=args.n_jobs, occupancy=args.occupancy, plot_width=args.width, plot_height=args.height,
          save_viz=not args.not_save_pics, unique_id=unique_id, pdb=pdb,
          verbose=args.verbose, binding_site_cutoff=args.binding_site_cutoff,
          parallel_strategy=None if args.parallel_strategy == 'auto' else args.parallel_strategy,
          ligand_sdf=args.ligand_sdf,
          water_bridge=args.water_bridge, water_selection=args.water_selection,
          water_bridge_order=args.water_bridge_order, water_cutoff=args.water_cutoff,
          show_percentage=not args.no_show_percentage)
    finally:
        logging.shutdown()


if __name__ == '__main__':
    main()
