"""Convert GROMACS gro/xtc trajectories to pdb/dcd pairs for OpenMMDL analysis.

OpenMMDL Analysis (PLIP-based) expects a **PDB topology** and a coordinate-only
**DCD trajectory**. A DCD carries no topology, so the PDB and DCD are written from
the *same* atom selection and atom order to keep them consistent.

Two ways to drive it:

* ``--i`` : batch mode. Each item is either a StreaMD working directory (its
  ``md_fit.xtc`` is converted to ``md_fit.pdb``/``md_fit.dcd`` in place) or a
  trajectory file (converted to ``<name>.pdb``/``<name>.dcd`` next to it).
* ``-f``/``-o`` : explicit single conversion with a chosen output prefix.

A DCD stores coordinates only, so a topology is always required. When ``-s`` is not
given it is auto-detected next to each trajectory (``md_out.tpr`` -> ``md_out.gro``,
or ``<name>.tpr``/``<name>.gro`` for a user file).

The input trajectory should already be PBC-treated (whole molecules, protein
centered/fitted) - in a StreaMD directory that is ``md_fit.xtc``. This script does
no imaging; it only changes the file format.
"""

import argparse
import logging
import os
import sys
import warnings
from functools import partial

import MDAnalysis as mda

from streamd.utils.utils import filepath_type

# Trajectory converted for a StreaMD working directory, and the topologies tried
# for it (deffnm is 'md_out'; md_fit.xtc is the full-System fitted trajectory).
STREAMD_TRAJ = 'md_fit.xtc'
STREAMD_TOPOLOGIES = ('md_out.tpr', 'md_out.gro')


def _load_universe(topology, trajectory):
    """Load a universe, falling back from an unparseable TPR to a matching GRO.

    Newer GROMACS ``.tpr`` versions are not always readable by MDAnalysis; when the
    TPR fails and a sibling ``.gro`` exists, retry with the GRO (same behaviour as
    the analysis module). Any other error is propagated unchanged.
    """
    try:
        return mda.Universe(topology, trajectory)
    except ValueError as exc:
        gro = f'{os.path.splitext(topology)[0]}.gro'
        if os.path.splitext(topology)[1] == '.tpr' and os.path.isfile(gro):
            logging.warning('Cannot parse %s with MDAnalysis (%s); using %s instead',
                            topology, exc, gro)
            return mda.Universe(gro, trajectory)
        raise


def _ensure_elements(universe):
    """Guess element names when the topology lacks them, so the PDB has an element column.

    GRO topologies carry no element information; PLIP/RDKit rely on the PDB element
    column, so guess it when missing. Falls back gracefully across MDAnalysis
    versions and leaves the writer to infer from atom names if guessing fails.
    """
    if hasattr(universe.atoms, 'elements'):
        return
    try:
        # MDAnalysis >= 2.8 guesser API
        universe.guess_TopologyAttrs(to_guess=['elements'])
    except (AttributeError, ValueError, TypeError):
        try:
            from MDAnalysis.topology.guessers import guess_types
            universe.add_TopologyAttr('elements', guess_types(universe.atoms.names))
        except Exception as exc:  # noqa: BLE001 - guessing is best-effort
            logging.warning('Could not guess element names; PDB element column may be '
                            'inferred from atom names by the writer. %s', exc)


def _resolve_topology(trajectory, out_pdb, explicit=None, extra_names=()):
    """Return a topology path for ``trajectory``, or ``None`` if none is found.

    ``explicit`` (``-s``) always wins. Otherwise candidates are tried in order:
    same-basename siblings (``<name>.tpr``/``.gro``/``.pdb``) then ``extra_names``
    (e.g. StreaMD's ``md_out.*``) in the trajectory's directory. The tool's own
    output PDB is never used as its topology.
    """
    if explicit:
        return explicit
    directory = os.path.dirname(trajectory)
    base = os.path.splitext(os.path.basename(trajectory))[0]
    candidates = [f'{base}.tpr', f'{base}.gro', f'{base}.pdb', *extra_names]
    for name in candidates:
        candidate = os.path.join(directory, name)
        if os.path.abspath(candidate) == os.path.abspath(out_pdb):
            continue  # never treat our own output as the topology
        if os.path.isfile(candidate):
            return candidate
    return None


def _resolve_target(item, topology_override):
    """Resolve one ``--i`` item to ``(trajectory, topology, out_prefix)``.

    ``item`` is a StreaMD working directory (uses ``md_fit.xtc`` -> ``md_fit.*``) or a
    trajectory file (uses ``<name>.*`` next to it). Raises ``FileNotFoundError`` with a
    clear message when the trajectory or a matching topology cannot be located.
    """
    item = os.path.abspath(item)
    if os.path.isdir(item):
        trajectory = os.path.join(item, STREAMD_TRAJ)
        if not os.path.isfile(trajectory):
            raise FileNotFoundError(
                f'{item}: no {STREAMD_TRAJ} found in this directory')
        extra_names = STREAMD_TOPOLOGIES
    elif os.path.isfile(item):
        trajectory = item
        extra_names = STREAMD_TOPOLOGIES
    else:
        raise FileNotFoundError(f'{item}: not an existing directory or file')

    out_prefix = os.path.splitext(trajectory)[0]
    out_pdb = f'{out_prefix}.pdb'
    topology = _resolve_topology(trajectory, out_pdb,
                                 explicit=topology_override, extra_names=extra_names)
    if topology is None:
        raise FileNotFoundError(
            f'{trajectory}: could not auto-detect a topology (looked for '
            f'<name>.tpr/.gro/.pdb and {", ".join(STREAMD_TOPOLOGIES)}). '
            f'Pass one with -s/--topology.')
    if os.path.abspath(topology) == os.path.abspath(out_pdb):
        raise FileNotFoundError(
            f'{trajectory}: the only topology found is the output path {out_pdb}; '
            f'pass a different topology with -s/--topology.')
    return trajectory, topology, out_prefix


def convert(topology, trajectory, out_pdb, out_dcd, selection='all',
            start=None, stop=None, step=None):
    """Write a topology PDB (first written frame) and a DCD trajectory from a gro/xtc pair.

    :param topology: Path to the topology/coordinate file (``.gro``/``.pdb``/``.tpr``).
    :param trajectory: Path to the trajectory file (``.xtc``/``.trr``).
    :param out_pdb: Output PDB path used as the OpenMMDL topology.
    :param out_dcd: Output DCD path used as the OpenMMDL trajectory.
    :param selection: MDAnalysis atom selection to keep (default ``'all'``).
    :param start: First frame index to write (inclusive). ``None`` means the first frame.
    :param stop: Stop frame index (exclusive). ``None`` means through the last frame.
    :param step: Frame stride. ``None`` means every frame.
    :return: Tuple ``(out_pdb, out_dcd, n_atoms, n_frames)``.
    """
    universe = _load_universe(topology, trajectory)
    _ensure_elements(universe)

    atomgroup = universe.select_atoms(selection)
    if atomgroup.n_atoms == 0:
        raise ValueError(f'Selection {selection!r} matched no atoms in {topology}')

    n_frames = 0
    wrote_topology = False
    with mda.Writer(out_dcd, n_atoms=atomgroup.n_atoms) as dcd:
        for _ in universe.trajectory[start:stop:step]:
            # The DCD stores coordinates only; the PDB (first written frame) defines
            # the atom order/identity the DCD is read against, so write both from the
            # same atomgroup on the same frame to guarantee they match.
            if not wrote_topology:
                atomgroup.write(out_pdb)
                wrote_topology = True
            dcd.write(atomgroup)
            n_frames += 1

    if not wrote_topology:
        raise ValueError('No frames were written; check --start/--stop/--step against '
                         f'the {universe.trajectory.n_frames}-frame trajectory')

    logging.info('Wrote %s (%d atoms) and %s (%d frames) from %s [selection %r]',
                 out_pdb, atomgroup.n_atoms, out_dcd, n_frames, trajectory, selection)
    return out_pdb, out_dcd, atomgroup.n_atoms, n_frames


def main():
    parser = argparse.ArgumentParser(
        description='Convert GROMACS gro/xtc to pdb/dcd for OpenMMDL analysis. PDB and '
                    'DCD are written from the same selection so their atom order matches. '
                    'Feed an already PBC-treated/fitted trajectory (StreaMD md_fit.xtc).')
    parser.add_argument('--i', nargs='*', metavar='DIR_OR_XTC', default=None,
                        help='Batch mode. StreaMD directories (converts md_fit.xtc -> '
                             'md_fit.pdb/.dcd in place) and/or trajectory files '
                             '(converts <name>.xtc -> <name>.pdb/.dcd next to it).')
    parser.add_argument('-f', '--trajectory', metavar='md_fit.xtc', default=None,
                        type=partial(filepath_type, ext=('xtc', 'trr')),
                        help='Explicit single trajectory (alternative to --i). '
                             'Preferably already centered/fitted.')
    parser.add_argument('-o', '--output', metavar='PREFIX', default=None,
                        help='Output prefix for -f mode; writes PREFIX.pdb and PREFIX.dcd. '
                             'Defaults to the trajectory path without extension.')
    parser.add_argument('-s', '--topology', metavar='md_out.tpr', default=None,
                        type=partial(filepath_type, ext=('gro', 'pdb', 'tpr')),
                        help='Topology override applied to all inputs. When omitted it is '
                             'auto-detected next to each trajectory.')
    parser.add_argument('--selection', metavar='SELECTION', default='all',
                        help="MDAnalysis atom selection to keep (default: 'all'). "
                             "E.g. 'not resname SOL NA CL' to drop water/ions, or "
                             "'protein or resname UNL'.")
    parser.add_argument('--start', metavar='INT', type=int, default=None,
                        help='First frame index to write (inclusive).')
    parser.add_argument('--stop', metavar='INT', type=int, default=None,
                        help='Stop frame index (exclusive).')
    parser.add_argument('--step', metavar='INT', type=int, default=None,
                        help='Write every STEP-th frame (subsample the trajectory).')
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    # Quiet MDAnalysis' own per-file chatter: attribute-guessing INFO lines and the
    # PDBWriter "Found no information for attr" / "missing chainIDs" warnings are all
    # expected when a topology comes from GROMACS (a gro has no altLocs/chainIDs/
    # record_types), and drown out this tool's messages in batch mode.
    logging.getLogger('MDAnalysis').setLevel(logging.WARNING)
    warnings.filterwarnings('ignore', message='Found no information for attr',
                            category=UserWarning)
    warnings.filterwarnings('ignore', message='Found missing chainIDs',
                            category=UserWarning)

    if (args.i is not None) == (args.trajectory is not None):
        parser.error('provide exactly one of --i (batch) or -f/--trajectory (single).')

    # Resolve every requested conversion to (trajectory, topology, out_prefix).
    targets = []
    failures = 0
    if args.i is not None:
        if not args.i:
            parser.error('--i was given without any directories or files.')
        for item in args.i:
            try:
                targets.append(_resolve_target(item, args.topology))
            except FileNotFoundError as exc:
                logging.error('Skipping %s', exc)
                failures += 1
    else:
        trajectory = args.trajectory
        out_prefix = os.path.splitext(args.output or trajectory)[0]
        topology = _resolve_topology(trajectory, f'{out_prefix}.pdb',
                                     explicit=args.topology, extra_names=STREAMD_TOPOLOGIES)
        if topology is None:
            parser.error(
                f'{trajectory}: could not auto-detect a topology; pass one with -s/--topology.')
        targets.append((trajectory, topology, out_prefix))

    for trajectory, topology, out_prefix in targets:
        try:
            convert(topology=topology, trajectory=trajectory,
                    out_pdb=f'{out_prefix}.pdb', out_dcd=f'{out_prefix}.dcd',
                    selection=args.selection,
                    start=args.start, stop=args.stop, step=args.step)
        except (ValueError, OSError) as exc:
            logging.error('Failed to convert %s: %s', trajectory, exc)
            failures += 1

    if failures:
        sys.exit(1)


if __name__ == '__main__':
    main()
