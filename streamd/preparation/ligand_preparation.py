"""Routines for preparing ligand parameters and inputs."""

import itertools
import logging
import os

from glob import glob
import re

import parmed as pmd

from rdkit import Chem
from rdkit.Chem import rdmolops

from streamd.utils.dask_init import init_dask_cluster, calc_dask
from streamd.utils.utils import run_check_subprocess

# Supported ligand force fields (AmberTools GAFF family). The selected value is
# applied consistently to Antechamber atom typing (``antechamber -at``),
# ``parmchk2 -s`` and the LEaP parameter library sourced in ``tleap.in``.
# Partial charges use AM1-BCC (``antechamber -c bcc``) for standard organic ligands;
# boron-containing ligands are parameterized via Gaussian using RESP charges. The
# charge method is independent of the GAFF/GAFF2 choice.
LIGAND_LEAPRC = {
    'gaff': 'leaprc.gaff',
    'gaff2': 'leaprc.gaff2',
}

# Name of the per-ligand file that records the force field a preparation was run with.
# It is written as soon as preparation commits (before the expensive Antechamber/tleap
# steps) and read before reusing existing files, so a GAFF topology is never silently
# reused for a GAFF2 request (and vice versa) and an interrupted run can be resumed with
# the same force field. Reuse of a *completed* preparation is separately gated on the
# presence of the .itp/posre outputs, so a partial run is never skipped-and-reused.
LIGAND_FORCEFIELD_METADATA = 'ligand_forcefield.txt'

# Generated files moved aside (not deleted) before re-parameterizing a ligand with
# a different force field, so stale artifacts cannot be reused (in particular the
# cached .mol2, which would otherwise make Antechamber atom typing be skipped) while
# the potentially expensive previous parameterization is preserved for inspection.
_STALE_LIGAND_SUFFIXES = ('.mol2', '.mol', '.frcmod', '.itp', '.top', '.gro',
                          '.prmtop', '.inpcrd', '.lib')


def read_prepared_ligand_forcefield(wdir):
    """Return the recorded ligand force field, or ``None``.

    The record is written when preparation commits (before the heavy steps), so a
    missing record means a legacy preparation created before ``--ligand_forcefield``
    was introduced (such topologies were always GAFF).
    """
    path = os.path.join(wdir, LIGAND_FORCEFIELD_METADATA)
    if os.path.isfile(path):
        with open(path) as inp:
            value = inp.read().strip()
        return value or None
    return None


def write_prepared_ligand_forcefield(wdir, ligand_forcefield):
    """Record the force field a ligand directory is (being) prepared with."""
    with open(os.path.join(wdir, LIGAND_FORCEFIELD_METADATA), 'w') as out:
        out.write(f'{ligand_forcefield}\n')


# Files whose presence in a run directory means an MD simulation/analysis was already
# produced there (energy minimization, equilibration, production, analysis), so a
# force-field switch must not silently reuse/mix it.
_RUN_OUTPUT_MARKERS = ('em.gro', 'nvt.gro', 'npt.gro',
                       'md_out.tpr', 'md_out.xtc', 'md_out.gro', 'md_out.cpt',
                       'md_fit.xtc')


def run_has_md_outputs(run_dir):
    """True if a run directory already contains MD simulation/analysis outputs."""
    return any(os.path.isfile(os.path.join(run_dir, m)) for m in _RUN_OUTPUT_MARKERS)


def ensure_run_forcefield_compatible(run_dir, ligand_forcefield):
    """Abort if a run directory already holds MD outputs built with a different force field.

    A run directory that already contains MD simulation/analysis files (energy
    minimization, equilibration, production or analysis outputs) prepared with a
    different force field must not be reused, as it would mix force fields. Run
    directories with no MD outputs yet are allowed (their inputs are refreshed by the
    caller). Legacy run directories (MD outputs but no record) are conservatively
    assumed to be GAFF.

    :raises ValueError: if the run has MD outputs and its recorded/assumed force field
        differs from the requested one.
    """
    if not run_has_md_outputs(run_dir):
        return
    recorded = read_prepared_ligand_forcefield(run_dir) or 'gaff'
    if recorded != ligand_forcefield:
        raise ValueError(
            f'Run directory {run_dir} already contains MD simulation files prepared with the '
            f'"{recorded}" ligand force field, but "{ligand_forcefield}" was requested. Rerun with '
            f'--ligand_forcefield {recorded} to match the existing run, or remove {run_dir} (or use a '
            f'fresh --wdir) to run with "{ligand_forcefield}". Use --wdir_to_continue to extend the existing run.')


def resolve_ligand_forcefield(search_dirs, requested, explicit):
    """Resolve the effective ligand force field for a run.

    When the user did not explicitly request one (``explicit`` is False), adopt the
    force field recorded by a previous run under ``search_dirs`` - so, e.g., a working
    directory prepared with GAFF2 is not silently reverted to the default. If the user
    was explicit, or nothing is recorded, ``requested`` is returned unchanged. Backup
    subdirectories (``backup_*``) are ignored.

    :raises ValueError: if records disagree (ambiguous) and no explicit choice was given.
    """
    if explicit:
        return requested
    recorded = set()
    for base in search_dirs:
        for path in glob(os.path.join(base, '**', LIGAND_FORCEFIELD_METADATA), recursive=True):
            if any(part.startswith('backup_') for part in os.path.normpath(path).split(os.sep)):
                continue  # ignore preserved previous parameterizations
            try:
                with open(path) as inp:
                    value = inp.read().strip()
            except OSError:
                continue
            if value:
                recorded.add(value)
    if not recorded:
        return requested
    if len(recorded) > 1:
        raise ValueError(
            f'Multiple ligand force fields are recorded in this working directory ({sorted(recorded)}). '
            f'Specify --ligand_forcefield explicitly to choose one.')
    adopted = recorded.pop()
    if adopted != requested:
        logging.info(f'Adopting previously used ligand force field "{adopted}" '
                     f'(pass --ligand_forcefield to override).')
    return adopted


def _ligand_artifact_paths(wdir, molid):
    """Regenerable ligand artifacts (excluding the force-field metadata record).

    Includes intermediate files (``.mol2``/``.frcmod``/``tleap.in`` ...) so a
    partial/interrupted preparation is detected, not only completed topologies.
    """
    paths = [os.path.join(wdir, f'{molid}{suffix}') for suffix in _STALE_LIGAND_SUFFIXES]
    paths.append(os.path.join(wdir, f'posre_{molid}.itp'))
    paths.append(os.path.join(wdir, 'tleap.in'))
    return paths


def has_existing_ligand_files(wdir, molid):
    """True if any artifact from a previous (complete or partial) preparation exists."""
    return any(os.path.isfile(p) for p in _ligand_artifact_paths(wdir, molid))


def backup_stale_ligand_files(wdir, molid, prepared_forcefield):
    """Move a previous ligand parameterization aside before re-parameterizing.

    Files are moved (never deleted) into a ``backup_<prepared_forcefield>``
    subdirectory of the ligand's own working directory, so an expensive previous
    parameterization (e.g. AM1-BCC or Gaussian/RESP derived) is preserved and can
    be inspected or restored. Returns the backup directory path (or ``None`` if
    there was nothing to move).
    """
    to_move = _ligand_artifact_paths(wdir, molid) + [os.path.join(wdir, LIGAND_FORCEFIELD_METADATA)]
    existing = [p for p in to_move if os.path.isfile(p)]
    if not existing:
        return None

    # pick a non-colliding backup directory name (repeated switches keep every snapshot)
    backup_dir = os.path.join(wdir, f'backup_{prepared_forcefield}')
    suffix_n = 1
    while os.path.exists(backup_dir):
        backup_dir = os.path.join(wdir, f'backup_{prepared_forcefield}_{suffix_n}')
        suffix_n += 1
    os.makedirs(backup_dir)

    for path in existing:
        try:
            os.rename(path, os.path.join(backup_dir, os.path.basename(path)))
        except OSError as e:
            # A stale file that cannot be moved would be left in place and reused,
            # silently combining the previous force field's atom types/charges with a
            # different parameter set. This is a correctness failure, so abort rather
            # than continue with a partially moved (inconsistent) directory.
            raise RuntimeError(
                f'Failed to move stale ligand file {path} to {backup_dir}: {e}. Aborting so the '
                f'previous "{prepared_forcefield}" parameterization of {molid} is not partially reused '
                f'with a different ligand force field. Free or remove {path} (or the whole {molid} '
                f'directory) and rerun.') from e
    logging.warning(f'Previous "{prepared_forcefield}" parameterization of {molid} was moved to {backup_dir}')
    return backup_dir


def reorder_hydrogens(mol):
    """Move hydrogens to follow their heavy atom in atom order."""
    new_order = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 1:  # Not a hydrogen
            new_order.append(atom.GetIdx())
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 1:  # Is a hydrogen
                    new_order.append(neighbor.GetIdx())
    mol = rdmolops.RenumberAtoms(mol, new_order)
    return mol

def supply_mols_tuple(fname, preset_resid=None, protein_resid_set=None):
    """Yield RDKit molecules with generated residue identifiers."""
    def generate_resid(protein_resid_set):
        ascii_uppercase_digits = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
        for i in itertools.product(ascii_uppercase_digits, repeat=3):
            resid = ''.join(i)
            if resid != 'UNL' and (protein_resid_set is None or resid not in protein_resid_set):
                yield resid

    def add_ids(mol, n, input_fname, resid):
        mol.SetProp('resid', resid)
        if not mol.HasProp('_Name') or not mol.GetProp('_Name'):
            mol.SetProp('_Name', f'{input_fname}_{n}')
        return mol

    resid_generator = generate_resid(protein_resid_set=protein_resid_set)

    if fname.endswith('.sdf'):
        for n, mol in enumerate(Chem.SDMolSupplier(fname, removeHs=False), 1):
            if mol:
                if preset_resid is None:
                    resid = next(resid_generator)
                else:
                    resid = preset_resid

                mol = add_ids(mol, n, input_fname=os.path.basename(fname).rstrip('.sdf'), resid=resid)
                yield (mol, mol.GetProp('_Name'), resid)

    if fname.endswith('.mol') or fname.endswith('.mol2'):
        if fname.endswith('.mol2'):
            mol = pmd.load_file(fname).to_structure()
            mol = mol.rdkit_mol
        else:
            mol = Chem.MolFromMolFile(fname, removeHs=False)

        if mol:
            if preset_resid is None:
                resid = next(resid_generator)
            else:
                resid = preset_resid
            mol = add_ids(mol, n=1, input_fname=os.path.basename(fname).rstrip('.mol2').rstrip('.mol'), resid=resid)
            yield (mol, mol.GetProp('_Name'), resid)


def check_mols(fname):
    """Return number of molecules and list of problematic ones.
    :param fname:
    :return: number of mols, list of problem_mols
    """
    def check_if_problem(mol, n):
        molid = mol.GetProp('_Name') if mol.HasProp('_Name') else n
        try:
            Chem.SanitizeMol(mol)
        except:
            return molid
        return None
        
    problem_mols = []
    number_of_mols = 0
    if fname.endswith('.sdf'):
        n = 0
        for n, mol in enumerate(Chem.SDMolSupplier(fname, removeHs=False, sanitize=False), 1):
            problem_molid = check_if_problem(mol, n)
            if problem_molid:
                problem_mols.append(problem_molid)
        number_of_mols = n

    if fname.endswith('.mol') or fname.endswith('.mol2'):
        if fname.endswith('.mol2'):
            mol = pmd.load_file(fname).to_structure()
            mol = mol.rdkit_mol
        else:
            mol = Chem.MolFromMolFile(fname, removeHs=False, sanitize=False)
        problem_molid = check_if_problem(mol, 1)
        if problem_molid:
            problem_mols.append(problem_molid)
        number_of_mols = 1

    return number_of_mols, problem_mols


def make_all_itp(fileitp_input_list, fileitp_output_list, out_file):
    """Merge ligand ITP files and collect unique atom types."""
    atom_type_list = []
    start_columns = None
    # '[ atomtypes ]\n; name    at.num    mass    charge ptype  sigma      epsilon\n'
    for itp_input, itp_output in zip(fileitp_input_list, fileitp_output_list):
        with open(itp_input) as input:
            data = input.read()
        start = data.find('[ atomtypes ]')
        end = data.find('[ moleculetype ]') - 1
        atom_type_list.extend(data[start:end].split('\n')[2:])
        if start_columns is None:
            start_columns = data[start:end].split('\n')[:2]
        new_data = data[:start] + data[end + 1:]
        with open(itp_output, 'w') as output:
            output.write(new_data)

    atom_type_uniq = [i for i in set(atom_type_list) if i]
    with open(out_file, 'w') as output:
        output.write('\n'.join(start_columns) + '\n')
        output.write('\n'.join(atom_type_uniq) + '\n')

def prepare_tleap(tleap_template, tleap, molid, conda_env_path, ligand_forcefield='gaff'):
    """Fill a tleap template with ligand identifiers, environment path and force field.

    ``ligand_forcefield`` selects the LEaP parameter library sourced by tleap.
    Only validated values from :data:`LIGAND_LEAPRC` are inserted into the
    template (no unrestricted text interpolation).
    """
    if ligand_forcefield not in LIGAND_LEAPRC:
        raise ValueError(f'Unsupported ligand force field: {ligand_forcefield}. '
                         f'Choose one of: {", ".join(LIGAND_LEAPRC)}')
    leaprc = LIGAND_LEAPRC[ligand_forcefield]
    with open(tleap_template) as inp:
        data = inp.read()
    # order matters: substitute the leaprc placeholder before the generic
    # 'ligand' -> molid replacement so the mapped value is left untouched
    new_data = (data.replace('env_path', conda_env_path)
                    .replace('leaprc_ligff', leaprc)
                    .replace('ligand', molid))
    with open(tleap, 'w') as output:
        output.write(new_data)


def prepare_gaussian_files(file_template, file_out, ncpu, opt_restart=False, gaussian_basis=r'B3LYP/6-31G*',
                           gaussian_memory='60GB'):
    """Create Gaussian input from template adjusting resources and options."""
    with open(file_template) as inp:
        data = inp.read()
    standard_basis = r'B3LYP/6-31G\*'
    new_data = re.sub('%NProcShared=[0-9]*', f'%NProcShared={ncpu}', data)
    new_data = re.sub('%Mem=[0-9a-zA-Z]*', f'%Mem={gaussian_memory}', new_data)
    if 'SCF=' not in new_data:
        new_data = re.sub(standard_basis, f'{gaussian_basis} SCF=XQC', new_data)
    else:
        new_data = re.sub(standard_basis, gaussian_basis, new_data)
    if opt_restart and 'Opt' in new_data and 'Opt=Restart' not in new_data:
        new_data = re.sub('Opt', 'Opt=Restart', new_data)

    with open(file_out, 'w') as output:
        output.write(new_data)


def prep_ligand(mol_tuple, script_path, project_dir, wdir_ligand,
                conda_env_path, bash_log, gaussian_exe=None,  no_dr=False,
                activate_gaussian=None, gaussian_basis='B3LYP/6-31G*', gaussian_memory='60GB', ncpu=1,
                mol2_file=None, env=None, ligand_forcefield='gaff'):
    """Prepare force-field parameters for a single ligand.

    :param ligand_forcefield: AmberTools parameter set (``gaff`` or ``gaff2``)
        applied to Antechamber atom typing, ``parmchk2`` and LEaP. Charges use
        AM1-BCC, except boron-containing ligands, which are parameterized via
        Gaussian using RESP. Default: ``gaff``.
    """
    if ligand_forcefield not in LIGAND_LEAPRC:
        raise ValueError(f'Unsupported ligand force field: {ligand_forcefield}. '
                         f'Choose one of: {", ".join(LIGAND_LEAPRC)}')
    mol, molid, resid = mol_tuple

    wdir_ligand_cur = os.path.abspath(os.path.join(wdir_ligand, molid))
    os.makedirs(wdir_ligand_cur, exist_ok=True)

    itp_file = os.path.join(wdir_ligand_cur, f'{molid}.itp')
    posre_file = os.path.join(wdir_ligand_cur, f'posre_{molid}.itp')

    # Determine the force field of any existing preparation - complete OR partial. The
    # requested force field is recorded (below) as soon as preparation commits, BEFORE
    # the heavy Antechamber/tleap steps, so an interrupted run is attributable to its
    # force field and can be resumed (or correctly rejected on a switch) rather than
    # guessed at. A missing record therefore means a legacy preparation created before
    # --ligand_forcefield existed, which was always GAFF.
    # The record is READ here, before it is (re)written, so a genuine force-field switch
    # is still detected as a mismatch. This check must run before the completed-output
    # reuse below, otherwise a stale {molid}.mol2 from an interrupted run would be reused
    # (Antechamber skipped) and silently combined with a different parmchk2 -s / leaprc set.
    recorded_forcefield = read_prepared_ligand_forcefield(wdir_ligand_cur)
    prepared_forcefield = recorded_forcefield or (
        'gaff' if has_existing_ligand_files(wdir_ligand_cur, molid) else None)

    if prepared_forcefield is not None and prepared_forcefield != ligand_forcefield:
        # Existing files were built with a different force field. Move them aside so
        # nothing stale is reused (in particular a mismatched {molid}.mol2 that would
        # make Antechamber atom typing be skipped), then regenerate from scratch.
        logging.warning(
            f'Existing {molid} files were prepared with the "{prepared_forcefield}" ligand force field, but '
            f'"{ligand_forcefield}" was requested. They will be regenerated with "{ligand_forcefield}" '
            f'(this re-runs Antechamber atom typing, parmchk2 and LEaP); the previous files are moved to a '
            f'backup subdirectory rather than deleted.')
        backup_stale_ligand_files(wdir_ligand_cur, molid, prepared_forcefield)
    elif os.path.isfile(itp_file) and os.path.isfile(posre_file):
        # A complete parameterization with a matching force field exists - reuse it.
        logging.warning(
            f'{molid}.itp and posre_{molid}.itp files already exist and were prepared with the '
            f'"{ligand_forcefield}" ligand force field. Mol preparation step will be skipped for such molecule')
        if not os.path.isfile(os.path.join(wdir_ligand_cur, 'resid.txt')):
            with open(os.path.join(wdir_ligand_cur, 'resid.txt'), 'w') as out:
                out.write(f'{molid}\t{resid}\n')
        return wdir_ligand_cur
    # else: nothing exists yet, or a partial preparation with the SAME force field - fall
    # through and (re)prepare. A consistent same-force-field {molid}.mol2 left by an
    # interrupted run may be reused below, which is a safe resume.

    # Record the requested force field now, before the (expensive) Antechamber/tleap
    # steps. If preparation is interrupted, a later run with the same --ligand_forcefield
    # resumes from the partial files, while a run with a different force field reads this
    # value first (above) and backs the partial work up instead of reusing it.
    write_prepared_ligand_forcefield(wdir_ligand_cur, ligand_forcefield)

    if not mol2_file or not os.path.isfile(mol2_file):
        mol2_file = os.path.join(wdir_ligand_cur, f'{molid}.mol2')
        mol_file = os.path.join(wdir_ligand_cur, f'{molid}.mol')
        # if addH:
        mol = Chem.AddHs(mol, addCoords=True)
        # reorder hydrogens for Gromacs GPU update functionality
        mol = reorder_hydrogens(mol)

        Chem.MolToMolFile(mol, mol_file)

        charge = rdmolops.GetFormalCharge(mol)

        # generate mol2
        if not os.path.isfile(mol2_file):
            # boron atom
            if mol.HasSubstructMatch(Chem.MolFromSmarts("[#5]")):
                if gaussian_exe:
                    for file in glob(os.path.join(script_path, 'com', '*.com')):
                        prepare_gaussian_files(file_template=file,
                                               file_out=os.path.join(wdir_ligand_cur, os.path.basename(file)),
                                               ncpu=ncpu,
                                               gaussian_basis=gaussian_basis, gaussian_memory=gaussian_memory)
                    cmd = f'script_path={script_path} lfile={mol_file} input_dirname={wdir_ligand_cur} ' \
                          f'resid={resid} molid={molid} charge={charge} gaussian_version={gaussian_exe} ' \
                          f'ligand_forcefield="{ligand_forcefield}" ' \
                          f'activate_gaussian="{activate_gaussian if activate_gaussian else ""}" ' \
                          f'bash {os.path.join(script_path, "script_sh", "ligand_mol2prep_by_gaussian.sh")} ' \
                          f' >> {os.path.join(wdir_ligand_cur, bash_log)} 2>&1'
                    if not run_check_subprocess(cmd, molid, log=os.path.join(wdir_ligand_cur, bash_log), env=env):
                        return None
                else:
                    return None
            else:
                cmd = f'script_path={script_path} lfile={mol_file} input_dirname={wdir_ligand_cur} ' \
                      f'resid={resid} molid={molid} charge={charge} ligand_forcefield="{ligand_forcefield}" dr=yes bash {os.path.join(script_path, "script_sh", "ligand_mol2prep.sh")} ' \
                      f' >> {os.path.join(wdir_ligand_cur, bash_log)} 2>&1',
                if not run_check_subprocess(cmd, molid, log=os.path.join(wdir_ligand_cur, bash_log), env=env,
                                            ignore_error=True if no_dr else False):
                    if not no_dr:
                        return None
                    else:
                        logging.warning(f'Ambertools structure checking returned an error for the {mol_file} file.'
                                        f'Check the input structure carefully. Continue with -dr no mode.')
                        cmd = f'script_path={script_path} lfile={mol_file} input_dirname={wdir_ligand_cur} ' \
                          f'resid={resid} molid={molid} charge={charge} ligand_forcefield="{ligand_forcefield}" dr=no bash {os.path.join(script_path, "script_sh","ligand_mol2prep.sh")} ' \
                          f' >> {os.path.join(wdir_ligand_cur, bash_log)} 2>&1',
                        if not run_check_subprocess(cmd, molid, log=os.path.join(wdir_ligand_cur, bash_log), env=env):
                            return None
    else:
        mol2 = pmd.load_file(mol2_file).to_structure()
        mol2.residues[0].name = 'UNL'
        mol2.save(os.path.join(wdir_ligand_cur, f'{molid}.mol2'))
        logging.info(f'No mol2 file will be generated. {mol2_file} will be used instead')
        # Antechamber is skipped for a pre-existing MOL2 input, so its atom types and
        # partial charges are NOT reassigned. Only parmchk2 (-s) and LEaP (leaprc.*)
        # use the selected force field; the MOL2 is trusted to already match it.
        logging.warning(
            f'The supplied MOL2 file {mol2_file} is used as-is: Antechamber is skipped for '
            f'pre-existing MOL2 inputs, so its atom types and partial charges are NOT reassigned '
            f'(AM1-BCC is not applied). They are trusted to already be consistent with the selected '
            f'ligand force field "{ligand_forcefield}". Only parmchk2 (-s {ligand_forcefield}) and '
            f'LEaP (source {LIGAND_LEAPRC[ligand_forcefield]}) are applied. If the MOL2 atom types do '
            f'not correspond to {ligand_forcefield}, provide a .mol/.sdf input instead so that '
            f'Antechamber can assign {ligand_forcefield} atom types.')

    prepare_tleap(os.path.join(script_path, 'tleap.in'), tleap=os.path.join(wdir_ligand_cur, 'tleap.in'),
                  molid=molid, conda_env_path=conda_env_path, ligand_forcefield=ligand_forcefield)
    cmd = f'script_path={script_path} input_dirname={wdir_ligand_cur} ' \
          f'molid={molid} ligand_forcefield="{ligand_forcefield}" bash {os.path.join(script_path, "script_sh","ligand_prep.sh")} ' \
          f' >> {os.path.join(wdir_ligand_cur, bash_log)} 2>&1'
    if not run_check_subprocess(cmd, molid, log=os.path.join(wdir_ligand_cur, bash_log), env=env):
        return None

    # create log for molid resid corresponding
    with open(os.path.join(wdir_ligand_cur, 'resid.txt'), 'w') as out:
        out.write(f'{molid}\t{resid}\n')
    return wdir_ligand_cur


def prepare_input_ligands(ligand_fname, preset_resid, protein_resid_set, script_path, project_dir, wdir_ligand,
                          no_dr, gaussian_exe, activate_gaussian, gaussian_basis, gaussian_memory,
                          hostfile, ncpu, bash_log, ligand_forcefield='gaff'):
    """Prepare parameterization inputs for multiple ligands.
     :param ligand_fname:
    :param preset_resid:
    :param script_path:
    :param project_dir:
    :param wdir_ligand:
    :param no_dr:
    :param gaussian_exe: str or None
    :param activate_gaussian: str or None
    :param hostfile:
    :param ncpu:
    :param bash_log:
    :param ligand_forcefield: AmberTools parameter set (``gaff`` or ``gaff2``)
        applied to Antechamber atom typing, ``parmchk2`` and LEaP. Default: ``gaff``.
    :return:
    """
    lig_wdirs = []

    if ligand_fname.endswith('.mol2'):
        mol_tuple = next(supply_mols_tuple(ligand_fname, preset_resid=preset_resid, protein_resid_set=protein_resid_set))
        res = prep_ligand(mol_tuple=mol_tuple, script_path=script_path,
            project_dir=project_dir, wdir_ligand=wdir_ligand,
            conda_env_path=os.environ["CONDA_PREFIX"],
            ncpu = ncpu, mol2_file = ligand_fname,
            bash_log=bash_log, env = os.environ.copy(),
            ligand_forcefield=ligand_forcefield)

        if res:
            lig_wdirs.append(res)

    else:
        standard_mols, boron_containing_mols = [], []
        for mol_tuple in supply_mols_tuple(ligand_fname, preset_resid=preset_resid, protein_resid_set=protein_resid_set):
            mol = mol_tuple[0]
            if mol.HasSubstructMatch(Chem.MolFromSmarts("[#5]")):
                boron_containing_mols.append(mol_tuple)
            else:
                standard_mols.append(mol_tuple)

        dask_client, cluster = None, None
        # prepare boron-containig mols
        if boron_containing_mols:
            if gaussian_exe:
                try:
                    dask_client, cluster = init_dask_cluster(hostfile=hostfile,
                                                             n_tasks_per_node=1,
                                                             use_multi_servers=True if len(boron_containing_mols) > 1 else False,
                                                             ncpu=ncpu)
                    for res in calc_dask(prep_ligand, boron_containing_mols, dask_client,
                                         script_path=script_path, project_dir=project_dir,
                                         wdir_ligand=wdir_ligand, conda_env_path=os.environ["CONDA_PREFIX"],
                                         gaussian_exe=gaussian_exe, activate_gaussian=activate_gaussian,
                                         gaussian_basis=gaussian_basis, gaussian_memory=gaussian_memory,
                                         ncpu=ncpu, bash_log=bash_log, no_dr=no_dr,
                                         env=os.environ.copy(), ligand_forcefield=ligand_forcefield):
                        if res:
                            lig_wdirs.append(res)
                finally:
                    if dask_client:
                        dask_client.shutdown()
                    if cluster:
                        cluster.close()
            else:
                logging.warning(
                    f'There are molecules from {ligand_fname} which have Boron atom and to prepare such molecules you need to set up Gaussian.'
                    f' Please restart the run again and use --gaussian_exe arguments')

        if standard_mols:
            try:
                dask_client, cluster = init_dask_cluster(hostfile=hostfile,
                                                         n_tasks_per_node=min(ncpu, len(standard_mols)),
                                                         use_multi_servers= True if len(standard_mols) > ncpu else False,
                                                         ncpu=ncpu)
                for res in calc_dask(prep_ligand, standard_mols, dask_client,
                                     script_path=script_path, project_dir=project_dir,
                                     wdir_ligand=wdir_ligand, no_dr=no_dr,
                                     conda_env_path=os.environ["CONDA_PREFIX"],
                                     ncpu=ncpu, bash_log=bash_log,
                                     env=os.environ.copy(), ligand_forcefield=ligand_forcefield):
                    if res:
                        lig_wdirs.append(res)
            finally:
                if dask_client:
                    dask_client.shutdown()
                if cluster:
                    cluster.close()

    return lig_wdirs
