"""Tests for the ProLIF issue #358 backbone-nitrogen H-bond acceptor fix.

See https://github.com/chemosim-lab/ProLIF/issues/358. ProLIF matches its H-bond
SMARTS on the protein split into single-residue fragments. Cutting the peptide bond
severs the backbone nitrogen from the previous residue's carbonyl and - because the
fragments are not re-sanitized - the amide nitrogen keeps its SP2 state while losing
the connectivity that made ProLIF's amide exclusion fire. It then wrongly matches as
an H-bond acceptor. The fix ANDs an extra exclusion into the aliphatic-nitrogen
branch of the acceptor SMARTS that is specific to that fragmented state.

The suite covers three levels:
  * string surgery on the SMARTS (fast, no chem deps);
  * the actual chemical (mis)classification, using an RDKit fragment reconstructed
    the same way ProLIF fragments residues - this is the real issue #358 regression;
  * ProLIF integration (the override is accepted and removes the false contact).
"""
import pytest

from streamd.prolif.run_prolif import (
    BACKBONE_N_ACCEPTOR_EXCLUSION,
    inject_backbone_n_acceptor_exclusion,
    fixed_hbond_acceptor_smarts,
)

# The stock ProLIF 2.1.0 H-bond acceptor SMARTS, used as a representative input so the
# string/chemistry tests do not depend on the installed ProLIF version.
DEFAULT_ACCEPTOR_2_1_0 = (
    "[$([N&!$([NX3]-*=[O,N,P,S])&!$([NX3]-[a])&!$([Nv4+1])&!$(N=C(-[C,N])-N)])"
    ",$([n+0&!X3&!$([n&r5]:[n+&r5])])"
    ",$([O&!$([OX2](C)C=O)&!$(O(~a)~a)&!$(O=N-*)&!$([O-]-N=O)])"
    ",$([o+0])"
    ",$([F&$(F-[#6])&!$(F-[#6][F,Cl,Br,I])])]"
)


def _split_top_level_branches(smarts):
    """Split the outer atom query into its top-level comma-separated branches."""
    body = smarts.strip()[1:-1]
    branches, depth, start = [], 0, 0
    for i, ch in enumerate(body):
        if ch in '([':
            depth += 1
        elif ch in ')]':
            depth -= 1
        elif ch == ',' and depth == 0:
            branches.append(body[start:i])
            start = i + 1
    branches.append(body[start:])
    return branches


def _fragmented_backbone_residue(embed_3d=False):
    """Reconstruct ProLIF's fragmented backbone residue (the core of issue #358).

    Builds N-acetyl-L-alanine-N'-methylamide (CH3-CO-NH-CH(CH3)-CO-NH-CH3) and cuts
    BOTH flanking peptide (amide C-N) bonds exactly the way
    ``prolif.utils.split_mol_by_residues`` does: ``FragmentOnBonds(addDummies=False)``
    followed by ``GetMolFrags(sanitizeFrags=False)``. Because the fragments are not
    re-sanitized, the central alanine's backbone nitrogen keeps the SP2 / degree-2 /
    valence-3 state (explicit backbone H, degree-2 carbonyl two bonds away) that the
    bug - and the fix - depend on. This reproduces the exact atom state observed on a
    real ProLIF-fragmented protein residue (LEU106 in the issue's 8C2Z test file).

    :param embed_3d: also generate a 3D conformer (needed for a Fingerprint run).
    :return: ``(residue_mol, nitrogen_index)``.
    """
    from rdkit import Chem
    mol = Chem.AddHs(Chem.MolFromSmiles("CC(=O)NC(C)C(=O)NC"))
    if embed_3d:
        from rdkit.Chem import AllChem
        AllChem.EmbedMolecule(mol, randomSeed=0xF00D)
        AllChem.MMFFOptimizeMolecule(mol)
    amide = Chem.MolFromSmarts("[CX3](=O)-[NX3]")
    bonds = [mol.GetBondBetweenAtoms(c, n).GetIdx()
             for c, _o, n in mol.GetSubstructMatches(amide)]
    fragmented = Chem.FragmentOnBonds(mol, bonds, addDummies=False)
    for frag in Chem.GetMolFrags(fragmented, asMols=True, sanitizeFrags=False):
        for atom in frag.GetAtoms():
            # the central residue's backbone N is the only degree-2 nitrogen
            if atom.GetAtomicNum() == 7 and atom.GetDegree() == 2:
                return frag, atom.GetIdx()
    raise AssertionError("failed to reconstruct the fragmented backbone residue")


def _matched_atoms(mol, smarts):
    from rdkit import Chem
    query = Chem.MolFromSmarts(smarts)
    return {idx for match in mol.GetSubstructMatches(query) for idx in match}


# --------------------------------------------------------------------------- #
# 1. SMARTS string surgery
# --------------------------------------------------------------------------- #

def test_exclusion_is_injected_into_nitrogen_branch():
    """The exclusion clause is added inside the aliphatic-nitrogen branch."""
    patched = inject_backbone_n_acceptor_exclusion(DEFAULT_ACCEPTOR_2_1_0)
    assert patched != DEFAULT_ACCEPTOR_2_1_0
    assert BACKBONE_N_ACCEPTOR_EXCLUSION in patched
    n_branch = next(b for b in _split_top_level_branches(patched) if b.startswith('$([N'))
    assert BACKBONE_N_ACCEPTOR_EXCLUSION in n_branch


def test_only_nitrogen_branch_is_modified():
    """Oxygen/fluorine/aromatic branches are untouched and the branch count is kept."""
    before = _split_top_level_branches(DEFAULT_ACCEPTOR_2_1_0)
    after = _split_top_level_branches(inject_backbone_n_acceptor_exclusion(DEFAULT_ACCEPTOR_2_1_0))
    assert len(before) == len(after)
    for b_old, b_new in zip(before, after):
        if b_old.startswith('$([N'):
            assert BACKBONE_N_ACCEPTOR_EXCLUSION in b_new
        else:
            assert b_new == b_old


def test_injection_is_idempotent():
    """Applying the fix twice does not duplicate the exclusion."""
    once = inject_backbone_n_acceptor_exclusion(DEFAULT_ACCEPTOR_2_1_0)
    twice = inject_backbone_n_acceptor_exclusion(once)
    assert once == twice
    assert once.count(BACKBONE_N_ACCEPTOR_EXCLUSION) == 1


def test_returns_unchanged_when_no_nitrogen_branch():
    """A pattern without an aliphatic-nitrogen branch is returned unchanged."""
    only_oxygen = "[$([O&!$([OX2](C)C=O)]),$([o+0])]"
    assert inject_backbone_n_acceptor_exclusion(only_oxygen) == only_oxygen


def test_returns_unchanged_for_empty_or_non_bracketed_input():
    """Empty / None / non-bracketed input is returned unchanged (no crash)."""
    assert inject_backbone_n_acceptor_exclusion("") == ""
    assert inject_backbone_n_acceptor_exclusion(None) is None
    assert inject_backbone_n_acceptor_exclusion("N") == "N"


@pytest.mark.parametrize("smarts", ["[", "[]", "[$([N)", "[$([N])", "[$([N])]", "[N,", "[O]"])
def test_malformed_bracketed_patterns_do_not_crash(smarts):
    """Structurally malformed bracketed inputs never raise (returned as a string)."""
    assert isinstance(inject_backbone_n_acceptor_exclusion(smarts), str)


def test_patched_smarts_is_valid_rdkit_pattern():
    """The patched SMARTS parses as a valid RDKit query."""
    Chem = pytest.importorskip("rdkit.Chem")
    patched = inject_backbone_n_acceptor_exclusion(DEFAULT_ACCEPTOR_2_1_0)
    assert Chem.MolFromSmarts(patched) is not None


# --------------------------------------------------------------------------- #
# 2. Chemical regression: the actual issue #358 misclassification
# --------------------------------------------------------------------------- #

def test_reconstructed_fragment_matches_real_backbone_n_state():
    """Guard: the reconstruction has the SP2 / degree-2 / valence-3 fragmented N state.

    If RDKit ever changes so this no longer holds, the regression assertions below
    could pass for the wrong reason - fail loudly here instead.
    """
    Chem = pytest.importorskip("rdkit.Chem")
    res, n_idx = _fragmented_backbone_residue()
    n = res.GetAtomWithIdx(n_idx)
    assert n.GetHybridization() == Chem.HybridizationType.SP2
    assert n.GetDegree() == 2
    assert n.GetTotalValence() == 3
    assert n.GetFormalCharge() == 0
    assert n.GetTotalNumHs() == 1


def test_fragmented_backbone_n_matches_default_but_not_patched():
    """The heart of issue #358: the stock SMARTS accepts the fragmented backbone N as
    an acceptor (the bug); the patched SMARTS rejects it (the fix)."""
    pytest.importorskip("rdkit.Chem")
    res, n_idx = _fragmented_backbone_residue()
    default_hits = _matched_atoms(res, DEFAULT_ACCEPTOR_2_1_0)
    patched_hits = _matched_atoms(res, inject_backbone_n_acceptor_exclusion(DEFAULT_ACCEPTOR_2_1_0))
    assert n_idx in default_hits       # bug present with the stock pattern
    assert n_idx not in patched_hits   # removed by the fix


@pytest.mark.parametrize(
    "smiles, atomic_num, label",
    [
        ("CC=O", 8, "carbonyl oxygen"),
        ("CC(=O)[O-]", 8, "carboxylate oxygen"),
        ("c1ccncc1", 7, "pyridine aromatic nitrogen"),
        ("Fc1ccccc1", 9, "aromatic organofluorine"),
        ("CF", 9, "aliphatic organofluorine"),
    ],
)
def test_common_acceptors_are_preserved(smiles, atomic_num, label):
    """The patch does not remove genuine, unrelated acceptors."""
    Chem = pytest.importorskip("rdkit.Chem")
    mol = Chem.MolFromSmiles(smiles)
    default_hits = _matched_atoms(mol, DEFAULT_ACCEPTOR_2_1_0)
    patched_hits = _matched_atoms(mol, inject_backbone_n_acceptor_exclusion(DEFAULT_ACCEPTOR_2_1_0))
    target = {a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == atomic_num}
    assert default_hits & target, f"{label} should be an acceptor originally"
    assert patched_hits & target, f"{label} should still be an acceptor after the patch"


@pytest.mark.parametrize("smiles", ["NCC=O", "CNCC=O"])
def test_intact_amine_acceptor_unaffected_by_patch(smiles):
    """Intact (sp3, non-fragmented) amine alpha to a carbonyl is unchanged.

    The exclusion is deliberately specific to the SP2 / degree-2 fragmented state, so
    it does not over-reach onto real ligand chemistry. This is what makes it safe to
    apply the same acceptor override to HBAcceptor (matched against the ligand) as
    well as HBDonor (matched against the protein).
    """
    Chem = pytest.importorskip("rdkit.Chem")
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    default_hits = mol.GetSubstructMatches(Chem.MolFromSmarts(DEFAULT_ACCEPTOR_2_1_0))
    patched_hits = mol.GetSubstructMatches(
        Chem.MolFromSmarts(inject_backbone_n_acceptor_exclusion(DEFAULT_ACCEPTOR_2_1_0)))
    assert default_hits == patched_hits


# --------------------------------------------------------------------------- #
# 3. fixed_hbond_acceptor_smarts() behaviour (deterministic, via explicit input)
# --------------------------------------------------------------------------- #

def test_fixed_helper_patches_explicit_pattern():
    """Given a stock pattern, the helper returns it patched."""
    patched = fixed_hbond_acceptor_smarts(DEFAULT_ACCEPTOR_2_1_0)
    assert patched is not None
    assert BACKBONE_N_ACCEPTOR_EXCLUSION in patched


def test_fixed_helper_returns_none_when_upstream_already_fixed():
    """When the base pattern already carries the exclusion (e.g. ProLIF fixes it
    upstream), the helper returns None so no redundant override is applied."""
    already_fixed = inject_backbone_n_acceptor_exclusion(DEFAULT_ACCEPTOR_2_1_0)
    assert fixed_hbond_acceptor_smarts(already_fixed) is None


def test_fixed_helper_returns_none_when_pattern_unrecognized():
    """When the base pattern has no aliphatic-nitrogen branch to patch, returns None."""
    assert fixed_hbond_acceptor_smarts("[$([O&!$([OX2](C)C=O)]),$([o+0])]") is None


def test_fixed_hbond_acceptor_smarts_from_installed_prolif():
    """From the installed ProLIF: either no override (already fixed / unreadable) or a
    valid patched pattern that carries the exclusion. Both are acceptable, so this
    stays green once ProLIF fixes the issue upstream."""
    pytest.importorskip("prolif")
    Chem = pytest.importorskip("rdkit.Chem")
    patched = fixed_hbond_acceptor_smarts()
    if patched is None:
        return
    assert BACKBONE_N_ACCEPTOR_EXCLUSION in patched
    assert Chem.MolFromSmarts(patched) is not None


# --------------------------------------------------------------------------- #
# 4. ProLIF integration
# --------------------------------------------------------------------------- #

def test_prolif_fingerprint_accepts_acceptor_override():
    """ProLIF accepts our acceptor override for HBDonor/HBAcceptor via ``parameters``
    (the wiring StreaMD relies on)."""
    plf = pytest.importorskip("prolif")
    patched = fixed_hbond_acceptor_smarts()
    if patched is None:
        pytest.skip("installed ProLIF default acceptor could not be patched")
    fp = plf.Fingerprint(
        ["HBDonor", "HBAcceptor"],
        parameters={"HBDonor": {"acceptor": patched}, "HBAcceptor": {"acceptor": patched}},
    )
    registered = set(getattr(fp, "interactions", {}) or {})
    assert {"HBDonor", "HBAcceptor"} <= registered


def test_prolif_fingerprint_removes_false_backbone_hbdonor():
    """End-to-end regression: a ligand O-H donating to a fragmented backbone nitrogen
    is reported as an HBDonor with ProLIF's default SMARTS, and the issue #358 patch
    removes it.

    Uses a reconstructed fragmented residue plus a methanol donor placed against its
    backbone N. If the 3D setup does not reproduce the false contact on a given RDKit
    build the test skips (so it never fails spuriously), but when the contact is
    reproduced the fix must remove it.
    """
    plf = pytest.importorskip("prolif")
    Chem = pytest.importorskip("rdkit.Chem")
    import numpy as np
    from rdkit.Chem import AllChem
    from rdkit.Geometry import Point3D

    res, n_idx = _fragmented_backbone_residue(embed_3d=True)
    n_pos = np.array(res.GetConformer().GetAtomPosition(n_idx))
    for atom in res.GetAtoms():
        info = Chem.AtomPDBResidueInfo()
        info.SetResidueName("ALA")
        info.SetResidueNumber(1)
        info.SetChainId("A")
        info.SetName(f" {atom.GetSymbol():<2} ")
        atom.SetMonomerInfo(info)
    protein = plf.Molecule.from_rdkit(res)

    ligand_rd = Chem.AddHs(Chem.MolFromSmiles("CO"))
    AllChem.EmbedMolecule(ligand_rd, randomSeed=1)
    o_idx = next(a.GetIdx() for a in ligand_rd.GetAtoms() if a.GetAtomicNum() == 8)
    oh_h = next(nb.GetIdx() for nb in ligand_rd.GetAtomWithIdx(o_idx).GetNeighbors()
                if nb.GetAtomicNum() == 1)
    conf = ligand_rd.GetConformer()
    shift = (n_pos + np.array([2.9, 0.0, 0.0])) - np.array(conf.GetAtomPosition(o_idx))
    for i in range(ligand_rd.GetNumAtoms()):
        conf.SetAtomPosition(i, Point3D(*(np.array(conf.GetAtomPosition(i)) + shift)))
    # point the hydroxyl hydrogen toward the backbone nitrogen
    o_now = np.array(conf.GetAtomPosition(o_idx))
    conf.SetAtomPosition(oh_h, Point3D(*(n_pos + (o_now - n_pos) * 0.35)))
    ligand = plf.Molecule.from_rdkit(ligand_rd)

    def contacts(parameters):
        fp = plf.Fingerprint(["HBDonor", "HBAcceptor"], parameters=parameters)
        ifp = fp.generate(ligand, protein, metadata=True)
        found = set()
        for interactions in ifp.values():
            found.update(interactions.keys())
        return found

    patched = inject_backbone_n_acceptor_exclusion(DEFAULT_ACCEPTOR_2_1_0)
    old = contacts(None)
    if "HBDonor" not in old:
        pytest.skip("3D setup did not reproduce the false backbone-N HBDonor on this build")
    new = contacts({"HBDonor": {"acceptor": patched}, "HBAcceptor": {"acceptor": patched}})
    assert "HBDonor" not in new
