"""Tests for the ProLIF issue #358 backbone-nitrogen H-bond acceptor fix.

See https://github.com/chemosim-lab/ProLIF/issues/358: ProLIF matches its H-bond
SMARTS on the protein split into single-residue fragments, which severs the
backbone nitrogen from the previous residue's carbonyl and makes the amide N wrongly
match as an H-bond acceptor. The fix ANDs an extra exclusion into the aliphatic
nitrogen branch of the acceptor pattern.
"""
import pytest

from streamd.prolif.run_prolif import (
    BACKBONE_N_ACCEPTOR_EXCLUSION,
    inject_backbone_n_acceptor_exclusion,
    fixed_hbond_acceptor_smarts,
)

# The stock ProLIF 2.1.0 H-bond acceptor SMARTS (used as a representative input so
# the string-level tests do not depend on the installed ProLIF version).
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


def test_exclusion_is_injected_into_nitrogen_branch():
    """The exclusion clause is added inside the aliphatic-nitrogen branch."""
    patched = inject_backbone_n_acceptor_exclusion(DEFAULT_ACCEPTOR_2_1_0)
    assert patched != DEFAULT_ACCEPTOR_2_1_0
    assert BACKBONE_N_ACCEPTOR_EXCLUSION in patched
    n_branch = next(b for b in _split_top_level_branches(patched) if b.startswith('$([N'))
    assert BACKBONE_N_ACCEPTOR_EXCLUSION in n_branch


def test_only_nitrogen_branch_is_modified():
    """Oxygen/fluorine/aromatic branches are untouched and the count is preserved."""
    before = _split_top_level_branches(DEFAULT_ACCEPTOR_2_1_0)
    after = _split_top_level_branches(inject_backbone_n_acceptor_exclusion(DEFAULT_ACCEPTOR_2_1_0))
    assert len(before) == len(after)
    for b_old, b_new in zip(before, after):
        if b_old.startswith('$([N'):
            assert BACKBONE_N_ACCEPTOR_EXCLUSION in b_new
        else:
            assert b_new == b_old


def test_exclusion_placed_before_atom_closing_bracket():
    """The clause is inserted inside the nitrogen atom query, not appended outside it."""
    patched = inject_backbone_n_acceptor_exclusion(DEFAULT_ACCEPTOR_2_1_0)
    # the clause must appear immediately before the closing ']' of the N atom query
    assert (BACKBONE_N_ACCEPTOR_EXCLUSION + "])") in patched


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


def test_returns_unchanged_for_empty_or_malformed_input():
    """Empty or non-bracketed input is returned unchanged (no crash)."""
    assert inject_backbone_n_acceptor_exclusion("") == ""
    assert inject_backbone_n_acceptor_exclusion(None) is None
    assert inject_backbone_n_acceptor_exclusion("N") == "N"


def test_patched_smarts_is_valid_rdkit_pattern():
    """The patched SMARTS parses as a valid RDKit query."""
    Chem = pytest.importorskip("rdkit.Chem")
    patched = inject_backbone_n_acceptor_exclusion(DEFAULT_ACCEPTOR_2_1_0)
    assert Chem.MolFromSmarts(patched) is not None


def test_fixed_hbond_acceptor_smarts_from_installed_prolif():
    """The runtime helper derives a patched, valid pattern from the installed ProLIF."""
    pytest.importorskip("prolif")
    Chem = pytest.importorskip("rdkit.Chem")
    patched = fixed_hbond_acceptor_smarts()
    assert patched is not None
    assert BACKBONE_N_ACCEPTOR_EXCLUSION in patched
    assert Chem.MolFromSmarts(patched) is not None
