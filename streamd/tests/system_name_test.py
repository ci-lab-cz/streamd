import numpy as np
import pandas as pd
import pytest

from streamd.analysis.md_system_analysis import _parse_system_name, _system_metadata
from streamd.analysis.run_analysis import _ensure_replica_cols


def test_parse_system_name_replica():
    """Parse explicit replica suffixes from system names."""
    protein, replica = _parse_system_name("protein_HIS_replica7")
    assert protein == "protein_HIS"
    assert replica == 7


@pytest.mark.parametrize(
    'system_name, ligand_name, expected',
    [
        ('protein_cmp17_replica2', 'cmp17', ('protein_replica2', 'protein', 2)),
        ('protein_cmp17', 'cmp17', ('protein', 'protein', 1)),
        ('protein_without_ligand_name', 'cmp17', ('protein_without_ligand_name', 'protein_without_ligand_name', 1)),
        ('protein_cmp17_domain_cmp17_replica2', 'cmp17', ('protein_cmp17_domain_replica2', 'protein_cmp17_domain', 2)),
        ('protein_only_replica3', None, ('protein_only_replica3', 'protein_only', 3)),
    ],
)
def test_system_metadata_removes_only_trailing_ligand_token(system_name, ligand_name, expected):
    """Remove only the final ligand token when deriving protein metadata."""
    assert _system_metadata(system_name, ligand_name) == expected


def test_ensure_replica_cols_infers_metadata_per_row():
    """Infer protein and replica metadata across proteins and replicas."""
    df = _ensure_replica_cols(pd.DataFrame({
        'system': [
            'protein_a_replica1',
            'protein_a_replica2',
            'protein_b_replica3',
        ],
    }))

    assert df['protein_name'].tolist() == ['protein_a', 'protein_a', 'protein_b']
    assert df['replica'].tolist() == [1, 2, 3]


def test_ensure_replica_cols_preserves_and_fills_metadata():
    """Fill missing metadata per row without replacing existing values."""
    df = _ensure_replica_cols(pd.DataFrame({
        'system': ['a_replica1', 'b_replica2', 'c_replica3'],
        'protein_name': ['custom_a', np.nan, 'custom_c'],
        'replica': [9, 10, np.nan],
    }))

    assert df['protein_name'].tolist() == ['custom_a', 'b', 'custom_c']
    assert df['replica'].tolist() == [9, 10, 3]


def test_ensure_replica_cols_requires_system():
    """Reject RMSD metadata repair when the system column is absent."""
    with pytest.raises(ValueError, match="missing the required 'system' column"):
        _ensure_replica_cols(pd.DataFrame({'backbone': [0.1]}))
