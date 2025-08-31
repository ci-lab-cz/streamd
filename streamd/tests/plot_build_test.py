"""Tests for plotting helpers."""

import pytest
import tempfile
import pandas as pd
from pathlib import Path

from streamd.analysis.plot_build import plot_rmsd, plot_rmsd_mean_std

analysis_test = pytest.mark.skipif(
    "not config.getoption('--run-analysis')",
    reason="Only run when --run-analysis is given",
)


@pytest.fixture
def rmsd_df():
    """Provide a small RMSD DataFrame."""
    return pd.DataFrame({
        'time(ns)': [0, 1, 2, 3, 4, 5],
        'replica1': [2.2, 2.1, 2.2, 2.0, 1.8, 1.75],
        'replica2': [0.55, 0.58, 0.65, 0.7, 0.68, 0.72]
    })


@pytest.fixture
def rmsd_summary_df():
    """Provide a summary RMSD DataFrame."""
    return pd.DataFrame({
        'RMSD_mean': [1.0, 3.0, 4.5, 2.5],
        'RMSD_std': [0.2, 0.4, 0.6, 0.3],
        'rmsd_system': ['A', 'A', 'B', 'B'],
        'time_range': ['0-10', '10-20', '20-30', '30-40'],
        'system': ['X', 'X', 'Y', 'Y']
    })


@analysis_test
def test_plot_rmsd_creates_file(rmsd_df):
    """Ensure RMSD plot is written to disk."""
    with tempfile.TemporaryDirectory() as tmpdir:
        output_path = Path(tmpdir) / "rmsd_plot.png"
        plot_rmsd(rmsd_df, "TestSystem", output_path)
        assert output_path.exists() # Check if output exists

@analysis_test
def test_plot_rmsd_mean_std_creates_file(rmsd_summary_df):
    """Ensure RMSD mean/std plot is written to disk."""
    with tempfile.TemporaryDirectory() as tmpdir:
        output_path = Path(tmpdir) / "rmsd_mean_std.html"
        plot_rmsd_mean_std(
            data=rmsd_summary_df,
            paint_by_col='system',
            show_legend=True,
            out_name=str(output_path),
            title="RMSD Plot"
        )
        assert output_path.exists() # check if output exists
