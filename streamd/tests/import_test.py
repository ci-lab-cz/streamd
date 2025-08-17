"""Basic test on import of streamd."""
import subprocess
import sys


def test_imports():
    """Ensure package imports and dependencies resolve."""
    assert "streamd" in sys.modules
    assert subprocess.check_call('gmx') == 0