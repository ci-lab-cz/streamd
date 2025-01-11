"""Basic test on import of streamd"""
import subprocess
import sys

def test_imports():
    assert "streamd" in sys.modules
    assert subprocess.check_call('gmx') == 0