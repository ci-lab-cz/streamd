"""Basic sanity check for package import."""

import importlib


def test_import_streamd() -> None:
    """The package should be importable without optional binaries."""
    importlib.import_module("streamd")

