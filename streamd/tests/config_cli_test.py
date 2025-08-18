"""Tests for merging YAML config with CLI arguments.

These tests focus solely on the argument parsing logic used by the
command‑line entry points.  Instead of importing the heavy StreaMD
modules, we build minimal ``argparse`` parsers that mirror the options
relevant to each script and apply the same precedence rules:

* extra keys in the YAML file are ignored;
* values provided on the command line override YAML defaults;
* if a value is missing from the CLI it is taken from the config file;
* CLI‑only options remain available even if no config is supplied.

This keeps the tests lightweight while validating the configuration
mechanism.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import yaml

from streamd.utils.utils import parse_with_config


def _write_config(tmp_path: Path, data: dict) -> Path:
    path = tmp_path / "config.yml"
    path.write_text(yaml.dump(data))
    return path


def _make_parser(options: Iterable[Tuple[str, dict]]) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config")
    for name, kwargs in options:
        parser.add_argument(name, **kwargs)
    return parser


def test_run_md_argument_parsing(tmp_path: Path) -> None:
    """CLI values override YAML and unused keys are ignored for MD parser."""

    parser = _make_parser(
        [
            ("--protein", {}),
            ("--ligand", {}),
            ("--ncpu", {"type": int, "default": 1}),
        ]
    )

    config = {"ligand": "LIG_A", "ncpu": 4, "unused": "x"}
    config_path = _write_config(tmp_path, config)

    cli = ["--config", str(config_path), "--protein", "PRO.pdb", "--ncpu", "8"]
    args, _ = parse_with_config(parser, cli)

    assert args.protein == "PRO.pdb"  # CLI only
    assert args.ligand == "LIG_A"  # from config
    assert args.ncpu == 8  # CLI overrides config
    assert not hasattr(args, "unused")


def test_run_gbsa_argument_parsing(tmp_path: Path) -> None:
    """Ensure GBSA parser merges config and CLI options correctly."""

    parser = _make_parser(
        [
            ("--tpr", {}),
            ("--xtc", {}),
            ("--ligand_id", {}),
            ("--ncpu", {"type": int, "default": 1}),
        ]
    )

    config = {"ligand_id": "ABC", "ncpu": 4, "unused": "x"}
    config_path = _write_config(tmp_path, config)

    cli = [
        "--config",
        str(config_path),
        "--tpr",
        "file.tpr",
        "--xtc",
        "file.xtc",
        "--ncpu",
        "2",
    ]
    args, _ = parse_with_config(parser, cli)

    assert args.tpr == "file.tpr"
    assert args.xtc == "file.xtc"
    assert args.ligand_id == "ABC"  # from config
    assert args.ncpu == 2  # CLI overrides config
    assert not hasattr(args, "unused")


def test_run_prolif_argument_parsing(tmp_path: Path) -> None:
    """Validate ProLIF parser config/CLI precedence."""

    parser = _make_parser(
        [
            ("--tpr", {}),
            ("--xtc", {}),
            ("--ligand", {}),
            ("--ncpu", {"type": int, "default": 1}),
        ]
    )

    config = {"ligand": "LIG1", "ncpu": 4, "unused": "x"}
    config_path = _write_config(tmp_path, config)

    cli = [
        "--config",
        str(config_path),
        "--tpr",
        "file.tpr",
        "--xtc",
        "file.xtc",
        "--ncpu",
        "2",
    ]
    args, _ = parse_with_config(parser, cli)

    assert args.tpr == "file.tpr"
    assert args.xtc == "file.xtc"
    assert args.ligand == "LIG1"  # from config
    assert args.ncpu == 2  # CLI overrides config
    assert not hasattr(args, "unused")

