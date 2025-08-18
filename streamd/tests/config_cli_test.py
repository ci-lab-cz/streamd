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
from typing import Iterable, Tuple

import yaml
import pytest

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


@pytest.mark.parametrize(
    "options, config, cli, expected",
    [
        (
            [
                ("--protein", {}),
                ("--ligand", {}),
                ("--ncpu", {"type": int, "default": 1}),
            ],
            {"ligand": "LIG_A", "ncpu": 4, "unused": "x"},
            ["--protein", "PRO.pdb", "--ncpu", "8"],
            {"protein": "PRO.pdb", "ligand": "LIG_A", "ncpu": 8},
        ),
        (
            [
                ("--tpr", {}),
                ("--xtc", {}),
                ("--ligand_id", {}),
                ("--ncpu", {"type": int, "default": 1}),
                ("--debug", {"action": "store_true", "default": False}),
            ],
            {"ligand_id": "ABC", "ncpu": 4, "unused": "x", "debug": False},
            [
                "--tpr",
                "file.tpr",
                "--xtc",
                "file.xtc",
                "--ncpu",
                "2",
                "--debug",
            ],
            {
                "tpr": "file.tpr",
                "xtc": "file.xtc",
                "ligand_id": "ABC",
                "ncpu": 2,
                "debug": True,
            },
        ),
        (
            [
                ("--tpr", {}),
                ("--xtc", {}),
                ("--ligand", {}),
                ("--ncpu", {"type": int, "default": 1}),
            ],
            {"ligand": "LIG1", "ncpu": 4, "unused": "x"},
            ["--tpr", "file.tpr", "--xtc", "file.xtc", "--ncpu", "2"],
            {
                "tpr": "file.tpr",
                "xtc": "file.xtc",
                "ligand": "LIG1",
                "ncpu": 2,
            },
        ),
    ],
    ids=["md", "gbsa", "prolif"],
)
def test_argument_parsing(
    tmp_path: Path,
    options: Iterable[Tuple[str, dict]],
    config: dict,
    cli: list[str],
    expected: dict,
) -> None:
    """Common test ensuring YAML/CLI merging works for all entry points."""

    parser = _make_parser(options)
    config_path = _write_config(tmp_path, config)
    args, _ = parse_with_config(parser, ["--config", str(config_path), *cli])

    for key, value in expected.items():
        assert getattr(args, key) == value
    assert not hasattr(args, "unused")


def test_list_argument_from_config(tmp_path: Path) -> None:
    """Lists provided via YAML are split and converted by ``nargs``/``type``."""

    parser = _make_parser([
        ("--steps", {"nargs": "*", "type": int, "default": []}),
    ])

    # Provided as a single whitespace-separated string in YAML
    config = {"steps": "1 2"}
    config_path = _write_config(tmp_path, config)

    args, _ = parse_with_config(parser, ["--config", str(config_path)])
    assert args.steps == [1, 2]

    # CLI values still override config-provided defaults
    args, _ = parse_with_config(
        parser, ["--config", str(config_path), "--steps", "3", "4"]
    )
    assert args.steps == [3, 4]

