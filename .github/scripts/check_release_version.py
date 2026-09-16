"""Validate that a version-file change is safe to publish to PyPI."""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import subprocess
import sys
from urllib.error import HTTPError
from urllib.parse import quote
from urllib.request import Request, urlopen

from packaging.version import InvalidVersion, Version


VERSION_PATTERN = re.compile(
    r"^\s*__version__\s*=\s*(['\"])([^'\"]+)\1\s*$",
    re.MULTILINE,
)


def extract_version(source: str) -> str:
    """Return the single PEP 440 version assigned in *source*."""
    matches = VERSION_PATTERN.findall(source)
    if len(matches) != 1:
        raise ValueError("version file must contain exactly one __version__ assignment")

    version = matches[0][1]
    try:
        parsed = Version(version)
    except InvalidVersion as exc:
        raise ValueError(f"invalid PEP 440 version: {version!r}") from exc
    if version != str(parsed):
        raise ValueError(
            f"version must use canonical PEP 440 form: {version!r} normalizes to {str(parsed)!r}"
        )
    return version


def ensure_version_increased(previous: str, current: str) -> None:
    """Reject unchanged versions and version downgrades."""
    if Version(current) <= Version(previous):
        raise ValueError(
            f"version must increase: previous={previous!r}, current={current!r}"
        )


def pypi_version_exists(project: str, version: str) -> bool:
    """Return whether PyPI already contains *version* of *project*."""
    url = (
        f"https://pypi.org/pypi/{quote(project, safe='')}/"
        f"{quote(version, safe='')}/json"
    )
    request = Request(url, headers={"User-Agent": "streamd-release-check/1"})
    try:
        with urlopen(request, timeout=15) as response:
            return response.status == 200
    except HTTPError as exc:
        if exc.code == 404:
            return False
        raise RuntimeError(f"PyPI version check failed with HTTP {exc.code}") from exc


def previous_version(before: str, version_file: str) -> str:
    """Read the version file from the commit preceding the push."""
    result = subprocess.run(
        ["git", "show", f"{before}:{version_file}"],
        check=True,
        capture_output=True,
        text=True,
    )
    return extract_version(result.stdout)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--before", required=True, help="Git SHA before the push")
    parser.add_argument("--project", default="streamd")
    parser.add_argument("--version-file", default="streamd/__init__.py")
    args = parser.parse_args()

    try:
        current = extract_version(Path(args.version_file).read_text(encoding="utf-8"))
        previous = previous_version(args.before, args.version_file)
        ensure_version_increased(previous, current)
        if pypi_version_exists(args.project, current):
            raise ValueError(
                f"{args.project} {current} already exists on PyPI; versions cannot be replaced"
            )
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as exc:
        print(f"release version validation failed: {exc}", file=sys.stderr)
        return 1

    print(current)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
