"""Tests for the release version guard used by GitHub Actions."""

import importlib.util
from pathlib import Path
from urllib.error import HTTPError

import pytest


SCRIPT = Path(__file__).parents[2] / ".github" / "scripts" / "check_release_version.py"
SPEC = importlib.util.spec_from_file_location("check_release_version", SCRIPT)
release_version = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(release_version)


@pytest.mark.parametrize(
    "source, expected",
    [
        ("__version__ = '0.6.1'", "0.6.1"),
        ('__version__ = "1.0rc1"\n', "1.0rc1"),
    ],
)
def test_extract_version(source, expected):
    assert release_version.extract_version(source) == expected


@pytest.mark.parametrize(
    "source",
    [
        "",
        "__version__ = 'not a version'",
        "__version__ = 'v1.0'",
        "__version__ = '1.0'\n__version__ = '1.1'",
    ],
)
def test_extract_version_rejects_invalid_source(source):
    with pytest.raises(ValueError):
        release_version.extract_version(source)


def test_version_must_increase():
    release_version.ensure_version_increased("0.6", "0.6.1")

    with pytest.raises(ValueError, match="must increase"):
        release_version.ensure_version_increased("0.6", "0.6")
    with pytest.raises(ValueError, match="must increase"):
        release_version.ensure_version_increased("0.6", "0.5.1")


class _Response:
    status = 200

    def __enter__(self):
        return self

    def __exit__(self, *_args):
        return None


def test_existing_pypi_version_is_detected(monkeypatch):
    monkeypatch.setattr(release_version, "urlopen", lambda *a, **k: _Response())
    assert release_version.pypi_version_exists("streamd", "0.6") is True


def test_missing_pypi_version_is_available(monkeypatch):
    def not_found(request, timeout):
        raise HTTPError(request.full_url, 404, "Not Found", {}, None)

    monkeypatch.setattr(release_version, "urlopen", not_found)
    assert release_version.pypi_version_exists("streamd", "999.0") is False


def test_pypi_server_error_fails_closed(monkeypatch):
    def server_error(request, timeout):
        raise HTTPError(request.full_url, 503, "Unavailable", {}, None)

    monkeypatch.setattr(release_version, "urlopen", server_error)
    with pytest.raises(RuntimeError, match="HTTP 503"):
        release_version.pypi_version_exists("streamd", "0.6.1")
