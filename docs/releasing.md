# Releasing StreaMD to PyPI

PyPI publication is triggered only when a push to `master` changes
`streamd/__init__.py`. Set `__version__` to a new, higher PEP 440 version and
commit the change. The release workflow then validates the version, runs the
lightweight tests, builds and checks the source and wheel distributions, and
publishes them to PyPI.

Before the first automated release, a StreaMD PyPI maintainer must add a GitHub
Trusted Publisher to the existing `streamd` PyPI project with these values:

- Owner: `ci-lab-cz`
- Repository: `streamd`
- Workflow: `release.yml`
- Environment: `pypi`

The workflow fails without publishing if the version did not increase, already
exists on PyPI, the tests fail, or either distribution is invalid. PyPI versions
are immutable, so a failed release must be corrected with a new version when
any artifact reached PyPI.
