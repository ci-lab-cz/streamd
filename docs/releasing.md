# Releasing StreaMD to PyPI

PyPI publication is triggered only when a push to `master` changes
`streamd/__init__.py`. Set `__version__` to a new, higher PEP 440 version and
commit the change. The release workflow then validates the version, runs the
lightweight tests, builds and checks the source and wheel distributions, and
publishes them to PyPI. After PyPI succeeds, it creates a `v<version>` Git tag
and GitHub Release at the same commit, generates release notes, and attaches the
source and wheel distributions.

Before the first automated release, a StreaMD PyPI maintainer must add a GitHub
Trusted Publisher to the existing `streamd` PyPI project with these values:

- Owner: `ci-lab-cz`
- Repository: `streamd`
- Workflow: `release.yml`
- Environment: `pypi`

The workflow fails without publishing if the version did not increase, already
exists on PyPI, the tests fail, or either distribution is invalid. PyPI versions
are immutable, so a failed release must be corrected with a new version when
any artifact reached PyPI. If PyPI succeeds but the final GitHub Release job
fails, use GitHub Actions' **Re-run failed jobs** operation; do not rerun the
whole workflow, because the version now exists on PyPI.
