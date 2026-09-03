<div align="center">
    <img src="docs/images/logo.png" height="200px">
    <h3>Medchem - Molecular filtering for drug discovery</h3>
</div>

---

[![PyPI](https://img.shields.io/pypi/v/medchem)](https://pypi.org/project/medchem/)
[![Conda](https://img.shields.io/conda/v/conda-forge/medchem?label=conda&color=success)](https://anaconda.org/conda-forge/medchem)
[![PyPI - Downloads](https://img.shields.io/pypi/dm/medchem)](https://pypi.org/project/medchem/)
[![Conda](https://img.shields.io/conda/dn/conda-forge/medchem)](https://anaconda.org/conda-forge/medchem)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/medchem)](https://pypi.org/project/medchem/)
[![license](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](LICENSE.md)
[![GitHub Repo stars](https://img.shields.io/github/stars/datamol-io/medchem)](https://github.com/datamol-io/medchem/stargazers)
[![GitHub Repo stars](https://img.shields.io/github/forks/datamol-io/medchem)](https://github.com/datamol-io/medchem/network/members)
[![test](https://github.com/datamol-io/medchem/actions/workflows/test.yml/badge.svg)](https://github.com/datamol-io/medchem/actions/workflows/test.yml)
[![release](https://github.com/datamol-io/medchem/actions/workflows/release.yml/badge.svg)](https://github.com/datamol-io/medchem/actions/workflows/release.yml)
[![code-check](https://github.com/datamol-io/medchem/actions/workflows/code-check.yml/badge.svg)](https://github.com/datamol-io/medchem/actions/workflows/code-check.yml)
[![doc](https://github.com/datamol-io/medchem/actions/workflows/doc.yml/badge.svg)](https://github.com/datamol-io/medchem/actions/workflows/doc.yml)

Medchem is a Python library that proposes multiple molecular medchem filters to a wide range of use cases relevant in a drug discovery context.

## Updates

Medchem 2.1.0 updates the supported Python and RDKit stack, adds RDKit
SpacialScore support, corrects the Toxicophore Michael-acceptor rule, and
refreshes the optional Lilly MedChem Rules integration against upstream 2.1.
The Lilly wrapper now preserves every input row, uses the reference thresholds
and query set, supports parallel batches, and streams its native stages.

See the [complete changelog](CHANGELOG.md) and the
[2.1.0 upgrade guide](docs/migration.md). These notes describe the upcoming
release; PyPI and conda-forge still provide the published stable versions.

Release maintainers: see the [manual release guide](docs/releasing.md).

## Installation

### Core install (without Lilly)

Every Medchem filter except the optional Lilly MedChem Rules works with a plain
install from PyPI or conda-forge — no extra steps and no native build:

```bash
# uv (recommended)
uv add medchem

# pip
pip install medchem

# conda-forge
micromamba install -c conda-forge medchem
```

Medchem 2.1.0 supports Python 3.11 through 3.14 and RDKit 2024.09 or newer. See
the [migration guide](docs/migration.md) for dependency details.

### With the optional Lilly MedChem Rules

`LillyDemeritsFilters` needs the upstream Lilly command-line tools, which are
**not** bundled with the package. Installing Medchem — through pip, conda-forge
or uv — never downloads or compiles Lilly. After installing Medchem, run the
explicit installer once; it behaves the same regardless of how Medchem was
installed, because it compiles the checksum-pinned upstream 2.1 tools from source
beside the active Python:

```bash
# pip or conda-forge environment
medchem install-lilly

# uv-managed project
uv run medchem install-lilly
```

The installer verifies the source-archive checksum, builds three upstream
executables (not a Python extension), and runs Lilly's complete 35,862-molecule
regression suite. Toolchain requirements: Linux needs a C++ compiler, GNU Make,
zlib, and Ruby (for the upstream test driver); macOS needs the Xcode
command-line tools and Ruby. Native Windows is not supported by upstream 2.1
because its query reader cannot resolve the bundled manifests — Windows users
should run it through WSL.

There is no `medchem[lilly]` extra: a Python extra can only pull other
distributions, not run a post-install compiler, and no platform-specific Lilly
wheels are published yet. Do not use the conda-forge `lilly-medchem-rules`
package either — it is still version 1.0.1 and is incompatible with the vendored
2.1 rule set.

## Documentation

Visit <https://medchem-docs.datamol.io/>.

## Development lifecycle

### Setup dev environment

```bash
uv sync --all-extras
uv run medchem install-lilly
```

`env.yml` remains available when a Conda development environment is required.

### Tests

You can run tests locally with:

```bash
uv run python -m pytest -m "not integration"
uv run python -m pytest -m integration --no-cov -n 0
```

The first command is the fast core suite. The second runs the available Lilly
2.1 checks and executable tutorials. GitHub Actions validates the Python core
on Linux, Windows, macOS Apple Silicon, and macOS Intel. The pinned Lilly
release is built and tested separately on Linux and both macOS architectures.

## License

Under the Apache-2.0 license. See [LICENSE.md](LICENSE.md).
Bundled Lilly query data retain their [upstream attribution](medchem/data/LILLY_NOTICE).
The optional native tools are a separate upstream distribution.

## Citation

[![DOI](https://zenodo.org/badge/653814138.svg)](https://doi.org/10.5281/zenodo.14588937)
