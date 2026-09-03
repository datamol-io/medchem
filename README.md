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
[2.1.0 upgrade guide](docs/migration.md).

## Installation

Install from PyPI or conda-forge:

```bash
# uv (recommended)
uv add medchem

# pip
pip install medchem

# conda-forge
micromamba install -c conda-forge medchem
```

Medchem 2.1.0 supports Python 3.11 through 3.14 and RDKit 2024.09 or newer. See
the [upgrade guide](docs/migration.md) for details.

### Optional Lilly MedChem Rules

`LillyDemeritsFilters` uses the upstream Lilly command-line tools. Install them
once after Medchem:

```bash
# pip or conda-forge environment
medchem install-lilly

# uv-managed project
uv run medchem install-lilly
```

The installer compiles the pinned upstream 2.1 tools from source, so a C++
toolchain is required: Linux needs a C++ compiler, GNU Make, zlib, and Ruby;
macOS needs the Xcode command-line tools and Ruby. On Windows, run it through
WSL.

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
