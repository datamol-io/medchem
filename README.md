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
[![license](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://github.com/datamol-io/medchem/blob/main/LICENSE.md)
[![GitHub Repo stars](https://img.shields.io/github/stars/datamol-io/medchem)](https://github.com/datamol-io/medchem/stargazers)
[![GitHub Repo stars](https://img.shields.io/github/forks/datamol-io/medchem)](https://github.com/datamol-io/medchem/network/members)
[![test](https://github.com/datamol-io/medchem/actions/workflows/test.yml/badge.svg)](https://github.com/datamol-io/medchem/actions/workflows/test.yml)
[![release](https://github.com/datamol-io/medchem/actions/workflows/release.yml/badge.svg)](https://github.com/datamol-io/medchem/actions/workflows/release.yml)
[![code-check](https://github.com/datamol-io/medchem/actions/workflows/code-check.yml/badge.svg)](https://github.com/datamol-io/medchem/actions/workflows/code-check.yml)
[![doc](https://github.com/datamol-io/medchem/actions/workflows/doc.yml/badge.svg)](https://github.com/datamol-io/medchem/actions/workflows/doc.yml)

Medchem is a Python library that proposes multiple molecular medchem filters to a wide range of use cases relevant in a drug discovery context.

## Installation

```bash
micromamba install -c conda-forge medchem

# or using pip
pip install medchem
```

## Documentation

Visit <https://medchem-docs.datamol.io/>.

## Development lifecycle

### Setup dev environment

```bash
micromamba create -n medchem -f env.yml
micromamba activate medchem

pip install --no-deps -e .
```

### Tests

You can run tests locally with:

```bash
pytest
```

## License

Under the Apache-2.0 license. See [LICENSE.md](LICENSE.md).
