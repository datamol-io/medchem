# Migrating to Medchem 3.x

Medchem 3.x refreshes the supported scientific Python stack and delivery
process without adding filtering features. Existing rules, thresholds and
return formats remain unchanged unless a dependency required a compatibility
correction.

## Supported runtime

- Python 3.11 through 3.14
- RDKit 2024.09 or newer
- NumPy 1.26 or newer
- pandas 2.2 or newer
- NetworkX 3.2 or newer
- Datamol 0.12.5 or newer

Create a fresh environment rather than upgrading a long-lived 2.x environment
in place:

```bash
python -m venv .venv
. .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install medchem
```

Recent RDKit releases can emit a different canonical SMILES for the same
chemical graph. Code and tests should compare chemical structure and filter
invariants rather than relying on a version-specific serialization.

## Lilly MedChem Rules

The Lilly MedChem Rules remain an optional integration because they require
external command-line tools. Importing Medchem or
`medchem.structural.lilly_demerits` no longer searches for those executables;
they are resolved only when a Lilly filter is executed.

Install the tools from conda-forge before using that integration:

```bash
mamba install -c conda-forge lilly-medchem-rules
```

## Development and releases

Create the complete development environment with:

```bash
mamba env create -n medchem -f env.yml
mamba activate medchem
```

CI tests the supported Python and RDKit series on Linux, macOS and Windows.
The Lilly integration and executable notebooks run in their own conda-backed
job. Documentation, formatting and package distributions are checked
separately. Published GitHub releases are built once and uploaded to PyPI with
short-lived OpenID Connect credentials.
