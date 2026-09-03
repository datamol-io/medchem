# Upgrading to Medchem 2.1.0

Medchem 2.1.0 refreshes the supported scientific Python stack and delivery
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

The Lilly MedChem Rules remain an optional integration because their exact
query language and demerit engine are implemented by the native upstream
tools. Medchem now tracks upstream 2.1.0, including its corrected and added
queries. The previous vendored rule set differed from the 2.1 reference result
on 528 of the 35,862 oracle records.

Install the checksum-pinned compatible tools beside the active Python:

```bash
medchem install-lilly
```

The command downloads the official source archive, verifies its SHA-256,
builds the three required executables and runs the upstream regression suite.
It is an explicit optional operation, separate from `pip install medchem`, and
does not create Python bindings. The upstream test driver also requires Ruby.
Native Windows is not supported because the upstream 2.1 query reader cannot
resolve the bundled manifests; use WSL for the Lilly integration on Windows.
The conda-forge build remains at 1.0.1 and cannot reproduce the 2.1 rule set,
so it was removed from the development environment. A fresh Python or Rust
rewrite was deliberately not substituted for the reference chemistry because
it would not preserve the richer LillyMol query semantics.

Importing Medchem still performs no executable lookup or installation. The
tools are resolved lazily, including beside an unactivated Python interpreter.
Raw SMILES are passed to LillyMol without an RDKit pre-validation step, so
structures supported by LillyMol's valence model are no longer rejected early.

The Python wrapper now matches upstream's default 7/25/40 atom thresholds,
retains rows suppressed by LillyMol as explicit failures, and connects the
native stages with portable OS pipes. This resolves the open Lilly correctness,
missing-row, parallel-batching, and avoidable-I/O issues. The vendored 2.1
queries reproduce all five upstream outcome files on the 35,862-molecule
reference corpus.

The open SpacialScore pull request was ported to the current RDKit API and
completed with public export, validation, tests, and documentation. The
Toxicophore Michael-acceptor SMARTS now require the imine or thiocarbonyl
double bond rather than treating ordinary amines and thiols as acceptors.

For uv-managed projects, use:

```bash
uv add medchem
uv run medchem install-lilly
```

A direct `medchem[lilly]` extra needs a separately published native-wheel
distribution for every supported platform. Extras cannot themselves run a
post-install compiler, so this release keeps the honest two-step installation
instead of exposing an extra that would not install Lilly.

## Development and releases

Create the complete development environment with:

```bash
uv sync --all-extras
```

`env.yml` remains a supported Conda alternative. CI uses uv and tests the
supported Python and RDKit series on Linux, Windows, macOS Apple
Silicon, and macOS Intel. The Lilly integration is built and checked on Linux
and both macOS architectures; executable notebooks run in a separate Linux job.
Documentation, formatting and package distributions are checked separately.
Publication remains a manual action using the `PYPI_API_TOKEN` secret.
The release action reruns tests, validates both installed distributions and
builds documentation before uploading. See the [release guide](releasing.md)
for dry runs, prereleases and the separate conda-forge recipe updates.
