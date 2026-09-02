# Medchem Changelog

This file records user-visible changes. See the [migration guide](docs/migration.md)
for upgrade instructions and [GitHub releases](https://github.com/datamol-io/medchem/releases)
for earlier release notes.

## Next major release (unreleased)

### Highlights

- Refresh Medchem for Python 3.11–3.14 and RDKit 2024.09 or newer.
- Update the optional Lilly MedChem Rules integration from the obsolete 1.0.1
  conda-forge build to the checksum-pinned upstream 2.1 reference release.
- Integrate the community SpacialScore contribution and resolve all previously
  open correctness and performance issues covered by this release.

### Added

- Add normalized and unnormalized RDKit SpacialScore calculations and expose
  SpacialScore through the complexity filter with a custom threshold file.
  The original contribution from
  PR #22 is preserved in the `dev` history.
- Add `medchem install-lilly`, an explicit optional installer that downloads
  the official source archive, verifies its SHA-256 checksum, compiles the
  three reference executables and runs the upstream 35,862-molecule regression
  suite.
- Add regression coverage for Lilly kill rules, suppressed input molecules,
  default thresholds, raw-SMILES valence behavior, batching and subprocess
  pipeline failures.

### Changed

- Require NumPy 1.26+, pandas 2.2+ and Datamol 0.12.5+ alongside the supported
  Python and RDKit versions.
- Update the bundled Lilly queries to upstream 2.1. The previous query set
  disagreed with 528 of the 35,862 upstream oracle records; the new set
  reproduces all five reference outcome files.
- Pass raw SMILES to LillyMol without RDKit pre-validation so the reference
  engine remains authoritative for structures with different valence support.
- Use the upstream 7/25/40 default atom thresholds.
- Stream the complete Lilly pipeline through OS subprocess pipes and forward
  batch parallelism options correctly.
- Keep Lilly entirely optional: importing or installing Medchem does not locate,
  download or compile the native tools. The integration is validated on Linux,
  macOS Apple Silicon and macOS Intel. Native Windows users should use WSL
  because the upstream query reader cannot resolve its manifests there.

### Fixed

- Preserve one output row per Lilly input when a native stage suppresses an
  invalid molecule, marking the missing result as an explicit failure instead
  of raising an indexing error (#35).
- Apply the upstream Lilly kill rules correctly (#36).
- Correct the Toxicophore `michael_acceptors_double` SMARTS so ordinary amines
  and thiols are not classified as Michael acceptors (#32).
- Make complete Lilly runs use a pipeline (#34) and make parallel batch options
  effective (#33).

### Compatibility and delivery

- Use uv for development and CI while retaining `env.yml` as a Conda
  alternative.
- Validate core Medchem on Linux, Windows, macOS Apple Silicon and macOS Intel;
  run Lilly, notebook, documentation, formatting and distribution checks in
  separate jobs.
- Keep publication manual through the `release` action and `PYPI_API_TOKEN`,
  with PEP 740 attestations. Release tests and isolated wheel/source checks
  gate publication; prereleases never replace the stable documentation.
- Add a non-publishing dry run and a [release guide](docs/releasing.md).
  Conda-forge remains a separate channel requiring recipe updates.
