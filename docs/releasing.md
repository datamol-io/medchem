# Releasing medchem

Publication is manual. Merging into `main`, pushing a tag or publishing a
GitHub Release does not upload a package to PyPI.

## Before publishing

1. Finalize the release notes in `CHANGELOG.md` with an exact version heading,
   for example `## 2.1.0 - YYYY-MM-DD`. Date the release and update README
   and website wording that still calls it unreleased or specific to `dev`.
2. Merge the release changes into `main`, preserving contributor history.
3. Confirm that the repository secret `PYPI_API_TOKEN` contains a valid PyPI
   token authorized for `medchem`. The workflow checks that the secret is
   present but cannot validate its scope without contacting PyPI.
4. Delete any stale tag that already claims the target version. A `2.1.0` tag
   was pushed to GitHub during earlier development but never published to PyPI;
   the publish job creates the release tag itself and fails if the name already
   exists. Remove the old tag first with `git push origin :refs/tags/2.1.0`
   (and `git tag -d 2.1.0` locally) so the action can tag the tested commit.

## Run the release action

Open **Actions → release → Run workflow**. Select `main`, enter
`2.1.0` in **release-version**, and leave **dry-run** checked for a
rehearsal. Dry runs may also run from `dev` and accept unfinished release
notes. They create a tag only inside the temporary runner checkout, never on
GitHub, and cannot publish packages or documentation.

For publication, launch the action again from `main` with **dry-run**
unchecked. Use canonical versions without a `v` prefix, such as `2.1.0`
or `2.1.0rc1`; `a`, `b` and `rc` suffixes identify prereleases.

The action validates the version and release notes, reruns the complete test
and quality workflows on the selected commit, and builds the distributions.
Both wheel and source installations are checked with Python's isolated mode,
including their version and import location. Documentation must also build
successfully before anything is uploaded.

The publish job sends the artifacts and PEP 740 attestations to PyPI using
`PYPI_API_TOKEN`. OpenID Connect is used to sign the attestations, not to
authenticate the upload; no PyPI Trusted Publisher registration is required.

Only after PyPI succeeds does the action create the GitHub tag and Release
at the tested commit, then deploy versioned documentation. A prerelease never
moves the `stable` documentation alias. A failed upload leaves that alias
unchanged. If a later step fails, rerun the failed jobs rather than rebuilding
an already published version.

## Conda-forge

The [feedstock](https://github.com/conda-forge/medchem-feedstock) is a separate
release channel. Its update bot proposes version changes after PyPI publication,
but maintainers must review dependencies and run the recipe tests there.
The package release action does not publish to conda-forge.

Raise Python to 3.11+, update the direct dependencies and keep Lilly optional.
The native Lilly tools are installed separately with `medchem install-lilly`;
do not add a compilation step to the Python package installation.

Do not merge a version-only feedstock update for this release. Its
dependency metadata and import/CLI tests must reflect the new installation
boundaries.
