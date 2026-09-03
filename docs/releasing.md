# Releasing medchem

Publication is manual and driven by the **release** GitHub Action. Merging,
pushing a tag or drafting a GitHub Release does not, on its own, upload anything
to PyPI.

## Before publishing

1. Finalize `CHANGELOG.md`: give the release an exact `## <version> - YYYY-MM-DD`
   heading (the workflow rejects any other form) and refresh README or docs
   wording that still calls the release unreleased.
2. Merge the release into `main`, preserving contributor history.
3. Confirm the project is registered on PyPI as a Trusted Publisher for this
   repository, the `release.yml` workflow and the `pypi` environment. No API
   token is used; publication authenticates over OpenID Connect.

## Run the release action

Open **Actions → release → Run workflow** from `main` and enter the version in
**release-version** — canonical, without a `v` prefix (`2.1.0`, or `2.1.0rc1`
where `a`, `b` and `rc` mark prereleases). Start with **dry-run** enabled. A dry
run may also be launched from `dev`, accepts unfinished release notes and tags
only the throwaway runner checkout, so it never touches GitHub or PyPI.

Rerun from `main` with **dry-run** disabled to publish. The action validates the
version and notes, reruns the full test and quality suites on the target commit,
builds the wheel and sdist, installs each in isolation to verify its version and
import path, and builds the documentation before anything is uploaded.

Upload uses PyPI Trusted Publishing over OpenID Connect — no API token — and the
PEP 740 attestations are accepted as part of that flow (PyPI rejects
attestations on token-based uploads). Only after PyPI accepts the upload does the
action tag the commit, create the GitHub Release and deploy versioned
documentation. A prerelease never moves the `stable` docs alias, and a failed
upload leaves it untouched — rerun the failed jobs rather than rebuilding an
already published version.

## Conda-forge

The [feedstock](https://github.com/conda-forge/medchem-feedstock) is a separate
channel. Its bot proposes a version bump after PyPI publication, but maintainers
must review the dependencies and run the recipe tests there; the release action
never publishes to conda-forge.

Keep Lilly optional in the recipe: the native tools install separately through
`medchem install-lilly`, so do not add a compilation step to the package
install. A version-only feedstock bump is not enough — its dependency metadata
and import/CLI tests must reflect the new installation boundaries.
