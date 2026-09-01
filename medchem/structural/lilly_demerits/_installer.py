"""Reproducible installer for the optional Lilly MedChem Rules executables."""

from __future__ import annotations

import hashlib
import os
import shutil
import subprocess
import sys
import tarfile
import tempfile
import urllib.request
from pathlib import Path

LILLY_VERSION = "2.1.0"
LILLY_ARCHIVE_URL = (
    "https://github.com/IanAWatson/Lilly-Medchem-Rules/" f"archive/refs/tags/v{LILLY_VERSION}.tar.gz"
)
LILLY_ARCHIVE_SHA256 = "57c2c3889bfe412e97d4e35c0a93e07b0870327476c1d99d73a20a7a8cb5788e"
LILLY_BINARIES = ("mc_first_pass", "tsubstructure", "iwdemerit")
LILLY_VERSION_MARKER = ".lilly-medchem-rules-version"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _extract_archive(archive_path: Path, destination: Path) -> Path:
    destination = destination.resolve()
    destination.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive_path, mode="r:gz") as archive:
        members = archive.getmembers()
        for member in members:
            target = (destination / member.name).resolve()
            if not target.is_relative_to(destination):
                raise RuntimeError(f"Unsafe path in Lilly archive: {member.name}")
            if member.issym() or member.islnk() or member.isdev():
                raise RuntimeError(f"Unsupported entry in Lilly archive: {member.name}")
        archive.extractall(destination, members=members)

    roots = [path for path in destination.iterdir() if path.is_dir()]
    if len(roots) != 1:
        raise RuntimeError("The Lilly archive does not contain one source directory")
    return roots[0]


def install_lilly_rules(
    prefix: Path | str | None = None,
    *,
    jobs: int | None = None,
    run_tests: bool = True,
    force: bool = False,
) -> dict[str, str]:
    """Build and install the pinned Lilly MedChem Rules release.

    The source archive is checksum-verified, compiled with its upstream
    Makefile, and copied beside the active Python interpreter. This keeps the
    optional integration out of Medchem's runtime dependencies while avoiding
    the obsolete conda-forge 1.0.1 build.

    Args:
        prefix: Installation prefix. Defaults to the active Python prefix.
        jobs: Parallel compiler jobs. Defaults to at most eight CPU cores.
        run_tests: Run the upstream 35,862-molecule regression suite.
        force: Rebuild even when the pinned version is already installed.

    Returns:
        A mapping from executable names to installed paths.
    """
    if sys.platform == "win32" and not os.environ.get("MSYSTEM"):
        raise RuntimeError(
            "Building Lilly MedChem Rules on Windows requires an MSYS2 MSYS shell "
            "with the gcc, make, and zlib-devel packages installed"
        )

    install_prefix = Path(prefix).expanduser().resolve() if prefix else Path(sys.prefix).resolve()
    executable_suffix = ".exe" if sys.platform == "win32" else ""
    bin_dir = install_prefix / ("Scripts" if sys.platform == "win32" else "bin")
    marker = bin_dir / LILLY_VERSION_MARKER
    installed = {name: str(bin_dir / f"{name}{executable_suffix}") for name in LILLY_BINARIES}

    if (
        not force
        and marker.is_file()
        and marker.read_text(encoding="utf-8").strip() == LILLY_VERSION
        and all(Path(path).is_file() for path in installed.values())
    ):
        return installed

    make = shutil.which("make")
    if make is None:
        platform_hint = " (use an MSYS2 MSYS shell on Windows)" if sys.platform == "win32" else ""
        raise RuntimeError(
            "Building Lilly MedChem Rules requires `make`, a C++ compiler, and zlib" + platform_hint
        )
    if run_tests and shutil.which("ruby") is None:
        raise RuntimeError(
            "Validating Lilly MedChem Rules requires Ruby; install Ruby or use "
            "`medchem install-lilly --no-test` to skip the upstream regression suite"
        )

    compiler_jobs = jobs if jobs is not None else min(os.cpu_count() or 1, 8)
    if compiler_jobs < 1:
        raise ValueError("jobs must be at least 1")

    with tempfile.TemporaryDirectory(prefix="medchem-lilly-") as temporary:
        temporary_dir = Path(temporary)
        archive_path = temporary_dir / f"lilly-medchem-rules-{LILLY_VERSION}.tar.gz"
        request = urllib.request.Request(
            LILLY_ARCHIVE_URL,
            headers={"User-Agent": f"medchem-lilly-installer/{LILLY_VERSION}"},
        )
        with urllib.request.urlopen(request, timeout=120) as response, archive_path.open("wb") as out:
            shutil.copyfileobj(response, out)

        checksum = _sha256(archive_path)
        if checksum != LILLY_ARCHIVE_SHA256:
            raise RuntimeError(
                "Lilly source checksum mismatch: " f"expected {LILLY_ARCHIVE_SHA256}, received {checksum}"
            )

        source_dir = _extract_archive(archive_path, temporary_dir / "source")
        subprocess.run([make, f"-j{compiler_jobs}"], cwd=source_dir, check=True)
        if run_tests:
            subprocess.run([make, "test"], cwd=source_dir, check=True)

        bin_dir.mkdir(parents=True, exist_ok=True)
        for name in LILLY_BINARIES:
            candidates = [
                source_dir / "bin" / name,
                source_dir / "bin" / f"{name}.exe",
                source_dir / "Molecule" / name,
                source_dir / "Molecule" / f"{name}.exe",
            ]
            source_binary = next((path for path in candidates if path.is_file()), None)
            if source_binary is None:
                raise RuntimeError(f"Lilly build did not produce the {name} executable")
            destination = Path(installed[name])
            shutil.copy2(source_binary, destination)
            destination.chmod(0o755)

        marker.write_text(f"{LILLY_VERSION}\n", encoding="utf-8")

    return installed
