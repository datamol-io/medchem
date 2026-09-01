import io
import re
import shutil
import subprocess
import sys
import tempfile
from collections.abc import Sequence
from pathlib import Path

from loguru import logger

import pandas as pd

from ._installer import LILLY_VERSION
from ._installer import LILLY_VERSION_MARKER


def find_lilly_binaries() -> dict[str, str]:
    """Find the command-line tools required by the Lilly filters.

    Returns:
        A mapping from each command name to its resolved executable path.

    Raises:
        ImportError: If any of the optional Lilly command-line tools is missing.
    """
    binaries_list = ["mc_first_pass", "tsubstructure", "iwdemerit"]
    prefix_dirs = [
        Path(sys.executable).resolve().parent,
        Path(sys.prefix) / "bin",
        Path(sys.prefix) / "Scripts",
        Path(sys.prefix) / "Library" / "bin",
    ]
    binary_paths: dict[str, str] = {}
    for binary_name in binaries_list:
        binary_path = shutil.which(binary_name)

        # Direct interpreter invocation does not necessarily activate its
        # environment, so the conda/venv binary directory may be absent from
        # PATH. Resolve executables installed beside the interpreter too.
        if binary_path is None:
            candidates = [directory / binary_name for directory in prefix_dirs]
            if sys.platform == "win32":  # pragma: no cover - exercised in CI
                candidates.extend(directory / f"{binary_name}.exe" for directory in prefix_dirs)
            binary_path = next((str(path) for path in candidates if path.is_file()), None)

        if binary_path is None:
            raise ImportError(
                "The Lilly command-line tools required by `medchem.structural.lilly_demerits` "
                "are not installed. Install the checksum-pinned compatible release with "
                "`medchem install-lilly`."
            )

        binary_paths[binary_name] = binary_path

    markers = {
        Path(binary_path).resolve().parent / LILLY_VERSION_MARKER for binary_path in binary_paths.values()
    }
    if len(markers) != 1:
        raise ImportError("The Lilly command-line tools resolve from different installations")

    marker = markers.pop()
    installed_version = marker.read_text(encoding="utf-8").strip() if marker.is_file() else None
    if installed_version != LILLY_VERSION:
        version = installed_version or "unverified"
        raise ImportError(
            f"Lilly MedChem Rules {version} is not compatible with Medchem's {LILLY_VERSION} "
            "rule set. Run `medchem install-lilly --force`."
        )

    return binary_paths


def rreplace(input_str: str, old: str, rep: str, occurrence: int) -> str:
    """Replace the last ``occurrence`` instances of ``old`` with ``rep``."""
    tmp = input_str.rsplit(old, occurrence)
    return rep.join(tmp)


def parse_output(rowlist: Sequence[str]) -> pd.DataFrame:
    """Parse content of `rowlist` to dataframe"""
    rows = [
        rreplace(
            re.sub(
                r"\s+\([0-9]*\s+(matches\sto\s)",
                ' "',
                line.strip(),
                count=1,
            ),
            "')",
            "'\"",
            1,
        ).strip("'")
        for line in rowlist
        if line.strip()
    ]
    if not rows:
        return pd.DataFrame(
            {
                "smiles": pd.Series(dtype="object"),
                "ID": pd.Series(dtype="int64"),
                "reasons": pd.Series(dtype="object"),
            }
        )

    content = "\n".join(rows)
    flux = io.StringIO(content)
    df = pd.read_csv(flux, sep=" ", doublequote=True, names=["smiles", "ID", "reasons"])
    df["ID"] = pd.to_numeric(df["ID"])
    df["reasons"] = df["reasons"].apply(lambda x: x.strip("'") if x and isinstance(x, str) else x)
    return df


def run_cmd(cmd: Sequence[str | Path], shell: bool = False) -> subprocess.CompletedProcess:
    """Run a command"""
    res = subprocess.run(cmd, capture_output=True, shell=shell, check=False)
    if res.returncode != 0:
        logger.error("".join(res.stderr.decode("utf-8")))
        logger.error(" ".join(map(str, cmd)))
        res.check_returncode()
    return res


def run_pipeline(commands: Sequence[Sequence[str | Path]], output_path: str | Path) -> None:
    """Run commands as a portable subprocess pipeline.

    This uses OS pipes directly rather than a shell command string, so the
    Lilly pipeline has the same execution model on POSIX and Windows.
    """
    if not commands:
        raise ValueError("A pipeline requires at least one command")

    processes: list[subprocess.Popen] = []
    error_streams = [tempfile.TemporaryFile() for _ in commands]
    try:
        with Path(output_path).open("wb") as output:
            previous_stdout = None
            for index, command in enumerate(commands):
                is_last = index == len(commands) - 1
                process = subprocess.Popen(
                    command,
                    stdin=previous_stdout,
                    stdout=output if is_last else subprocess.PIPE,
                    stderr=error_streams[index],
                )
                if previous_stdout is not None:
                    previous_stdout.close()
                previous_stdout = process.stdout
                processes.append(process)

            for process in reversed(processes):
                process.wait()

        for command, process, error_stream in zip(commands, processes, error_streams):
            if process.returncode:
                error_stream.seek(0)
                stderr = error_stream.read().decode("utf-8", errors="replace")
                logger.error(stderr)
                logger.error(" ".join(map(str, command)))
                raise subprocess.CalledProcessError(
                    process.returncode,
                    command,
                    stderr=stderr,
                )
    finally:
        for process in processes:
            if process.poll() is None:
                process.kill()
                process.wait()
        for error_stream in error_streams:
            error_stream.close()
