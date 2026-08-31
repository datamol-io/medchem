import io
import re
import shutil
import subprocess
from collections.abc import Sequence
from pathlib import Path

from loguru import logger

import pandas as pd


def find_lilly_binaries() -> dict[str, str]:
    """Find the command-line tools required by the Lilly filters.

    Returns:
        A mapping from each command name to its resolved executable path.

    Raises:
        ImportError: If any of the optional Lilly command-line tools is missing.
    """
    binaries_list = ["mc_first_pass", "tsubstructure", "iwdemerit"]
    binary_paths: dict[str, str] = {}
    for binary_name in binaries_list:
        binary_path = shutil.which(binary_name)

        if binary_path is None:
            raise ImportError(
                "The Lilly command-line tools required by `medchem.structural.lilly_demerits` "
                "are not installed. Install them with "
                "`mamba install -c conda-forge lilly-medchem-rules`."
            )

        binary_paths[binary_name] = binary_path

    return binary_paths


def rreplace(input_str: str, old: str, rep: str, occurrence: int) -> str:
    """Replace the last ``occurrence`` instances of ``old`` with ``rep``."""
    tmp = input_str.rsplit(old, occurrence)
    return rep.join(tmp)


def parse_output(rowlist: Sequence[str]) -> pd.DataFrame:
    """Parse content of `rowlist` to dataframe"""
    content = "\n".join(
        [
            rreplace(
                re.sub(r"\s+\([0-9]*\s+(matches\sto\s)", ' "', line.strip(), 1),
                "')",
                "'\"",
                1,
            ).strip("'")
            for line in rowlist
        ]
    )
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
