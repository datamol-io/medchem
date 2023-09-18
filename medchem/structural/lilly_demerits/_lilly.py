import io
import re
import shutil
import subprocess

from loguru import logger

import pandas as pd


def find_lilly_binaries():
    """Find the required binaries to run the Lilly medchem filters"""
    binaries_list = ["mc_first_pass", "tsubstructure", "iwdemerit"]
    binary_paths = {}
    for binary_name in binaries_list:
        binary_path = shutil.which(binary_name)

        if binary_path is None:
            raise ImportError(
                "The Lilly binaries required to use the `medchem.structural.demerits` module seems to be missing. "
                "Install with `mamba install lilly-medchem-rules`."
            )

        binary_paths[binary_name] = binary_path

    return binary_paths


def rreplace(input_str, old, rep, occurrence):
    """Replace last occurence of 'old' in 'input_str' by 'rep'"""
    tmp = input_str.rsplit(old, occurrence)
    return rep.join(tmp)


def parse_output(rowlist):
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


def run_cmd(cmd, shell=False):
    """Run a command"""
    res = subprocess.run(cmd, capture_output=True, shell=shell, check=False)
    if res.returncode != 0:
        logger.error("".join(res.stderr.decode("utf-8")))
        logger.error(" ".join(cmd))
        res.check_returncode()
    return res
