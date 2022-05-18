import sys
import os

import typer
import numpy as np

from loguru import logger
from tqdm.auto import tqdm

import pandas as pd
import datamol as dm
import medchem
from rdkit.Chem.GraphDescriptors import BertzCT

from medchem.complexity import _complexity_calc as calc

app = typer.Typer()


def _silent(func):
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except:
            return None

    return wrapper


def compute_props(row):
    row = row.copy()
    with dm.without_rdkit_log():
        row["mol"] = row["smiles"].apply(dm.to_mol)
        row["bertz"] = row["mol"].apply(BertzCT)
        row["sas"] = row["mol"].apply(_silent(dm.descriptors.sas))
        row["qed"] = row["mol"].apply(_silent(dm.descriptors.qed))
        row["whitlock"] = row["mol"].apply(_silent(calc.WhitlockCT))
        row["barone"] = row["mol"].apply(_silent(calc.BaroneCT))
        row["smcm"] = row["mol"].apply(_silent(calc.SMCM))
        row["twc"] = row["mol"].apply(_silent(calc.TWC))
        row = row.drop(columns=["mol"]).reset_index(drop=True)
    return row


@app.command()
def process(
    input_path: str = typer.Option(..., help="Path to the train the csv file to use."),
    output_path: str = typer.Option(..., help="Output path"),
    partition_name: int = typer.Option(..., help="Partition to use"),
):

    # Parameters
    smiles_col = "smiles"
    verbose = True
    batch_size = 256
    columns = ["smiles", "zinc_id", "parquet_partition"]

    # Load data
    logger.info(f"Loading partition: {partition_name}")
    data = pd.read_parquet(
        input_path,
        filters=[("parquet_partition", "=", partition_name)],
        columns=columns,
    )

    # Load data
    logger.info("Build molecule objects")
    data["mol"] = dm.parallelized(
        dm.to_mol, data[smiles_col].values, n_jobs=-1, progress=verbose
    )
    data = data.reset_index(drop=True)
    # Sanity reset here
    logger.info("Processing dataset")
    n_splits = len(data) // batch_size
    processed_data = dm.parallelized(
        compute_props, np.array_split(data, n_splits), n_jobs=-1, progress=verbose
    )
    processed_data = pd.concat(processed_data, ignore_index=True)
    processed_data.to_csv(output_path, index=False)

    processed_data = data.drop(columns=["mol"])
    processed_data.to_csv(
        output_path.format(partition_name=partition_name), index=False
    )
    logger.info("Done filtering")


if __name__ == "__main__":
    app()
