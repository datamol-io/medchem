import typer

from loguru import logger

import medchem as mc
import datamol as dm
import pandas as pd

app = typer.Typer(help="The Medchem CLI", add_completion=False)


@app.command(help="Filtering for common structural alerts")
def common_alerts(
    input: str = typer.Argument(..., help="Input file path (CSV, Excel, Parquet, JSON, SDF)"),
    output: str = typer.Argument(..., help="Output file path (CSV, Excel, Parquet, JSON, SDF)"),
    smiles_column: str = typer.Option("smiles", help="The name of the column containing the SMILES"),
    n_jobs: int = typer.Option(-1, help="Number of jobs to use"),
    scheduler: str = typer.Option("auto", help="The scheduler to use"),
    keep_details: bool = typer.Option(False, help="Whether to keep the details or not"),
    progress: bool = typer.Option(True, help="Whether to show progress or not"),
    progress_leave: bool = typer.Option(False, help="Whether to leave the progress bar or not"),
):
    logger.info(f"Loading data from {input}")
    data = dm.io.open_df(input)

    logger.info(f"Filtering {data.shape[0]} molecules")
    filters_obj = mc.structural.CommonAlertsFilters()

    results = filters_obj(
        mols=data[smiles_column].tolist(),
        n_jobs=n_jobs,
        scheduler=scheduler,
        keep_details=keep_details,
        progress=progress,
        progress_leave=progress_leave,
    )

    results = pd.concat([data[[smiles_column]], results.drop(columns="mol")], axis=1)

    logger.info(f"Saving results to {output}")
    dm.io.save_df(results, output)


@app.command(help="Filtering with NIBR structural alerts")
def nibr_filters(
    input: str = typer.Argument(..., help="Input file path (CSV, Excel, Parquet, JSON, SDF)"),
    output: str = typer.Argument(..., help="Output file path (CSV, Excel, Parquet, JSON, SDF)"),
    smiles_column: str = typer.Option("smiles", help="The name of the column containing the SMILES"),
    n_jobs: int = typer.Option(-1, help="Number of jobs to use"),
    scheduler: str = typer.Option("threads", help="The scheduler to use"),
    keep_details: bool = typer.Option(False, help="Whether to keep the details or not"),
    progress: bool = typer.Option(True, help="Whether to show progress or not"),
    progress_leave: bool = typer.Option(False, help="Whether to leave the progress bar or not"),
):
    logger.info(f"Loading data from {input}")
    data = dm.io.open_df(input)

    logger.info(f"Filtering {data.shape[0]} molecules")
    filters_obj = mc.structural.NIBRFilters()

    results = filters_obj(
        mols=data[smiles_column].tolist(),
        n_jobs=n_jobs,
        scheduler=scheduler,
        keep_details=keep_details,
        progress=progress,
        progress_leave=progress_leave,
    )

    results = pd.concat([data[[smiles_column]], results.drop(columns="mol")], axis=1)

    logger.info(f"Saving results to {output}")
    dm.io.save_df(results, output)


if __name__ == "__main__":
    app()  # pragma: no cover
