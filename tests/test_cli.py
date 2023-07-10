from typer.testing import CliRunner

import datamol as dm

from medchem.cli import app

runner = CliRunner()


def test_cli():
    result = runner.invoke(app, ["--help"])
    assert result.exit_code == 0


def test_common_alerts(tmp_path):
    input_path = str(tmp_path / "input.csv")
    output_path = str(tmp_path / "output.csv")

    # Load a dataset
    data = dm.data.solubility()
    data = data.sample(50, random_state=20)

    data.to_csv(input_path, index=False)

    result = runner.invoke(
        app,
        [
            "common-alerts",
            str(input_path),
            str(output_path),
            "--smiles-column",
            "smiles",
            "--keep-details",
        ],
    )
    assert result.exit_code == 0

    # Check the output
    results = dm.io.open_df(output_path)
    assert set(results.columns) == {"smiles", "pass_filter", "status", "reasons", "details"}
    assert results.shape[0] == 50


def test_nibr_filters(tmp_path):
    input_path = str(tmp_path / "input.csv")
    output_path = str(tmp_path / "output.csv")

    # Load a dataset
    data = dm.data.solubility()
    data = data.sample(50, random_state=20)

    data.to_csv(input_path, index=False)

    result = runner.invoke(
        app,
        [
            "nibr-filters",
            str(input_path),
            str(output_path),
            "--smiles-column",
            "smiles",
            "--keep-details",
        ],
    )
    assert result.exit_code == 0

    # Check the output
    results = dm.io.open_df(output_path)
    assert set(results.columns) == {
        "smiles",
        "pass_filter",
        "reasons",
        "severity",
        "status",
        "n_covalent_motif",
        "special_mol",
        "details",
    }
    assert results.shape[0] == 50
