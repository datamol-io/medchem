from typer.testing import CliRunner

from medchem.cli import app

runner = CliRunner()


def test_cli():
    result = runner.invoke(app, ["--help"])
    assert result.exit_code == 0


def test_dummy():
    result = runner.invoke(app, ["dummy1"])
    assert result.exit_code == 0

    result = runner.invoke(app, ["dummy2"])
    assert result.exit_code == 0
