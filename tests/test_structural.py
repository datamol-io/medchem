import sys
import subprocess

import pandas as pd
import pytest

import medchem as mc
import datamol as dm

from medchem.structural.lilly_demerits import LillyDemeritsFilters
from medchem.structural.lilly_demerits._lilly import find_lilly_binaries
from medchem.structural.lilly_demerits._lilly import materialize_query_manifest
from medchem.structural.lilly_demerits._lilly import run_pipeline

try:
    find_lilly_binaries()
    HAS_LILLY_TOOLS = True
except ImportError:
    HAS_LILLY_TOOLS = False
requires_lilly = pytest.mark.skipif(
    not HAS_LILLY_TOOLS,
    reason="requires the optional Lilly MedChem Rules command-line tools",
)


def test_common_alerts():
    alerts = mc.structural.CommonAlertsFilters()

    data = dm.data.freesolv()
    data = data.iloc[:50]

    data["mol"] = data["smiles"].apply(dm.to_mol)

    results = alerts(
        mols=data["mol"].tolist(),
        n_jobs=-1,
        scheduler="processes",
        keep_details=True,
    )

    assert results["pass_filter"].sum() == 44
    assert results["reasons"].unique().tolist() == [
        None,
        "halogen_heteroatom;sulfonyl_halide",
        "primary_halide_sulfate",
        "non_ring_CH2O_acetal;phosphorus_sulfur_bond",
        "aldehyde",
        "gte_10_carbon_sb_chain;gte_8_CF2_or_CH2",
    ]

    assert set(results.columns.tolist()) == {"mol", "pass_filter", "status", "reasons", "details"}


def test_common_alerts_invalid():
    alerts = mc.structural.CommonAlertsFilters()

    results = alerts(mols=[None, "CC9888", "CCCCO"])

    assert results["mol"].isnull().sum() == 2
    assert results.shape == (3, 4)
    assert results["reasons"].tolist() == ["invalid", "invalid", None]


def test_common_alerts_list():
    ll = mc.structural.CommonAlertsFilters.list_default_available_alerts()
    assert ll.columns.tolist() == ["rule_set_name", "smarts", "catalog_description", "rule_set", "source"]


def test_nibr():
    nibr_filters = mc.structural.NIBRFilters()

    data = dm.data.solubility()
    data = data.iloc[:50]

    results = nibr_filters(
        mols=data["mol"].tolist(),
        n_jobs=-1,
        scheduler="processes",
        keep_details=True,
    )

    assert results["pass_filter"].sum() == 49
    assert set(results.columns.tolist()) == {
        "mol",
        "reasons",
        "severity",
        "status",
        "n_covalent_motif",
        "special_mol",
        "pass_filter",
        "details",
    }


def test_nibr_invalid():
    nibr_filters = mc.structural.NIBRFilters()

    results = nibr_filters(mols=[None, "CC9888", "CCCCO"])

    assert results["mol"].isnull().sum() == 2
    assert results.shape == (3, 7)
    assert results["reasons"].tolist() == ["invalid", "invalid", None]
    assert set(results.columns.tolist()) == {
        "mol",
        "reasons",
        "severity",
        "status",
        "n_covalent_motif",
        "special_mol",
        "pass_filter",
    }


@pytest.mark.integration
@pytest.mark.lilly
@requires_lilly
def test_lilly_demerits():
    dfilters = LillyDemeritsFilters()

    data = dm.data.solubility()
    data = data.sample(50, random_state=20)

    results = dfilters(mols=data["mol"].tolist())

    assert results["pass_filter"].sum() == 29
    assert set(results.columns.tolist()) == {
        "smiles",
        "reasons",
        "step",
        "demerit_score",
        "status",
        "pass_filter",
        "mol",
    }


@pytest.mark.integration
@pytest.mark.lilly
@requires_lilly
def test_lilly_demerits_config():
    test_config = {
        "output": "test",
        "min_atoms": 7,
        "soft_max_atoms": 30,
        "hard_max_atoms": 50,
        "smarts": [],
        "nodemerit": False,
        "dthresh": 160,
        "odm": [],
        "okiso": False,
        "noapdm": False,
    }

    dfilters = LillyDemeritsFilters(**test_config)

    data = dm.data.solubility()
    data = data.sample(50, random_state=20)

    results = dfilters(mols=data["mol"].tolist())

    assert results["pass_filter"].sum() == 30
    assert set(results.columns.tolist()) == {
        "smiles",
        "reasons",
        "step",
        "demerit_score",
        "status",
        "pass_filter",
        "mol",
    }


@pytest.mark.integration
@pytest.mark.lilly
@requires_lilly
def test_lilly_21_regression_examples():
    """Pin rule-query and executable behavior that changed after 1.0.1."""
    smiles = [
        "C1=CC=CC(=C1C(C(=O)OCC)C[N+](=O)[O-])Cl",
        "C1=CC=CC=C1C(C(=O)OCC)C[N+](=O)[O-]",
        "C1=CC=C(Cl)C=C1C(C(=O)OCC)C[N+](=O)[O-]",
        "C1=CC(=CC=C1C(C(=O)OCC)C[N+](=O)[O-])C",
        "C1=CC(=CC=C1C(C(=O)OCC)C[N+](=O)[O-])C(F)(F)F",
        "C(=O)(O)CCCC=CCC=CCC=CCC(O)C(O)CCCCC",
        "C(=O)(O)CCCC=CCC=CCC=CCC=CCCCCC",
        "C(=O)(O)CCCC=CCC=CCC(O)C(O)CC=CCCCCC",
    ]

    results = LillyDemeritsFilters(
        min_atoms=7,
        soft_max_atoms=25,
        hard_max_atoms=40,
    )(smiles)

    assert results["pass_filter"].tolist() == [False] * 5 + [True] * 3
    assert results["demerit_score"].tolist() == [110, 110, 110, 110, 110, 65, 90, 90]


@pytest.mark.integration
@pytest.mark.lilly
@requires_lilly
def test_lilly_uses_reference_valence_model_for_raw_smiles():
    smiles = "P(F)(F)(F)(F)(F)F.C1=CC(=CC=C1N(C)C)N(C)C"

    result = LillyDemeritsFilters()([smiles])

    assert result.loc[0, "smiles"] == smiles
    assert result.loc[0, "reasons"] == "phenylenediamine"


@pytest.mark.integration
@pytest.mark.lilly
@requires_lilly
def test_lilly_defaults_match_upstream_kill_thresholds():
    smiles = "Brc1cc(N=c2[nH]c(=N)[n]c(CSc3c(O)cc4c(=O)[n](ccc4c3)C)[nH]2)c(C)cc1"

    result = LillyDemeritsFilters()([smiles])

    assert not result.loc[0, "pass_filter"]
    assert result.loc[0, "demerit_score"] == 134
    assert "too_many_atoms:D40" in result.loc[0, "reasons"]


@pytest.mark.integration
@pytest.mark.lilly
@requires_lilly
def test_lilly_handles_smiles_suppressed_by_the_native_tools():
    result = LillyDemeritsFilters()(["C1CC"])

    assert result.loc[0, "smiles"] == "C1CC"
    assert result.loc[0, "reasons"] == "lillymol_invalid"
    assert result.loc[0, "status"] == "exclude"
    assert not result.loc[0, "pass_filter"]


@pytest.mark.integration
@pytest.mark.lilly
@requires_lilly
@pytest.mark.parametrize("stop_after_step", [0, 1, 2])
def test_lilly_can_stop_after_intermediate_steps(stop_after_step):
    smiles = ["CCO", "[Na+].[Cl-]", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"]

    result = LillyDemeritsFilters(stop_after_step=stop_after_step)(smiles)

    assert result["mol"].tolist() == smiles
    assert result["pass_filter"].dtype == bool
    assert result["step"].max() <= stop_after_step + 1


def test_lilly_rejects_invalid_stop_step():
    with pytest.raises(ValueError, match="between 0 and 3"):
        LillyDemeritsFilters(stop_after_step=4)


def test_lilly_forwards_parallel_batch_configuration(monkeypatch):
    captured = {}

    def fake_parallelized(_fn, batches, **kwargs):
        captured.update(kwargs)
        return [pd.DataFrame({"batch_size": [len(batch)]}) for batch in batches]

    monkeypatch.setattr(dm, "parallelized", fake_parallelized)

    result = LillyDemeritsFilters()(["CC"] * 11, n_jobs=3, batch_size=5)

    assert result["batch_size"].tolist() == [4, 4, 3]
    assert captured["n_jobs"] == 3
    assert captured["scheduler"] == "threads"


def test_lilly_pipeline_connects_native_stages(tmp_path):
    output = tmp_path / "pipeline-output.txt"
    commands = [
        [sys.executable, "-c", "import sys; sys.stdout.write('alpha\\nbeta\\n')"],
        [sys.executable, "-c", "import sys; sys.stdout.write(sys.stdin.read().upper())"],
    ]

    run_pipeline(commands, output)

    assert output.read_text(encoding="utf-8") == "ALPHA\nBETA\n"


def test_lilly_pipeline_propagates_stage_failures(tmp_path):
    output = tmp_path / "pipeline-output.txt"
    commands = [[sys.executable, "-c", "import sys; sys.stderr.write('failed stage'); sys.exit(7)"]]

    with pytest.raises(subprocess.CalledProcessError) as error:
        run_pipeline(commands, output)

    assert error.value.returncode == 7
    assert error.value.stderr == "failed stage"


def test_lilly_pipeline_requires_a_stage(tmp_path):
    with pytest.raises(ValueError, match="at least one command"):
        run_pipeline([], tmp_path / "pipeline-output.txt")


def test_lilly_query_manifests_use_resolvable_relative_paths(tmp_path):
    query_dir = tmp_path / "queries"
    query_dir.mkdir()
    first_query = query_dir / "first.qry"
    second_query = query_dir / "second.qry"
    first_query.touch()
    second_query.touch()
    source = query_dir / "manifest"
    source.write_bytes(f"first.qry\r\n\r\n{second_query.resolve()}\r\n".encode())
    destination = tmp_path / "resolved-manifest"

    result = materialize_query_manifest(source, destination)

    assert result == str(destination.resolve())
    entries = destination.read_text(encoding="utf-8").splitlines()
    assert [(destination.parent / entry).resolve() for entry in entries] == [
        first_query.resolve(),
        second_query.resolve(),
    ]
    assert b"\r" not in destination.read_bytes()


def test_lilly_query_manifest_rejects_missing_files(tmp_path):
    source = tmp_path / "manifest"
    source.write_text("missing.qry\n", encoding="utf-8")

    with pytest.raises(FileNotFoundError, match="missing.qry"):
        materialize_query_manifest(source, tmp_path / "resolved-manifest")


@pytest.mark.integration
@pytest.mark.lilly
@requires_lilly
def test_demerits_invalid():
    dfilters = LillyDemeritsFilters()

    with pytest.raises(ValueError):
        dfilters(mols=[None, "CC9888", "CCCCO"])


def test_lilly_tools_are_resolved_lazily(monkeypatch):
    import medchem.structural.lilly_demerits._demerits as demerits_module

    def missing_tools():
        raise ImportError("missing Lilly tools")

    monkeypatch.setattr(demerits_module, "find_lilly_binaries", missing_tools)

    dfilters = LillyDemeritsFilters()
    with pytest.raises(ImportError, match="missing Lilly tools"):
        dfilters(mols=["CC"])


def test_lilly_tools_are_found_beside_interpreter(monkeypatch, tmp_path):
    import medchem.structural.lilly_demerits._lilly as lilly_module

    interpreter = tmp_path / "bin" / "python"
    interpreter.parent.mkdir()
    interpreter.touch()
    for name in ("mc_first_pass", "tsubstructure", "iwdemerit"):
        (interpreter.parent / name).touch()
    (interpreter.parent / ".lilly-medchem-rules-version").write_text("2.1.0\n")

    monkeypatch.setattr(lilly_module.shutil, "which", lambda name: None)
    monkeypatch.setattr(sys, "executable", str(interpreter))
    monkeypatch.setattr(sys, "prefix", str(tmp_path))

    binaries = lilly_module.find_lilly_binaries()

    assert binaries == {
        name: str(interpreter.parent / name) for name in ("mc_first_pass", "tsubstructure", "iwdemerit")
    }


def test_lilly_windows_layout_uses_scripts_and_exe(monkeypatch, tmp_path):
    import medchem.structural.lilly_demerits._installer as installer_module
    import medchem.structural.lilly_demerits._lilly as lilly_module

    scripts = tmp_path / "Scripts"
    scripts.mkdir()
    expected = {}
    for name in ("mc_first_pass", "tsubstructure", "iwdemerit"):
        binary = scripts / f"{name}.exe"
        binary.touch()
        expected[name] = str(binary)
    (scripts / ".lilly-medchem-rules-version").write_text("2.1.0\n")

    monkeypatch.setattr(sys, "platform", "win32")
    monkeypatch.setattr(sys, "prefix", str(tmp_path))
    monkeypatch.setattr(sys, "executable", str(tmp_path / "python.exe"))
    monkeypatch.setattr(lilly_module.shutil, "which", lambda name: None)
    monkeypatch.setenv("MSYSTEM", "MSYS")

    assert installer_module.install_lilly_rules(prefix=tmp_path) == expected
    assert lilly_module.find_lilly_binaries() == expected


def test_lilly_windows_build_requires_msys2(monkeypatch, tmp_path):
    import medchem.structural.lilly_demerits._installer as installer_module

    monkeypatch.setattr(sys, "platform", "win32")
    monkeypatch.delenv("MSYSTEM", raising=False)

    with pytest.raises(RuntimeError, match="MSYS2 MSYS shell"):
        installer_module.install_lilly_rules(prefix=tmp_path)


def test_lilly_rejects_unverified_binary_version(monkeypatch, tmp_path):
    import medchem.structural.lilly_demerits._lilly as lilly_module

    binaries = {}
    for name in ("mc_first_pass", "tsubstructure", "iwdemerit"):
        binary = tmp_path / name
        binary.touch()
        binaries[name] = str(binary)

    monkeypatch.setattr(lilly_module.shutil, "which", binaries.get)

    with pytest.raises(ImportError, match="unverified.*medchem install-lilly"):
        lilly_module.find_lilly_binaries()
