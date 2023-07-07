import pytest

import medchem as mc
import datamol as dm

from medchem.structural.demerits import DemeritsFilters


def test_common_alerts():
    alerts = mc.structural.CommonAlertsFilters()

    data = dm.data.solubility()
    data = data.iloc[:50]

    results = alerts(mols=data["mol"].tolist(), n_jobs=-1, scheduler="auto")

    assert results["pass_filter"].sum() == 29
    assert results["reasons"].unique().tolist() == [
        None,
        "Aliphatic long chain",
        "isolated alkene",
        "Aliphatic long chain;isolated alkene",
        "I1 Aliphatic methylene chains 7 or more long;Aliphatic long chain;isolated alkene",
        "polyene",
        "triple bond",
        "Aliphatic long chain;triple bond",
        "I1 Aliphatic methylene chains 7 or more long;Aliphatic long chain;triple bond",
    ]

    assert set(results.columns.tolist()) == {"mol", "pass_filter", "status", "reasons"}


def test_common_alerts_invalid():
    alerts = mc.structural.CommonAlertsFilters()

    results = alerts(mols=[None, "CC9888", "CCCCO"])

    assert results["mol"].isnull().sum() == 2
    assert results.shape == (3, 4)
    assert results["reasons"].tolist() == ["invalid", "invalid", None]


def test_common_alerts_list():
    l = mc.structural.CommonAlertsFilters.list_default_available_alerts()
    assert l.columns.tolist() == ["rule_set_name", "smarts", "catalog_description", "rule_set", "source"]


def test_nibr():
    nibr_filters = mc.structural.NIBRFilters()

    data = dm.data.solubility()
    data = data.iloc[:50]

    results = nibr_filters(mols=data["mol"].tolist(), n_jobs=-1, scheduler="auto", keep_details=True)

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


def test_demerits():
    dfilters = DemeritsFilters()

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


def test_demerits_config():
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

    dfilters = DemeritsFilters(**test_config)

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


def test_demerits_invalid():
    dfilters = DemeritsFilters()

    with pytest.raises(ValueError):
        dfilters(mols=[None, "CC9888", "CCCCO"])
