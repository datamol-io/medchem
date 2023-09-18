import pytest

import medchem as mc
import datamol as dm

from medchem.structural.lilly_demerits import LillyDemeritsFilters


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


def test_lilly_demerits():
    dfilters = LillyDemeritsFilters()

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


def test_demerits_invalid():
    dfilters = LillyDemeritsFilters()

    with pytest.raises(ValueError):
        dfilters(mols=[None, "CC9888", "CCCCO"])
