import medchem as mc
import datamol as dm


def test_common_alerts():
    alerts = mc.structural.CommonAlerts()

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

    assert results.columns.tolist() == ["mol", "pass_filter", "status", "reasons"]


def test_common_alerts_invalid():
    alerts = mc.structural.CommonAlerts()

    results = alerts(mols=[None, "CC9888", "CCCCO"])

    assert results["mol"].isnull().sum() == 2
    assert results.shape == (3, 4)
    assert results["reasons"].tolist() == ["invalid", "invalid", None]


def test_common_alerts_list():
    l = mc.structural.CommonAlerts.list_default_available_alerts()
    assert l.columns.tolist() == ["rule_set_name", "smarts", "catalog_description", "rule_set", "source"]
