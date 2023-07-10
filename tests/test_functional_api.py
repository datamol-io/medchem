import datamol as dm
import medchem as mc


def test_alert_filter():
    data = dm.data.freesolv()
    data = data.iloc[:10]

    results = mc.functional.alert_filter(
        mols=data["smiles"].tolist(),
        alerts=["BMS"],
        n_jobs=-1,
        progress=True,
        return_idx=False,
    )

    assert results.tolist() == [True, False, True, True, True, True, True, True, True, True]

    results = mc.functional.alert_filter(
        mols=data["smiles"].tolist(),
        alerts=["BMS"],
        n_jobs=-1,
        progress=True,
        return_idx=True,
    )

    assert results.tolist() == [0, 2, 3, 4, 5, 6, 7, 8, 9]


def test_nibr_filter():
    data = dm.data.freesolv()
    data = data.iloc[:10]

    results = mc.functional.nibr_filter(
        mols=data["smiles"].tolist(),
        n_jobs=-1,
        progress=True,
        return_idx=False,
    )

    assert results.tolist() == [True, False, True, True, True, True, True, True, True, True]

    results = mc.functional.nibr_filter(
        mols=data["smiles"].tolist(),
        n_jobs=-1,
        progress=True,
        return_idx=True,
    )

    assert results.tolist() == [0, 2, 3, 4, 5, 6, 7, 8, 9]


def test_catalog_filter():
    data = dm.data.freesolv()
    data = data.iloc[:10]

    results = mc.functional.catalog_filter(
        mols=data["smiles"].tolist(),
        catalogs=["tox", "pains"],
        n_jobs=-1,
        progress=True,
        progress_leave=True,
        return_idx=False,
        batch_size=100,
        scheduler="threads",
    )

    assert results.tolist() == [True, True, False, True, False, True, True, True, True, True]

    results = mc.functional.catalog_filter(
        mols=data["smiles"].tolist(),
        catalogs=["tox", "pains"],
        n_jobs=-1,
        progress=True,
        progress_leave=True,
        return_idx=True,
        batch_size=100,
        scheduler="threads",
    )

    assert results.tolist() == [0, 1, 3, 5, 6, 7, 8, 9]


def test_chemical_group_filter():
    data = dm.data.freesolv()
    data = data.iloc[:10]

    cg = mc.groups.ChemicalGroup("common_organic_solvents")

    results = mc.functional.chemical_group_filter(
        mols=data["smiles"].tolist(),
        chemical_group=cg,
        n_jobs=-1,
        progress=True,
        progress_leave=True,
        return_idx=False,
        scheduler="threads",
    )

    assert results.tolist() == [False, True, False, True, False, False, True, False, False, False]

    results = mc.functional.chemical_group_filter(
        mols=data["smiles"].tolist(),
        chemical_group=cg,
        n_jobs=-1,
        progress=True,
        progress_leave=True,
        return_idx=True,
        scheduler="threads",
    )

    assert results.tolist() == [1, 3, 6]


def test_rules_filter():
    rfilter = mc.rules.RuleFilters(rule_list=["rule_of_five", "rule_of_oprea", "rule_of_cns"])

    data = dm.data.cdk2()
    data = data.iloc[:10]

    results = mc.functional.rules_filter(
        mols=data["smiles"].tolist(),
        rules=rfilter,
        return_idx=False,
    )

    assert results.tolist() == [False, False, False, True, True, False, False, False, False, True]

    results = mc.functional.rules_filter(
        mols=data["smiles"].tolist(),
        rules=rfilter,
        return_idx=True,
    )
    assert results.tolist() == [3, 4, 9]


def test_complexity_filter():
    data = dm.data.cdk2()
    data = data.iloc[:10]

    results = mc.functional.complexity_filter(
        mols=data["smiles"].tolist(),
        return_idx=False,
    )

    assert results.tolist() == [True, True, True, True, True, True, True, True, True, True]

    results = mc.functional.complexity_filter(
        mols=data["smiles"].tolist(),
        return_idx=True,
    )

    assert results.tolist() == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]


def test_bredt_filter():
    bredt_test_set = [
        "C1C2=C1C2",
        "C1CC=C=CC1",
        "C1C2CCCC=C12",
        "C1C2=C1CCCC2",  # is ok
        "C1CC2=CCC1C2",
        "C1CC2=CC1CC2",
        "C1CC2CCC1C=C2",  # is ok
        "C1CC2=CCC1CC2",
    ]

    results = mc.functional.bredt_filter(mols=bredt_test_set, return_idx=False)

    assert results.tolist() == [False, False, False, True, False, False, True, False]

    results = mc.functional.bredt_filter(mols=bredt_test_set, return_idx=True)

    assert results.tolist() == [3, 6]


def test_molecular_graph_filter():
    data = dm.data.freesolv()
    data = data.iloc[:10]

    results = mc.functional.molecular_graph_filter(mols=data["smiles"].tolist(), return_idx=False)

    assert results.tolist() == [True, True, True, True, True, True, True, True, True, True]

    results = mc.functional.molecular_graph_filter(mols=data["smiles"].tolist(), return_idx=True)

    assert results.tolist() == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]


def test_lilly_demerit_filter():
    data = dm.data.cdk2()
    data = data.iloc[:10]

    results = mc.functional.lilly_demerit_filter(mols=data["smiles"].tolist(), return_idx=True)

    assert results.tolist() == [0, 1, 2, 3, 4, 5, 6, 7]

    results = mc.functional.lilly_demerit_filter(mols=data["smiles"].tolist(), return_idx=False)
    assert results.tolist() == [True, True, True, True, True, True, True, True, False, False]


def test_protecting_groups_filter():
    data = dm.data.freesolv()
    data = data.iloc[:10]

    results = mc.functional.protecting_groups_filter(mols=data["smiles"].tolist(), return_idx=False)

    assert results.tolist() == [True, True, True, True, True, True, True, True, True, True]

    results = mc.functional.protecting_groups_filter(mols=data["smiles"].tolist(), return_idx=True)

    assert results.tolist() == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
