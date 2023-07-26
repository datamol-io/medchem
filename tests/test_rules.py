from typing import Callable

import pytest

import pandas as pd
import datamol as dm
import numpy as np

import medchem as mc


def test_list_descriptors():
    list_descs = mc.rules.list_descriptors()
    assert isinstance(list_descs, list)
    assert len(list_descs) > 0
    assert all(isinstance(x, str) for x in list_descs)


def test_has_spider_chains():
    mols_with_appendages = [
        "CCNCC1CCC(CNC(C)C)C1",
        "CCCCC1=CN(C=C1CCC)C1=CC(=CC=C1)C(C)CCO",
        "CCCCC1=CC2=C(OC(CN(C)CC)O2)C=C1",
    ]
    mols_without_appendages = [
        "CC(C)C(O)C(O)C(N)=O",
        "OC1=CC=CC(CC2CC2)=C1",
        "CCN(C)C1=CC=CC(CC(O)CC(C)O)=C1",
        "CN1CCC(CC(=O)NCC2=CC=CC=C2)C1",
        "CCCCC1=CC(=CC=C1)C(CC)C(C)C",
    ]

    assert all(
        mc.rules.has_spider_chains(dm.to_mol(x)) for x in mols_with_appendages
    ), "Fail appendage test for mols with appendages"
    assert all(
        not mc.rules.has_spider_chains(dm.to_mol(x)) for x in mols_without_appendages
    ), "Fail appendage test for mols without appendages"


def test_fraction_ring_system():
    mols = [
        "CCNCC1CCC(CNC(C)C)C1",
        "CCCCC1=CN(C=C1CCC)C1=CC(=CC=C1)C(C)CCO",
        "CCCCC1=CN(C=C1CCC)C1=CC(=CC=C1)C(C)CCCCC1CCO1",
        "CC(C)C(O)C(O)C(N)=O",
        "C1CCCCCC1",
    ]
    mols = [dm.to_mol(x) for x in mols]
    expected_results = [5 / 14, 11 / 23, 20 / 28, 0, 1]  # n_atoms_rings / n_atoms
    np.testing.assert_allclose(expected_results, [mc.rules.fraction_atom_in_scaff(x) for x in mols])


def test_n_fused_aromatic_rings():
    smiles = [
        "C1CNC2CCCNC2C1",
        "C1=CC2=C(C=C1)C1=C(C=C2)N=CC=C1",
        "C1CC1C1=CN=CC(=C1)C1=CC2=C(OCO2)C=C1",
        "C1OC2=C(O1)C=C(C=C2)C1=CN=C2C=CC=CC2=C1",
        "O=C(NC1=CC=C2N=CNC2=C1)C1=CC=C2C=CC=CC2=C1",
    ]
    mols = [dm.to_mol(x) for x in smiles]
    expected_results_pairwise = [0, 2, 0, 1, 2]
    expected_results = [0, 1, 0, 1, 2]
    expected_results_no_all_aro = [0, 1, 1, 2, 2]
    results = [mc.rules.n_fused_aromatic_rings(x) for x in mols]
    results_pairwise = [mc.rules.n_fused_aromatic_rings(x, pairwise=True) for x in mols]
    results_no_all_aro = [mc.rules.n_fused_aromatic_rings(x, require_all_aromatic=False) for x in mols]
    assert results == expected_results, "Fail fused aromatic ring test"
    assert expected_results_pairwise == results_pairwise, "Fail pairwise fused aromatic ring test"
    assert results_no_all_aro == expected_results_no_all_aro, "Fail not fully aromatic ring test"


def test_n_heavy_metals():
    smiles = [
        r"[Pd++].[H][C@]1(C)\C2=C\C3=C(C(C)=O)C(C)=C([N-]3)\C=C3/N=C(/C(/CC(=O)OC)=C4\[N-]\C(=C/C(=N2)[C@]1([H])CC)C(C)=C4C(=O)NCCS(O)(=O)=O)[C@@]([H])(CCC(O)=O)[C@]3([H])C",
        "[Na+].[Na+].[Na+].OC(CC([O-])=O)(CC([O-])=O)C([O-])=O",
        "CC(=O)OC1=CC=CC=C1C(O)=O",
    ]
    mols = [dm.to_mol(x) for x in smiles]
    output = [mc.rules.n_heavy_metals(m) for m in mols]
    expected_output = [1, 0, 0]
    assert output == expected_output, "Fail heavy atom test"


def test_rule_of_generative_stereo_case():
    smiles = [
        "CCC(O)C1=CC=C(C=C1)S(=O)(=O)C1CCNC(C)C1",
        "CC(F)C(O)C1=CC=C(C=C1)S(=O)(=O)C1CCNC(C)C1",
    ]
    mols = [dm.to_mol(x) for x in smiles]
    output = [mc.rules.basic_rules.rule_of_generative_design(x) for x in mols]
    excepted_output = [True, True]
    assert output == excepted_output, "Fail simple rule of generative"

    output = [mc.rules.basic_rules.rule_of_generative_design_strict(x) for x in mols]
    excepted_output = [True, False]  # second molecule fails due to n_stereo_center
    assert output == excepted_output, "Fail rule of generative stereo center test"


def test_basic_rules():
    data = dm.data.freesolv()
    data = data.iloc[:50]

    all_basic_rules = mc.rules.RuleFilters.list_available_rules()["name"].values
    for rule_name in all_basic_rules:
        rule_fn = getattr(mc.rules.basic_rules, rule_name)
        rule_fn(data["smiles"].values[10])


def test_rfilter():
    data = dm.data.freesolv()
    data = data.iloc[:50]

    def my_rules(mol, **kwargs):
        mol = dm.to_mol(mol)

        if mol is None:
            raise ValueError("Molecule is None")

        return mol.GetNumHeavyAtoms() / dm.descriptors.clogp(mol) > 2

    rfilter = mc.rules.RuleFilters(
        rule_list=[my_rules, "rule_of_five"],
        rule_list_names=["my_rules", "rule_of_five"],
    )
    out = rfilter(data["smiles"].tolist())
    assert out.shape[-1] == 5
    assert out.shape[0] == 50
    assert out["my_rules"].sum() == 45
    assert out["rule_of_five"].sum() == 48
    assert out["pass_all"].sum() == 44
    assert out["pass_any"].sum() == 49


def test_rule_basic_raise_error_if_none():
    all_basic_rules = mc.rules.RuleFilters.list_available_rules()["name"].values
    for rule_name in all_basic_rules:
        rule_fn = getattr(mc.rules.basic_rules, rule_name)

        with pytest.raises(ValueError):
            rule_fn(None)


def test_rule_filters_no_names():
    def my_rules(mol, **kwargs):
        return True

    rfilter = mc.rules.RuleFilters(rule_list=[my_rules, "rule_of_five"])

    assert list(rfilter.rules.keys()) == ["my_rules", "rule_of_five"]


def test_rule_filters_class():
    rfilter = mc.rules.RuleFilters(rule_list=["rule_of_three", "rule_of_five"])
    assert len(rfilter) == 2
    assert len(rfilter.rules) == 2
    assert isinstance(rfilter.rules["rule_of_three"], Callable)


def test_rule_filters_invalid():
    rfilter = mc.rules.RuleFilters(rule_list=["rule_of_five"])

    with pytest.raises(ValueError):
        rfilter(mols=[None, "CCCCO"], fail_if_invalid=True)

    results = rfilter(mols=[None, "CCCCO"], fail_if_invalid=False)
    assert results["mol"].isnull().sum() == 1

    assert pd.isna(results.iloc[0]["pass_all"])
    assert pd.isna(results.iloc[0]["pass_any"])
    assert pd.isna(results.iloc[0]["rule_of_five"])

    assert results.iloc[1]["pass_all"] is True
    assert results.iloc[1]["pass_any"] is True
    assert results.iloc[1]["rule_of_five"] is True


def test_rule_filters_list():
    ll = mc.rules.RuleFilters.list_available_rules()
    assert isinstance(ll, pd.DataFrame)
    assert "name" in ll.columns

    ll = mc.rules.RuleFilters.list_available_rules_names()
    assert isinstance(ll, list)
    assert all(isinstance(x, str) for x in ll)
    assert len(ll) > 0

    ll = mc.rules.RuleFilters.list_available_rules("building block")
    assert isinstance(ll, pd.DataFrame)
    assert "name" in ll.columns

    ll = mc.rules.RuleFilters.list_available_rules_names("building block", "fragment")
    assert isinstance(ll, list)
    assert all(isinstance(x, str) for x in ll)
    assert len(ll) > 0
