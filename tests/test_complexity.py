import pytest

import datamol as dm

from medchem.complexity import WhitlockCT
from medchem.complexity import BaroneCT
from medchem.complexity import SMCM
from medchem.complexity import TWC
from medchem.complexity import ComplexityFilter


def test_available_list_cf():
    listed = set(ComplexityFilter.list_default_available_filters())

    assert listed == set(ComplexityFilter.COMPLEXITY_FNS.keys())

    with pytest.raises(ValueError):
        ComplexityFilter(limit="fake")

    with pytest.raises(ValueError):
        # sas is not defined for the zinc_filter
        ComplexityFilter(limit="99", complexity_metric="sas", threshold_stats_file="zinc_12")

    assert len(ComplexityFilter.list_default_percentile()) > 0


def test_whitlockct():
    """Test WhitlockCT"""

    smiles = [
        "CC(=O)O[C@@H]1C2=C(C)[C@H](C[C@@](O)([C@@H](OC(=O)C)[C@@H]3[C@@]4(CO[C@@H]4C[C@H](O)[C@@]3(C)C1=O)OC(C)=O)C2(C)C)OC(=O)[C@H](O)[C@@H](NC(=O)C)C",
        "C1(C)=CC(C)(C)C2(CCC3)CCC(C)C123",
        "O=C1C=CC(=O)C=C1",
    ]

    expected_outputs = [67, 20, 14]
    mols = [dm.to_mol(x) for x in smiles]

    assert [WhitlockCT(m) for m in mols] == expected_outputs

    with pytest.raises(ValueError):
        WhitlockCT(None)


def test_baronect():
    """Test BaroneCT"""
    smiles_non_chiral = [
        "CC1=CCCC1=O",
        "O=C1CCC2(C)C(=O)CCCC2=C1",
        "CC1(C)CCCC2(C)C3CCC(C2O)C13",
    ]
    smiles_chiral = [
        "O=C1CCC2(C)C(=O)CCCC2=C1",
        "O=NCNCCCCc1cccc(OC2=COC=C2)c1",
        "c1ccccc1C(=O)C1=NC=C(C=N1)CCCCCCNC(=O)C(N)N=O",
        "N1CNCC1SCC=CC=C(CC(C)=O)Cc1cccc(C)c1",
    ]
    expected_outputs = [135, 270, 288] + [290, 363, 512, 419]
    mols_non_chiral = [dm.to_mol(x) for x in smiles_non_chiral]
    mols_chiral = [dm.to_mol(x) for x in smiles_chiral]
    outputs = [BaroneCT(mol) for mol in mols_non_chiral] + [BaroneCT(mol, chiral=True) for mol in mols_chiral]

    assert outputs == expected_outputs


def test_smcm():
    """Test SMCM"""
    mols = [
        dm.to_mol(x)
        for x in [
            "CC1=C[C@@H]2[C@]([C@@H](C1=O)O)([C@]3([C@H]([C@H]([C@H]([C@@]34CO4)O2)O)OC(=O)C)C)CO",
            "CN(C)CCOCCOCCN(C)C",
            "OC(=O)CC(C)(O)CC(O)=O",
        ]
    ]
    expected_outputs = [59.376, 20.034, 20.485]
    outputs = [SMCM(mol) for mol in mols]
    assert outputs == expected_outputs


def test_twc():
    """Test TWC"""
    smiles = [
        "CC(C)C(C)C(C)CC(C)C",
        "C[Pt](C)(C)(C)(C)CCCC",
        "CC(C)C(C)CC(C)CCCC",
        "CC(C)CCC(C)C(C)CCC",
        "CCC(C)C(C(C)C)C(C)CC",
        "CCC(C)CC(C)CCCCC",
        "CCC(C)CCCC(CC)CC",
    ]

    mols = [dm.to_mol(x) for x in smiles]
    expected_outputs = [21784, 21784, 40145, 40145, 69926, 31474, 31474]
    outputs = [TWC(mol, log10=False) for mol in mols]
    assert outputs == expected_outputs

    TWC(mols[0], log10=True)


def test_complexity_filter():
    cf = ComplexityFilter(threshold_stats_file="zinc_12")

    smiles = dm.data.cdk2().smiles.iloc[:10].values

    expected_outputs = [False, True, False, False, False, False, True, True, True, False]
    outputs = [cf(dm.to_mol(x)) for x in smiles]

    cf2 = ComplexityFilter(
        limit="99",
        complexity_metric="sas",
        threshold_stats_file="zinc_15_available",
    )
    _ = [cf2(dm.to_mol(x)) for x in smiles]

    cf3 = ComplexityFilter(
        limit="max",
        complexity_metric="clogp",
        threshold_stats_file="zinc_15_available",
    )
    _ = [cf3(dm.to_mol(x)) for x in smiles]

    cf4 = ComplexityFilter(
        limit="999",
        complexity_metric="whitlock",
        threshold_stats_file="zinc_15_available",
    )
    _ = [cf4(dm.to_mol(x)) for x in smiles]

    assert outputs == expected_outputs
