import unittest as ut
import datamol as dm
from medchem.rules import RuleFilters
from medchem.rules import basic_rules
from medchem.rules._utils import n_fused_aromatic_rings, n_heavy_metals, has_flagels


def test_has_flagels():
    mols_with_flagels = [
        "CCNCC1CCC(CNC(C)C)C1",
        "CCCCC1=CN(C=C1CCC)C1=CC(=CC=C1)C(C)CCO",
        "CCCCC1=CC2=C(OC(CN(C)CC)O2)C=C1",
    ]
    mols_without_flagels = [
        "CC(C)C(O)C(O)C(N)=O",
        "OC1=CC=CC(CC2CC2)=C1",
        "CCN(C)C1=CC=CC(CC(O)CC(C)O)=C1",
        "CN1CCC(CC(=O)NCC2=CC=CC=C2)C1",
        "CCCCC1=CC(=CC=C1)C(CC)C(C)C",
    ]

    assert all(
        has_flagels(dm.to_mol(x)) for x in mols_with_flagels
    ), "Fail flagel test for mols with flagels"
    assert all(
        not has_flagels(dm.to_mol(x)) for x in mols_without_flagels
    ), "Fail flagel test for mols without flagels"


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
    results = [n_fused_aromatic_rings(x) for x in mols]
    results_pairwise = [n_fused_aromatic_rings(x, pairwise=True) for x in mols]
    results_no_all_aro = [
        n_fused_aromatic_rings(x, require_all_aromatic=False) for x in mols
    ]
    assert results == expected_results, "Fail fused aromatic ring test"
    assert (
        expected_results_pairwise == results_pairwise
    ), "Fail pairwise fused aromatic ring test"
    assert (
        results_no_all_aro == expected_results_no_all_aro
    ), "Fail not fully aromatic ring test"


def test_n_heavy_metals():
    smiles = [
        "[Pd++].[H][C@]1(C)\C2=C\C3=C(C(C)=O)C(C)=C([N-]3)\C=C3/N=C(/C(/CC(=O)OC)=C4\[N-]\C(=C/C(=N2)[C@]1([H])CC)C(C)=C4C(=O)NCCS(O)(=O)=O)[C@@]([H])(CCC(O)=O)[C@]3([H])C",
        "[Na+].[Na+].[Na+].OC(CC([O-])=O)(CC([O-])=O)C([O-])=O",
        "CC(=O)OC1=CC=CC=C1C(O)=O",
    ]
    mols = [dm.to_mol(x) for x in smiles]
    output = [n_heavy_metals(m) for m in mols]
    expected_output = [1, 0, 0]
    assert output == expected_output, "Fail heavy atom test"


class Test_Rules(ut.TestCase):

    data = dm.data.freesolv()

    def test_basic_rules(self):
        all_basic_rules = RuleFilters.list_available_rules()["name"].values
        for rule_name in all_basic_rules:
            rule_fn = getattr(basic_rules, rule_name)
            rule_fn(self.data.smiles.values[10])

    def test_rule_obj(self):
        def my_rules(mol, **kwargs):
            mol = dm.to_mol(mol)
            return mol.GetNumHeavyAtoms() / dm.descriptors.clogp(mol) > 2

        rule_obj = RuleFilters(
            rule_list=[my_rules, "rule_of_five"],
            rule_list_names=["my_rules", "rule_of_five"],
        )
        out = rule_obj(self.data.smiles.values[:100])
        self.assertEqual(out.shape[-1], 2)
        self.assertEqual(out.shape[0], 100)
        self.assertEqual(out["my_rules"].sum(), 88)
        self.assertEqual(out["rule_of_five"].sum(), 96)


if __name__ == "__main__":
    ut.main()
