import unittest as ut
import datamol as dm
from medchem.rules import RuleFilters
from medchem.rules import basic_rules


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
