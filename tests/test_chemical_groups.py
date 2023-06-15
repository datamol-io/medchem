import unittest as ut
import datamol as dm
from medchem.groups import ChemicalGroup
from medchem.utils.loader import get_data_path


class Test_ChemicalGroup(ut.TestCase):
    def test_chemical_group(self):
        c_group = ChemicalGroup(groups="rings_in_drugs")
        rings_numbers = 92
        self.assertEqual(len(c_group), rings_numbers)
        self.assertEqual(len(c_group.smarts), rings_numbers)
        out_false = c_group.has_match("CCCCCCCCCC")
        self.assertFalse(out_false)
        out_true = c_group.has_match("C1CCCCC1")
        self.assertTrue(out_true)

    def test_chemical_group_query(self):
        c_group = ChemicalGroup(groups="rings_in_drugs")
        mol = dm.to_mol("CCS(=O)(=O)N1CC(C1)(CC#N)N2C=C(C=N2)C3=C4C=CNC4=NC=N3")
        # we have two rings of interest for smarts
        out = c_group.get_matches(mol, use_smiles=False)
        self.assertEqual(out.shape[0], 2)
        self.assertSetEqual(set(out.group.unique()), {"rings_in_drugs"})
        self.assertSetEqual(set(out.name.unique()), {"diazine", "1H-pyrrole"})
        # however, if we use smiles, we would have 3
        # this is a bug that needs to be fixed
        out_smiles = c_group.get_matches(mol, use_smiles=True)
        self.assertEqual(out_smiles.shape[0], 3)
        self.assertSetEqual(set(out_smiles.group.unique()), {"rings_in_drugs"})
        self.assertSetEqual(set(out_smiles.name.unique()), {"diazine", "1H-pyrrole", "1H-pyrazole"})

    def test_external_bank(self):
        c_group = ChemicalGroup(groups_db=get_data_path("smarts_bank.csv"))
        mol = dm.to_mol("CCS(=O)(=O)N1CC(C1)(CC#N)N2C=C(C=N2)C3=C4C=CNC4=NC=N3")
        out = c_group.get_matches(mol, use_smiles=False)
        expected_match = set(["HBA", "HBD", "Hydrogen", "SP3 Nitrogen", "SP2 Carbon"])
        self.assertTrue(expected_match.issubset(set(out["name"].values)))


if __name__ == "__main__":
    ut.main()
