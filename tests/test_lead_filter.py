import unittest as ut
import datamol as dm
import pandas as pd
import numpy as np
from medchem import catalog
from medchem.filter import lead
from medchem.alerts import AlertFilters, NovartisFilters
from medchem.groups import ChemicalGroup

from rdkit.Chem.Descriptors import MolWt


class Test_LeadFilter(ut.TestCase):

    data = dm.data.freesolv()
    screening_smiles_set = [
        # smiles, covalent, special, severity
        ("Nc1nc(Sc2cccc(Cl)c2)c(C#N)c(c3ccc4OCCOc4c3)c1C#N", 1, 0, 10),
        ("CCOC(=O)NC1(NCc2ccccc2)Oc3ccccc3O1", 0, 0, 10),
        ("Oc1nnc(SCC(=O)N2CCOCC2)n1c3ccccc3", 0, 0, 1),
        ("Fc1ccc(C(=O)NNC(=O)c2ccc(cc2)n3cnnn3)c(F)c1", 0, 0, 0),
        (
            "COc1cc(Cl)ccc1N=C(S)N(CCN2CCOCC2)C3CCN(CC3)C(=O)C",
            float("nan"),
            float("nan"),
            0,
        ),
        ("CCCCCCCCCCCCCNC(=O)[C@H](CO)\\N=C\\c1ccccc1", 0, 0, 2),
    ]
    bredt_test_set = [
        "C1C2=C1C2",
        "C1CC=C=CC1",
        "C1C2CCCC=C12",
        "C1C2=C1CCCC2",  # is ok
        "C1CC2=CCC1C2",
        "C1CC2=CC1CC2",
        "C1CC2CCC1C=C2",  # is ok
        "C1CC2=CCC1CC2",
    ] + list(data.smiles.values[:10])

    def test_alert_filter(self):
        ok_mols = lead.alert_filter(
            self.data.smiles.values, alerts=["Glaxo", "BMS"], n_jobs=2, return_idx=False
        )
        ok_index_2 = lead.alert_filter(
            self.data.smiles.values, alerts=["Glaxo", "bms"], n_jobs=2, return_idx=True
        )
        self.assertEqual(sum(ok_mols), len(ok_index_2))
        self.assertTrue(sum(ok_mols) == 503)

        # test filter object
        rule_dict = dict(MW=[100, 200])
        ok_index_3 = lead.alert_filter(
            self.data.smiles.values,
            alerts=["Glaxo", "bms"],
            n_jobs=2,
            rule_dict=rule_dict,
            return_idx=True,
        )
        self.assertTrue(
            all(
                (x >= 100 and x <= 200)
                for x in self.data.iloc[ok_index_3].smiles.apply(dm.to_mol).apply(MolWt)
            )
        )
        al_filter = AlertFilters(alerts_set=["BMS"])
        df = al_filter(self.data.smiles.values[:10])
        expected_cols = set(
            ["_smiles", "status", "reasons", "MW", "LogP", "HBD", "HBA", "TPSA"]
        )
        self.assertTrue(
            len(expected_cols.intersection(set(df.columns))) == len(expected_cols)
        )

    def test_catalog_filter(self):
        idx = lead.catalog_filter(
            self.data.smiles.values,
            catalogs=["pains"],
            return_idx=False,
        )
        self.assertTrue(sum(idx) == 632)
        idx = lead.catalog_filter(
            self.data.smiles.values,
            catalogs=["brenk"],
            return_idx=True,
        )
        # brenk should filter alkyl halide
        self.assertTrue(640 not in idx)
        # ensure not error is raised on this
        for catalog_name in catalog.list_named_catalogs():
            lead.catalog_filter(self.data.smiles.values, catalogs=[catalog_name])

    def test_screening_filter(self):
        df = pd.DataFrame.from_records(
            self.screening_smiles_set,
            columns=["smiles", "covalent", "special_mol", "severity"],
        )
        idx = lead.screening_filter(
            df.smiles, n_jobs=2, max_severity=5, return_idx=True
        )
        expected_idx = list(df[df.severity < 5].index)
        self.assertEqual(len(set(idx).difference(set(expected_idx))), 0)

        nov_filters = NovartisFilters()
        out = nov_filters(df.smiles)
        expected_cols = set(
            [
                "_smiles",
                "status",
                "reasons",
                "severity",
                "covalent",
                "special_mol",
            ]
        )
        self.assertTrue(
            len(expected_cols.intersection(set(out.columns))) == len(expected_cols)
        )
        np.testing.assert_array_equal(df.severity, out.severity)
        np.testing.assert_array_equal(df.covalent, out.covalent)
        np.testing.assert_array_equal(df.special_mol, out.special_mol)

    def test_bredt_filter(self):
        """Test whether the input molecules pass all bredt filters"""

        output = lead.bredt_filter(self.bredt_test_set, n_jobs=2)
        expected_results = [False, False, False, True, False, False, True, False] + [
            True
        ] * 10
        np.testing.assert_array_equal(output, expected_results)

    def test_chemical_groups(self):
        """Test whether the input molecules contains any privileged scaffold"""

        c_group = ChemicalGroup(groups=["privileged_scaffolds"])
        mols = self.data.smiles.values
        output = lead.chemical_group_filter(mols, c_group)
        self.assertEqual(len(output), len(mols))
        self.assertEqual(sum(output), 622)

    def test_rules_filter(self):
        """Test rule filtering"""
        out = lead.rules_filter(
            self.data.smiles.values,
            rules=["rule_of_five", "rule_of_gsk_4_400"],
            n_jobs=2,
        )
        self.assertEqual(sum(out), 598)


if __name__ == "__main__":
    ut.main()
