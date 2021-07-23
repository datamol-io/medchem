import unittest as ut
import datamol as dm
import pandas as pd
from medchem.filter import lead
from rdkit.Chem.Descriptors import MolWt


class Test_LeadFilter(ut.TestCase):

    data = dm.data.freesolv()
    screening_smiles_set = [
        # smiles, covalent, special, severity
        ("Nc1nc(Sc2cccc(Cl)c2)c(C#N)c(c3ccc4OCCOc4c3)c1C#N", 1, 0, 10),
        ("CCOC(=O)NC1(NCc2ccccc2)Oc3ccccc3O1", 0, 0, 10),
        ("Oc1nnc(SCC(=O)N2CCOCC2)n1c3ccccc3", 0, 0, 1),
        ("Fc1ccc(C(=O)NNC(=O)c2ccc(cc2)n3cnnn3)c(F)c1", 0, 0, 0),
        ("COc1cc(Cl)ccc1N=C(S)N(CCN2CCOCC2)C3CCN(CC3)C(=O)C", 0, 0, 0),
        ("CCCCCCCCCCCCCNC(=O)[C@H](CO)\N=C\c1ccccc1", 0, 0, 2),
    ]

    def test_alert_filter(self):
        ok_mols = lead.alert_filter(
            self.data.smiles.values, alerts=["Glaxo", "BMS"], n_jobs=2, return_idx=False
        )
        ok_index_2 = lead.alert_filter(
            self.data.smiles.values, alerts=["Glaxo", "bms"], n_jobs=2, return_idx=False
        )
        self.assertEquals(sum(ok_mols) == len(ok_index_2))
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
                for x in data.iloc[ok_index_3].smiles.apply(dm.to_mol).apply(MolWt)
            )
        )
        al_filter = lead.AlertFilters(alerts_set=["BMS"])
        df = al_filter(self.data.smiles.values[:10])
        expected_cols = set(
            ["_smiles", "status", "reasons", "MW", "LogP", "HBD", "HBA", "TPSA"]
        )
        self.assertTrue(
            len(expected_cols.intersection(set(df.columns)) == len(expected_cols))
        )

    def test_common_filter(self):
        idx = lead.common_filter(
            self.data.smiles.values,
            pains=True,
            brenk=False,
            nih=False,
            zinc=False,
            return_idx=False,
        )
        self.assertTrue(sum(idx) == 632)
        idx = lead.common_filter(
            self.data.smiles.values,
            pains=False,
            brenk=True,
            nih=False,
            zinc=False,
            return_idx=True,
        )
        # brenk should filter alkyl halide
        self.assertTrue(640 not in idx)

    def test_screening_filter(self):
        df = pd.DataFrame.from_records(
            screening_smiles_set, columns=["smiles", "covalent", "special", "severity"]
        )
        idx = lead.screening_filter(
            df.smiles, n_jobs=2, max_severity=2, return_idx=True
        )
        expected_idx = list(df[df.severity < 5].index)
        self.assertListEqual(len(set(idx).difference(set(expected_idx))), 0)

        nov_filters = lead.NovartisFilters()


if __name__ == "__main__":
    ut.main()
