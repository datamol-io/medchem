import unittest as ut
import datamol as dm
from medchem.complexity import _complexity_calc as calc
from medchem.complexity import ComplexityFilter
from medchem.filter.lead import complexity_filter


class Test_Complexity(ut.TestCase):
    def test_whitlockct(self):
        """Test WhitlockCT"""

        smiles = [
            "CC(=O)O[C@@H]1C2=C(C)[C@H](C[C@@](O)([C@@H](OC(=O)C)[C@@H]3[C@@]4(CO[C@@H]4C[C@H](O)[C@@]3(C)C1=O)OC(C)=O)C2(C)C)OC(=O)[C@H](O)[C@@H](NC(=O)C)C",
            "C1(C)=CC(C)(C)C2(CCC3)CCC(C)C123",
            "O=C1C=CC(=O)C=C1",
        ]

        expected_outputs = [67, 20, 14]
        mols = [dm.to_mol(x) for x in smiles]
        self.assertListEqual([calc.WhitlockCT(m) for m in mols], expected_outputs)

    def test_baronect(self):
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
        outputs = [calc.BaroneCT(mol) for mol in mols_non_chiral] + [
            calc.BaroneCT(mol, chiral=True) for mol in mols_chiral
        ]
        self.assertListEqual(outputs, expected_outputs)

    def test_smcm(self):
        """Tsest SMCM"""
        mols = [
            dm.to_mol(x)
            for x in [
                "CC1=C[C@@H]2[C@]([C@@H](C1=O)O)([C@]3([C@H]([C@H]([C@H]([C@@]34CO4)O2)O)OC(=O)C)C)CO",
                "CN(C)CCOCCOCCN(C)C",
                "OC(=O)CC(C)(O)CC(O)=O",
            ]
        ]
        expected_outputs = [59.376, 20.034, 20.485]
        outputs = [calc.SMCM(mol) for mol in mols]
        self.assertListEqual(outputs, expected_outputs)

    def test_twc(self):
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
        outputs = [calc.TWC(mol, log10=False) for mol in mols]
        self.assertListEqual(outputs, expected_outputs)

    def test_complexity_filter(self):
        cf = ComplexityFilter()
        smiles = dm.data.cdk2().smiles.iloc[:10].values
        expected_outputs = [
            False,
            True,
            False,
            False,
            False,
            False,
            True,
            True,
            True,
            False,
        ]
        outputs = [cf(dm.to_mol(x)) for x in smiles]
        output_filters = complexity_filter(smiles, complexity_threshold="zinc_12")
        self.assertListEqual(outputs, expected_outputs)

    def test_available_list_cf(self):
        listed = set(ComplexityFilter.list_default_available_filters())

        self.assertSetEqual(listed, set(ComplexityFilter.COMPLEXITY_FNS.keys()))
        with self.assertRaises(ValueError):
            cf = ComplexityFilter(limit="fake")
        with self.assertRaises(ValueError):
            # sas is not defined for the default
            cf = ComplexityFilter(limit="99", complexity_metric="sas")


if __name__ == "__main__":
    ut.main()
