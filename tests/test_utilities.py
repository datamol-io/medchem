import unittest as ut
import datamol as dm
import datamol as dm
from medchem.filter import lead
from medchem.catalog import merge_catalogs, NamedCatalogs
from medchem.utils.smarts import SMARTSUtils


class Test_Utils(ut.TestCase):

    data = dm.data.freesolv()

    def test_catalog_merge(self):
        mols = self.data["smiles"].apply(dm.to_mol).values
        tox = NamedCatalogs.tox(
            pains_a=False, pains_b=False, pains_c=False, nih=True, brenk=True
        )
        nih = NamedCatalogs.nih()
        brenk = NamedCatalogs.brenk()
        merged_tox = merge_catalogs(nih, brenk)
        toxic1 = [tox.HasMatch(m) for m in mols]
        toxic2 = [merged_tox.HasMatch(m) for m in mols]
        self.assertEqual(toxic1, toxic2)


class Test_SMARTSUtils(ut.TestCase):
    def test_ortho_meta_para(self):
        smiles = [
            "C1CC(C)C(O)CC1",  # ortho
            "c1c(C)cc(O)cc1",  # meta + aromatic
            "c1(C)ccc(O)cc1",  # para + aromat
        ]
        mols = [dm.to_mol(x) for x in smiles]
        sm1 = "[#6;!R]"
        sm2 = "[#8]"
        ortho_query = dm.from_smarts(SMARTSUtils.ortho(sm1, sm2))
        ortho_aro_query = dm.from_smarts(
            SMARTSUtils.ortho(sm1, sm2, aromatic_only=True)
        )
        meta_query = dm.from_smarts(SMARTSUtils.meta(sm1, sm2))
        para_query = dm.from_smarts(SMARTSUtils.para(sm1, sm2))
        expected_output = [True, True, True, False, False]
        output = [
            mols[0].HasSubstructMatch(ortho_query),  # True
            mols[1].HasSubstructMatch(meta_query),  # True
            mols[2].HasSubstructMatch(para_query),  # True
            mols[0].HasSubstructMatch(para_query),  # False
            mols[0].HasSubstructMatch(ortho_aro_query),  # False
        ]

        self.assertEqual(output, expected_output)

    def test_long_chain(self):
        smiles = [
            "CCCCCCCC",  # octane
            "CCCC=CCCN",  #  hept-3-en-1-amine
            "CCOCCC",  #  1-ethoxypropane
            "CC(C)(C)CN",  #  2,2-dimethylpropan-1-amine
            "CCC(C)(C)CN",  # 2,2-dimethylbutan-1-amine
            "C1CCCC1CC",  # ethylcyclopentane
        ]
        mols = [dm.to_mol(x) for x in smiles]

        chain1 = dm.from_smarts(
            SMARTSUtils.aliphatic_chain(min_size=5, unsaturated_bondtype=None)
        )
        chain2 = dm.from_smarts(
            SMARTSUtils.aliphatic_chain(
                min_size=6, unbranched=True, unsaturated_bondtype=dm.SINGLE_BOND
            )
        )
        chain3 = dm.from_smarts(
            SMARTSUtils.aliphatic_chain(
                min_size=6,
                unbranched=True,
                unsaturated_bondtype=dm.DOUBLE_BOND,
                allow_hetero_atoms=False,
            )
        )
        chain1_out = [m.HasSubstructMatch(chain1) for m in mols]
        chain2_out = [m.HasSubstructMatch(chain2) for m in mols]
        chain3_out = [m.HasSubstructMatch(chain3) for m in mols]
        self.assertEqual(chain1_out, [True, True, True, False, True, False])
        self.assertEqual(chain2_out, [True, False, True, False, False, False])
        self.assertEqual(chain3_out, [True, True, False, False, False, False])

    def test_smarts_fragment(self):
        smiles = ["COCC(O)=O", "OCCO", "O.CCO.O"]
        mols = [dm.to_mol(x) for x in smiles]
        query1 = SMARTSUtils.same_fragment("[8]", "[8]", "[8]")
        query2 = SMARTSUtils.different_fragment("[8]", "[8]", "[8]")
        self.assertEqual(query1, "([8]).([8]).([8])")
        self.assertEqual(query2, "([8].[8].[8])")
        # query1_out = [m.HasSubstructMatch(dm.from_smarts(query1)) for m in mols]
        # query2_out = [m.HasSubstructMatch(dm.from_smarts(query2)) for m in mols]
        # self.assertEqual(query1_out, [True, False, False])
        # self.assertEqual(query2_out, [False, False, True])

    def test_atom_in_env(self):
        mol = dm.to_mol("c1[c:1](OC)cc(F)cc1")
        # we are trying to match the carbon with atom map 1
        # these queries should be equivalent
        atom_id = tuple(
            [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomMapNum() == 1]
        )
        query1 = SMARTSUtils.atom_in_env(
            "[#6;r6]", "[*][OD2][C&D1]", "[c]aa[F]", union=False
        )
        query2 = SMARTSUtils.atom_in_env("[#6;r6][OD2][C&D1]", "[c]aa[F]", union=False)
        res1 = mol.GetSubstructMatch(dm.from_smarts(query1))
        res2 = mol.GetSubstructMatch(dm.from_smarts(query2))
        self.assertEqual(res1, res2)
        self.assertEqual(res1, atom_id)


if __name__ == "__main__":
    ut.main()
