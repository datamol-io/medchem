import unittest as ut
import datamol as dm
import datamol as dm
from medchem.catalog import merge_catalogs, NamedCatalogs
from medchem.utils.smarts import SMARTSUtils
from medchem.utils.matches import Constraints
from medchem.utils.graph import score_symmetry


class Test_Utils(ut.TestCase):
    data = dm.data.freesolv()

    def test_catalog_merge(self):
        mols = self.data["smiles"].apply(dm.to_mol).values
        tox = NamedCatalogs.tox(pains_a=False, pains_b=False, pains_c=False, nih=True, brenk=True)
        nih = NamedCatalogs.nih()
        brenk = NamedCatalogs.brenk()
        merged_tox = merge_catalogs(nih, brenk)
        toxic1 = [tox.HasMatch(m) for m in mols]
        toxic2 = [merged_tox.HasMatch(m) for m in mols]
        self.assertEqual(toxic1, toxic2)


class TestGraph(ut.TestCase):
    test_mols = [
        ("O=C(O)c1cc(-n2ccnc2)cc(-n2ccnc2)c1", False),  # high symmetry ~0.857
        # fully symmetric
        ("OC(O)c1cc(-n2ccnc2)cc(-n2ccnc2)c1", True),
        (  # symmetry should be 1, modif of above
            "CC(C)(C)[C@@H]1COC(C2(C3=N[C@H](C(C)(C)C)CO3)Cc3ccccc3C2)=N1",
            True,
        ),
        ("c1ccc2oc(-c3ccc(-c4nc5ccccc5o4)s3)nc2c1", True),
        # symmetrical star graph
        ("CC(C)(C)", True),
        ("NC(C)(N)C", True),
        # non symmetrical
        ("NC(C)(C)C", False),
        ("CCCC1=CNC(CC)=C1C", False),  # modification of the above  # not symmetrical
        # disconnection does not prevent symmetric detection
        ("Cl.[O-]C(=O)CC([O-])=O.CCCC1=CNC=C1CCC", True),
        ("Cl.CC([O-])=O.CCCC1=CNC=C1CCC", True),  # regular salt, should be symmetrical
        (
            "[O-]C(=O)C1=CC=CC=C1.CCCC1=CNC=C1CCC",
            False,
        ),  # not a salt, should not be symmetrical
    ]

    def test_symmetry(self):
        mols, expected_vals = zip(*self.test_mols)
        scores = [score_symmetry(x, exclude_self_mapped_edged=False) for x in mols]
        fully_symmetrical = [x == 1 for x in scores]
        self.assertListEqual(list(expected_vals), fully_symmetrical)

        # if we don't exclude self_mapped edges, then we will have non symmetrical here
        expected_mol2 = score_symmetry(mols[1], exclude_self_mapped_edged=True)
        self.assertLess(expected_mol2, 1)
        self.assertGreater(expected_mol2, 0.8)


class TestConstraints(ut.TestCase):
    def test_constraints(self):
        def my_constraints(mol):
            # we want to either (have phenol) OR (have less than 7 atoms and not ring)
            return mol.HasSubstructMatch(dm.to_mol("Oc1ccccc1")) or (
                mol.GetNumAtoms() < 7 and dm.descriptors.n_rings(mol) < 1
            )

        smiles = [
            "CN(C)C(=O)c1cncc(C)c1",  # match, n_atoms < 7 and no ring
            "Cc1cncc(CC2CCCCC2)c1",  # not match, n_atoms < 7 but ring
            "Cc1cncc(c1)-c1ccc(O)cc1",  # match phenol
            "Cc1cncc(c1)-c1cccc2nc[nH]c12",  # no match n_atoms >= 7
        ]
        expected_results = []
        mols = [dm.to_mol(x) for x in smiles]
        core = dm.from_smarts("[C;H3]c1cncc([*:1])c1")

        # now let's set the constraints query
        for atom in core.GetAtoms():
            if (
                atom.GetAtomMapNum() == 1
            ):  # we add a recursive query to check again on any match that starts with this atom position
                atom.SetProp("query", "my_constraints")
        constraint_fns = dict(my_constraints=my_constraints)
        constraint = Constraints(core, constraint_fns)
        matches = [constraint(mol) for mol in mols]
        expected_results = [True, False, True, False]
        self.assertListEqual(expected_results, matches)


class Test_SMARTSUtils(ut.TestCase):
    def test_standard_attachment(self):
        smiles = [
            "C1CC(C)C(O)CC1[1*]",  # ortho
            "c1c(C)cc(O)cc1[*:2]",  # meta + aromatic
            "c1(C)ccc(O)cc1[*]",  # para + aromat
        ]

        out = [SMARTSUtils.standardize_attachment(sm, "[*:1]") for sm in smiles]
        self.assertTrue(all("[*:1]" in x for x in out))

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
        ortho_aro_query = dm.from_smarts(SMARTSUtils.ortho(sm1, sm2, aromatic_only=True))
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

        chain1 = dm.from_smarts(SMARTSUtils.aliphatic_chain(min_size=5, unsaturated_bondtype=None))
        chain2 = dm.from_smarts(
            SMARTSUtils.aliphatic_chain(min_size=6, unbranched=True, unsaturated_bondtype=dm.SINGLE_BOND)
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
        atom_id = tuple([atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomMapNum() == 1])
        query1 = SMARTSUtils.atom_in_env("[#6;r6]", "[*][OD2][C&D1]", "[c]aa[F]", union=False)
        query2 = SMARTSUtils.atom_in_env("[#6;r6][OD2][C&D1]", "[c]aa[F]", union=False)
        res1 = mol.GetSubstructMatch(dm.from_smarts(query1))
        res2 = mol.GetSubstructMatch(dm.from_smarts(query2))
        self.assertEqual(res1, res2)
        self.assertEqual(res1, atom_id)


if __name__ == "__main__":
    ut.main()
