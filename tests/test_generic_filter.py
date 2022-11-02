import unittest as ut
import numpy as np
import datamol as dm

from medchem.filter import generic


class Test_GenericFilter(ut.TestCase):

    data = dm.data.freesolv()
    smiles_with_issues = [
        # none
        "Nc1noc2cc(-c3noc(C(F)(F)F)n3)ccc12",
        # hologenicity, atom_list_filter, ring_infraction
        "Nc1noc2cc(-c3noc(n3)C(Cl)(Cl)Cl)c(Cl)c(B3C=C3)c12",
        # num_atom_filter
        "C[C@H](Nc1ncc(-c2noc(C(F)(F)F)n2)cc1Cl)c1ccncc1",
        # macrocycle, num_atom_filter
        "NC1=NOC2=CC(C3=NOC(=N3)C(F)(F)F)=C3CCCCCCCC3=C12",
        # ring infraction, macrocycle at 10, num_atom_filter at 20
        "Nc1noc2c(cc(cc12)C1NCO1)C1CCCCCCCCC1",
    ]

    def test_num_atom_filter(self):
        mols = [dm.to_mol(x) for x in self.data.smiles.values]
        passing_idx = [i for i, m in enumerate(mols) if m.GetNumAtoms() < 20]
        idx = generic.num_atom_filter(
            self.data.smiles.values, max_atoms=20, return_idx=True
        )
        np.testing.assert_array_equal(idx, passing_idx)

    def test_ring_infraction_filter(self):
        out = generic.ring_infraction_filter(self.smiles_with_issues, return_idx=True)
        self.assertListEqual([0, 2, 3], out)

    def test_macrocycle_filter(self):
        out = generic.macrocycle_filter(
            self.smiles_with_issues, max_cycle_size=7, return_idx=False
        )
        out = list(out)
        self.assertListEqual(out, [True, True, True, False, False])

    def test_atom_list_filter(self):
        out = generic.atom_list_filter(
            self.smiles_with_issues, unwanted_atom_list=["B"], return_idx=True
        )
        out = list(out)
        self.assertListEqual([0, 2, 3, 4], out)

    def test_halogenicity_filter(self):
        out = generic.halogenicity_filter(self.smiles_with_issues, return_idx=True)
        out = list(out)
        self.assertListEqual([0, 2, 3, 4], out)

    def test_n_stereo_center(self):

        smiles = [
            "CC(C)(O)C(O)[C@H](O)C(O)CO",  # fail, too many undefined stereocenter
            "C1[C@H]([C@@H]([C@H]([C@H](O1)O)O)O)O",  # fail, too many stereocenter
            "C[C@@H]1CC[C@H]([C@@H](C1)O)C(C)C",  # pass
        ]

        out = generic.num_stereo_center_filter(smiles, return_idx=False)
        out = list(out)
        self.assertListEqual([False, False, True], out)

    def test_symmetry_filter(self):
        symmetric_mols = [
            "O=C(O)c1cc(-n2ccnc2)cc(-n2ccnc2)c1",
            "CC(C)(C)[C@@H]1COC(C2(C3=N[C@H](C(C)(C)C)CO3)Cc3ccccc3C2)=N1",
            "c1ccc2oc(-c3ccc(-c4nc5ccccc5o4)s3)nc2c1",
            "Cc1cc(O)c(C(c2ccc(Cl)cc2)c2c(O)cc(C)[nH]c2=O)c(=O)[nH]1",
        ]

        random_mols = [
            "CCn1cc(S(=O)(=O)n2cc(Cl)cn2)cn1",
            "Cc1cc(C)c(S(=O)(=O)Nc2cc(C)ccc2C)c(C)c1",
            "c1ccc(CCC2CCN(CCC3COCCO3)CC2)cc1",
            "CCCC1CCC([C@H]2CC[C@H](C(=O)O)CC2)CC1",
        ]
        symmetric = generic.symmetry_filter(symmetric_mols, symmetry_threshold=0.8)
        np.testing.assert_array_equal(symmetric, [False] * len(symmetric))

        not_symmetric = generic.symmetry_filter(random_mols, symmetry_threshold=0.8)
        np.testing.assert_array_equal(not_symmetric, [True] * len(not_symmetric))


if __name__ == "__main__":
    ut.main()
