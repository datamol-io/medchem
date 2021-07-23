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
    ]

    def test_num_atom_filter(self):
        mols = [dm.to_mol(x) for x in self.data.smiles.values]
        passing_idx = [i for i, m in enumerate(mols) if m.GetNumAtoms() <= 20]
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
        self.assertListEqual(out, [True, True, True, False])

    def test_atom_list_filter(self):
        out = generic.atom_list_filter(
            self.smiles_with_issues, unwanted_atom_list=["B"], return_idx=True
        )
        self.assertListEqual([0, 2, 3], out)

    def test_halogenicity_filter(self):
        out = generic.halogenicity_filter(self.smiles_with_issues, return_idx=True)
        self.assertListEqual([0, 2, 3], out)


if __name__ == "__main__":
    ut.main()
