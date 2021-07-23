import unittest as ut
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
        passing_idx = [i for i, m in enumerate(mols) if m.GetNumAtoms() < 20]
        idx = generic.num_atom_filter(
            self.data.smiles.values, max_atoms=20, return_idx=True
        )
        self.assertListEqual(passing_idx, idx)

    def test_ring_infraction_filter(self):
        idx = generic.ring_infraction_filter(self.smiles_with_issues, return_idx=True)
        self.assertListEqual([0, 2, 3], idx)

    def test_macrocycle_filter(self):
        idx = generic.macrocycle_filter(
            self.smiles_with_issues, max_cycle_size=7, return_idx=False
        )
        self.assertListEqual([0, 1, 2], [False, False, False, True])

    def test_atom_list_filter(self):
        idx = generic.atom_list_filter(self.smiles_with_issues, ["B"], return_idx=True)
        self.assertListEqual([0, 2, 3], idx)

    def test_halogenicity_filter(self):
        idx, _ = generic.halogenicity_filter(self.smiles_with_issues, True)
        self.assertListEqual([0, 2, 3], idx)


if __name__ == "__main__":
    ut.main()
