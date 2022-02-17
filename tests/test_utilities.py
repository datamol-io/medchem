import unittest as ut
import datamol as dm
import pandas as pd
import numpy as np
import datamol as dm
from medchem.filter import lead
from medchem.catalog import merge_catalogs, NamedCatalogs


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


if __name__ == "__main__":
    ut.main()
