from typing import List
from typing import Optional
from typing import Union

import os
import functools
import pandas as pd
import datamol as dm

from medchem.utils import get_data
from medchem import catalog


class ChemicalGroup:
    """Build a library of chemical groups using a list of structures parsed from a file

    The default library of structure has been curated from https://github.com/Sulstice/global-chem.

    !!! warning
        The SMILES and SMARTS used in the default list of substructures do not result in the same matches.
        Unless specified otherwise, the SMILES will be used in the matching done by this class,
        whereas due to RDKit's limitation, the SMARTS will be used in the matching done by the generated catalog.
        For more information see this discussion: https://github.com/valence-platform/medchem/pull/19,

    """

    def __init__(
        self,
        groups: Union[str, List[str]] = None,
        medchem_only: bool = True,
        n_jobs: Optional[int] = None,
        groups_db: Optional[os.PathLike] = None,
    ):
        """Build a chemical group library

        Args:
            groups: List of groups to use. Defaults to None where all functional groups are used
            medchem_only: Whether to filter out any subset that is not classified as `medicinal_chemistry. Defaults to True.
            n_jobs: Optional number of jobs to run in parallel for internally building the data. Defaults to None.
            groups_db: Path to a file containing the dump of the chemical groups. Defaults is internal dataset
        """

        if isinstance(groups, str):
            groups = [groups]
        if groups is None:
            groups = []
        self.groups = groups
        self.medchem_only = medchem_only
        self.n_jobs = n_jobs or 0
        if groups_db is None:
            groups_db = get_data("chemical_groups.csv")
        data = pd.read_csv(groups_db)
        self.data = data[data.hierarchy.str.contains("|".join(self.groups))]
        if self.medchem_only:
            self.data = self.data[self.data.hierarchy.str.contains("medicinal")]

        self._initialize_data()

    def _initialize_data(self):
        """Initialize the data by precomputing some features"""
        self.data["mol_smarts"] = dm.parallelized(
            dm.from_smarts,
            self.data["smarts"].values,
            n_jobs=self.n_jobs,
            progress=False,
        )
        self.data["mol"] = dm.parallelized(
            dm.to_mol, self.data["smiles"].values, n_jobs=self.n_jobs, progress=False
        )

    def __len__(self):
        return len(self.data)

    @property
    def iupac(self):
        """Get the IUPAC name of the chemical groups in this instance"""
        return self.data.iupac.tolist()

    @property
    def smiles(self):
        """Get the SMILES of the chemical groups in this instance"""
        return self.data.smiles.tolist()

    @property
    def mols(self):
        """Get the Molecule object of the SMILES for the chemical groups in this instance"""
        return self.data.mol.tolist()

    @property
    def smarts(self):
        """Get the SMARTS of the chemical groups in this instance"""
        return self.data.smarts.tolist()

    @property
    def mol_smarts(self):
        """Get the SMARTS of the chemical groups in this instance"""
        return self.data.mol_smarts.tolist()

    @property
    def dataframe(self):
        """Get the dataframe of the chemical groups"""
        return self.data

    @functools.lru_cache(maxsize=32)
    def get_catalog(self):
        """Build an rdkit catalog from the current chemical group data"""
        return catalog.from_smarts(
            self.mol_smarts,
            self.iupac,
            entry_as_inds=False,
        )

    @staticmethod
    def list_groups():
        """List all the chemical groups available"""
        return catalog.list_chemical_groups(hierachy=False)

    @staticmethod
    def list_hierarchy_groups():
        """List all the hierarchy in chemical groups available.
        To get the full hierarchy on each path, split by the `.` character.
        """
        return catalog.list_chemical_groups(hierachy=True)

    def get_matches(self, mol: Union[dm.Mol, str], use_smiles: bool = True):
        """Get all the functional groups in this instance that matches the input molecule

        Args:
            mol: input molecule
            use_smiles: whether to use the smiles representation of the catalog or the smarts
        """
        if isinstance(mol, str):
            mol = dm.to_mol(mol)
        if mol is None:
            return None

        def matcher(query):
            return mol.GetSubstructMatches(query) or None

        if use_smiles:
            matches = dm.parallelized(
                matcher, self.data.mol.values, n_jobs=self.n_jobs, progress=False
            )
        else:
            matches = dm.parallelized(
                matcher, self.data.mol_smarts.values, n_jobs=self.n_jobs, progress=False
            )
        out = self.data[["iupac", "smiles", "smarts", "group"]].copy()
        out["matches"] = matches
        out = out.dropna(subset=["matches"])
        return out

    def has_match(self, mol: Union[dm.Mol, str]):
        """Check whether the input molecule has any functional group in this instance

        Args:
            mol: input molecule
        """
        matches = self.get_matches(mol)
        return matches is not None and len(matches) > 0
