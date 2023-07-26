from typing import List
from typing import Optional
from typing import Union

import os
import functools
import pandas as pd
import datamol as dm

from rdkit.Chem import GetMolFrags
from rdkit.Chem.rdmolops import ReplaceCore
from rdkit.Chem.rdmolops import AdjustQueryParameters
from rdkit.Chem.rdmolops import AdjustQueryProperties
from rdkit.Chem.rdmolops import AdjustQueryWhichFlags

from medchem.catalogs import catalog_from_smarts
from medchem.utils.loader import get_data_path


def list_default_chemical_groups(hierarchy: bool = False) -> list:
    """List all the available chemical groups.

    !!! note
        chemical groups defines how a collection of patterns are organized.
        They **do not** correspond to individual pattern name.

    Args:
        hierarchy: whether to return the full hierarchy or the group name only

    Returns:
        List of chemical groups
    """
    data = pd.read_csv(get_data_path("chemical_groups.csv"))
    if hierarchy:
        return list(data.hierarchy.unique())
    return list(data.group.unique())


def list_functional_group_names(unique: bool = True) -> list:
    """
    List common functional group names

    Args:
        unique: whether to return only unique names

    Returns:
        List of functional group names
    """
    data = pd.read_csv(get_data_path("chemical_groups.csv"))
    data = data[data.hierarchy.str.contains("functional_groups")]
    return list(data.name.unique())


@functools.lru_cache(maxsize=None)
def get_functional_group_map() -> dict:
    """
    Map functional groups to their corresponding SMARTS string.

    Returns:
        List of functional group names
    """
    data = pd.read_csv(get_data_path("chemical_groups.csv"))
    data = data[data.hierarchy.str.contains("functional_groups")]
    # EN: any group that is not unique should be dropped
    # this is because of the `basic_groups` hierarchy that is loosely defined
    data = data.drop_duplicates(subset=["name"], keep=False)
    data = data.sort_values("name", ascending=True)
    return dict(zip(data["name"], data["smarts"]))


class ChemicalGroup:
    """Build a library of chemical groups using a list of structures parsed from a file

    The default library of structure has been curated from https://github.com/Sulstice/global-chem and additional open source data.

    !!! note
        For new chemical groups, please minimally provide the 'smiles'/'smarts', 'name' and "group" and optional 'hierarchy' columns

    !!! warning
        The SMILES and SMARTS used in the default list of substructures do not result in the same matches.
        Unless specified otherwise, the SMILES will be used in the matching done by this class,
        whereas due to RDKit's limitation, the SMARTS will be used in the matching done by the generated catalog.

    """

    def __init__(
        self,
        groups: Optional[Union[str, List[str]]] = None,
        n_jobs: Optional[int] = None,
        groups_db: Optional[Union[os.PathLike, str]] = None,
    ):
        """Build a chemical group library

        Args:
            groups: List of groups to use. Defaults to None where all functional groups are used
            n_jobs: Optional number of jobs to run in parallel for internally building the data. Defaults to None.
            groups_db: Path to a file containing the dump of the chemical groups. Default is internal dataset
        """

        if isinstance(groups, str):
            groups = [groups]
        if groups is None:
            groups = []
        self.groups = groups
        self.n_jobs = n_jobs or 0
        if groups_db is None:
            groups_db = get_data_path("chemical_groups.csv")
        self.data = pd.read_csv(groups_db)
        if "hierarchy" not in self.data.columns:
            self.data["hierarchy"] = self.data["group"]
        if self.groups:
            self.data = self.data[self.data.hierarchy.str.contains("|".join(self.groups))]
        # EN: fill smiles and smarts with empty string when they are missing
        # this prevent error when applying dm.to_mol, dm.from_smarts and building the FilterCatalog
        self.data["smiles"] = self.data["smiles"].fillna("")
        self.data["smarts"] = self.data["smarts"].fillna("")
        self._initialize_data()

    def filter(self, names: List[str], fuzzy: bool = False):
        """Filter the group to restrict to only the name in input

        Args:
            names: list of names to use for filters
            fuzzy: whether to use exact or fuzzy matching
        """
        if names is None or len(names) == 0:
            return self

        if fuzzy:
            names = [f".*{name}.*" for name in names]
            self.data = self.data[self.data["name"].str.match(r"|".join(names))]
        else:
            self.data = self.data[self.data["name"].isin(names)]

        return self

    def _initialize_data(self):
        """Initialize the data by precomputing some features"""
        with dm.without_rdkit_log():
            self.data["mol_smarts"] = dm.parallelized(
                dm.from_smarts,
                self.data["smarts"].values,
                n_jobs=self.n_jobs,
                progress=False,
            )
            self.data["mol"] = dm.parallelized(
                dm.to_mol,
                self.data["smiles"].values,
                n_jobs=self.n_jobs,
                progress=False,
            )

    def __len__(self):
        return len(self.data)

    @property
    def name(self):
        """Get the Name of the chemical groups in this instance"""
        return self.data.name.tolist()

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
    def mol_adjusted(self):
        """Get the Molecules object of the SMILES, adjusted for stricter match for the chemical groups in this instance"""
        # EN: beware of pickling issues with AdjustQueryParameters
        query_adjust = AdjustQueryParameters()
        query_adjust.adjustRingChain = True
        query_adjust.adjustHeavyDegree = True
        query_adjust.adjustDegree = True
        query_adjust.adjustDegreeFlags = (
            AdjustQueryWhichFlags.ADJUST_IGNORERINGS | AdjustQueryWhichFlags.ADJUST_IGNORECHAINS
        )
        mols = self.data.mol.apply(lambda x: AdjustQueryProperties(x, query_adjust))
        return mols.tolist()

    @property
    def dataframe(self):
        """Get the dataframe of the chemical groups"""
        return self.data

    # @functools.lru_cache(maxsize=32)
    def get_catalog(self, exact_match: bool = True):
        """Build an rdkit catalog from the current chemical group data

        Args:
            exact_match: whether to adjust the queries for a more stringent match
        """

        if exact_match:
            return catalog_from_smarts(
                self.mol_adjusted,
                self.name,
                entry_as_inds=False,
            )
        return catalog_from_smarts(
            self.mol_smarts,
            self.name,
            entry_as_inds=False,
        )

    def list_groups(self):
        """List all the chemical groups available"""
        return list(self.data.group.unique())

    def list_hierarchy_groups(self):
        """List all the hierarchy in chemical groups available.
        To get the full hierarchy on each path, split by the `.` character.
        """
        return list(self.data.hierarchy.unique())

    def get_matches(
        self,
        mol: Union[dm.Mol, str],
        use_smiles: bool = True,
        exact_match: bool = False,
        terminal_only: bool = False,
    ):
        """Get all the functional groups in this instance that matches the input molecule

        Args:
            mol: input molecule
            use_smiles: whether to use the smiles representation of the catalog or the smarts
            exact_match: whether to use exact matching by adjusting the query
            terminal_only: ensure whether the matches to the functional group are terminal,
                meaning that any subgraph matching should not be in the middle of the molecules.
        """
        if isinstance(mol, str):
            mol = dm.to_mol(mol)

        if mol is None:
            return None

        # EN: beware of pickling issues with AdjustQueryParameters
        query_adjust = AdjustQueryParameters()
        query_adjust.adjustRingChain = True
        query_adjust.adjustDegree = True

        def matcher(query):
            if exact_match:
                query = AdjustQueryProperties(query, query_adjust)
            matches = mol.GetSubstructMatches(query, uniquify=True) or None
            if matches is not None and terminal_only:
                side_chains = [ReplaceCore(mol, query, match) for match in matches]
                side_chains = [GetMolFrags(x, asMols=True) for x in side_chains]
                matches = [x for x, y in zip(matches, side_chains) if len(y) == 1]
            return matches or None

        if use_smiles:
            matches = dm.parallelized(matcher, self.data.mol.values, n_jobs=self.n_jobs, progress=False)
        else:
            matches = dm.parallelized(
                matcher, self.data.mol_smarts.values, n_jobs=self.n_jobs, progress=False
            )
        out = self.data[["name", "smiles", "smarts", "group"]].copy()
        out["matches"] = matches
        out = out.dropna(subset=["matches"])
        return out

    def has_match(self, mol: Union[dm.Mol, str], exact_match: bool = False, terminal_only: bool = False):
        """Check whether the input molecule has any functional group in this instance

        Args:
            mol: input molecule
            exact_match: whether to use exact matching by adjusting the query
            terminal_only: ensure the matches to the functional group are terminal
        """
        matches = self.get_matches(mol, exact_match=exact_match, terminal_only=terminal_only)
        return matches is not None and len(matches) > 0
