from typing import List
from typing import Union
from typing import Optional
import copy
import functools
import pandas as pd
import numpy as np
from rdkit.Chem import FilterCatalog
from medchem.utils import get_data


def merge(*catalogs):
    """Merge several catalogs into a single one

    Returns:
        catalog (FilterCatalog): merged catalog
    """
    params = FilterCatalog.FilterCatalogParams()
    for catlg in catalogs:
        params.AddCatalog(catlg)
    return FilterCatalog.FilterCatalog(params)


def from_smarts(
    smarts: List[str],
    labels: Optional[List[str]] = None,
    mincounts: Optional[List[int]] = None,
    maxcounts: Optional[List[int]] = None,
    entry_as_inds: bool = False,
):
    """Load catalog from a list of smarts

    Args:
        smarts: list of input smarts to add to the catalog
        labels: list of label for each smarts
        mincounts: minimum count before a match is recognized
        maxcounts: maximum count for a match to be valid
        entry_as_inds: whether to use index for entry id or the label

    Returns:
        catalog (FilterCatalog): merged catalog
    """
    catalog = FilterCatalog.FilterCatalog()
    if labels is None:
        labels = smarts
    if mincounts is None:
        mincounts = [1] * len(smarts)

    for i, (sm, lb, count) in enumerate(zip(smarts, labels, mincounts)):
        if maxcounts is None:
            fil = FilterCatalog.SmartsMatcher(lb, sm, count)
        else:
            fil = FilterCatalog.SmartsMatcher(lb, sm, count, maxcounts[i])
        entry_name = lb
        if entry_as_inds:
            entry_name = i
        catalog.AddEntry(FilterCatalog.FilterCatalogEntry(entry_name, fil))
    return catalog


class NamedCatalogs:
    """
    Holder for substructure matching catalogs
    """

    @staticmethod
    @functools.lru_cache(maxsize=32)
    def tox(
        pains_a: bool = True,
        pains_b: bool = True,
        pains_c: bool = False,
        brenk: bool = True,
        nih: bool = False,
        zinc: bool = False,
    ):
        """Common toxicity and interference catalog

        Args:
            pains_a: whether to include PAINS filters from assay A
            pains_b: whether to include PAINS filters from assay B
            pains_c: whether to include PAINS filters from assay C
            brenk: whether to include BRENK filters
            nih: whether to include NIH filters
            zinc: whether to include ZINC filters
        """
        params = FilterCatalog.FilterCatalogParams()
        if pains_a:
            params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_A)
        if pains_b:
            params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_B)
        if pains_c:
            params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_C)
        if brenk:
            params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK)
        if nih:
            params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.NIH)
        if zinc:
            params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.ZINC)
        catalog = FilterCatalog.FilterCatalog(params)
        return catalog

    @staticmethod
    @functools.lru_cache(maxsize=32)
    def pains():
        return FilterCatalog.FilterCatalog(
            FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS
        )

    @staticmethod
    @functools.lru_cache(maxsize=32)
    def nih():
        return FilterCatalog.FilterCatalog(
            FilterCatalog.FilterCatalogParams.FilterCatalogs.NIH
        )

    @staticmethod
    @functools.lru_cache(maxsize=32)
    def zinc():
        return FilterCatalog.FilterCatalog(
            FilterCatalog.FilterCatalogParams.FilterCatalogs.ZINC
        )

    @staticmethod
    @functools.lru_cache(maxsize=32)
    def brenk():
        return FilterCatalog.FilterCatalog(
            FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK
        )

    @staticmethod
    @functools.lru_cache(maxsize=32)
    def alerts(subset: Optional[Union[List[str], str]] = None):
        """Alerts filter catalogs used in most pharma
        Args:
            subset: subset of providers to consider

        Returns:
            catalog (FilterCatalog): filter catalog
        """
        rd_filters = pd.read_csv(get_data("rd_alerts.csv"))
        if subset is not None:
            if isinstance(subset, str):
                subset = [subset]
            subset = [x.lower() for x in subset]
            rd_filters = rd_filters[rd_filters.rule_set_name.str.lower().isin(subset)]
        mincount = np.maximum(rd_filters["mincount"], 1).astype(int)
        labels = nibr_filters.apply(
            lambda x: "{0}||{1}_min({2})||{3}".format(
                x["rule_set_name"],
                x["description"],
                x["mincount"],
                x["priority"],
            ),
            axis=1,
        )
        return from_smarts(
            rd_filters["smarts"].values, labels, mincount, entry_as_inds=True
        )

    @staticmethod
    def dundee():
        return NamedCatalogs.alerts(subset=["Dundee"])

    @staticmethod
    def bms():
        return NamedCatalogs.alerts(subset=["BMS"])

    @staticmethod
    def glaxo():
        return NamedCatalogs.alerts(subset=["Glaxo"])

    @staticmethod
    def schembl():
        return NamedCatalogs.alerts(subset=["SureChEMBL"])

    @staticmethod
    def mlsmr():
        return NamedCatalogs.alerts(subset=["MLSMR"])

    @staticmethod
    def inpharmatica():
        return NamedCatalogs.alerts(subset=["Inpharmatica"])

    @staticmethod
    def lint():
        return NamedCatalogs.alerts(subset=["LINT"])

    @staticmethod
    @functools.lru_cache(maxsize=32)
    def nibr():
        """Catalog from NIBR"""
        nibr_filters = pd.read_csv(get_data("nibr.csv"))
        mincount = np.maximum(nibr_filters["mincount"], 1).astype(int)
        labels = nibr_filters.apply(
            lambda x: "{0}||{1}_min({2})||{3}||{4}||{5}".format(
                x["rule_set_name"],
                x["description"],
                x["mincount"],
                x["severity"],
                x["covalent"],
                x["special_mol"],
            ),
            axis=1,
        )
        return from_smarts(
            nibr_filters["smarts"].values, labels, mincount, entry_as_inds=False
        )
