from typing import List
from typing import Union
from typing import Optional

import warnings
import functools

from loguru import logger

import pandas as pd
import numpy as np
import datamol as dm

with warnings.catch_warnings():

    # Remove annoying `RuntimeWarnings`:
    # `<frozen importlib._bootstrap>:219: RuntimeWarning: to-Python converter for boost::shared_ptr<RDKit::FilterCatalogEntry const> already registered; second conversion method ignored.`
    warnings.simplefilter("ignore")

    from rdkit.Chem import FilterCatalog

from medchem.utils import get_data


def list_named_catalogs():
    """
    List all available named catalogs. This list will ignore all chemical groups
    For a list of chemical group to be queried using NamedCatalog.chemical_groups, use `medchem.group.list_default_chemical_groups`
    """
    return [
        x
        for x in NamedCatalogs.__dict__.keys()
        if (not x.startswith("_") and x not in ["alerts", "chemical_groups"])
    ]


def merge_catalogs(*catalogs):
    """Merge several catalogs into a single one

    Returns:
        catalog (FilterCatalog): merged catalog
    """
    if len(catalogs) == 1:
        return catalogs[0]
    params = FilterCatalog.FilterCatalogParams()
    missing_catalogs = []
    for catlg in catalogs:
        if isinstance(catlg, FilterCatalog.FilterCatalogParams.FilterCatalogs):
            params.AddCatalog(catlg)
        else:
            missing_catalogs.append(catlg)
    parameterized_catalogs = FilterCatalog.FilterCatalog(params)
    for catlg in missing_catalogs:
        for entry_nums in range(catlg.GetNumEntries()):
            entry = catlg.GetEntryWithIdx(entry_nums)
            parameterized_catalogs.AddEntry(entry)
    return parameterized_catalogs


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
        catalog (FilterCatalog): merged catalogs
    """

    with dm.without_rdkit_log():
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
            entry_name = str(lb)
            if entry_as_inds:
                entry_name = str(i)

            if fil.IsValid():
                catalog.AddEntry(FilterCatalog.FilterCatalogEntry(entry_name, fil))
            else:
                logger.warning(f"It seem like the SMARTS {sm} is invalid")
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
        catalogs = []
        if pains_a:
            catalogs.append(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_A)
        if pains_b:
            catalogs.append(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_B)
        if pains_c:
            catalogs.append(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_C)
        if brenk:
            catalogs.append(FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK)
        if nih:
            catalogs.append(FilterCatalog.FilterCatalogParams.FilterCatalogs.NIH)
        if zinc:
            catalogs.append(FilterCatalog.FilterCatalogParams.FilterCatalogs.ZINC)
        catalog = merge_catalogs(*catalogs)
        return catalog

    @staticmethod
    @functools.lru_cache(maxsize=32)
    def pains():
        return FilterCatalog.FilterCatalog(
            FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS
        )

    @staticmethod
    @functools.lru_cache(maxsize=32)
    def pains_a():
        return FilterCatalog.FilterCatalog(
            FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_A
        )

    @staticmethod
    @functools.lru_cache(maxsize=32)
    def pains_b():
        return FilterCatalog.FilterCatalog(
            FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_B
        )

    @staticmethod
    @functools.lru_cache(maxsize=32)
    def pains_c():
        return FilterCatalog.FilterCatalog(
            FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_C
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
    def alerts(subset: Optional[Union[List[str], str]] = None):
        """Alerts filter catalogs commonly used in molecule filtering

        Args:
            subset: subset of providers to consider

        Returns:
            catalog (FilterCatalog): filter catalog
        """
        rd_filters = pd.read_csv(get_data("common_alert_collection.csv"))
        if subset is not None:
            if isinstance(subset, str):
                subset = [subset]
            subset = [x.lower() for x in subset]
            rd_filters = rd_filters[rd_filters.rule_set_name.str.lower().isin(subset)]
        mincount = np.maximum(rd_filters["mincount"], 1).astype(int)

        labels = rd_filters.apply(
            lambda x: "{0}||{1}_min({2})||{3}".format(
                x["rule_set_name"],
                x["description"],
                x["mincount"],
                x["priority"],
            ),
            axis=1,
        )
        return from_smarts(
            rd_filters["smarts"].values,
            labels,
            mincount,
            entry_as_inds=False,
        )

    @staticmethod
    @functools.lru_cache(maxsize=32)
    def dundee():
        return NamedCatalogs.alerts(subset=["Dundee"])

    @staticmethod
    @functools.lru_cache(maxsize=32)
    def bms():
        return NamedCatalogs.alerts(subset=["BMS"])

    @staticmethod
    @functools.lru_cache(maxsize=32)
    def glaxo():
        return NamedCatalogs.alerts(subset=["Glaxo"])

    @staticmethod
    @functools.lru_cache(maxsize=32)
    def schembl():
        return NamedCatalogs.alerts(subset=["SureChEMBL"])

    @staticmethod
    @functools.lru_cache(maxsize=32)
    def mlsmr():
        return NamedCatalogs.alerts(subset=["MLSMR"])

    @staticmethod
    @functools.lru_cache(maxsize=32)
    def inpharmatica():
        return NamedCatalogs.alerts(subset=["Inpharmatica"])

    @staticmethod
    @functools.lru_cache(maxsize=32)
    def lint():
        return NamedCatalogs.alerts(subset=["LINT"])

    @staticmethod
    @functools.lru_cache(maxsize=32)
    def alarm_nmr():
        return NamedCatalogs.alerts(subset=["Alarm-NMR"])

    @staticmethod
    @functools.lru_cache(maxsize=32)
    def alphascreen():
        return NamedCatalogs.alerts(subset=["AlphaScreen-Hitters"])

    @staticmethod
    @functools.lru_cache(maxsize=32)
    def gst_hitters():
        return NamedCatalogs.alerts(subset=["GST-Hitters"])

    @staticmethod
    @functools.lru_cache(maxsize=32)
    def his_hitters():
        return NamedCatalogs.alerts(subset=["HIS-Hitters"])

    @staticmethod
    @functools.lru_cache(maxsize=32)
    def luciferase():
        return NamedCatalogs.alerts(subset=["LuciferaseInhibitor"])

    @staticmethod
    @functools.lru_cache(maxsize=32)
    def dnabinder():
        return NamedCatalogs.alerts(subset=["DNABinder"])

    @staticmethod
    @functools.lru_cache(maxsize=32)
    def chelator():
        return NamedCatalogs.alerts(subset=["Chelator"])

    @staticmethod
    @functools.lru_cache(maxsize=32)
    def hitters():
        return NamedCatalogs.alerts(subset=["Frequent-Hitter"])

    @staticmethod
    @functools.lru_cache(maxsize=32)
    def electrophilic():
        return NamedCatalogs.alerts(subset=["Electrophilic"])

    @staticmethod
    @functools.lru_cache(maxsize=32)
    def carcinogen(include_non_genotoxic: bool = True):
        catalogs = ["Genotoxic-Carcinogenicity"]
        if include_non_genotoxic:
            catalogs += ["Non-Genotoxic-Carcinogenicity"]
        return NamedCatalogs.alerts(subset=catalogs)

    @staticmethod
    @functools.lru_cache(maxsize=32)
    def ld50_oral():
        return NamedCatalogs.alerts(subset=["LD50-Oral"])

    @staticmethod
    @functools.lru_cache(maxsize=32)
    def reactive_unstable_toxic():
        return NamedCatalogs.alerts(subset=["Reactive-Unstable-Toxic"])

    @staticmethod
    @functools.lru_cache(maxsize=32)
    def skin():
        return NamedCatalogs.alerts(subset=["Skin"])

    @staticmethod
    @functools.lru_cache(maxsize=32)
    def toxicophore():
        return NamedCatalogs.alerts(subset=["Toxicophore"])

    @staticmethod
    @functools.lru_cache(maxsize=32)
    def nibr():
        """Catalog from NIBR

        !!! warning
            This includes all the compounds in the catalog, regardless of severity (FLAG, EXCLUDE, ANNOTATION)
            You likely don't want to use this for blind prioritization
        """
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

    @staticmethod
    @functools.lru_cache(maxsize=32)
    def bredt():
        """Bredt fitler rules
        Also see example of usage by surge's
        https://github.com/StructureGenerator/SURGE/blob/main/doc/surge1_0.pdf

        """
        bredt_df = pd.read_csv(get_data("bredt.csv"))
        return from_smarts(
            bredt_df["smarts"].values,
            bredt_df["labels"].values,
            entry_as_inds=True,
        )

    @staticmethod
    @functools.lru_cache(maxsize=32)
    def unstable_graph(max_severity: int = 5):
        """Unstable molecular graph to filter out especially for generative models
        Args:
            max_severity: maximum severity to consider for graph rules to be acceptable
        """
        graph_df = pd.read_csv(get_data("graph.csv"))
        # only apply rules with severity >= max_severity
        graph_df = graph_df[graph_df["severity"] >= max_severity]
        return from_smarts(
            graph_df["smarts"].values,
            graph_df["labels"].values,
            entry_as_inds=True,
        )

    @staticmethod
    @functools.lru_cache(maxsize=32)
    def chemical_groups(filters: Union[str, List[str]] = "medicinal"):
        """Chemical group filter catalogs

        Args:
            filters: list of tag to filter the catalog on.
        """
        chemical_groups = pd.read_csv(get_data("chemical_groups.csv"))
        # we cannot keep the nan values in the dataframe
        chemical_groups = chemical_groups.dropna(subset=["smarts"])
        if isinstance(filters, str):
            filters = [filters]
            chemical_groups = chemical_groups[
                chemical_groups.hierarchy.str.contains("|".join(filters))
            ]

        return from_smarts(
            chemical_groups["smarts"].values,
            chemical_groups["name"].values,
            entry_as_inds=True,
        )
