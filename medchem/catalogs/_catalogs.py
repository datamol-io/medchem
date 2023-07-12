from typing import List
from typing import Union
from typing import Optional
from typing import Sequence

import functools

from loguru import logger

import pandas as pd
import numpy as np
import datamol as dm

from rdkit.Chem import rdfiltercatalog  # type: ignore

from medchem.utils.loader import get_data_path


def list_named_catalogs():
    """
    List all available named catalogs. This list will not report chemical groups !

    !!! tip
        For a list of chemical groups that can be queried using `NamedCatalog.chemical_groups`,
        use `medchem.groups.list_default_chemical_groups`
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
    params = rdfiltercatalog.FilterCatalogParams()
    missing_catalogs = []
    for catlg in catalogs:
        if isinstance(catlg, rdfiltercatalog.FilterCatalogParams.FilterCatalogs):
            params.AddCatalog(catlg)
        else:
            missing_catalogs.append(catlg)
    parameterized_catalogs = rdfiltercatalog.FilterCatalog(params)
    for catlg in missing_catalogs:
        for entry_nums in range(catlg.GetNumEntries()):
            entry = catlg.GetEntryWithIdx(entry_nums)
            parameterized_catalogs.AddEntry(entry)
    return parameterized_catalogs


def catalog_from_smarts(
    smarts: Union[Sequence[str], np.ndarray, pd.Series],
    labels: Optional[Union[Sequence[str], np.ndarray, pd.Series]] = None,
    mincounts: Optional[Union[Sequence[int], np.ndarray, pd.Series]] = None,
    maxcounts: Optional[Union[Sequence[int], np.ndarray, pd.Series]] = None,
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
        catalog (FilterCatalog): filter catalog built from the input smarts
    """

    with dm.without_rdkit_log():
        catalog = rdfiltercatalog.FilterCatalog()
        if labels is None:
            labels = smarts
        if mincounts is None:
            mincounts = [1] * len(smarts)

        for i, (sm, lb, count) in enumerate(zip(smarts, labels, mincounts)):
            if maxcounts is None:
                fil = rdfiltercatalog.SmartsMatcher(lb, sm, count)
            else:
                fil = rdfiltercatalog.SmartsMatcher(lb, sm, count, maxcounts[i])
            entry_name = str(lb)
            if entry_as_inds:
                entry_name = str(i)

            if fil.IsValid():
                catalog.AddEntry(rdfiltercatalog.FilterCatalogEntry(entry_name, fil))
            else:
                logger.warning(f"It seem like the SMARTS {sm} is invalid")
    return catalog


class NamedCatalogs:
    """
    Holder for substructure matching catalogs. This class provides several popular and custom catalog
    that can be used for substructure matching.

    All the catalogs are cached using `functools.lru_cache` to avoid reloading them every time.

    !!! note
        A filter catalog is a collection of substructures and molecular patterns (SMARTS) used to flag molecules with **(un)desirable structural properties**.
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
            pains_a: whether to include PAINS filters (assay A)
            pains_b: whether to include PAINS filters (assay B)
            pains_c: whether to include PAINS filters (assay C)
            brenk: whether to include BRENK filters, also known as Dundee filters
            nih: whether to include NIH filters
            zinc: whether to include ZINC filters
        """
        catalogs = []
        if pains_a:
            catalogs.append(rdfiltercatalog.FilterCatalogParams.FilterCatalogs.PAINS_A)
        if pains_b:
            catalogs.append(rdfiltercatalog.FilterCatalogParams.FilterCatalogs.PAINS_B)
        if pains_c:
            catalogs.append(rdfiltercatalog.FilterCatalogParams.FilterCatalogs.PAINS_C)
        if brenk:
            catalogs.append(rdfiltercatalog.FilterCatalogParams.FilterCatalogs.BRENK)
        if nih:
            catalogs.append(rdfiltercatalog.FilterCatalogParams.FilterCatalogs.NIH)
        if zinc:
            catalogs.append(rdfiltercatalog.FilterCatalogParams.FilterCatalogs.ZINC)
        catalog = merge_catalogs(*catalogs)
        return catalog

    @staticmethod
    @functools.lru_cache(maxsize=32)
    def pains():
        return rdfiltercatalog.FilterCatalog(rdfiltercatalog.FilterCatalogParams.FilterCatalogs.PAINS)

    @staticmethod
    @functools.lru_cache(maxsize=32)
    def pains_a():
        return rdfiltercatalog.FilterCatalog(rdfiltercatalog.FilterCatalogParams.FilterCatalogs.PAINS_A)

    @staticmethod
    @functools.lru_cache(maxsize=32)
    def pains_b():
        return rdfiltercatalog.FilterCatalog(rdfiltercatalog.FilterCatalogParams.FilterCatalogs.PAINS_B)

    @staticmethod
    @functools.lru_cache(maxsize=32)
    def pains_c():
        return rdfiltercatalog.FilterCatalog(rdfiltercatalog.FilterCatalogParams.FilterCatalogs.PAINS_C)

    @staticmethod
    @functools.lru_cache(maxsize=32)
    def nih():
        return rdfiltercatalog.FilterCatalog(rdfiltercatalog.FilterCatalogParams.FilterCatalogs.NIH)

    @staticmethod
    @functools.lru_cache(maxsize=32)
    def zinc():
        return rdfiltercatalog.FilterCatalog(rdfiltercatalog.FilterCatalogParams.FilterCatalogs.ZINC)

    @staticmethod
    @functools.lru_cache(maxsize=32)
    def brenk():
        return rdfiltercatalog.FilterCatalog(rdfiltercatalog.FilterCatalogParams.FilterCatalogs.BRENK)

    @staticmethod
    def alerts(subset: Optional[Union[List[str], str]] = None):
        """Common alerts filter catalogs commonly used in molecule filtering

        Args:
            subset: subset of alert providers to consider

        Returns:
            catalog (FilterCatalog): filter catalog
        """
        rd_filters = pd.read_csv(get_data_path("common_alerts_collection.csv"))
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
        return catalog_from_smarts(
            rd_filters["smarts"],
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
            This will return all the compounds in the catalog, regardless of their severity (`FLAG`, `EXCLUDE`, `ANNOTATION`).
            You likely don't want to use this for blind prioritization.
        """
        nibr_filters = pd.read_csv(get_data_path("nibr.csv"))
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
        return catalog_from_smarts(
            nibr_filters["smarts"],
            labels,
            mincount,
            entry_as_inds=False,
        )

    @staticmethod
    @functools.lru_cache(maxsize=32)
    def bredt():
        """Bredt Filter Rules: a catalog for filtering unstable molecules, ideal for molecules generated by
        deep learning models or chemical space enumeration.
        See example of usage by [surge](https://github.com/StructureGenerator/SURGE/blob/main/doc/surge1_0.pdf)
        """
        bredt_df = pd.read_csv(get_data_path("bredt.csv"))
        return catalog_from_smarts(
            bredt_df["smarts"],
            bredt_df["labels"],
            entry_as_inds=True,
        )

    @staticmethod
    @functools.lru_cache(maxsize=32)
    def unstable_graph(severity_threshold: int = 5):
        """Unstable molecular graph to filter out, ideal for generated *de novo* molecules.

        !!! warning
            This method returns problematic patterns and thus patterns with higher severity than the threshold.

        Args:
            severity_threshold: minimum severity for a pattern to **be returned**.
        """
        graph_df = pd.read_csv(get_data_path("graph.csv"))
        # only apply rules with severity >= severity_threshold
        graph_df = graph_df[graph_df["severity"] >= severity_threshold]
        return catalog_from_smarts(
            graph_df["smarts"],
            graph_df["labels"],
            entry_as_inds=True,
        )

    @staticmethod
    @functools.lru_cache(maxsize=32)
    def chemical_groups(filters: Union[str, List[str]] = "medicinal"):
        """Chemical group filter catalogs

        Args:
            filters: list of tag to filter the catalog on.
        """
        chemical_groups = pd.read_csv(get_data_path("chemical_groups.csv"))
        # we cannot keep the nan values in the dataframe
        chemical_groups = chemical_groups.dropna(subset=["smarts"])
        if isinstance(filters, str):
            filters = [filters]
            chemical_groups = chemical_groups[chemical_groups.hierarchy.str.contains("|".join(filters))]

        return catalog_from_smarts(
            chemical_groups["smarts"],
            chemical_groups["name"],
            entry_as_inds=True,
        )
