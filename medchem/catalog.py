import functools
import pandas as pd
from rdkit.Chem import FilterCatalog
from medchem.utils import get_data


class NamedCatalogs:
    """
    Holder for substructure matching catalogs
    """

    @staticmethod
    @functools.lru_cache(maxsize=32)
    def tox_catalog(
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
    def nibr_catalog():
        """Catalog from NIBR"""
        nibr_filters = pd.read_csv(get_data("nibr.csv"))
        catalog = FilterCatalog.FilterCatalog()
        for i in range(nibr_filters.shape[0]):
            mincount = int(max(nibr_filters["mincount"].values[i], 1))
            pname = nibr_filters["description"].values[i]
            sname = nibr_filters["rule_set_name"].values[i]
            severity = nibr_filters["severity"].values[i]
            covalent = nibr_filters["covalent"].values[i]
            special_mol = nibr_filters["special_mol"].values[i]
            pname_final = "{0}_min({1})||{2}||{3}||{4}".format(
                pname,
                mincount,
                severity,
                covalent,
                special_mol,
            )
            fil = FilterCatalog.SmartsMatcher(
                pname_final, nibr_filters["smarts"].values[i], mincount
            )
            catalog.AddEntry(FilterCatalog.FilterCatalogEntry(pname_final, fil))
            catalog.GetEntry(i).SetProp("Scope", sname)

        return catalog
