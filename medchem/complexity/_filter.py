from typing import Optional
from typing import Union

import functools

import datamol as dm
import numpy as np
import pandas as pd

from rdkit.Chem import GraphDescriptors

from medchem.utils.loader import get_data_path

from ._calc import WhitlockCT
from ._calc import BaroneCT
from ._calc import SMCM
from ._calc import TWC


class ComplexityFilter:
    """Complexity-based molecular filters which compare some molecular complexity metrics to
    precomputed thresholds on large datasets of commercially available compounds.

    The default threshold have been calculated using all commercially available compounds of the zinc-15 dataset.

    !!! example "Usage"
        ```python
        from medchem.complexity import ComplexityFilter
        complexity_filter = ComplexityFilter(limit="99", complexity_metric="bertz")
        complexity_filter("CC(=O)Nc1ccc(cc1)O") # True
        ```
    """

    COMPLEXITY_FNS = {
        "bertz": GraphDescriptors.BertzCT,  # bertz complexity index
        "sas": dm.descriptors.sas,  # sas score
        "qed": dm.descriptors.qed,  # qed score
        "clogp": dm.descriptors.clogp,  # clogp score
        "whitlock": WhitlockCT,  # whitlock complexity index
        "barone": BaroneCT,  # barone complexity index
        "smcm": SMCM,  # synthetic and molecular complexity
        "twc": TWC,  # Total walk count complexity
    }

    def __init__(
        self,
        limit: str = "99",
        complexity_metric: str = "bertz",
        threshold_stats_file: Optional[str] = "zinc_15_available",
    ):
        """
        Default complexity limit is set on at least 1 exceeding metric on the 999th permille level

        Args:
            limit: The complexity percentile outlier limit to be used (should be a string expressed as an integer)
            complexity_metric: The complexity filter name to be used.
                Use `ComplexityFilter.list_default_available_filters` to list default filters.
                The following complexity metrics are supported by default
                * "bertz": bertz complexity index
                * "sas": synthetic accessibility score  (`zinc_15_available` only)
                * "qed": qed score  (`zinc_15_available` only)
                * "clogp": clogp for how greasy a molecule is compared to other in the same mw range  (`zinc_15_available` only)
                * "whitlock": whitlock complexity index
                * "barone": barone complexity index
                * "smcm": synthetic and molecular complexity
                * "twc":  total walk count complexity  (`zinc_15_available` only)
            threshold_stats_file: The path to or type the threshold file to be used.
                The default available threshold stats files are
                * "zinc_12"
                * "zinc_15_available"

        """
        self.threshold_df = self.load_threshold_stats_file(threshold_stats_file)

        if limit not in self.threshold_df.percentile.unique():
            raise ValueError(f"Invalid value {limit}")

        available_filters = set(self.list_default_available_filters()).intersection(self.threshold_df.columns)

        if complexity_metric not in available_filters:
            raise ValueError(f"Invalid value {complexity_metric}. Should be one of {available_filters}")

        self.limit_index = limit
        self.complexity_metric = complexity_metric
        cur_df = self.threshold_df[self.threshold_df.percentile == self.limit_index]
        self.filter_selection_df = cur_df[[self.complexity_metric, "mw_bins"]].sort_values(
            "mw_bins",
            ascending=True,
        )

    @classmethod
    def list_default_available_filters(cls):
        """
        Return a list of unique filter names
        """
        return list(ComplexityFilter.COMPLEXITY_FNS.keys())

    @classmethod
    @functools.lru_cache()
    def list_default_percentile(cls, threshold_stats_file: Optional[str] = None):
        """
        Return the default percentile list for the threshold file
        """
        df = cls.load_threshold_stats_file(threshold_stats_file)
        return df["percentile"].unique().tolist()

    @classmethod
    def load_threshold_stats_file(cls, path: Optional[str] = None):
        """
        Load threshold file to compute the percentille depending on the MW for each complexity_metric
        Args:
            path: path to the threshold file
        """
        if path is None or path == "zinc_12":
            path = "zinc_12_stats.csv"
            path = get_data_path(path, module="medchem.data.complexity")

        elif path == "zinc_15_available":
            path = "zinc_15_available_stats.csv"
            path = get_data_path(path, module="medchem.data.complexity")

        data = pd.read_csv(path)
        return data

    def __call__(self, mol: Union[dm.Mol, str]):
        """
        Check whether the input structure is too complex given this instance of the complexity filter
        Return False is the molecule is too complex, else True

        Args:
            mol: input molecule
        """

        if isinstance(mol, str):
            mol = dm.to_mol(mol)

        mw = dm.descriptors.mw(mol)
        ind = np.digitize(mw, self.filter_selection_df["mw_bins"], right=True)

        fn = ComplexityFilter.COMPLEXITY_FNS[self.complexity_metric]
        threshold = self.filter_selection_df[self.complexity_metric].values[ind]

        with dm.without_rdkit_log():
            pred = fn(mol)

        return np.isnan(pred) or pred < threshold
