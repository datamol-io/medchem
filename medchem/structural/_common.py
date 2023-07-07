from typing import List
from typing import Optional
from typing import Union
from typing import Sequence

import os
import functools

import pandas as pd
import datamol as dm

from loguru import logger

from medchem.utils.loader import get_data_path


class CommonAlertsFilters:
    """Filtering class for building a library based on a list of structural alerts

    To list the available alerts, use the `list_default_available_alerts` method.
    """

    def __init__(
        self,
        alerts_set: Optional[Union[str, List[str]]] = None,
        alerts_db_path: Optional[Union[os.PathLike, str]] = None,
    ):
        """Filtering molecules based on chemical alerts

        Args:
            alerts_set: Filter set to use. Default is BMS+Dundee+Glaxo.
            alerts_db_path: Alerts file to use. Default is internal.
        """
        if alerts_db_path is None:
            alerts_db_path = get_data_path(filename="common_alerts_collection.csv")

        if alerts_set is None:
            alerts_set = ["BMS", "Dundee", "Glaxo"]
        elif isinstance(alerts_set, str):
            alerts_set = [alerts_set]

        alerts_set = [x.lower() for x in set(alerts_set)]

        self.alerts_df = self._build_alerts(alerts_set=alerts_set, alerts_db_path=alerts_db_path)

    def _build_alerts(
        self,
        alerts_set: Union[str, List[str]],
        alerts_db_path: Union[os.PathLike, str],
    ):
        # Load the db
        alerts_df = pd.read_csv(alerts_db_path)

        # Select only the filters we want
        mask = alerts_df["rule_set_name"].str.lower().isin(alerts_set)
        alerts_df = alerts_df[mask]
        alerts_df = alerts_df.reset_index(drop=True)

        with dm.without_rdkit_log():
            alerts_df["smarts_mol"] = alerts_df["smarts"].apply(dm.from_smarts)

        # Check for invalid SMARTS
        mask_invalid = alerts_df["smarts_mol"].isnull()
        if mask_invalid.sum() > 0:
            invalid_smarts = alerts_df.loc[mask_invalid, "smarts"].tolist()  # pragma: no cover
            logger.warning(
                f"The following SMARTS are invalid and will be ignored: {invalid_smarts}"
            )  # pragma: no cover

        return alerts_df

    @staticmethod
    @functools.lru_cache()
    def list_default_available_alerts():
        """
        Return a list of unique rule set names
        """
        alerts_db = get_data_path(filename="common_alerts_collection.csv")
        rule_df = pd.read_csv(alerts_db)
        rule_list = (
            rule_df.groupby("rule_set_name")
            .agg(
                {
                    "smarts": "count",
                    "catalog_description": "first",
                    "rule_set": "first",
                    "source": "first",
                }
            )
            .sort_values("rule_set")
            .reset_index()
        )

        return rule_list

    def evaluate(self, mol: Union[str, dm.Mol]):
        """
        Evaluate structure alerts on a molecule

        Args:
            mol: input molecule
        """
        if isinstance(mol, str):
            with dm.without_rdkit_log():
                _mol = dm.to_mol(mol)
        else:
            _mol = mol

        if _mol is None:
            return pd.Series(
                {
                    "mol": _mol,
                    "pass_filter": False,
                    "status": "exclude",
                    "reasons": "invalid",
                }
            )

        # Match the molecule against the alerts
        n_matches = self.alerts_df["smarts_mol"].apply(lambda x: len(_mol.GetSubstructMatches(x)))
        has_alert = n_matches > 0

        # Prepare the result
        result = pd.Series()
        result["mol"] = _mol

        if any(has_alert):
            result["pass_filter"] = False
            result["status"] = "exclude"
            result["reasons"] = ";".join(self.alerts_df[has_alert]["description"].dropna().unique().tolist())

            # NOTE(hadim): we could eventually backpropagate the full alert definition of the matching alerts +
            # match atom indices and also the molecule object with the highlighted atoms.

        else:
            result["pass_filter"] = True
            result["status"] = "ok"
            result["reasons"] = None

        return result

    def __call__(
        self,
        mols: Sequence[Union[str, dm.Mol]],
        n_jobs: Optional[int] = -1,
        progress: bool = False,
        progress_leave: bool = False,
        scheduler: str = "auto",
    ):
        """Run alert evaluation on this list of molecule and return the full dataframe

        Args:
            mols: list of input molecule object.
            n_jobs: number of jobs to run in parallel.
            progress: whether to show progress or not.
            progress_leave: whether to leave the progress bar or not.
            scheduler: which scheduler to use. If "auto", will use "processes" if `len(mols) > 500` else "threads".
        """

        if scheduler == "auto":
            if len(mols) > 500:
                scheduler = "processes"  # pragma: no cover
            else:
                scheduler = "threads"

        results = dm.parallelized(
            self.evaluate,
            mols,
            progress=progress,
            n_jobs=n_jobs,
            scheduler=scheduler,
            tqdm_kwargs=dict(
                desc="Filter by alerts",
                leave=progress_leave,
            ),
        )
        results = pd.DataFrame(results)

        return results
