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

    The output of the filter are explained below:
    - **status**: one of `["exclude", "flag", "annotations", "ok"]` (ordered by quality).
        Generally, you can keep anything without the "exclude" label, as long as you also apply
        a maximum severity score for compounds that collects too many flags.
    - **reasons**: list of reasons why the compound was flagged.
    - **pass_filter**: whether the compound passed the filter or not.
    - **details**: optional additional details of the evaluation, including matching patterns if `keep_details` is True.
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
            alerts_set = ["BMS"]
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
        Return a list of unique alert set names
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

    def _evaluate(
        self,
        mol: Union[str, dm.Mol],
        keep_details: bool = False,
    ):
        """
        Evaluate structural alerts on a molecule

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
        matches = self.alerts_df["smarts_mol"].apply(lambda x: _mol.GetSubstructMatches(x))
        n_matches = matches.apply(len)
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

        if keep_details:
            details = self.alerts_df[has_alert].copy()
            details = details.drop(columns=["smarts_mol"])
            details["matches"] = matches[has_alert].tolist()
            result["details"] = details.to_dict(orient="records")

        return result

    def __call__(
        self,
        mols: Sequence[Union[str, dm.Mol]],
        n_jobs: Optional[int] = -1,
        progress: bool = False,
        progress_leave: bool = False,
        scheduler: str = "auto",
        batch_size: Optional[int] = None,
        keep_details: bool = False,
    ) -> pd.DataFrame:
        """Run alert evaluation on this list of molecule and return the full dataframe

        Args:
            mols: list of input molecule object.
            n_jobs: number of jobs to run in parallel.
            progress: whether to show progress or not.
            progress_leave: whether to leave the progress bar or not.
            scheduler: which scheduler to use. If "auto", will use "processes" if `len(mols) > 500` else "threads".
            batch_size: batch size to use for parallelization.
            keep_details: whether to keep the details of the evaluation or not.
        """

        if scheduler == "auto":
            if len(mols) > 500:
                scheduler = "processes"  # pragma: no cover
            else:
                scheduler = "threads"
        if batch_size:
            results = dm.parallelized_with_batches(
                lambda batch: [
                    functools.partial(self._evaluate, keep_details=keep_details)(mol) for mol in batch
                ],
                mols,
                progress=progress,
                n_jobs=n_jobs,
                scheduler=scheduler,
                batch_size=batch_size,
                tqdm_kwargs=dict(
                    desc="Common alerts filtering",
                    leave=progress_leave,
                ),
            )
        else:
            results = dm.parallelized(
                functools.partial(self._evaluate, keep_details=keep_details),
                mols,
                progress=progress,
                n_jobs=n_jobs,
                scheduler=scheduler,
                tqdm_kwargs=dict(
                    desc="Common alerts filtering",
                    leave=progress_leave,
                ),
            )
        results = pd.DataFrame(results)

        return results
