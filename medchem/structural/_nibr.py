from typing import Sequence
from typing import Optional
from typing import Union

import functools

import pandas as pd
import datamol as dm

from medchem.catalogs import NamedCatalogs


class NIBRFilters:
    """Filtering class for building a screening deck following the novartis filtering process
    published in https://dx.doi.org/10.1021/acs.jmedchem.0c01332.

    The output of the filter are explained below:
    - **status**: one of `["exclude", "flag", "annotations", "ok"]` (ordered by quality).
        Generally, you can keep anything without the "exclude" label, as long as you also apply
        a maximum severity score for compounds that collects too many flags.
    - **n_covalent_motif**: number of potentially covalent motifs contained in the compound
    - **severity**: how severe are the issues with the molecules:
        - `0`: compound has no flags, might have annotations;
        - `1-9`:  number of flags the compound raises;
        - `>= 10`:  default exclusion criterion used in the paper
    - **special_mol**: whether the compound/parts of the compound belongs to a special class of molecules
        (e.g peptides, glycosides, fatty acid). In that case, you should review the rejection reasons.
    - **pass_filter**: whether the compound passed the filter or not.
    - **details**: optional additional details of the evaluation, including matching patterns if `keep_details` is True.
    """

    def __init__(self):
        self.catalog = NamedCatalogs.nibr()

    def _evaluate(
        self,
        mol: Union[str, dm.Mol],
        keep_details: bool = False,
    ):
        """Evaluate one molecule."""

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
                    "reasons": "invalid",
                    "severity": 0,
                    "status": "exclude",
                    "n_covalent_motif": 0,
                    "special_mol": 0,
                }
            )

        # Get the matches
        entries = self.catalog.GetMatches(_mol)

        if len(entries) == 0:
            return pd.Series(
                {
                    "mol": _mol,
                    "pass_filter": True,
                    "reasons": None,
                    "severity": 0,
                    "status": "ok",
                    "n_covalent_motif": 0,
                    "special_mol": 0,
                }
            )

        result = pd.Series()
        result["mol"] = _mol

        # Iterate over all the matchign entries
        result_entries = []
        for entry in entries:
            pname = entry.GetDescription()
            _, name, severity, n_covalent_motif, special_mol = pname.split("||")

            result_entry = {}
            result_entry["name"] = name
            result_entry["severity"] = int(severity)
            result_entry["n_covalent_motif"] = int(n_covalent_motif)
            result_entry["special_mol"] = bool(special_mol)
            result_entries.append(result_entry)

        result_entries = pd.DataFrame(result_entries)

        # Now build a flat result from the detailed results per entries

        result["reasons"] = "; ".join(result_entries["name"].tolist())

        # severity of 2 means EXCLUDE
        if 2 in result_entries["severity"].values:
            result["severity"] = 10
            result["status"] = "exclude"
        else:
            result["severity"] = result_entries["severity"].sum()

            if 1 in result_entries["severity"].values:
                result["status"] = "flag"
            elif 0 in result_entries["severity"].values:
                result["status"] = "annotations"

        # get number of covalent flags and special molecule flags
        result["n_covalent_motif"] = result_entries["n_covalent_motif"].sum()
        result["special_mol"] = result_entries["special_mol"].sum()

        if result["status"] == "exclude":
            result["pass_filter"] = False
        else:
            result["pass_filter"] = True

        if keep_details:
            result["details"] = result_entries.to_dict(orient="records")

        return result

    def __call__(
        self,
        mols: Sequence[Union[str, dm.Mol]],
        n_jobs: Optional[int] = -1,
        progress: bool = False,
        progress_leave: bool = False,
        scheduler: str = "threads",
        keep_details: bool = False,
    ):
        """Run alert evaluation on this list of molecule and return the full dataframe

        Args:
            mols: list of input molecule object.
            n_jobs: number of jobs to run in parallel.
            progress: whether to show progress or not.
            progress_leave: whether to leave the progress bar or not.
            scheduler: which scheduler to use. The `processes` scheduler works but is very
                inefficient due to RDKit Catalog serialization which tends to be very slow.
            keep_details: whether to keep the details of the evaluation or not.
        """

        results = dm.parallelized(
            functools.partial(self._evaluate, keep_details=keep_details),
            mols,
            progress=progress,
            n_jobs=n_jobs,
            scheduler=scheduler,
            tqdm_kwargs=dict(
                desc="NIBR filtering",
                leave=progress_leave,
            ),
        )
        results = pd.DataFrame(results)

        return results
