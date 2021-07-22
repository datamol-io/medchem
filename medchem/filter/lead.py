from typing import Iterable
from typing import List
from typing import Optional
from typing import Dict
from typing import Union

import os
import sys
import copy
import functools
import numpy as np
import multiprocessing as mp
import pandas as pd
import datamol as dm

from tqdm.auto import tqdm
from loguru import logger
from rdkit.Chem import rdchem
from rdkit.Chem import MolFromSmarts
from rdkit.Chem.Descriptors import MolWt, MolLogP, NumHDonors, NumHAcceptors, TPSA
from medchem.utils import get_data
from medchem.catalog import NamedCatalogs


class AlertsFilters:
    """
    Class for managing filters
    """

    def __init__(
        self,
        alerts_set: Union[str, List[str]] = ["BMS"],
        alerts_db: Optional[os.PathLike] = None,
    ):
        """Filtering molecules based on chemical alerts

        Args:
            alerts_set: Alerts catalog to use. Default is BMS
            alerts_db: Alerts file to use. Default is internal
        """
        if alerts_db is None:
            alerts_db = get_data(file="rd_alerts.csv")
        self.rule_df = pd.read_csv(alerts_db)
        self.rule_list = []

        if isinstance(alerts_set, str):
            alerts_set = [alerts_set]
        self.alerts_set = list(set(alerts_set))
        self._build_rule_list()

    def _build_rule_list(self):
        """
        Build the rule sets defined in alerts_set for this object
        """
        self.rule_df = self.rule_df[self.rule_df.rule_set_name.isin(self.alerts_set)]
        tmp_rule_list = self.rule_df[
            ["rule_id", "smarts", "mincount", "description"]
        ].values.tolist()
        for rule_id, smarts, mincount, desc in tmp_rule_list:
            smarts_mol = MolFromSmarts(smarts)
            if smarts_mol:
                self.rule_list.append([smarts_mol, mincount, desc])
            else:
                logger.warning(f"Error parsing SMARTS for rule {rule_id}")

    def get_alert_sets(self):
        """
        Return a list of unique rule set names
        """
        return self.rule_df.rule_set_name.unique()

    def evaluate(self, mol: Union[str, rdchem.Mol]):
        """
        Evaluate structure alerts on a molecule

        Args:
            mol: input molecule

        Returns:
            list of alerts matched
        """
        mol = dm.to_mol(mol)
        if mol is None:
            return [mol, "INVALID", -999, -999, -999, -999, -999] + [1] * len(
                self.rule_list
            )

        desc_list = [
            MolWt(mol),
            MolLogP(mol),
            NumHDonors(mol),
            NumHAcceptors(mol),
            TPSA(mol),
        ]
        alerts = [
            int(len(mol.GetSubstructMatches(patt)) >= mincount)
            for patt, mincount, desc in self.rule_list
        ]
        status = "Ok"
        reasons = None
        if any(alerts):
            status = "Exclude"
            reasons = "; ".join(
                [x[-1] for i, x in enumerate(self.rule_list) if alerts[i]]
            )

        return [dm.to_smiles(mol), status, reasons] + desc_list + alerts

    def __call__(
        self,
        mols: Iterable[Union[str, rdchem.Mol]],
        n_jobs: Optional[int] = None,
        progress: bool = False,
        include_all_alerts: bool = False,
    ):
        """Run alert evaluation on this list of molecule and return the full dataframe

        Args:
            mols: input list of molecules
            n_jobs: number of jobs
            progress: whether to show progress or not
            include_all_alerts: whether to include all of the alerts that match as columns
        """
        if n_jobs is not None:
            alert_filter = copy.deepcopy(self)
            results = dm.parallelized(
                alert_filter.evaluate, mols, n_jobs=n_jobs, progress=progress
            )
        else:
            iter_mols = mols
            if progress:
                iter_mols = tqdm(mols)
            results = [self.evaluate(mol) for mol in iter_mols]

        df = pd.DataFrame(
            results,
            columns=[
                "_smiles",
                "status",
                "reasons",
                "MW",
                "LogP",
                "HBD",
                "HBA",
                "TPSA",
            ]
            + [str(x[-1]) for x in self.rule_list],
        )
        if not include_all_alerts:
            df = df[
                ["_smiles", "status", "reasons", "MW", "LogP", "HBD", "HBA", "TPSA"]
            ]
        return df


class NovartisFilters:
    """
    Filtering class for building a screening deck following the novartis filtering process
    published in https://dx.doi.org/10.1021/acs.jmedchem.0c01332.

    This filters also provide a severity score:
        - 0 -> compound has no flags, might have annotations;
        - 1-9 number of flags the compound raises;
        - >= 10 exclusion criterion for our newly designed screening deck
    """

    def __call__(
        self,
        mols: Iterable[Union[str, rdchem.Mol]],
        n_jobs: Optional[int] = None,
        progress: bool = False,
    ):
        """Run alert evaluation on this list of molecule and return the full dataframe

        Args:
            mols: input list of molecules
            n_jobs: number of jobs
            progress: whether to show progress or not
        """

        catalog = FilterCatalog.nibr_catalog()
        if n_jobs is not None:
            mols = dm.parallelized(dm.to_mol, mols)
            matches = dm.parallelized(
                catalog.HasMatch, mols, n_jobs=n_jobs, progress=progress
            )
        else:
            mols = [dm.to_mol(x) for x in mols]
            iter_mols = mols
            if progress:
                iter_mols = tqdm(mols)
            matches = [catalog.HasMatch(mol) for mol in iter_mols]

        results = []
        for i, (mol, entries) in enumerate(zip(mols, matches)):
            status = "Ok"
            smiles = None
            reasons = None
            co = None
            sm = None
            sc = 0
            try:
                smiles = dm.to_smiles(mol)
                if len(list(entries)):
                    # initialize empty lists
                    names, severity, covalent, special_mol = ([] for _ in range(4))
                    # get the matches
                    for entry in entries:
                        pname = entry.GetDescription()
                        name, sev, cov, m = pname.split("||")
                        names.append(name)
                        severity.append(int(sev))
                        covalent.append(int(cov))
                        special_mol.append(int(m))
                    # concatenate all matching filters
                    reasons = "; ".join(names)
                    # severity of 2 means EXCLUDE
                    if severity.count(2):
                        sc = 10
                        status = "Fail"
                    else:
                        sc = sum(severity)
                        if severity.count(1):
                            status = "Flag"
                        elif severity.count(0):
                            status = "Annotations"
                    # get number of covalent flags and special molecule flags
                    co = sum(covalent)
                    sm = sum(special_mol)
            except Exception as e:
                logger.warning(f"Fail on molecule at index {i}")

            results.append([smiles, status, reasons, sc, co, sm])
        df = pd.DataFrame(
            results,
            columns=[
                "_smiles",
                "status",
                "reasons",
                "severity",
                "covalent",
                "special_mol",
            ],
        )
        return df


def lead_filter(
    mols: Iterable[Union[str, rdchem.Mol]],
    alerts: List[str],
    alerts_db: Optional[os.PathLike] = None,
    n_jobs: Optional[int] = 1,
    rule_dict: Dict = None,
    return_idx=False,
):
    r"""Filter a dataset of molecules, based on structural alerts and specific rules.

    Arguments:
        mols: List of molecules to filter
        alerts: List of alert collections to screen for.
            Supported collections are:
            * 'Glaxo'
            * 'Dundee'
            * 'BMS'
            * 'PAINS'
            * 'SureChEMBL'
            * 'MLSMR'
            * 'Inpharmatica'
            * 'LINT
        alerts_db: Path to the alert file name.
            The internal default file (alerts.csv) will be used if not provided
        n_jobs: Number of cpu to use
        rule_dict: Dictionary with additional rules to apply during the filtering.
            For example, such dictionary for drug-like compounds would look like this:
            >>> rule_dict
                {"MW": [0, 500], "LogP": [-0.5, 5], "HBD": [0, 5], "HBA": [0, 10], "TPSA": [0, 150]}
        return_idx: Whether to return the filtered index

    Returns:
        filtered_mask: boolean array (or index array) where true means the molecule is ok.
    """

    custom_filters = AlertsFilters(alerts_set=alerts, alerts_db=alerts_db)
    df = custom_filters(mols, n_jobs=n_jobs, progress=False)
    df = df[df.status == "Ok"]
    if rule_dict is not None and len(rule_dict) > 0:
        df = df[
            (df.MW.between(*rule_dict["MW"]))
            & (df.LogP.between(*rule_dict["LogP"]))
            & (df.HBD.between(*rule_dict["HBD"]))
            & (df.HBA.between(*rule_dict["HBA"]))
            & (df.TPSA.between(*rule_dict["TPSA"]))
        ]

    filtered_idx = df.index.values
    filtered_mask = np.zeros(len(mols), dtype=bool)
    filtered_mask[filtered_idx] = True
    if return_idx:
        return filtered_idx
    return filtered_mask


def screening_filter(
    mols: Iterable[Union[str, rdchem.Mol]],
    n_jobs: Optional[int] = None,
    max_severity: int = 10,
    return_idx: bool = False,
):
    """
    Filter a set of molecules based on novartis screening deck curation process
    Schuffenhauer, A. et al. Evolution of Novartis' small molecule screening deck design, J. Med. Chem. (2020)
    DOI. https://dx.doi.org/10.1021/acs.jmedchem.0c01332

    Args:
        mols: list of input molecules
        n_jobs: number of parallel job to run. Sequential by default
        max_severity: maximum severity allowed. Default is <10
        return_idx: Whether to return the filtered index

    Returns:
        filtered_mask: boolean array (or index array) where true means the molecule is ok.

    """

    filt_obj = NovartisFilter()
    df = filt_obj(mols, n_jobs=n_jobs)
    df = df[(df.status == "Ok") & (df.severity < max_severity)]
    filtered_idx = df.index.values
    filtered_mask = np.zeros(len(mols), dtype=bool)
    filtered_mask[filtered_idx] = True
    if return_idx:
        return filtered_idx
    return filtered_mask


def common_filter(
    mols: Iterable[Union[str, rdchem.Mol]],
    pains_a: bool = True,
    pains_b: bool = True,
    pains_c: bool = True,
    brenk: bool = True,
    nih: bool = False,
    zinc: bool = False,
    pains: bool = False,
    n_jobs: Optional[int] = None,
    return_idx: bool = False,
):
    """Filter a list of compounds according to common toxicity alerts

    Args:
        mols: list of input molecules
        pains_a: whether to include PAINS filters from assay A
        pains_b: whether to include PAINS filters from assay B
        pains_c: whether to include PAINS filters from assay C
        brenk: whether to include BRENK filters
        nih: whether to include NIH filters
        zinc: whether to include ZINC filters
        pains: whether to include all PAINS filters
        n_jobs: number of parallel job to run. Sequential by default
        return_idx: whether to return index of a boolean mask

    Returns:
        filtered_mask: boolean array (or index array) where true means the molecule is not toxic.
    """

    if pains:
        pains_a = pains_b = pains_c = True

    catalog = NamedCatalogs.tox_catalog(
        pains_a=pains_a,
        pains_b=pains_b,
        pains_c=pains_c,
        brenk=brenk,
        zinc=zinc,
        nih=nih,
    )
    toxic = [False] * len(mols)
    if n_jobs is not None:
        mols = dm.parallelized(dm.to_mol, mols)
        toxic = dm.parallelized(catalog.HasMatch, mols, n_jobs=n_jobs)
    else:
        mols = [dm.to_mol(x) for x in mols]
        toxic = [catalog.HasMatch(mol) for mol in mols]

    filtered_idx = [i for i, bad in enumerate(toxic) if not bad]
    if return_idx:
        return filtered_idx
    return np.bitwise_not(toxic)
