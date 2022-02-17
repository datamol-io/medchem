from typing import Iterable
from typing import List
from typing import Optional
from typing import Union

import os
import copy
import pandas as pd
import datamol as dm

from tqdm.auto import tqdm
from loguru import logger
from rdkit.Chem import rdchem
from rdkit.Chem.Descriptors import MolWt, MolLogP, NumHDonors, NumHAcceptors, TPSA
from medchem.utils import get_data
from medchem.catalog import NamedCatalogs


class NovartisFilters:
    """
    Filtering class for building a screening deck following the novartis filtering process
    published in https://dx.doi.org/10.1021/acs.jmedchem.0c01332.

    The output of the filter are explained below:
    - **status**: one of `["Exclude", "Flag", "Annotations", "Ok"]` (ordered by quality).
        Generally, you can keep anything without the "Exclude" label, as long as you also apply
        a maximum severity score for compounds that collects too many flags.
    - **covalent**: number of potentially covalent motifs contained in the compound
    - **severity**: how severe are the issues with the molecules:
        - `0`: compound has no flags, might have annotations;
        - `1-9`:  number of flags the compound raises;
        - `>= 10`:  default exclusion criterion used in the paper
    - **special_mol**: whether the compound/parts of the compound belongs to a special class of molecules
        (e.g peptides, glycosides, fatty acid). In that case, you should review the rejection reasons.
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

        catalog = NamedCatalogs.nibr()
        if n_jobs is not None:
            if isinstance(mols[0], str):
                mols = dm.parallelized(
                    dm.to_mol, mols, n_jobs=n_jobs, progress=progress
                )
            matches = dm.parallelized(
                catalog.GetMatches,
                mols,
                n_jobs=n_jobs,
                progress=progress,
                scheduler="threads",
            )
        else:
            mols = [dm.to_mol(x) if isinstance(x, str) else x for x in mols]
            iter_mols = mols
            if progress:
                iter_mols = tqdm(mols)
            matches = [catalog.GetMatches(mol) for mol in iter_mols]

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
                        _, name, sev, cov, m = pname.split("||")
                        names.append(name)
                        severity.append(int(sev))
                        covalent.append(int(cov))
                        special_mol.append(int(m))
                    # concatenate all matching filters
                    reasons = "; ".join(names)
                    # severity of 2 means EXCLUDE
                    if severity.count(2):
                        sc = 10
                        status = "Exclude"
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


class ChEMBLFilters:
    """
    Filtering class for building a library based on structural alerts provided by the ChEMBL database

    The following set of alerts is supported:
        * 'Glaxo'
        * 'Dundee'
        * 'BMS'
        * 'PAINS'
        * 'SureChEMBL'
        * 'MLSMR'
        * 'Inpharmatica'
        * 'LINT
    """

    def __init__(
        self,
        alerts_set: Union[str, List[str]] = None,
        alerts_db: Optional[os.PathLike] = None,
    ):
        """Filtering molecules based on chemical alerts

        Args:
            alerts_set: Alerts catalog to use. Default is BMS+Dundee+Glaxo
            alerts_db: Alerts file to use. Default is internal
        """
        if alerts_db is None:
            alerts_db = get_data(file="rd_alerts.csv")
        self.rule_df = pd.read_csv(alerts_db)
        self.rule_list = []
        if alerts_set is None:
            alerts_set = ["BMS", "Dundee", "Glaxo"]
        if isinstance(alerts_set, str):
            alerts_set = [alerts_set]
        self.alerts_set = [x.lower() for x in set(alerts_set)]
        self._build_rule_list()

    def _build_rule_list(self):
        """
        Build the rule sets defined in alerts_set for this object
        """
        self.rule_df = self.rule_df[
            self.rule_df.rule_set_name.str.lower().isin(self.alerts_set)
        ]
        tmp_rule_list = self.rule_df[
            ["rule_id", "smarts", "mincount", "description"]
        ].values.tolist()
        for rule_id, smarts, mincount, desc in tmp_rule_list:
            smarts_mol = dm.from_smarts(smarts)
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
        if isinstance(mol, str):
            mol = dm.to_mol(mol)
        if mol is None:
            return [mol, "Exclude", -999, -999, -999, -999, -999] + [1] * len(
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
