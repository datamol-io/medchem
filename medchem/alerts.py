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
from rdkit.Chem import MolFromSmarts
from rdkit.Chem.Descriptors import MolWt, MolLogP, NumHDonors, NumHAcceptors, TPSA
from medchem.utils import get_data


class AlertFilters:
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
