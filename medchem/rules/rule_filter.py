from sched import scheduler
from typing import Union
from typing import Callable
from typing import Optional
from typing import List

import functools
import pandas as pd
import datamol as dm
from tqdm.auto import tqdm
from medchem.utils import loader
from medchem.rules import basic_rules


def _compute_batch_props(mols):
    """Wrapper for making batch computing property under parallelization"""
    properties_fn = {
        "mw": dm.descriptors.mw,
        "n_rings": dm.descriptors.n_rings,
        "n_heavy_atoms": dm.descriptors.n_heavy_atoms,
        "n_rotatable_bonds": dm.descriptors.n_rotatable_bonds,
        "tpsa": dm.descriptors.tpsa,
        "clogp": dm.descriptors.clogp,
        "n_hba": dm.descriptors.n_hba,
        "n_hbd": dm.descriptors.n_hbd,
        "n_lipinski_hba": dm.descriptors.n_lipinski_hba,
        "n_lipinski_hbd": dm.descriptors.n_lipinski_hbd,
    }
    return dm.descriptors.compute_many_descriptors(
        mols, properties_fn, add_properties=False
    )


class RuleFilters:
    """
    Build a filter based on a compound phychem properties. For a list of default rules, use `RuleFilters.list_available_rules()`.
    Most of these rules have been collected from the litterature including https://fafdrugs4.rpbs.univ-paris-diderot.fr/descriptors.html
    """

    def __init__(
        self,
        rule_list: List[Union[str, Callable]],
        rule_list_names: Optional[List[str]] = None,
        precompute_props: bool = True,
    ):
        """Build a rule filtering object

        Args:
            rule_list: list of rules to apply. Either a callable that takes a molecule as input (with kwargs) or a string
                of the name of a pre-defined rule as defined in the basic_rules module
            rule_list_names: Name of the rules passed as inputs. Defaults to None.
            precompute_props: Whether to precompute the properties for all molecules to speed up redundant calculation. Defaults to True.
        """
        self.precompute_props = precompute_props
        self.rules = self._build_rules(rule_list)
        self.rule_names = rule_list_names or []
        if self.rule_names and len(self.rule_names) != len(self.rules):
            raise ValueError("rule_list_names must be the same length as rule_list")

    def _build_rules(self, rule_list: List[Union[str, Callable]]):
        """Build the list of rules to apply

        Args:
            rule_list: list of rules to apply. Either a callable that takes `mol` as first input (and additional kwargs) or a string
                of the name of a pre-defined rule as defined in the basic_rules module
        """
        rules = []
        for rule_name in rule_list:
            if isinstance(rule_name, str):
                rule = getattr(basic_rules, rule_name, None)
                if rule is None:
                    raise ValueError(f"Rule {rule_name} not found")
            elif callable(rule_name):
                rule = rule_name
            else:
                raise ValueError(
                    f"Unsupported rule {rule_name} of type {type(rule_name)}!"
                )
            rules.append(rule)
        return rules

    def __len__(self):
        """Return the number of rules inside this filter"""
        return len(self.rules)

    def __getitems__(self, ind):
        """Return a specific rule"""
        return self.rules[ind]

    def __call__(
        self,
        mols: List[Union[str, dm.Mol]],
        n_jobs: Optional[int] = None,
        progress: bool = False,
        scheduler: str = "processes",
    ):
        """Compute the rules for a list of molecules

        Args:
            mols: list of input molecule object.
            n_jobs: number of jobs to run in parallel. Defaults to None.
            progress: whether to show progress or not. Defaults to False.
            scheduler: which scheduler to use. Defaults to "processes".

        Returns:
            df: Dataframe where each row is a molecule and each column is a the outcomes of applying self.rules[column].
        """

        mols = dm.parallelized(
            dm.to_mol,
            mols,
            n_jobs=n_jobs,
            progress=progress,
            scheduler=scheduler,
            tqdm_kwargs={"desc": "Mol Convert", "leave": False},
        )
        if self.precompute_props:
            inputs = dm.parallelized(
                _compute_batch_props,
                mols,
                progress=progress,
                n_jobs=n_jobs,
                scheduler=scheduler,
                tqdm_kwargs=dict(leave=False),
            )
        else:
            inputs = [{} for _ in mols]

        inputs = list(inputs)
        for i, mol in enumerate(mols):
            inputs[i]["mol"] = mol

        computed_rules = {}
        for i, rule in enumerate(tqdm(self.rules, disable=not progress)):
            computed_rules[i] = dm.parallelized(
                rule,
                inputs,
                arg_type="kwargs",
                progress=progress,
                n_jobs=n_jobs,
                scheduler=scheduler,
                tqdm_kwargs=dict(desc="Props"),
            )

        df = pd.DataFrame.from_dict(computed_rules, orient="index").transpose()
        if self.rule_names is not None and len(self.rule_names) > 0:
            df.columns = self.rule_names
        return df

    @staticmethod
    @functools.lru_cache(maxsize=32)
    def list_available_rules(query: Union[str, List[str]] = None):
        """List all the available rules and they properties"""
        df = pd.read_csv(loader.get_data("medchem_rule_list.csv"))
        if query is not None:
            if isinstance(query, (list, tuple)):
                query = "|".join(query)
            df = df[df["description"].str.contains(query)]
        return df
