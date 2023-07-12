from typing import Union
from typing import Callable
from typing import Optional
from typing import List
from typing import Sequence
from typing import Dict

import functools

import pandas as pd
import datamol as dm

from medchem.utils import loader

from ..utils.func import get_function_name
from . import basic_rules


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
    return dm.descriptors.compute_many_descriptors(mols, properties_fn, add_properties=False)


class RuleFilters:
    """
    Build a filter based on a compound phychem properties. For a list of default rules, use `RuleFilters.list_available_rules()`.
    Most of these rules have been collected from the litterature including https://fafdrugs4.rpbs.univ-paris-diderot.fr/descriptors.html
    """

    def __init__(
        self,
        rule_list: List[Union[str, Callable]],
        rule_list_names: Optional[List[Optional[str]]] = None,
    ):
        """Build a rule filtering object

        Args:
            rule_list: list of rules to apply. Either a callable that takes a molecule as input (with kwargs) or a string
                of the name of a pre-defined rule as defined in the basic_rules module
            rule_list_names: Name of the rules passed as inputs. Defaults to None.
        """
        self.rules = self._build_rules(rule_list, rule_list_names)

    def _build_rules(
        self,
        rule_list: List[Union[str, Callable]],
        rule_list_names: Optional[List[Optional[str]]] = None,
    ) -> Dict[str, Callable]:
        """Build the list of rules to apply."""

        if rule_list_names is not None and len(rule_list_names) != len(rule_list):
            raise ValueError("rule_list_names must be the same length as rule_list")

        if rule_list_names is None:
            _rule_list_names = [None] * len(rule_list)
        else:
            _rule_list_names = rule_list_names

            # If rule_list_names is provided, we check whether the list is unique
            if len(set(_rule_list_names)) != len(_rule_list_names):
                raise ValueError("rule_list_names must be unique")

        rules = {}
        for rule_name, rule_fn in zip(_rule_list_names, rule_list):
            if isinstance(rule_fn, str):
                rule = getattr(basic_rules, rule_fn, None)

                if rule is None:
                    raise ValueError(f"Rule {rule_fn} not found")

            elif callable(rule_fn):
                rule = rule_fn

            else:
                raise ValueError(f"Unsupported rule {rule_fn} of type {type(rule_fn)}!")

            if rule_name is None:
                if isinstance(rule_fn, str):
                    rule_name = rule_fn
                else:
                    rule_name = get_function_name(rule_fn)

            rules[rule_name] = rule

        return rules

    def __len__(self):
        """Return the number of rules inside this filter"""
        return len(self.rules)

    def __getitems__(self, name: str):
        """Return a specific rule"""
        return self.rules[name]

    def __call__(
        self,
        mols: Sequence[Union[str, dm.Mol]],
        n_jobs: Optional[int] = -1,
        progress: bool = False,
        progress_leave: bool = False,
        scheduler: str = "auto",
        keep_props: bool = False,
        fail_if_invalid: bool = True,
    ) -> pd.DataFrame:
        """Compute the rules for a list of molecules

        Args:
            mols: list of input molecule object.
            n_jobs: number of jobs to run in parallel.
            progress: whether to show progress or not.
            progress_leave: whether to leave the progress bar or not.
            scheduler: which scheduler to use. If "auto", will use "processes" if `len(mols) > 500` else "threads".
            keep_props: whether to keep the properties columns computed by the rules.
            fail_if_invalid: whether to fail if a rule fails or not.

        Returns:
            df: Dataframe where each row is a molecule and each column is a the outcomes of applying self.rules[column].
        """

        if scheduler == "auto":
            if len(mols) > 500:
                scheduler = "processes"
            else:
                scheduler = "threads"

        def _rule_fn(mol: Union[str, dm.Mol]):
            # Convert to mol object if needed
            if isinstance(mol, str):
                mol = dm.to_mol(mol)

            if mol is None:
                if fail_if_invalid:
                    raise ValueError("Molecule is None")
                else:
                    return pd.Series()

            datum = pd.Series()
            datum["mol"] = mol

            # Precompute properties
            props = _compute_batch_props(mol)

            if keep_props:
                datum = pd.concat([datum, pd.Series(props)])

            # Filter with the rules
            datum["pass_all"] = True
            datum["pass_any"] = False
            for rule_name, rule_fn in self.rules.items():
                # Do the computation
                datum[rule_name] = rule_fn(mol, **props)

                datum["pass_all"] &= datum[rule_name]
                datum["pass_any"] |= datum[rule_name]

            return datum

        results = dm.parallelized(
            _rule_fn,
            mols,
            progress=progress,
            n_jobs=n_jobs,
            scheduler=scheduler,
            tqdm_kwargs=dict(
                desc="Filter by rules",
                leave=progress_leave,
            ),
        )
        results = pd.DataFrame(results)

        return results

    @staticmethod
    @functools.lru_cache(maxsize=32)
    def list_available_rules(*query: str):
        """List all the available rules and they properties"""
        df = pd.read_csv(loader.get_data_path("medchem_rule_list.csv"))
        if len(query) > 0:
            query_str = "|".join(query)
            df = df[df["description"].str.contains(query_str)]
        return df

    @staticmethod
    @functools.lru_cache(maxsize=32)
    def list_available_rules_names(*query: str):
        """List only the names of the available rules"""
        df = RuleFilters.list_available_rules(*query)
        return df["name"].tolist()
