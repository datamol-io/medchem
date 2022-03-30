from typing import Union
from typing import Callable
from typing import Optional
from typing import List

import ast
import re
import datamol as dm

from functools import partial
from lark import Lark
from loguru import logger
from medchem.catalog import NamedCatalogs
from medchem.catalog import list_named_catalogs
from medchem.groups import list_functional_group_names
from medchem.groups import _get_functional_group_map
from medchem.rules._utils import list_descriptors
from medchem.rules.rule_filter import RuleFilters
from medchem.rules import basic_rules
from medchem.utils.loader import get_grammar
from medchem.query.parser import QueryParser


class QueryOperator:
    """A class to hold all the operators that can be used in queries"""

    AVAILABLES_PROPERTIES = list_descriptors()  # list of available molecular properties
    AVAILABLE_CATALOGS = list_named_catalogs()  # list of available catalogs
    AVAILABLE_RULES = (
        RuleFilters.list_available_rules_names()
    )  # list of available rules
    AVAILABLE_FUNCTIONAL_GROUPS = (
        list_functional_group_names()
    )  # list of available functional groups

    @staticmethod
    def hassubstructure(
        mol: Union[dm.Mol, str],
        query: str,
        is_smarts: bool = False,
        operator: str = "min",
        limit: int = 1,
    ):
        """Check if a molecule has substructure provided by a query

        Args:
            mol: input molecule
            query: input smarts query
            is_smarts: whether this is a smarts query or not
            operator: one of min or max to specify the min or max limit
            limit: limit of substructures to be found

        Returns:
            has_substructure (bool): whether the query is a subgraph of the molecule
        """
        if is_smarts:
            query = dm.from_smarts(query)
        else:
            query = dm.to_mol(query)
        if limit is None:
            limit = 1
        if operator is None:
            operator = "min"
        mol = dm.to_mol(mol)
        matches = mol.GetSubstructMatches(query, uniquify=True)
        coeff = -1 if operator == "min" else 1
        return len(matches) * coeff <= limit * coeff

    @staticmethod
    def hassuperstructure(mol: Union[dm.Mol, str], query: str):
        """Check if a molecule has a superstructure defined by a query.
        Note that a superstructure cannot be a query (smarts)

        Args:
            mol: input molecule
            query: input smarts query

        Returns:
            has_superstructure (bool): whether the molecule is a subgraph of the query
        """
        query = dm.to_mol(query)
        mol = dm.to_mol(mol)
        return query.HasSubstructMatch(mol)

    @staticmethod
    def hasalert(mol: Union[dm.Mol, str], alert: str):
        """Check if a molecule match a named alert catalog.
        The alert catalog needs to be one supported by the medchem package.

        Args:
            mol: input molecule
            alert: named catalog to apply as filter on the molecule

        Returns:
            has_alert (bool): whether the molecule has a given alert
        """

        mol = dm.to_mol(mol)
        if isinstance(alert, str):
            if alert not in QueryOperator.AVAILABLE_CATALOGS:
                raise ValueError(
                    f"Alert {alert} is not supported. Available alerts are: {QueryOperator.AVAILABLE_CATALOGS}"
                )
            alert = getattr(NamedCatalogs, alert, lambda: None)()
        return alert.HasMatch(mol)

    @staticmethod
    def matchrule(mol: Union[dm.Mol, str], rule: str):
        """Check if a molecule match a druglikeness rule

        Args:
            mol: input molecule
            rule: druglikeness rule check on the molecule.

        Returns:
            match_rule (bool): whether the molecule match the given rule
        """

        mol = dm.to_mol(mol)
        if isinstance(rule, str):
            if rule not in QueryOperator.AVAILABLE_RULES:
                raise ValueError(
                    f"Rule {rule} is not supported. Available rules are: {QueryOperator.AVAILABLE_RULES}"
                )
            rule = getattr(basic_rules, rule, None)
        return rule(mol)

    @staticmethod
    def hasgroup(mol: Union[dm.Mol, str], group: str):
        """Check if a molecule has a specific functional group.
        Internally, this is done fetching the smarts corresponding to the group
        then calling `QueryOperator.hassubstructure`

        Args:
            mol: input molecule
            group: functional group to check on the molecule.

        Returns:
            has_group (bool): whether the molecule has the given functional group
        """

        mol = dm.to_mol(mol)
        if group not in QueryOperator.AVAILABLE_FUNCTIONAL_GROUPS:
            raise ValueError(
                f"Functional Group {group} is not supported. Available functional group are: {QueryOperator.AVAILABLE_FUNCTIONAL_GROUPS}"
            )
        group_smarts = _get_functional_group_map()[group]
        return QueryOperator.hassubstructure(mol, group_smarts, is_smarts=True)

    @staticmethod
    def hasprop(mol: Union[dm.Mol, str], prop: str, comparator: Callable, limit: float):
        """Check if a molecule has a molecule property within desired range

        Args:
            mol: input molecule
            prop: molecular property to apply as filter on the molecule
            comparator: operator function to apply to check whether the molecule property matches the expected value
            limit: limit value for determining whether the molecule property is within desired range

        Returns:
            has_property (bool): whether the molecule has a given property within a desired range
        """

        mol = dm.to_mol(mol)
        if isinstance(prop, str):
            if prop not in QueryOperator.AVAILABLES_PROPERTIES:
                raise ValueError(
                    f"Property {prop} is not supported. Available properties are: {QueryOperator.AVAILABLES_PROPERTIES}"
                )
            prop = getattr(dm.descriptors, prop, None)

        computed = prop(mol)
        return comparator(computed, limit)

    @staticmethod
    def getprop(mol: Union[dm.Mol, str], prop: str):
        """Compute the molecular property if a molecule.
        This is an alternative to the hasprop function, that does not enforce any comparison.

        Args:
            mol: input molecule
            prop: molecular property to apply as filter on the molecule

        Returns:
            property (float):  computed property value
        """

        mol = dm.to_mol(mol)
        prop_fn = None
        if isinstance(prop, str):
            if prop not in QueryOperator.AVAILABLES_PROPERTIES:
                raise ValueError(
                    f"Property {prop} is not supported. Available properties are: {QueryOperator.AVAILABLES_PROPERTIES}"
                )
            prop_fn = getattr(dm.descriptors, prop, None)

        return prop_fn(mol)

    @staticmethod
    def like(
        mol: Union[dm.Mol, str],
        query: Union[dm.Mol, str],
        comparator: Callable[[float, float], bool],
        limit: float,
    ):
        """Check if a molecule is similar or distant enough from another molecule using tanimoto ECFP distance.
        and is useful for letting python handles the binary comparison operators.

        Args:
            mol: input molecule
            query: input molecule to compare with
            comparator: operator function to apply to check whether the molecule property matches the expected value.
                Takes computed_similarity and `limit` as arguments and returns a boolean.
            limit: limit value for determining whether the molecule property is within desired range

        Returns:
            is_similar (bool): whether the molecule is similar or distant enough from the query
        """

        mol = dm.to_mol(mol)
        query = dm.to_mol(query)
        sim_val = 1 - float(dm.cdist([mol], [query]))
        return comparator(sim_val, limit)

    @staticmethod
    def similarity(
        mol: Union[dm.Mol, str],
        query: Union[dm.Mol, str],
    ):
        """Compute the ECFP tanimoto similarity between two molecules.
        This is an alternative to the like function, that does not enforce any comparison,
        and is useful for letting python handles the binary comparison operators.

        Args:
            mol: input molecule
            query: input query molecule to compute similarity against

        Returns:
            similarity (float): computed similarity value between mol and query
        """

        mol = dm.to_mol(mol)
        query = dm.to_mol(query)
        sim_val = 1 - float(dm.cdist([mol], [query]))
        return sim_val


class _NodeEvaluator:
    """Representation and evaluation of a single node"""

    def __init__(self, node_expr: str):
        self.node_expr = node_expr
        self.node_fn = None
        if self.node_expr.startswith("fn("):
            # EN: revisit this with a regexp eventually for robustness
            node_expr = node_expr[3:-1]  # remove fn( and )
            node_arg_list = node_expr.split(", ")
            _fn = getattr(QueryOperator, node_arg_list[0], None)
            if _fn is None:
                raise ValueError("Unknown function {}".format(node_arg_list[0]))
            _kwargs = dict(
                (k, ast.literal_eval(v))
                for k, v in [x.split("=") for x in node_arg_list[1:]]
            )
            self.node_fn = partial(_fn, **_kwargs)

    def __call__(self, *args, **kwargs):
        if self.node_fn is None:
            return self.node_expr
        return self.node_fn(*args, **kwargs)


class _EvaluableQuery:
    """Parser of a query into a list of evaluable fonction node"""

    FN_PATTERN = re.compile("`(.*?)`")

    def __init__(self, query: str, verbose: bool = False):
        """Constructor for query evaluation"""
        self.query_nodes = [_NodeEvaluator(x) for x in self.FN_PATTERN.split(query)]
        self.verbose = verbose

    def __call__(self, mol: Union[dm.Mol, str], exec: bool = True):
        """Evaluate a query on an input molecule

        Args:
            mol: input molecule
            exec: whether to interpret the resulting query or not

        Returns:
            query string or boolean value corresponding to the query result

        """
        query_eval = " ".join([f"{node(mol)}" for node in self.query_nodes])
        if self.verbose:
            logger.debug(query_eval)
        if exec:
            # EN: eval is not safe, but we are using it
            # because ast.literal_eval cannot parse some tree structures
            return eval(query_eval)
        return query_eval


class QueryFilter:
    """
    Query filtering system based on a custom query grammar
    """

    def __init__(self, query: str, grammar: Optional[str] = None, parser: str = "lalr"):
        if parser not in ["earley", "lalr"]:
            raise AttributeError("parser must be either 'earley' or 'lalr'")
        # if existing file, then load and parse it
        self.grammar = get_grammar(grammar, as_string=True)
        self.query_parser = Lark(self.grammar, parser=parser)
        self.transformer = QueryParser()
        self._query_str = query
        self.query = self.transformer.transform(
            self.query_parser.parse(self._query_str)
        )
        self._evaluable_query = _EvaluableQuery(self.query)

    def __repr__(self) -> str:
        return self.query

    def __call__(
        self,
        mols: List[Union[str, dm.Mol]],
        scheduler="processes",
        n_jobs: int = -1,
        progress: bool = True,
    ):
        """Call the internal chemical filter that has been build

        Args:
            mols: list of input molecules to filter
            n_jobs: whether to run job in parallel and number of jobs to consider. Defaults to -1.
            scheduler: scheduler to use. Defaults to 'processes'.
            progress: whether to show job progress. Defaults to True.
        """
        tqdm_kwargs = {}
        if progress:
            tqdm_kwargs = {"desc": "Loading Mols", "leave": False}

        results = dm.parallelized(
            self._evaluable_query,
            mols,
            n_jobs=n_jobs,
            scheduler=scheduler,
            progress=progress,
            tqdm_kwargs=tqdm_kwargs,
        )
        return list(results)
