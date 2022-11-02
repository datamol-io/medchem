from typing import Iterable
from typing import List
from typing import Optional
from typing import Dict
from typing import Union
from typing import Any
from typing import Sequence

import os
import numpy as np
import datamol as dm

from functools import partial
from loguru import logger
from rdkit.Chem.rdfiltercatalog import FilterCatalog
from medchem import demerits
from medchem.alerts import AlertFilters
from medchem.alerts import NovartisFilters
from medchem.catalog import NamedCatalogs
from medchem.catalog import merge_catalogs
from medchem.complexity.complexity_filter import ComplexityFilter
from medchem.groups import ChemicalGroup
from medchem.rules import RuleFilters


def alert_filter(
    mols: Iterable[Union[str, dm.Mol]],
    alerts: List[str],
    alerts_db: Optional[os.PathLike] = None,
    n_jobs: Optional[int] = 1,
    rule_dict: Dict = None,
    return_idx: bool = False,
):
    r"""Filter a dataset of molecules, based on common structural alerts and specific rules.

    Arguments:
        mols: List of molecules to filter
        alerts: List of alert collections to screen for. See AlertFilters.list_default_available_alerts()
        alerts_db: Path to the alert file name.
            The internal default file (alerts.csv) will be used if not provided
        n_jobs: Number of cpu to use
        rule_dict: Dictionary with additional rules to apply during the filtering.
            For example, such dictionary for drug-like compounds would look like this:
            >>> rule_dict
                {"MW": [0, 500], "LogP": [-0.5, 5], "HBD": [0, 5], "HBA": [0, 10], "TPSA": [0, 150]}
        return_idx: Whether to return the filtered index

    Returns:
        filtered_mask: boolean array (or index array) where true means
            the molecule IS OK (not found in the alert catalog).
    """

    custom_filters = AlertFilters(alerts_set=alerts, alerts_db=alerts_db)
    df = custom_filters(mols, n_jobs=n_jobs, progress=False)
    df = df[df.status != "Exclude"]
    if rule_dict is not None and len(rule_dict) > 0:
        df = df[
            (df.MW.between(*rule_dict.get("MW", [-np.inf, np.inf])))
            & (df.LogP.between(*rule_dict.get("LogP", [-np.inf, np.inf])))
            & (df.HBD.between(*rule_dict.get("HBD", [-np.inf, np.inf])))
            & (df.HBA.between(*rule_dict.get("HBA", [-np.inf, np.inf])))
            & (df.TPSA.between(*rule_dict.get("TPSA", [-np.inf, np.inf])))
        ]

    filtered_idx = df.index.values.astype(int)
    filtered_mask = np.zeros(len(mols), dtype=bool)
    filtered_mask[filtered_idx] = True
    if return_idx:
        return filtered_idx
    return filtered_mask


def screening_filter(
    mols: Iterable[Union[str, dm.Mol]],
    n_jobs: Optional[int] = None,
    max_severity: int = 10,
    return_idx: bool = False,
):
    """
    Filter a set of molecules based on novartis screening deck curation process
    Schuffenhauer, A. et al. Evolution of Novartis' small molecule screening deck design, J. Med. Chem. (2020)
    DOI. https://dx.doi.org/10.1021/acs.jmedchem.0c01332

    !!! note
        The severity argument corresponds to the accumulated severity for a compounds accross all pattern in the
        catalog.
    Args:
        mols: list of input molecules
        n_jobs: number of parallel job to run. Sequential by default
        max_severity: maximum severity allowed. Default is <10
        return_idx: Whether to return the filtered index

    Returns:
        filtered_mask: boolean array (or index array) where true means the molecule
            IS NOT REJECTED (i.e not found in the alert catalog).

    """

    filt_obj = NovartisFilters()
    df = filt_obj(mols, n_jobs=n_jobs)
    df = df[(df.status != "Exclude") & (df.severity < max_severity)]
    filtered_idx = df.index.values
    filtered_mask = np.zeros(len(mols), dtype=bool)
    filtered_mask[filtered_idx] = True
    if return_idx:
        return filtered_idx
    return filtered_mask


def catalog_filter(
    mols: Sequence[Union[str, dm.Mol]],
    catalogs: List[Union[str, FilterCatalog]],
    return_idx: bool = False,
    n_jobs: Optional[int] = None,
    progress: bool = False,
    scheduler: str = "processes",
    batch_size: int = 100,
):
    """Filter a list of compounds according to catalog of structures alerts and patterns

    Args:
        mols: list of input molecules
        catalogs: list of catalogs (name or FilterCatalog)
        return_idx: whether to return index or a boolean mask
        n_jobs: number of parallel job to run. Sequential by default
        progress: whether to show progress bar
        scheduler: joblib scheduler to use
        batch_size: batch size for parallel processing. Note that `batch_size` should be
            increased if the number of used CPUs gets very large.

    Returns:
        filtered_mask: boolean array (or index array) where true means the molecule is not found in the catalog.
    """

    # Build and merge the catalogs
    named_catalogs = []
    for catalog in catalogs:
        if catalog == "bredt":
            logger.warning(
                "It is not recommended to use the 'bredt' catalog here. Please use the `bredt_filter` function instead or be sure to use kekulized molecules as inputs."
            )
        if catalog == "nibr":
            raise ValueError(
                "You shouldn't use the nibr catalog here. Please use the `screening_filter` function instead."
            )
        elif catalog == "bredt-kekulized":
            catalog = "bredt"
        if isinstance(catalog, str):
            catalog_fn = getattr(NamedCatalogs, catalog, None)
            if catalog_fn is None:
                logger.warning(f"Catalog {catalog} not found. Ignoring.")
            else:
                named_catalogs.append(catalog_fn())
        else:
            named_catalogs.append(catalog)
    if len(named_catalogs) < 1:
        raise ValueError("Please provide at least one catalog !")

    catalog = merge_catalogs(*named_catalogs)

    # Serialize the catalog
    catalog_state = catalog.Serialize()

    def _fn(mols_chunk):
        # Init the catalog from the serialized state
        catalog = FilterCatalog(catalog_state)

        # To mols
        mols_chunk = [dm.to_mol(m) for m in mols_chunk]

        # Match the mols
        return [catalog.HasMatch(m) for m in mols_chunk]

    # Batch the inputs
    n_batches = len(mols) // batch_size
    n_batches = max(n_batches, 1)
    mols_batches = np.array_split(mols, n_batches)

    # Run the matching
    toxic = dm.parallelized(
        _fn,
        mols_batches,
        n_jobs=n_jobs,
        scheduler=scheduler,
        progress=progress,
        tqdm_kwargs=dict(desc="Match", leave=True),
    )

    # Flatten the results
    toxic = [item for sublist in toxic for item in sublist]

    filtered_idx = [i for i, bad in enumerate(toxic) if not bad]
    if return_idx:
        return np.asarray(filtered_idx)
    return np.bitwise_not(toxic)


def chemical_group_filter(
    mols: Iterable[Union[str, dm.Mol]],
    chemical_group: ChemicalGroup,
    return_idx: bool = False,
    n_jobs: Optional[int] = None,
    progress: bool = False,
    scheduler: str = "threads",
):
    """Filter a list of compounds according to a chemical group instance.

    !!! note
        This function will return the list of molecules that DO NOT match the chemical group

    Args:
        mols: list of input molecules
        chemical_group: a chemical group instance with the required functional groups to use.
        return_idx: whether to return index or a boolean mask
        n_jobs: number of parallel job to run. Sequential by default
        progress: whether to show progress bar
        scheduler: joblib scheduler to use

    Returns:
        filtered_mask: boolean array (or index array) where true means the molecule DOES NOT MATCH the groups.
    """

    if isinstance(chemical_group, ChemicalGroup):
        chemical_group = chemical_group.get_catalog()
    return catalog_filter(
        mols,
        [chemical_group],
        return_idx=return_idx,
        n_jobs=n_jobs,
        progress=progress,
        scheduler=scheduler,
    )


def rules_filter(
    mols: Iterable[Union[str, dm.Mol]],
    rules: Union[List[Any], RuleFilters],
    return_idx: bool = False,
    n_jobs: Optional[int] = None,
    progress: bool = False,
    scheduler: str = "processes",
):
    """Filter a list of compounds according to a predefined set of rules

    Args:
        mols: list of input molecules
        rules: list of rules to apply to the input molecules.
        return_idx: whether to return index or a boolean mask
        n_jobs: number of parallel job to run. Sequential by default
        progress: whether to show progress bar
        scheduler: joblib scheduler to use

    Returns:
        filtered_mask: boolean array (or index array) where true means the molecule MATCH the rules.
    """

    if not isinstance(rules, RuleFilters):
        rules = RuleFilters(rules, precompute_props=True)
    df = rules(mols, n_jobs=n_jobs, progress=progress, scheduler=scheduler)
    filtered_df = df.all(axis=1, bool_only=True)
    if return_idx:
        return filtered_df.index.values[filtered_df.values]
    return filtered_df.values


def complexity_filter(
    mols: Iterable[Union[str, dm.Mol]],
    complexity_metric: str = "bertz",
    threshold_stats_file: str = "zinc_15_available",
    limit: str = "99",
    return_idx: bool = False,
    n_jobs: Optional[int] = None,
    progress: bool = False,
    scheduler: str = "processes",
):
    """Filter a list of compounds according to a chemical group instance

    Args:
        mols: list of input molecules
        complexity_metric: complexity metric to use
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
        threshold_stats_file: complexity threshold statistic origin to use
        limit: complexity outlier percentile to use
        return_idx: whether to return index or a boolean mask
        n_jobs: number of parallel job to run. Sequential by default
        progress: whether to show progress bar
        scheduler: joblib scheduler to use

    Also see:
        medchem.complexity.ComplexityFilter
    Returns:
        filtered_mask: boolean array (or index array) where true means the molecule MATCH the rules.
    """

    cf = ComplexityFilter(
        limit=limit,
        complexity_metric=complexity_metric,
        threshold_stats_file=threshold_stats_file,
    )

    mols = dm.parallelized(
        dm.to_mol,
        mols,
        n_jobs=n_jobs,
        progress=progress,
        tqdm_kwargs=dict(desc="To mol", leave=False),
    )

    not_complex = dm.parallelized(
        cf,
        mols,
        n_jobs=n_jobs,
        scheduler=scheduler,
        progress=progress,
        tqdm_kwargs=dict(desc="Complexity Eval", leave=False),
    )
    not_complex = np.asarray(not_complex)
    filtered_idx = [i for i, good in enumerate(not_complex) if good]
    if return_idx:
        return np.asarray(filtered_idx)
    return not_complex


def bredt_filter(
    mols: Sequence[Union[str, dm.Mol]],
    return_idx: bool = False,
    n_jobs: Optional[int] = None,
    progress: bool = False,
    scheduler: str = "threads",
    batch_size: int = 100,
):
    """Filter a list of compounds according to Bredt's rules
    https://en.wikipedia.org/wiki/Bredt%27s_rule

    Args:
        mols: list of input molecules
        return_idx: whether to return index or a boolean mask
        n_jobs: number of parallel job to run. Sequential by default
        progress: whether to show progress bar
        scheduler: joblib scheduler to use
        batch_size: batch size for parallel processing. Note that `batch_size` should be
            increased if the number of used CPUs gets very large.

    Returns:
        filtered_mask: boolean array (or index array) where true means the molecule is not toxic.
    """

    mols = dm.parallelized(
        partial(dm.to_mol, kekulize=True),
        mols,
        n_jobs=n_jobs,
        progress=progress,
        tqdm_kwargs=dict(desc="To mol", leave=False),
    )

    return catalog_filter(
        mols=mols,
        catalogs=["bredt-kekulized"],  # already kekulized mols
        return_idx=return_idx,
        n_jobs=n_jobs,
        progress=progress,
        scheduler=scheduler,
        batch_size=batch_size,
    )


def molecular_graph_filter(
    mols: Iterable[Union[str, dm.Mol]],
    max_severity: int = 5,
    return_idx: bool = False,
    n_jobs: Optional[int] = None,
    progress: bool = False,
    scheduler: str = "threads",
):
    """Filter a list of compounds according to unstable molecular graph filter list.

    This list was obtained from observation around The disallowed graphs are:

    * K3,3 or K2,4 structure
    * Cone of P4 or K4 with 3-ear
    * Node in more than one ring of length 3 or 4

    Args:
        mols: list of input molecules
        max_severity: maximum acceptable severity (1-10). Default is <5
        return_idx: whether to return index or a boolean mask
        n_jobs: number of parallel job to run. Sequential by default
        progress: whether to show progress bar
        scheduler: joblib scheduler to use

    Returns:
        filtered_mask: boolean array (or index array) where true means the molecule is not toxic.
    """
    if max_severity is None:
        max_severity = 5
    catalog = NamedCatalogs.unstable_graph(max_severity=max_severity)
    if isinstance(mols[0], str):
        mols = dm.parallelized(
            dm.to_mol,
            mols,
            n_jobs=n_jobs,
            progress=progress,
            tqdm_kwargs=dict(desc="To mol", leave=False),
        )
    toxic = dm.parallelized(
        catalog.HasMatch,
        mols,
        n_jobs=n_jobs,
        scheduler=scheduler,
        progress=progress,
        tqdm_kwargs=dict(desc="Match", leave=False),
    )
    filtered_idx = [i for i, bad in enumerate(toxic) if not bad]
    if return_idx:
        return filtered_idx
    return np.bitwise_not(toxic)


def lilly_demerit_filter(
    smiles: Iterable[str],
    max_demerits: Optional[int] = 160,
    return_idx: bool = False,
    n_jobs: Optional[int] = None,
    progress: bool = False,
    **kwargs,
):
    """Run Lilly demerit filtering on current list of molecules

    Args:
        smiles: list of input molecules as smiles preferably
        max_demerits: Cutoff to reject molecules Defaults to 160.
        return_idx: whether to return a mask or a list of valid indexes
        progress: whether to show progress bar
        kwargs: parameters specific to the `demerits.score` function

    Returns:
        filtered_mask: boolean array (or index array) where true means the molecule is ok.
    """

    if not isinstance(smiles[0], str):
        # canonicalize the smiles and ensure the input is smiles and not Chem.Mol
        def canonical_smi(smi):
            mol = dm.to_mol(smi)
            if mol is not None:
                return dm.to_smiles(mol)
            return mol

        smiles = dm.parallelized(
            canonical_smi,
            smiles,
            n_jobs=n_jobs,
            progress=progress,
            tqdm_kwargs=dict(desc="Canonical Smiles", leave=False),
        )

    df = demerits.batch_score(smiles, n_jobs=n_jobs, progress=progress, **kwargs)
    df = df[
        (~df.rejected) & ((df.demerit_score.isna()) | (df.demerit_score < max_demerits))
    ]

    filtered_idx = df["ID"].values.astype(int)
    filtered_mask = np.zeros(len(smiles), dtype=bool)
    filtered_mask[filtered_idx] = True
    if return_idx:
        return filtered_idx
    return filtered_mask


def protecting_groups_filter(
    mols: Iterable[Union[str, dm.Mol]],
    return_idx: bool = False,
    protecting_groups: str = [
        "fmoc",
        "tert-butoxymethyl",
        "tert-butyl carbamate",
        "tert-butyloxycarbonyl",
    ],
    n_jobs: Optional[int] = None,
    progress: bool = False,
    scheduler: str = "threads",
):
    """Filter a list of compounds according to match to  known protecting groups.
    Note that is a syntaxic sugar for calling chemical_group_filter with the protecting groups subset

    Args:
        mols: list of input molecules
        protecting_groups: type of protection group to consider if not provided, will use all (not advised)
        return_idx: whether to return index or a boolean mask
        n_jobs: number of parallel job to run. Sequential by default
        progress: whether to show progress bar
        scheduler: joblib scheduler to use

    Returns:
        filtered_mask: boolean array (or index array) where true means the molecule DOES NOT MATCH the groups.
    """

    chemical_group = ChemicalGroup("protecting_groups")
    chemical_group = chemical_group.filter(protecting_groups)
    return chemical_group_filter(
        mols,
        chemical_group,
        return_idx=return_idx,
        n_jobs=n_jobs,
        progress=progress,
        scheduler=scheduler,
    )
