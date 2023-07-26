from typing import List
from typing import Optional
from typing import Union
from typing import Any
from typing import Sequence

import os
import numpy as np
import datamol as dm

from functools import partial
from loguru import logger

from rdkit.Chem import rdfiltercatalog  # type: ignore

from medchem.structural import CommonAlertsFilters
from medchem.structural import NIBRFilters
from medchem.catalogs import NamedCatalogs
from medchem.catalogs import merge_catalogs
from medchem.complexity import ComplexityFilter
from medchem.groups import ChemicalGroup
from medchem.rules import RuleFilters


def alert_filter(
    mols: Sequence[Union[str, dm.Mol]],
    alerts: List[str],
    alerts_db: Optional[Union[os.PathLike, str]] = None,
    n_jobs: Optional[int] = 1,
    progress: bool = False,
    return_idx: bool = False,
) -> np.ndarray:
    """Filter a dataset of molecules, based on common structural alerts and specific rules.

    ??? tip "True is good"
        Returning `True` means the molecule does not match any of the structural alerts.

    !!! info "See Also"
        `alert_filter` is a convenient functional API for the `medchem.structural.CommonAlertsFilters` class.

    Args:
        mols: List of molecules to filter
        alerts: List of alert collections to screen for. See `CommonAlertsFilters.list_default_available_alerts()`
        alerts_db: Path to the alert file name.
            The internal default file (`alerts.csv`) will be used if not provided
        n_jobs: Number of workers to use
        progress: Whether to show progress bar
        return_idx: Whether to return the filtered index

    Returns:
        filtered_mask: boolean array (or index array) where true means
            the molecule IS OK (not found in the alert catalog).

    """

    custom_filters = CommonAlertsFilters(alerts_set=alerts, alerts_db_path=alerts_db)
    results = custom_filters(mols, n_jobs=n_jobs, progress=progress)

    if return_idx:
        return results.query("pass_filter == True").index.values.astype(int)

    return results["pass_filter"].to_numpy()


def nibr_filter(
    mols: Sequence[Union[str, dm.Mol]],
    n_jobs: Optional[int] = None,
    max_severity: int = 10,
    progress: bool = False,
    return_idx: bool = False,
) -> np.ndarray:
    """
    Filter a set of molecules based on the Novartis Institutes for BioMedical Research screening deck curation process
    [Schuffenhauer, A. et al. Evolution of Novartis' small molecule screening deck design, J. Med. Chem. (2020)](https://dx.doi.org/10.1021/acs.jmedchem.0c01332)


    The severity argument corresponds to the accumulated severity for a compounds accross all pattern in the catalog.

    ??? tip "True is good"
        Returning `True` means the molecule does not match any of the structural alerts.

    !!! info "See Also"
        `nibr_filter` is a convenient functional API for the `medchem.structural.NIBRFilters` class.

    Args:
        mols: list of input molecules
        n_jobs: number of parallel job to run. Sequential by default
        max_severity: maximum severity allowed. Default is <10
        progress: whether to show progress bar
        return_idx: Whether to return the filtered index

    Returns:
        filtered_mask: boolean array (or index array) where true means the molecule
            IS NOT REJECTED (i.e not found in the alert catalog).

    """

    filt_obj = NIBRFilters()
    results = filt_obj(mols, n_jobs=n_jobs, progress=progress)

    mask = results["severity"] < max_severity

    if return_idx:
        return results[mask].index.values.astype(int)

    return mask.to_numpy()


def catalog_filter(
    mols: Sequence[Union[str, dm.Mol]],
    catalogs: List[Union[str, rdfiltercatalog.FilterCatalog]],
    return_idx: bool = False,
    n_jobs: Optional[int] = -1,
    progress: bool = False,
    progress_leave: bool = False,
    scheduler: str = "processes",
    batch_size: int = 100,
) -> np.ndarray:
    """Filter a list of compounds according to a catalog of structural alerts and patterns


    ??? tip "True is good"
        Returning `True` means the molecule does not match any of the structural alerts.

    Args:
        mols: list of input molecules
        catalogs: list of catalogs (name or FilterCatalog)
        return_idx: whether to return index or a boolean mask
        n_jobs: number of parallel job to run. Sequential by default
        progress: whether to show progress bar
        progress_leave: whether to leave the progress bar after completion
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
            raise ValueError(
                "It is not recommended to use the 'bredt' catalog here. Please use the `bredt_filter`."
            )
        if catalog == "nibr":
            raise ValueError(
                "You shouldn't use the nibr catalog here. Please use the `screening_filter` function instead."
            )
        if catalog == "bredt-kekulized":
            # we know the compounds have been kekulized
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
        catalog = rdfiltercatalog.FilterCatalog(catalog_state)

        # To mols
        mols_chunk = [dm.to_mol(m) for m in mols_chunk]

        # Match the mols
        return [catalog.HasMatch(m) for m in mols_chunk]

    # Run the matching
    toxic = dm.parallelized_with_batches(
        _fn,
        mols,
        batch_size=batch_size,
        n_jobs=n_jobs,
        scheduler=scheduler,
        progress=progress,
        tqdm_kwargs=dict(
            desc="Filtering with catalogs",
            leave=progress_leave,
        ),
        flatten_results=True,
    )

    toxic = np.asarray(toxic)
    filtered_idx = [i for i, bad in enumerate(toxic) if not bad]

    if return_idx:
        return np.asarray(filtered_idx)

    return np.bitwise_not(toxic)


def chemical_group_filter(
    mols: Sequence[Union[str, dm.Mol]],
    chemical_group: ChemicalGroup,
    exact_match: bool = False,
    return_idx: bool = False,
    n_jobs: Optional[int] = None,
    progress: bool = False,
    progress_leave: bool = False,
    scheduler: str = "threads",
) -> np.ndarray:
    """Filter a list of compounds according to a chemical group instance.

    !!! warning
        This function will return the list of molecules that **DO NOT** match the chemical group.

    !!! info "See Also"
        Consider exploring the `medchem.groups.ChemicalGroup` class.

    Args:
        mols: list of input molecules
        chemical_group: a chemical group instance with the required functional groups to use.
        exact_match: whether to use an exact match of the chemical group patterns (will switch to smiles )
        return_idx: whether to return index or a boolean mask
        n_jobs: number of parallel job to run. Sequential by default
        progress: whether to show progress bar
        progress_leave: whether to leave the progress bar after completion
        scheduler: joblib scheduler to use

    Returns:
        filtered_mask: boolean array (or index array) where true means the molecule DOES NOT MATCH the groups.
    """

    if isinstance(chemical_group, ChemicalGroup):
        chemical_group = chemical_group.get_catalog(exact_match=exact_match)

    return catalog_filter(
        mols,
        [chemical_group],
        return_idx=return_idx,
        n_jobs=n_jobs,
        progress=progress,
        progress_leave=progress_leave,
        scheduler=scheduler,
    )


def rules_filter(
    mols: Sequence[Union[str, dm.Mol]],
    rules: Union[List[Any], RuleFilters],
    return_idx: bool = False,
    n_jobs: Optional[int] = None,
    progress: bool = False,
    progress_leave: bool = False,
    scheduler: str = "processes",
) -> np.ndarray:
    """Filter a list of compounds according to a predefined set of rules

    ??? tip "True is good"
        Returning `True` means the molecule passes all the rules.

    !!! info "See Also"
        Consider exploring the `medchem.rules.RuleFilters` class.

    Args:
        mols: list of input molecules
        rules: list of rules to apply to the input molecules.
        return_idx: whether to return index or a boolean mask
        n_jobs: number of parallel job to run. Sequential by default
        progress: whether to show progress bar
        scheduler: joblib scheduler to use

    Returns:
        filtered_mask: boolean array (or index array) where true means the molecule MATCH the rule constraints.
    """

    if not isinstance(rules, RuleFilters):
        rules = RuleFilters(rules)

    results = rules(
        mols,
        n_jobs=n_jobs,
        progress=progress,
        progress_leave=progress_leave,
        scheduler=scheduler,
    )

    mask = results["pass_all"].to_numpy()

    if return_idx:
        return results[mask].index.to_numpy()

    return mask


def complexity_filter(
    mols: Sequence[Union[str, dm.Mol]],
    complexity_metric: str = "bertz",
    threshold_stats_file: str = "zinc_15_available",
    limit: str = "99",
    return_idx: bool = False,
    n_jobs: Optional[int] = None,
    progress: bool = False,
    progress_leave: bool = False,
    scheduler: str = "processes",
) -> np.ndarray:
    """Filter a list of compounds according to a complexity metric


    ??? tip "True is good"
        Returning `True` means the molecule passes the complexity filters.

    !!! info "See Also"
        Consider exploring the `medchem.complexity.ComplexityFilter` class.

    Args:
        mols: list of input molecules
        complexity_metric: complexity metric to use
            Use `ComplexityFilter.list_default_available_filters` to list default filters.
            The following complexity metrics are supported by default:

            - `bertz`: bertz complexity index
            - `sas`: synthetic accessibility score  (`zinc_15_available` only)
            - `qed`: qed score  (`zinc_15_available` only)
            - `clogp`: clogp for how greasy a molecule is compared to other in the same mw range  (`zinc_15_available` only)
            - `whitlock`: whitlock complexity index
            - `barone`: barone complexity index
            - `smcm`: synthetic and molecular complexity
            - `twc`:  total walk count complexity  (`zinc_15_available` only)

        threshold_stats_file: complexity threshold statistics file to use
        limit: complexity outlier percentile to use
        return_idx: whether to return index or a boolean mask
        n_jobs: number of parallel job to run. Sequential by default
        progress: whether to show progress bar
        progress_leave: whether to leave the progress bar after completion
        scheduler: joblib scheduler to use

    Returns:
        filtered_mask: boolean array (or index array) where true means the molecule MATCH the rules.
    """

    cf = ComplexityFilter(
        limit=limit,
        complexity_metric=complexity_metric,
        threshold_stats_file=threshold_stats_file,
    )

    not_complex = dm.parallelized(
        cf,
        mols,
        n_jobs=n_jobs,
        scheduler=scheduler,
        progress=progress,
        tqdm_kwargs=dict(
            desc="Complexity Evaluation",
            leave=progress_leave,
        ),
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
    progress_leave: bool = False,
    scheduler: str = "threads",
    batch_size: int = 100,
) -> np.ndarray:
    """Filter a list of compounds according to [Bredt's rules](https://en.wikipedia.org/wiki/Bredt%27s_rule)

    ??? tip "True is good"
        Returning `True` means the molecule does not violate the Bredt's rules.

    Args:
        mols: list of input molecules
        return_idx: whether to return index or a boolean mask
        n_jobs: number of parallel job to run. Sequential by default
        progress: whether to show progress bar
        progress_leave: whether to leave the progress bar after completion
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
        progress_leave=progress_leave,
        scheduler=scheduler,
        batch_size=batch_size,
    )


def molecular_graph_filter(
    mols: Sequence[Union[str, dm.Mol]],
    max_severity: Optional[int] = 5,
    return_idx: bool = False,
    n_jobs: Optional[int] = None,
    progress: bool = False,
    progress_leave: bool = False,
    scheduler: str = "threads",
) -> np.ndarray:
    """Filter a list of compounds according to unstable molecular graph patterns.
    This list was obtained from observation around technically valid molecular graphs from
    deep generative models that are not stable.

    The disallowed graphs are:

    -  K3,3 or K2,4 structures
    -  Cone of P4 or K4 with 3-ear
    -  Node in more than one ring of length 3 or 4

    ??? tip "True is good"
        Returning `True` means the molecule does not violate the molecular graph instability rules.

    Args:
        mols: list of input molecules
        max_severity: maximum acceptable severity (1-10). Default is <5
        return_idx: whether to return index or a boolean mask
        n_jobs: number of parallel job to run. Sequential by default
        progress: whether to show progress bar
        progress_leave: whether to leave the progress bar after completion
        scheduler: joblib scheduler to use

    Returns:
        filtered_mask: boolean array (or index array) where true means the molecule is not toxic.
    """

    if max_severity is None:
        max_severity = 5

    catalog = NamedCatalogs.unstable_graph(severity_threshold=max_severity)

    if isinstance(mols[0], str):
        mols = dm.parallelized(
            dm.to_mol,
            mols,
            n_jobs=n_jobs,
            progress=progress,
            tqdm_kwargs=dict(desc="To mol", leave=progress_leave),
        )

    toxic = dm.parallelized(
        catalog.HasMatch,
        mols,
        n_jobs=n_jobs,
        scheduler=scheduler,
        progress=progress,
        tqdm_kwargs=dict(desc="Match", leave=progress_leave),
    )
    toxic = np.asarray(toxic)
    filtered_idx = np.where(~toxic)[0]

    if return_idx:
        return filtered_idx

    return np.bitwise_not(toxic)


def lilly_demerit_filter(
    mols: Sequence[Union[str, dm.Mol]],
    max_demerits: Optional[int] = 160,
    return_idx: bool = False,
    n_jobs: Optional[int] = None,
    progress: bool = False,
    progress_leave: bool = False,
    scheduler: str = "threads",
    batch_size: int = 5_000,
    **kwargs: Any,
) -> np.ndarray:
    """Run the Eli Lilly's demerit filter on current list of molecules

    ??? tip "True is good"
        Returning `True` means the molecule does not violate the demerit rules.

    !!! info "See Also"
        Consider exploring the `LillyDemeritsFilters` class in `medchem.structural.lilly_demerits`

    Args:
        mols: list of input molecules as smiles preferably
        max_demerits: Cutoff to reject molecules Defaults to 160.
        return_idx: whether to return a mask or a list of valid indexes
        progress: whether to show progress bar
        progress_leave: whether to leave the progress bar after completion
        scheduler: joblib scheduler to usescheduler
        batch_size: batch size for parallel processing.
        kwargs: parameters specific to the `demerits.score` function

    Returns:
        filtered_mask: boolean array (or index array) where true means the molecule is ok.
    """

    from medchem.structural.lilly_demerits import LillyDemeritsFilters

    dfilter = LillyDemeritsFilters()

    results = dfilter(
        mols,
        n_jobs=n_jobs,
        progress=progress,
        progress_leave=progress_leave,
        scheduler=scheduler,
        batch_size=batch_size,
        **kwargs,
    )

    results = results[
        (results["pass_filter"] == True)  # noqa: E712
        & ((results["demerit_score"].isna()) | (results["demerit_score"] < max_demerits))
    ]

    filtered_idx = results.index.values.astype(int)
    filtered_mask = np.zeros(len(mols), dtype=bool)
    filtered_mask[filtered_idx] = True

    if return_idx:
        return filtered_idx

    return filtered_mask


def protecting_groups_filter(
    mols: Sequence[Union[str, dm.Mol]],
    return_idx: bool = False,
    protecting_groups: List[str] = [
        "fmoc",
        "tert-butoxymethyl",
        "tert-butyl carbamate",
        "tert-butyloxycarbonyl",
    ],
    n_jobs: Optional[int] = None,
    progress: bool = False,
    progress_leave: bool = False,
    scheduler: str = "threads",
) -> np.ndarray:
    """Filter a list of compounds according to match to known protecting groups.


    !!! warning
         This function will return the list of molecules that **DO NOT** have the protecting groups.

     !!! info "See Also"
         This is a syntaxic sugar for calling chemical_group_filter with the protecting groups subset.

     Args:
         mols: list of input molecules
         protecting_groups: type of protection group to consider if not provided, will use all **(not advised)**
         return_idx: whether to return index or a boolean mask
         n_jobs: number of parallel job to run. Sequential by default
         progress: whether to show progress bar
         progress_leave: whether to leave the progress bar after completion
         scheduler: joblib scheduler to use

     Returns:
         filtered_mask: boolean array (or index array) where true means the molecule DOES NOT MATCH the groups.
    """

    chemical_group = ChemicalGroup("protecting_groups")
    if protecting_groups and len(protecting_groups) > 0:
        chemical_group = chemical_group.filter(protecting_groups)

    protected = dm.parallelized(
        partial(chemical_group.has_match, exact_match=True, terminal_only=True),
        mols,
        n_jobs=n_jobs,
        progress=progress,
        scheduler=scheduler,
        tqdm_kwargs={"desc": "Checking protecting groups", "leave": progress_leave},
    )
    filtered_idx = [i for i, bad in enumerate(protected) if not bad]

    if return_idx:
        return np.asarray(filtered_idx)

    return np.bitwise_not(protected)
