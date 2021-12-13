from functools import partial
from typing import Callable, Iterable
from typing import List
from typing import Optional
from typing import Dict
from typing import Union

import os
import numpy as np
import datamol as dm

from rdkit.Chem import rdchem
from medchem.catalog import NamedCatalogs
from medchem import demerits
from medchem.alerts import AlertFilters
from medchem.novartis import NovartisFilters
from joblib.externals.loky import set_loky_pickler
from joblib import parallel_backend
from joblib import Parallel, delayed
from joblib import wrap_non_picklable_objects


def alert_filter(
    mols: Iterable[Union[str, rdchem.Mol]],
    alerts: List[str],
    alerts_db: Optional[os.PathLike] = None,
    n_jobs: Optional[int] = 1,
    rule_dict: Dict = None,
    return_idx: bool = False,
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

    filt_obj = NovartisFilters()
    df = filt_obj(mols, n_jobs=n_jobs)
    df = df[(df.status != "Exclude") & (df.severity < max_severity)]
    filtered_idx = df.index.values
    filtered_mask = np.zeros(len(mols), dtype=bool)
    filtered_mask[filtered_idx] = True
    if return_idx:
        return filtered_idx
    return filtered_mask


def _apply_catalog(
    mols: Iterable[rdchem.Mol],
    match_fn: Callable,
    n_jobs: Optional[int] = None,
    scheduler: str = "threads",
):
    """Apply a catalog to a list of molecules and return a boolean mask of the results

    Args:
        mols: list of input molecules, needs to be a list of rdkit.Chem.rdchem.Mol
        catalog: catalog to apply
        n_jobs: number of parallel job to run. Sequential by default
        scheduler: joblib scheduler to use. Default is threads
    """

    if n_jobs is not None:
        toxic = dm.parallelized(match_fn, mols, n_jobs=n_jobs, scheduler=scheduler)
    else:
        toxic = [match_fn(mol) for mol in mols]

    return toxic


def common_filter(
    mols: Iterable[Union[str, rdchem.Mol]],
    pains_a: bool = True,
    pains_b: bool = True,
    pains_c: bool = False,
    brenk: bool = False,
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
        return_idx: whether to return index or a boolean mask

    Returns:
        filtered_mask: boolean array (or index array) where true means the molecule is not toxic.
    """

    if pains:
        pains_a = pains_b = pains_c = True

    catalog = NamedCatalogs.tox(
        pains_a=pains_a,
        pains_b=pains_b,
        pains_c=pains_c,
        brenk=brenk,
        zinc=zinc,
        nih=nih,
    )

    mols = dm.parallelized(dm.to_mol, mols, n_jobs=n_jobs)
    toxic = _apply_catalog(mols, catalog.HasMatch, n_jobs=n_jobs)

    filtered_idx = [i for i, bad in enumerate(toxic) if not bad]
    if return_idx:
        return filtered_idx
    return np.bitwise_not(toxic)


def bredt_filter(
    mols,
    n_jobs: Optional[int] = None,
    return_idx: bool = False,
):
    """Filter a list of compounds according to Bredt's rules
    https://en.wikipedia.org/wiki/Bredt%27s_rule

    Args:
        mols: list of input molecules
        n_jobs: number of parallel job to run. Sequential by default
        return_idx: whether to return index or a boolean mask

    Returns:
        filtered_mask: boolean array (or index array) where true means the molecule is not toxic.
    """

    catalog = NamedCatalogs.bredt()
    # molecules needs to be kekulized for bredt
    mols = dm.parallelized(partial(dm.to_mol, kekulize=True), mols, n_jobs=n_jobs)

    toxic = _apply_catalog(mols, catalog.HasMatch, n_jobs=n_jobs)
    filtered_idx = [i for i, bad in enumerate(toxic) if not bad]
    if return_idx:
        return filtered_idx
    return np.bitwise_not(toxic)


def molecular_graph_filter(
    mols,
    n_jobs: Optional[int] = None,
    max_severity: int = 5,
    return_idx: bool = False,
):
    """Filter a list of compounds according to Unstable molecular graph filter list

    Args:
        mols: list of input molecules
        max_severity: maximum severity allowed. Default is <=5
        n_jobs: number of parallel job to run. Sequential by default
        return_idx: whether to return index or a boolean mask

    Returns:
        filtered_mask: boolean array (or index array) where true means the molecule is not toxic.
    """
    if max_severity is None:
        max_severity = 5
    catalog = NamedCatalogs.unstable_graph(max_severity=max_severity)
    mols = dm.parallelized(dm.to_mol, mols, n_jobs=n_jobs)
    toxic = _apply_catalog(mols, catalog, n_jobs=n_jobs)
    filtered_idx = [i for i, bad in enumerate(toxic) if not bad]
    if return_idx:
        return filtered_idx
    return np.bitwise_not(toxic)


def demerit_filter(
    mols,
    max_demerits: Optional[int] = 160,
    return_idx: bool = False,
    n_jobs: Optional[int] = None,
    **kwargs
):
    """Run demerit filtering on current list of molecules

    Args:
        mols: list of input molecules
        max_demerits: Cutoff to reject molecules Defaults to 160.
        return_idx: whether to return a mask or a list of valid indexes
        kwargs: parameters specific to the `demerits.score` function

    Returns:
        filtered_mask: boolean array (or index array) where true means the molecule is ok.
    """

    if not isinstance(mols[0], str):
        # canonicalize the smiles and ensure the input is smiles and not Chem.Mol
        mols = dm.parallelized(dm.to_mol, mols, n_jobs=n_jobs)
        mols = dm.parallelized(dm.to_smiles, mols, n_jobs=n_jobs)
    df = demerits.batch_score(mols, n_jobs=n_jobs, **kwargs)
    df = df[
        (~df.rejected) & ((df.demerit_score.isna()) | (df.demerit_score < max_demerits))
    ]

    filtered_idx = df["ID"].values.astype(int)
    filtered_mask = np.zeros(len(mols), dtype=bool)
    filtered_mask[filtered_idx] = True
    if return_idx:
        return filtered_idx
    return filtered_mask
