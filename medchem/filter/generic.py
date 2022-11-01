from typing import Iterable
from typing import Optional
from typing import Callable
from typing import Union

import itertools
import numpy as np
import datamol as dm
from medchem.utils.graph import score_symmetry


def _generic_filter(
    mols: Iterable[Union[str, dm.Mol]],
    rejection_fn: Callable,
    return_idx: bool = False,
    n_jobs: Optional[int] = None,
    progress: bool = False,
    scheduler: Optional[str] = None,
):
    """Apply any generic filter to a molecule
    Args:
        mols: list of input molecules
        rejection_fn: function to reject a molecule. Return True is the molecule is rejected
        return_idx: whether to return index or a boolean mask
        n_jobs: number of parallel job to run. Sequential by default
        progress: whether to show progress bar
        scheduler: joblib scheduler to use

    Returns:
        filtered_mask: boolean array (or index array) where true means the molecule is ok.
    """
    if isinstance(mols[0], str):
        mols = dm.parallelized(
            dm.to_mol,
            mols,
            n_jobs=n_jobs,
            progress=progress,
            tqdm_kwargs=dict(desc="Macrocyle match", leave=False),
        )
    toxic = dm.parallelized(
        rejection_fn,
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


def macrocycle_filter(
    mols: Iterable[Union[str, dm.Mol]],
    max_cycle_size: int = 10,
    return_idx: bool = False,
    n_jobs: Optional[int] = None,
    progress: bool = False,
    scheduler: Optional[str] = None,
):
    """Find molecules that do not infringe the strict maximum cycle size.

    Args:
        mols: list of input molecules
        max_cycle_size: strict maximum macrocycle size
        return_idx: whether to return index or a boolean mask
        n_jobs: number of parallel job to run. Sequential by default
        progress: whether to show progress bar
        scheduler: joblib scheduler to use

    Returns:
        filtered_mask: boolean array (or index array) where true means the molecule is ok.

    """

    def reject_fn(mol):
        if mol is None:
            return True
        rinfo = mol.GetRingInfo().AtomRings()
        return any(len(r) >= max_cycle_size for r in rinfo)

    return _generic_filter(
        mols,
        reject_fn,
        return_idx=return_idx,
        n_jobs=n_jobs,
        progress=progress,
        scheduler=scheduler,
    )


def atom_list_filter(
    mols: Iterable[Union[str, dm.Mol]],
    unwanted_atom_list: Optional[Iterable] = None,
    wanted_atom_list: Optional[Iterable] = None,
    return_idx: bool = False,
    n_jobs: Optional[int] = None,
    progress: bool = False,
    scheduler: Optional[str] = None,
):
    """Find molecule without any atom from a set of unwanted atom symbols
    and with all atoms in the set of desirable atom list

    Args:
        mols: list of input molecules
        unwanted_atom_list: list of undesirable atom symbol
        wanted_atom_list: list of desirable atom symbol
        return_idx: whether to return index or a boolean mask
        n_jobs: number of parallel job to run. Sequential by default
        progress: whether to show progress bar
        scheduler: joblib scheduler to use

    Returns:
        filtered_mask: boolean array (or index array) where true means the molecule is ok.
    """
    unwanted_atom_list = unwanted_atom_list or []
    wanted_atom_list = wanted_atom_list or []

    def reject_fn(mol):
        if mol is None:
            return True
        for atom in mol.GetAtoms():
            cur_symb = atom.GetSymbol()
            if cur_symb in unwanted_atom_list or (
                len(wanted_atom_list) > 0 and cur_symb not in wanted_atom_list
            ):
                return True
        return False

    return _generic_filter(
        mols,
        reject_fn,
        return_idx=return_idx,
        n_jobs=n_jobs,
        progress=progress,
        scheduler=scheduler,
    )


def ring_infraction_filter(
    mols: Iterable[Union[str, dm.Mol]],
    hetcycle_min_size: int = 4,
    return_idx: bool = False,
    n_jobs: Optional[int] = None,
    progress: bool = False,
    scheduler: Optional[str] = None,
):
    """
    Find molecules that have a ring infraction filter.
    Returning True means the molecule is fine

    Args:
        mols: list of input molecules
        hetcycle_min_size: Minimum ring size before more than 1 hetero atom or any non single bond is allowed.
            This is a *strict threshold (>)*
        return_idx: whether to return index or a boolean mask
        n_jobs: number of parallel job to run. Sequential by default
        progress: whether to show progress bar
        scheduler: joblib scheduler to use

    Returns:
        filtered_mask: boolean array (or index array) where true means the molecule is ok.
    """

    ring_allene = dm.from_smarts("[R]=[R]=[R]")
    double_bond_small_ring = dm.from_smarts("[r3,r4]=[r3,r4]")

    def reject_fn(mol):
        if mol is None:
            return True
        rejected = mol.HasSubstructMatch(ring_allene) or mol.HasSubstructMatch(
            double_bond_small_ring
        )
        if rejected:
            return True
        rinfo = mol.GetRingInfo()
        bond_rings = rinfo.BondRings()
        for r in bond_rings:
            r_bonds = [mol.GetBondWithIdx(b) for b in r]
            r_bond_types = [b.GetBondType() for b in r_bonds]
            r_atom_content = set(
                itertools.chain(
                    *[
                        (b.GetBeginAtom().GetSymbol(), b.GetEndAtom().GetSymbol())
                        for b in r_bonds
                    ]
                )
            )
            n_ring_heteroatoms = sum([at not in ["C", "H"] for at in r_atom_content])
            if len(r) <= hetcycle_min_size and (
                n_ring_heteroatoms > 1
                or any(btype != dm.SINGLE_BOND for btype in r_bond_types)
            ):
                # too many heteroatoms in low ring size OR
                # alkene, allen or aromatic in small rings
                return True
        return False

    return _generic_filter(
        mols,
        reject_fn,
        return_idx=return_idx,
        n_jobs=n_jobs,
        progress=progress,
        scheduler=scheduler,
    )


def num_atom_filter(
    mols: Iterable[Union[str, dm.Mol]],
    min_atoms: Optional[int] = None,
    max_atoms: Optional[int] = None,
    return_idx: bool = False,
    n_jobs: Optional[int] = None,
    progress: bool = False,
    scheduler: Optional[str] = None,
):
    """
    Find a molecule that match the atom number constraints
    Returning True means the molecule is fine

    Args:
        mols: list of input molecules
        min_atoms: strict minimum number of atoms (atoms > min_atoms)
        max_atoms: strict maximum number of atoms (atoms < max_atoms)
        return_idx: whether to return index or a boolean mask
        n_jobs: number of parallel job to run. Sequential by default
        progress: whether to show progress bar
        scheduler: joblib scheduler to use

    Returns:
        filtered_mask: boolean array (or index array) where true means the molecule is ok.
    """

    def reject_fn(mol):
        if mol is None:
            return True
        num_atoms = mol.GetNumAtoms()
        return (min_atoms is not None and num_atoms <= min_atoms) or (
            max_atoms is not None and num_atoms >= max_atoms
        )

    return _generic_filter(
        mols,
        reject_fn,
        return_idx=return_idx,
        n_jobs=n_jobs,
        progress=progress,
        scheduler=scheduler,
    )


def halogenicity_filter(
    mols: Iterable[Union[str, dm.Mol]],
    thresh_F: int = 6,
    thresh_Br: int = 3,
    thresh_Cl: int = 3,
    return_idx: bool = False,
    n_jobs: Optional[int] = None,
    progress: bool = False,
    scheduler: Optional[str] = None,
):
    """Find molecule that do not exceed halogen threshold. These thresholds are:

    Args:
        mols: list of input molecules
        thresh_F: maximum number of fluorine
        thresh_Br: maximum number of bromine
        thresh_Cl: maximum number of chlorine
        return_idx: whether to return index or a boolean mask
        n_jobs: number of parallel job to run. Sequential by default
        progress: whether to show progress bar
        scheduler: joblib scheduler to use

    Returns:
        filtered_mask: boolean array (or index array) where true means the molecule is ok.

    """
    fluorine_smarts = dm.from_smarts("[F]")
    bromine_smarts = dm.from_smarts("[Br]")
    chlorine_smarts = dm.from_smarts("[Cl]")

    def reject_fn(mol):
        if mol is None:
            return True
        fluorine_saturation = (
            len(mol.GetSubstructMatches(fluorine_smarts, uniquify=True)) > thresh_F
        )
        bromine_saturation = (
            len(mol.GetSubstructMatches(bromine_smarts, uniquify=True)) > thresh_Br
        )
        chlorine_saturation = (
            len(mol.GetSubstructMatches(chlorine_smarts, uniquify=True)) > thresh_Cl
        )
        return fluorine_saturation or bromine_saturation or chlorine_saturation

    return _generic_filter(
        mols,
        reject_fn,
        return_idx=return_idx,
        n_jobs=n_jobs,
        progress=progress,
        scheduler=scheduler,
    )


def symmetry_filter(
    mols: Iterable[Union[str, dm.Mol]],
    symmetry_threshold: float = 0.8,
    return_idx: bool = False,
    n_jobs: Optional[int] = None,
    progress: bool = False,
    scheduler: Optional[str] = None,
):
    """Find molecules that are not symmetrical, given a symmetry threshold

    Args:
        mols: list of input molecules
        symmetry_threshold: threshold to consider a molecule highly symmetrical
        return_idx: whether to return index or a boolean mask
        n_jobs: number of parallel job to run. Sequential by default
        progress: whether to show progress bar
        scheduler: joblib scheduler to use

    Returns:
        filtered_mask: boolean array (or index array) where true means the molecule is ok.
    """

    def reject_fn(mol):
        return score_symmetry(mol) > symmetry_threshold

    return _generic_filter(
        mols,
        reject_fn,
        return_idx=return_idx,
        n_jobs=n_jobs,
        progress=progress,
        scheduler=scheduler,
    )
