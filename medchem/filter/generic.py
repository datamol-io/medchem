from typing import Iterable
import datamol as dm
from rdkit.Chem import MolFromSmarts


def macrocycle_filter(
    mols: Iterable, max_cycle_size: int = 10, return_idx: bool = False
):
    """Find molecules that do not infringe the maximum cycle size

    Args:
        mols: list of input molecules
        max_cycle_size: maximum macrocycle size
        return_idx: whether to return index or a boolean mask

    Returns:
        filtered_mask: boolean array (or index array) where true means the molecule is ok.

    """
    filtered_mask = []
    filtered_idx = []
    for i, mol in enumerate(mols):
        mol = dm.to_mol(mol)
        rinfo = mol.GetRingInfo().AtomRings()
        can_add = True
        for r in rinfo:
            if len(r) >= max_cycle_size:
                can_add = False
                break
        if can_add:
            filtered_idx.append(i)
        filtered_mask.append(can_add)
    if return_idx:
        return filtered_idx
    return filtered_mask


def atom_list_filter(mols: Iterable, atom_list: Iterable, return_idx: bool = False):
    """Find molecule without any atom from a set of unwanted atom symbols

    Args:
        mols: list of input molecules
        atom_list: list of undesirable atom symbol
        return_idx: whether to return index or a boolean mask

    Returns:
        filtered_mask: boolean array (or index array) where true means the molecule is ok.
    """
    filtered_mask = []
    filtered_idx = []
    for i, mol in enumerate(mols):
        mol = dm.to_mol(mol)
        add_mol = True
        for atom in mol.GetAtoms():
            if atom.GetSymbol() in atom_list:
                add_mol = False
                break
        filtered_mask.append(add_mol)
        if add_mol:
            filtered_idx.append(i)
    if return_idx:
        return filtered_idx
    return filtered_mask


def ring_infraction_filter(mols: Iterable, return_idx: bool = False):
    """
    Find molecules that have a ring infraction filter.
    Returning True means the molecule is fine

    Args:
        mols: list of input molecules
        return_idx: whether to return index or a boolean mask

    Returns:
        filtered_mask: boolean array (or index array) where true means the molecule is ok.
    """
    filtered_mask = []
    filtered_idx = []
    for i, molecule in enumerate(mols):
        molecule = dm.to_mol(molecule)
        ring_allene = molecule.HasSubstructMatch(MolFromSmarts("[R]=[R]=[R]"))
        double_bond_in_small_ring = molecule.HasSubstructMatch(
            MolFromSmarts("[r3,r4]=[r3,r4]")
        )
        if not any([ring_allene, double_bond_in_small_ring]):
            filtered_mask.append(True)
            filtered_idx.append(i)
        else:
            filtered_mask.append(False)
    if return_idx:
        return filtered_idx
    return filtered_mask


def num_atom_filter(
    mols: Iterable,
    min_atoms: int = None,
    max_atoms: int = None,
    return_idx: bool = False,
):
    """
    Find a molecule that match the atom number constraints
    Returning True means the molecule is fine

    Args:
        mols: list of input molecules
        min_atoms: minimum number of atoms
        max_atoms: maximum number of atoms
        return_idx: whether to return index or a boolean mask

    Returns:
        filtered_mask: boolean array (or index array) where true means the molecule is ok.
    """
    filtered_idx = []
    filtered_mask = []
    for i, mol in enumerate(mols):
        mol = dm.to_mol(mol)
        num_atoms = mol.GetNumAtoms()
        if (not min_atoms or min_atoms <= num_atoms) and (
            not max_atoms or num_atoms <= max_atoms
        ):
            filtered_mask.append(True)
            filtered_idx.append(i)
        else:
            filtered_mask.append(False)
    if return_idx:
        return filtered_idx
    return filtered_mask


def halogenicity_filter(
    mols: Iterable,
    thresh_F: int = 6,
    thresh_Br: int = 3,
    thresh_Cl: int = 3,
    return_idx: bool = False,
):
    """Find molecule that do not exceed halogen threshold. These thresholds are:

    Args:
        mols: list of input molecules
        thresh_F: maximum number of fluorine
        thresh_Br: maximum number of bromine
        thresh_Cl: maximum number of chlorine
        return_idx: whether to return index or a boolean mask

    Returns:
        filtered_mask: boolean array (or index array) where true means the molecule is ok.

    """
    filtered_mask = []
    filtered_idx = []
    for i, molecule in enumerate(mols):
        molecule = dm.to_mol(molecule)
        fluorine_saturation = (
            len(molecule.GetSubstructMatches(MolFromSmarts("[F]"))) > thresh_F
        )
        bromide_saturation = (
            len(molecule.GetSubstructMatches(MolFromSmarts("[Br]"))) > thresh_Br
        )
        chlorine_saturation = (
            len(molecule.GetSubstructMatches(MolFromSmarts("[Cl]"))) > thresh_Cl
        )
        if not any([chlorine_saturation, bromide_saturation, fluorine_saturation]):
            filtered_mask.append(True)
            filtered_idx.append(i)
        else:
            filtered_mask.append(False)
    if return_idx:
        return filtered_idx
    return filtered_mask
