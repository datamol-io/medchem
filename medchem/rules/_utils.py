import datamol as dm
from rdkit import Chem


def _in_range(x, min_val: float = -float("inf"), max_val: float = float("inf")):
    """Check if a value is in a range
    Args:
        x: value to check
        min_val: minimum value
        max_val: maximum value
    """
    return min_val <= x <= max_val


def _compute_ring_system(mol: dm.Mol, include_spiro: bool = True):
    """
    Compute the list of ring system in a molecule. This is based on RDKit's cookbook:
    https://www.rdkit.org/docs/Cookbook.html#rings-aromaticity-and-kekulization

    # EN: move to datamol

    Args:
        mol: input molecule
        include_spiro: whether to include spiro rings. Defaults to False.

    Returns:
        ring_system: list of ring system
    """
    ri = mol.GetRingInfo()
    systems = []
    for ring in ri.AtomRings():
        ringAts = set(ring)
        nSystems = []
        for system in systems:
            nInCommon = len(ringAts.intersection(system))
            if nInCommon and (include_spiro or nInCommon > 1):
                ringAts = ringAts.union(system)
            else:
                nSystems.append(system)
        nSystems.append(ringAts)
        systems = nSystems
    return systems


def _compute_charge(mol: dm.Mol):
    """
    Compute the charge of a molecule.

    Args:
        mol: input molecule

    Returns:
        charge: charge of the molecule
    """
    return Chem.rdmolops.GetFormalCharge(mol)


def _compute_refractivity(mol: dm.Mol):
    """
    Compute the refractivity of a molecule.

    Args:
        mol: input molecule

    Returns:
        mr: molecular refractivity of the molecule
    """
    return Chem.Crippen.MolMR(mol)


def _compute_rigid_bonds(mol: dm.Mol):
    """
    Compute the number of rigid bonds in a molecule.

    Args:
        mol: input molecule

    Returns:
        n_rigid_bonds: number of rigid bonds in the molecule
    """
    # rigid bonds are bonds that are not (single and not in rings)
    # I don't think this is the same as the number of non rotatable bonds ?
    # EN: move this to datamol ?
    non_rigid_bonds_count = dm.from_smarts("*-&!@*")
    n_rigid_bonds = mol.GetNumBonds() - len(
        mol.GetSubstructMatches(non_rigid_bonds_count)
    )
    return n_rigid_bonds


def _compute_n_stereo_center(mol: dm.Mol):
    """
    Compute the number of stereocenters in a molecule.

    Args:
        mol: input molecule

    Returns:
        n_stero_center: number of stereocenters in the molecule
    """
    n_stereo_center = 0
    try:
        Chem.FindPotentialStereo(mol, cleanIt=False)  # type: ignore
        n_stereo_center = Chem.rdMolDescriptors.CalcNumAtomStereoCenters(mol)
    except:
        pass
    return n_stereo_center


def _compute_n_charged_atoms(mol: dm.Mol):
    """
    Compute the number of charged atoms in a molecule.

    Args:
        mol: input molecule

    Returns:
        n_charged_atoms: number of charged atoms in the molecule
    """
    return sum([at.GetFormalCharge() != 0 for at in mol.GetAtoms()])
