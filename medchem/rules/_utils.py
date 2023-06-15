from typing import List
import itertools
import datamol as dm
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from loguru import logger
from medchem.utils.smarts import SMARTSUtils

_DESCRIPTOR_LIST = [
    "mw",
    "fsp3",
    "n_lipinski_hba",
    "n_lipinski_hbd",
    "n_rings",
    "n_hetero_atoms",
    "n_heavy_atoms",
    "n_rotatable_bonds",
    "n_radical_electrons",
    "n_NHOH",
    "n_NO",
    "tpsa",
    "qed",
    "clogp",
    "sas",
    "formal_charge",
    "refractivity",
    "n_aliphatic_carbocycles",
    "n_aliphatic_heterocycles",
    "n_aliphatic_rings",
    "n_aromatic_carbocycles",
    "n_aromatic_heterocycles",
    "n_aromatic_rings",
    "n_saturated_carbocycles",
    "n_saturated_heterocycles",
    "n_saturated_rings",
    "n_charged_atoms",
    "n_stereo_center",
    "n_aromatic_atoms",
    "n_rigid_bonds",
    "n_aromatic_atoms_proportion",
]


def _in_range(x, min_val: float = -float("inf"), max_val: float = float("inf")):
    """Check if a value is in a range
    Args:
        x: value to check
        min_val: minimum value
        max_val: maximum value
    """
    return min_val <= x <= max_val


def n_heavy_metals(mol: dm.Mol, allowed_metals: List[str] = ["Li", "Be", "K", "Na", "Ca", "Mg"]):
    """Count the number of heavy metals in a molecule

    Metal are defined using the M notation in marvinjs. It's quicker to exclude atoms than to list all metals

    Args:
        mol: input molecule
        allowed_metals: list of metals not counted as heavy metals. Default is ["Li", "Be", "K", "Na", "Ca", "Mg"]
    """
    non_metals = [
        2,
        5,
        6,
        7,
        8,
        9,
        10,
        14,
        15,
        16,
        17,
        18,
        33,
        34,
        35,
        36,
        52,
        53,
        54,
        85,
        86,
    ]

    heavy_metals = [
        x.GetAtomicNum()
        for x in mol.GetAtoms()
        if x.GetAtomicNum() not in non_metals and x.GetSymbol() not in allowed_metals
    ]
    return len(heavy_metals)


def has_spider_chains(mol: dm.Mol, min_flagel: int = 2, min_flagel_len: int = 4):
    """Check whether a molecule has multiple flagel like
    Args:
        mol: input molecule
        min_flagel: minimum number of flagels (>=)
        min_flagel_len: minimum length of a flagel (>=)
    """
    # we first need to remove all rings systems (so scaffold from the molecules)
    # then we can count the both condition for spider/flagel like look
    side_chains = [mol]
    flagel_query = SMARTSUtils.aliphatic_chain(min_size=min_flagel_len)
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    try:
        scaffold = dm.sanitize_mol(scaffold)
        side_chains = Chem.ReplaceCore(mol, scaffold, labelByIndex=False)
        if side_chains is not None:
            side_chains = list(Chem.GetMolFrags(side_chains, asMols=True))
            side_chains = [dm.to_smiles(x) for x in side_chains]
            side_chains = [SMARTSUtils.standardize_attachment(x, "[1*]") for x in side_chains]
            side_chains = [dm.to_mol(x) for x in side_chains]
            flagel_query = "[1*]-,=" + flagel_query
        else:
            side_chains = [mol]
    except Exception as e:
        logger.error(e)

    # extract side chains from the scaffold
    flagel_query = dm.from_smarts(flagel_query)
    matches = [x.HasSubstructMatch(flagel_query) for x in side_chains]
    return sum(matches) >= min_flagel


def n_fused_aromatic_rings(mol: dm.Mol, require_all_aromatic: bool = True, pairwise: bool = False):
    """Count the number of fused aromatic rings in a molecule

    !!! warning
        There is no such thing as a spiroaoaromatic ring in this implementation

    Args:
        mol: input molecule
        require_all_aromatic: whether to require all simple rings in the fused system to be aromatic
        pairwise: whether to compute the number of fused aromatic rings pairwise.
            meaning phenanthrene and anthracene would count for 2 fused aromatic rings each
    """
    # EN: might make sense to move this to datamol
    # This code can be spedt up by sacrificing readability, will revisit eventually

    ring_systems = mol.GetRingInfo()
    # we use bond since we are focusing on fused rings
    simple_rings = list(ring_systems.BondRings())
    rings = [set(r) for r in simple_rings]
    ring_map = [set([x]) for x in range(len(rings))]
    go_next = True
    while go_next:
        go_next = False
        for i, j in itertools.combinations(range(len(rings)), 2):
            if rings[i] & rings[j]:
                new_map = set().union(ring_map[i], ring_map[j])
                q = rings[i] | rings[j]
                min_ind, max_ind = min(i, j), max(i, j)
                del rings[max_ind], rings[min_ind]
                del ring_map[max_ind], ring_map[min_ind]
                rings.append(q)
                ring_map.append(new_map)
                go_next = True
                break
    # simple_rings: is the list of simple rings from ring info
    # rings: the list of rings after mergin fused rings
    # ring_map: the mapping between fused rings and the basic rings their contains
    if pairwise:
        # we need to count all pair in any fused rings with more than 2 fused rings
        fused_rings = []
        for _, ring_ids in enumerate(ring_map):
            for ring_1, ring_2 in itertools.combinations(ring_ids, 2):
                if set(simple_rings[ring_1]) & set(simple_rings[ring_2]):
                    fused_ring = set.union(set(simple_rings[ring_1]), (simple_rings[ring_2]))
                    fused_rings.append(fused_ring)
    else:
        fused_rings = [r for i, r in enumerate(rings) if len(ring_map[i]) >= 2]
    n_aromatic_fused_rings = 0
    for i, fused_ring_bonds in enumerate(fused_rings):
        aromatic_system = [mol.GetBondWithIdx(bond_id).GetIsAromatic() for bond_id in fused_ring_bonds]
        if require_all_aromatic:
            n_aromatic_fused_rings += all(aromatic_system)
        else:
            n_aromatic_fused_rings += any(aromatic_system)

    return n_aromatic_fused_rings


def fraction_atom_in_scaff(
    mol: dm.Mol,
):
    """Compute the fraction of atoms that belong to any ring system of the molecule
    as defined by the scaffold

    Args:
        mol: input molecule
    """
    n_heavy_atoms = mol.GetNumHeavyAtoms()
    if n_heavy_atoms < 1:
        return 0
    n_heavy_scaffold_atoms = 0
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    try:
        scaffold = dm.sanitize_mol(scaffold)
        n_heavy_scaffold_atoms = scaffold.GetNumHeavyAtoms()
    except:
        pass
    return n_heavy_scaffold_atoms / n_heavy_atoms


def list_descriptors():
    """List all descriptors available for computation"""
    # EN: this is until datamol release
    return _DESCRIPTOR_LIST
