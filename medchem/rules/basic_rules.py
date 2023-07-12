from typing import Union
from typing import Optional
from typing import Any

import datamol as dm

from ._utils import in_range
from ._utils import n_fused_aromatic_rings
from ._utils import n_heavy_metals
from ._utils import has_spider_chains
from ._utils import fraction_atom_in_scaff


def rule_of_five(
    mol: Union[dm.Mol, str],
    mw: Optional[float] = None,
    clogp: Optional[float] = None,
    n_lipinski_hbd: Optional[float] = None,
    n_lipinski_hba: Optional[float] = None,
    **kwargs: Any,
) -> bool:
    """Compute the Lipinski's rule-of-5 for a molecule. Also known as Pfizer's rule of five or RO5.
    This rule is a rule of thumb to evaluate the druglikeness of a chemical compounds.

    It computes: `MW <= 500 & logP <= 5  & HBD <= 5 & HBA <= 10`

    Args:
        mol: input molecule
        mw: precomputed molecular weight.
        clogp: precomputed cLogP.
        n_lipinski_hbd: precomputed number of HBD.
        n_lipinski_hba: precomputed number of HBA.
        **kwargs: Allow extra arguments for descriptors pre-computation.

    Returns:
        ro5: True if molecule is compliant, False otherwise
    """
    if isinstance(mol, str):
        mol = dm.to_mol(mol)

    if mol is None:
        raise ValueError("Molecule is None")

    mw = mw if mw is not None else dm.descriptors.mw(mol)
    clogp = clogp if clogp is not None else dm.descriptors.clogp(mol)
    n_lipinski_hbd = n_lipinski_hbd if n_lipinski_hbd is not None else dm.descriptors.n_lipinski_hbd(mol)
    n_lipinski_hba = n_lipinski_hba if n_lipinski_hba is not None else dm.descriptors.n_lipinski_hba(mol)
    return mw <= 500 and clogp <= 5 and n_lipinski_hbd <= 5 and n_lipinski_hba <= 10


def rule_of_five_beyond(
    mol: Union[dm.Mol, str],
    mw: Optional[float] = None,
    clogp: Optional[float] = None,
    n_hbd: Optional[float] = None,
    n_hba: Optional[float] = None,
    tpsa: Optional[float] = None,
    n_rotatable_bonds: Optional[int] = None,
    **kwargs: Any,
) -> bool:
    """Compute the Beyond rule-of-5 rule for a molecule. This rule illustrates the potential of compounds far beyond rule-of-5 space to
    modulate novel and difficult target classes that have large, flat, and groove-shaped binding sites and has been described in:

    [Doak, Bradley C., et al. (2015) How Beyond Rule of 5 Drugs and Clinical Candidates Bind to Their Targets.](https://pubmed.ncbi.nlm.nih.gov/26457449/)

    It computes: `MW <= 1000 & logP in [-2, 10] & HBD <= 6 & HBA <= 15 & TPSA <=250 & ROTBONDS <= 20`

    !!! note
        This is a very permissive rule and is likely to not be a good predictor for druglikeness as known for small molecules.

    Args:
        mol: input molecule
        mw: precomputed molecular weight.
        clogp: precomputed cLogP.
        n_hbd: precomputed number of HBD.
        n_hba: precomputed number of HBA.
        tpsa: precomputed TPSA.
        n_rotatable_bonds: precomputed number of rotatable bonds.
        **kwargs: Allow extra arguments for descriptors pre-computation.

    Returns:
        ro5: True if molecule is compliant, False otherwise
    """
    if isinstance(mol, str):
        mol = dm.to_mol(mol)

    if mol is None:
        raise ValueError("Molecule is None")

    mw = mw if mw is not None else dm.descriptors.mw(mol)
    clogp = clogp if clogp is not None else dm.descriptors.clogp(mol)
    n_hbd = n_hbd if n_hbd is not None else dm.descriptors.n_hbd(mol)
    n_hba = n_hba if n_hba is not None else dm.descriptors.n_hba(mol)
    tpsa = tpsa if tpsa is not None else dm.descriptors.tpsa(mol)
    n_rotatable_bonds = (
        n_rotatable_bonds if n_rotatable_bonds is not None else dm.descriptors.n_rotatable_bonds(mol)
    )
    return (
        mw <= 1000
        and in_range(clogp, -2, 10)
        and n_hbd <= 6
        and n_hba <= 15
        and tpsa <= 250
        and n_rotatable_bonds <= 20
    )


def rule_of_zinc(
    mol: Union[dm.Mol, str],
    mw: Optional[float] = None,
    clogp: Optional[float] = None,
    n_hba: Optional[float] = None,
    n_hbd: Optional[float] = None,
    tpsa: Optional[float] = None,
    n_rotatable_bonds: Optional[int] = None,
    n_rings: Optional[int] = None,
    charge: Optional[float] = None,
    **kwargs: Any,
) -> bool:
    """Compute the Zinc rule for a molecule. This rule is a rule of thumb to evaluate the druglikeness of a chemical compounds, based on:

    [Irwin & Schoichet (2005) ZINC - A Free Database of Commercially Available Compounds for Virtual Screening](https://pubmed.ncbi.nlm.nih.gov/15667143/)

     Also see: https://fafdrugs4.rpbs.univ-paris-diderot.fr/filters.html

     It computes: `MW in [60, 600] & logP < in [-4, 6] & HBD <= 6 & HBA <= 11 & TPSA <=150 & ROTBONDS <= 12 & RIGBONDS <= 50 & N_RINGS <= 7 & MAX_SIZE_RING <= 12 & N_CARBONS >=3 & HC_RATIO <= 2.0 & CHARGE in [-4, 4]`

     Args:
         mol: input molecule
         mw: precomputed molecular weight.
         clogp: precomputed cLogP.
         n_hba: precomputed number of HBA.
         n_hbd: precomputed number of HBD.
         tpsa: precomputed TPSA.
         n_rotatable_bonds: precomputed number of rotatable bonds.
         n_rings: precomputed number of rings in the molecules.
         charge: precomputed charge.
         **kwargs: Allow extra arguments for descriptors pre-computation.

    """
    if isinstance(mol, str):
        mol = dm.to_mol(mol)

    if mol is None:
        raise ValueError("Molecule is None")

    mw = mw if mw is not None else dm.descriptors.mw(mol)
    clogp = clogp if clogp is not None else dm.descriptors.clogp(mol)
    n_hba = n_hba if n_hba is not None else dm.descriptors.n_hba(mol)
    n_hbd = n_hbd if n_hbd is not None else dm.descriptors.n_hbd(mol)
    tpsa = tpsa if tpsa is not None else dm.descriptors.tpsa(mol)
    n_rotatable_bonds = (
        n_rotatable_bonds if n_rotatable_bonds is not None else dm.descriptors.n_rotatable_bonds(mol)
    )
    n_rings = n_rings if n_rings is not None else dm.descriptors.n_rings(mol)
    n_rigid_bonds = dm.descriptors.n_rigid_bonds(mol)
    ring_system = dm.compute_ring_system(mol, include_spiro=False)
    max_size_ring = 0 if len(ring_system) == 0 else max([len(x) for x in ring_system])
    n_carbons = len([at for at in mol.GetAtoms() if at.GetSymbol() == "C"])
    if n_carbons == 0:
        het_carb_ratio = float("inf")
    else:
        het_carb_ratio = dm.descriptors.n_hetero_atoms(mol) / n_carbons
    charge = charge if charge is not None else dm.descriptors.formal_charge(mol)
    return (
        in_range(mw, 60, 600)
        and in_range(clogp, -4, 6)
        and n_hbd <= 6
        and n_hba <= 11
        and tpsa <= 150
        and n_rotatable_bonds <= 12
        and n_rigid_bonds <= 50
        and n_rings <= 7
        and max_size_ring <= 12
        and n_carbons >= 3
        and in_range(het_carb_ratio, 0, 2.0)
        and in_range(charge, -4, 4)
    )


def rule_of_leadlike_soft(
    mol: Union[dm.Mol, str],
    mw: Optional[float] = None,
    clogp: Optional[float] = None,
    n_hba: Optional[float] = None,
    n_hbd: Optional[float] = None,
    tpsa: Optional[float] = None,
    n_rotatable_bonds: Optional[int] = None,
    n_rings: Optional[int] = None,
    n_hetero_atoms: Optional[int] = None,
    charge: Optional[float] = None,
    **kwargs: Any,
) -> bool:
    """
    Compute the Lead-Like Soft rule available in [FAF-Drugs4](https://fafdrugs4.rpbs.univ-paris-diderot.fr/filters.html).

    It computes:
    ```
    MW in [150, 400] & logP < in [-3, 4] & HBD <= 4 & HBA <= 7 & TPSA <=160 & ROTBONDS <= 9 &
    RIGBONDS <= 30 & N_RINGS <= 4 & MAX_SIZE_RING <= 18 & N_CARBONS in [3, 35] &  N_HETEROATOMS in [1, 15] &
    HC_RATIO in [0.1, 1.1] & CHARGE in [-4, 4] & N_ATOM_CHARGE <= 4 & N_STEREO_CENTER <= 2
    ```

    Args:
        mol: input molecule
        mw: precomputed molecular weight.
        clogp: precomputed cLogP.
        n_hba: precomputed number of HBA.
        n_hbd: precomputed number of HBD.
        tpsa: precomputed TPSA.
        n_rotatable_bonds: precomputed number of rotatable bonds.
        n_rings: precomputed number of rings in the molecules.
        n_hetero_atoms: precomputed number of heteroatoms.
        charge: precomputed charge.
        **kwargs: Allow extra arguments for descriptors pre-computation.

    """
    if isinstance(mol, str):
        mol = dm.to_mol(mol)

    if mol is None:
        raise ValueError("Molecule is None")

    mw = mw if mw is not None else dm.descriptors.mw(mol)
    clogp = clogp if clogp is not None else dm.descriptors.clogp(mol)
    n_hba = n_hba if n_hba is not None else dm.descriptors.n_hba(mol)
    n_hbd = n_hbd if n_hbd is not None else dm.descriptors.n_hbd(mol)
    tpsa = tpsa if tpsa is not None else dm.descriptors.tpsa(mol)
    n_rotatable_bonds = (
        n_rotatable_bonds if n_rotatable_bonds is not None else dm.descriptors.n_rotatable_bonds(mol)
    )
    n_rings = n_rings if n_rings is not None else dm.descriptors.n_rings(mol)
    n_hetero_atoms = n_hetero_atoms if n_hetero_atoms is not None else dm.descriptors.n_hetero_atoms(mol)
    charge = charge if charge is not None else dm.descriptors.formal_charge(mol)
    num_charged_atom = dm.descriptors.n_charged_atoms(mol)
    n_stereo_center = dm.descriptors.n_stereo_centers(mol)
    n_rigid_bonds = dm.descriptors.n_rigid_bonds(mol)
    ring_system = dm.compute_ring_system(mol, include_spiro=False)
    max_size_ring = 0 if len(ring_system) == 0 else max([len(x) for x in ring_system])
    n_carbons = len([at for at in mol.GetAtoms() if at.GetSymbol() == "C"])
    if n_carbons == 0:
        het_carb_ratio = float("inf")
    else:
        het_carb_ratio = n_hetero_atoms / n_carbons
    return (
        in_range(mw, 150, 400)
        and in_range(clogp, -3, 4)
        and n_hbd <= 4
        and n_hba <= 7
        and tpsa <= 160
        and n_rotatable_bonds <= 9
        and n_rigid_bonds <= 30
        and n_rings <= 4
        and max_size_ring <= 18
        and in_range(n_carbons, 3, 35)
        and in_range(n_hetero_atoms, 1, 15)
        and in_range(het_carb_ratio, 0.1, 1.1)
        and in_range(charge, -4, 4)
        and num_charged_atom <= 4
        and n_stereo_center <= 2
    )


def rule_of_druglike_soft(
    mol: Union[dm.Mol, str],
    mw: Optional[float] = None,
    clogp: Optional[float] = None,
    n_hba: Optional[float] = None,
    n_hbd: Optional[float] = None,
    tpsa: Optional[float] = None,
    n_rotatable_bonds: Optional[int] = None,
    n_rings: Optional[int] = None,
    n_hetero_atoms: Optional[int] = None,
    charge: Optional[float] = None,
    **kwargs: Any,
) -> bool:
    """
    Compute the DrugLike Soft rule available in [FAF-Drugs4](https://fafdrugs4.rpbs.univ-paris-diderot.fr/filters.html).

    It computes:
    ```
    MW in [100, 600] & logP < in [-3, 6] & HBD <= 7 & HBA <= 12 & TPSA <=180 & ROTBONDS <= 11 &
    RIGBONDS <= 30 & N_RINGS <= 6 & MAX_SIZE_RING <= 18 & N_CARBONS in [3, 35] &  N_HETEROATOMS in [1, 15] &
    HC_RATIO in [0.1, 1.1] & CHARGE in [-4, 4] & N_ATOM_CHARGE <= 4
    ```

    Args:
        mol: input molecule
        mw: precomputed molecular weight.
        clogp: precomputed cLogP.
        n_hba: precomputed number of HBA.
        n_hbd: precomputed number of HBD.
        tpsa: precomputed TPSA.
        n_rotatable_bonds: precomputed number of rotatable bonds.
        n_rings: precomputed number of rings in the molecules.
        n_hetero_atoms: precomputed number of heteroatoms.
        charge: precomputed charge.
        **kwargs: Allow extra arguments for descriptors pre-computation.
    """
    if isinstance(mol, str):
        mol = dm.to_mol(mol)

    if mol is None:
        raise ValueError("Molecule is None")

    mw = mw if mw is not None else dm.descriptors.mw(mol)
    clogp = clogp if clogp is not None else dm.descriptors.clogp(mol)
    n_hba = n_hba if n_hba is not None else dm.descriptors.n_hba(mol)
    n_hbd = n_hbd if n_hbd is not None else dm.descriptors.n_hbd(mol)
    tpsa = tpsa if tpsa is not None else dm.descriptors.tpsa(mol)
    n_rotatable_bonds = (
        n_rotatable_bonds if n_rotatable_bonds is not None else dm.descriptors.n_rotatable_bonds(mol)
    )
    n_rings = n_rings if n_rings is not None else dm.descriptors.n_rings(mol)
    n_hetero_atoms = n_hetero_atoms if n_hetero_atoms is not None else dm.descriptors.n_hetero_atoms(mol)
    charge = charge if charge is not None else dm.descriptors.formal_charge(mol)
    num_charged_atom = dm.descriptors.n_charged_atoms(mol)
    n_rigid_bonds = dm.descriptors.n_rigid_bonds(mol)
    ring_system = dm.compute_ring_system(mol, include_spiro=False)
    max_size_ring = 0 if len(ring_system) == 0 else max([len(x) for x in ring_system])
    n_carbons = len([at for at in mol.GetAtoms() if at.GetSymbol() == "C"])
    n_hydrogens = sum([at.GetTotalNumHs() for at in mol.GetAtoms()])

    if n_carbons == 0:
        het_carb_ratio = float("inf")
    else:
        het_carb_ratio = n_hydrogens / n_carbons

    return (
        in_range(mw, 100, 600)
        and in_range(clogp, -3, 6)
        and n_hbd <= 7
        and n_hba <= 12
        and tpsa <= 180
        and n_rotatable_bonds <= 11
        and n_rigid_bonds <= 30
        and n_rings <= 6
        and max_size_ring <= 18
        and in_range(n_carbons, 3, 35)
        and in_range(n_hetero_atoms, 1, 15)
        and in_range(het_carb_ratio, 0.1, 1.1)
        and in_range(charge, -4, 4)
        and num_charged_atom <= 4
    )


def rule_of_four(
    mol: Union[dm.Mol, str],
    mw: Optional[float] = None,
    clogp: Optional[float] = None,
    n_hba: Optional[float] = None,
    n_rings: Optional[int] = None,
    **kwargs: Any,
) -> bool:
    """Compute the rule-of-4 for a molecule. The rule-of-4 define a rule of thumb for PPI inhibitors,
    which are typically larger and more lipophilic than inhibitors of standard binding sites. It has been published in:
    [Morelli X, Bourgeas R, Roche P. (2011) Chemical and structural lessons from recent successes in protein–protein interaction inhibition](https://pubmed.ncbi.nlm.nih.gov/21684802/)

    Also see: [Shin et al. (2020) Current Challenges and Opportunities in Designing Protein–Protein Interaction Targeted Drugs](https://pubmed.ncbi.nlm.nih.gov/33209039/)

    It computes: `MW >= 400 & logP >= 4  & RINGS >=4 & HBA >= 4`


    !!! warning
        Do not use this for small molecules that are not PPI inhibitors !

    Args:
        mol: input molecule
        mw: precomputed molecular weight.
        clogp: precomputed cLogP.
        n_hba: precomputed number of HBA.
        n_rings: precomputed number of rings in the molecules.
        **kwargs: Allow extra arguments for descriptors pre-computation.

    Returns:
        ro4: True if molecule is compliant, False otherwise
    """
    if isinstance(mol, str):
        mol = dm.to_mol(mol)

    if mol is None:
        raise ValueError("Molecule is None")

    mw = mw if mw is not None else dm.descriptors.mw(mol)
    clogp = clogp if clogp is not None else dm.descriptors.clogp(mol)
    n_hba = n_hba if n_hba is not None else dm.descriptors.n_hba(mol)
    n_rings = n_rings if n_rings is not None else dm.descriptors.n_rings(mol)

    return mw >= 400 and clogp >= 4 and n_rings >= 4 and n_hba >= 4


def rule_of_three(
    mol: Union[dm.Mol, str],
    mw: Optional[float] = None,
    clogp: Optional[float] = None,
    n_hba: Optional[float] = None,
    n_hbd: Optional[float] = None,
    n_rotatable_bonds: Optional[int] = None,
    **kwargs: Any,
) -> bool:
    """Compute the rule-of-3. The rule-of-three is a rule of thumb for molecular fragments (and not small molecules) published in:
    [Congreve M, Carr R, Murray C, Jhoti H. (2003) `A "rule of three" for fragment-based lead discovery?`](https://pubmed.ncbi.nlm.nih.gov/14554012/)

    It computes: `MW <= 300 & logP <= 3 & HBA <= 3 & HBD <= 3 & ROTBONDS <= 3`

    !!! note
        TPSA is not used in this version of the rule of three. Other version uses `TPSA <= 60 AND logP in [-3, 3]` in addition

    Args:
        mol: input molecule
        mw: precomputed molecular weight.
        clogp: precomputed cLogP.
        n_hba: precomputed number of HBA.
        n_hbd: precomputed number of HBD.
        n_rotatable_bonds: precomputed number of rotatable bonds in the molecule.
        **kwargs: Allow extra arguments for descriptors pre-computation.

    Returns:
        ro3: True if molecule is compliant, False otherwise
    """

    if isinstance(mol, str):
        mol = dm.to_mol(mol)

    if mol is None:
        raise ValueError("Molecule is None")

    mw = mw if mw is not None else dm.descriptors.mw(mol)
    clogp = clogp if clogp is not None else dm.descriptors.clogp(mol)
    n_hba = n_hba if n_hba is not None else dm.descriptors.n_hba(mol)
    n_hbd = n_hbd if n_hbd is not None else dm.descriptors.n_hbd(mol)
    n_rotatable_bonds = (
        n_rotatable_bonds if n_rotatable_bonds is not None else dm.descriptors.n_rotatable_bonds(mol)
    )

    return mw <= 300 and clogp <= 3 and n_hbd <= 3 and n_hba <= 3 and n_rotatable_bonds <= 3


def rule_of_three_extended(
    mol: Union[dm.Mol, str],
    mw: Optional[float] = None,
    clogp: Optional[float] = None,
    n_hba: Optional[float] = None,
    n_hbd: Optional[float] = None,
    tpsa: Optional[float] = None,
    n_rotatable_bonds: Optional[int] = None,
    **kwargs: Any,
) -> bool:
    """Compute the extended rule-of-3 which is an extension of the rule of three that
    includes and `TPSA` and relaxes `HBA` constraints.

    It computes: `MW <= 300 & logP in [-3, 3]  & HBA <= 6 & HBD <= 3 & ROTBONDS <= 3 & TPSA <= 60`

    Args:
        mol: input molecule
        mw: precomputed molecular weight.
        clogp: precomputed cLogP.
        n_hba: precomputed number of HBA.
        n_hbd: precomputed number of HBD.
        tpsa: precomputed TPSA.
        n_rotatable_bonds: precomputed number of rotatable bonds in the molecule.
        **kwargs: Allow extra arguments for descriptors pre-computation.

    Returns:
        ro3: True if molecule is compliant, False otherwise
    """

    if isinstance(mol, str):
        mol = dm.to_mol(mol)

    if mol is None:
        raise ValueError("Molecule is None")

    mw = mw if mw is not None else dm.descriptors.mw(mol)
    clogp = clogp if clogp is not None else dm.descriptors.clogp(mol)
    n_hba = n_hba if n_hba is not None else dm.descriptors.n_hba(mol)
    n_hbd = n_hbd if n_hbd is not None else dm.descriptors.n_hbd(mol)
    tpsa = tpsa if tpsa is not None else dm.descriptors.tpsa(mol)
    n_rotatable_bonds = (
        n_rotatable_bonds if n_rotatable_bonds is not None else dm.descriptors.n_rotatable_bonds(mol)
    )

    return (
        mw <= 300
        and in_range(clogp, -3, 3)
        and n_hba <= 6
        and n_hbd <= 3
        and n_rotatable_bonds <= 3
        and tpsa <= 60
    )


def rule_of_two(
    mol: Union[dm.Mol, str],
    mw: Optional[float] = None,
    clogp: Optional[float] = None,
    n_hba: Optional[float] = None,
    n_hbd: Optional[float] = None,
    **kwargs: Any,
) -> bool:
    """
    Computes rules-of-2 for reagent (building block design). It aims for prioritization of reagents that typically
    do not add more than 200 Da in MW or 2 units of clogP. The rule of two has been described in:

    [Goldberg et al. (2015) Designing novel building blocks is an overlooked strategy to improve compound quality](http://csmres.co.uk/cs.public.upd/article-downloads/Designing-novel-building-blocks.pdf)

    !!! note
        Their analysis showed that molecular weight (MW) and clogP were important factors in the frequency of use of reagents.
        Other parameters, such as TPSA, HBA, HBD and ROTBONDS count, were less important.

    It computes `MW <= 200 & logP <= 2 & HBA <= 4 & HBD <= 2`

    Args:
        mol: input molecule
        mw: precomputed molecular weight.
        clogp: precomputed cLogP.
        n_hba: precomputed number of HBA.
        n_hbd: precomputed number of HBD.
        **kwargs: Allow extra arguments for descriptors pre-computation.

    Returns:
        ro2: True if molecule is compliant, False otherwise
    """
    if isinstance(mol, str):
        mol = dm.to_mol(mol)

    if mol is None:
        raise ValueError("Molecule is None")

    mw = mw if mw is not None else dm.descriptors.mw(mol)
    clogp = clogp if clogp is not None else dm.descriptors.clogp(mol)
    n_hba = n_hba if n_hba is not None else dm.descriptors.n_hba(mol)
    n_hbd = n_hbd if n_hbd is not None else dm.descriptors.n_hbd(mol)

    return mw <= 200 and clogp <= 2 and n_hba <= 4 and n_hbd <= 2


def rule_of_ghose(
    mol: Union[dm.Mol, str],
    mw: Optional[float] = None,
    clogp: Optional[float] = None,
    mr: Optional[float] = None,
    **kwargs: Any,
) -> bool:
    """
    Compute the Ghose filter. The Ghose filter is a drug-like filter described in:
    [Ghose, AK.; Viswanadhan, VN.; Wendoloski JJ. (1999) A knowledge-based approach in designing combinatorial or medicinal
    chemistry libraries for drug discovery.1. A qualitative and quantitative characterization of known drug databases.](https://pubs.acs.org/doi/10.1021/cc9800071)

    It computes: `MW in [160, 480] & logP in [-0.4, 5.6] & Natoms in [20, 70] & refractivity in [40, 130]`

    Args:
        mol: input molecule
        mw: precomputed molecular weight.
        clogp: precomputed cLogP.
        mr: precomputed molecule refractivity.
        **kwargs: Allow extra arguments for descriptors pre-computation.

    Returns:
        rog: True if molecule is compliant, False otherwise
    """

    if isinstance(mol, str):
        mol = dm.to_mol(mol)

    if mol is None:
        raise ValueError("Molecule is None")

    mw = mw if mw is not None else dm.descriptors.mw(mol)
    clogp = clogp if clogp is not None else dm.descriptors.clogp(mol)
    num_atoms = mol.GetNumAtoms()  # ghose seems to use total number of atoms not heavy
    mr = mr if mr is not None else dm.descriptors.refractivity(mol)
    return (
        in_range(mw, 160, 480)
        and in_range(clogp, -0.4, 5.6)
        and in_range(num_atoms, 20, 70)
        and in_range(mr, 40, 130)
    )


def rule_of_veber(
    mol: Union[dm.Mol, str],
    tpsa: Optional[float] = None,
    n_rotatable_bonds: Optional[int] = None,
    **kwargs: Any,
) -> bool:
    """
    Compute the Veber filter. The Veber filter is a druglike filter for orally active drugs described in:

    [Veber et. al. (2002) Molecular Properties That Influence the Oral Bioavailability of Drug Candidates](https://pubmed.ncbi.nlm.nih.gov/12036371/)

    It computes: `ROTBONDS <= 10 & TPSA < 140`

    Args:
        mol: input molecule
        tpsa: precomputed TPSA.
        n_rotatable_bonds: precomputed number of rotatable bonds.
        **kwargs: Allow extra arguments for descriptors pre-computation.

    Returns:
        rov: True if molecule is compliant, False otherwise
    """

    if isinstance(mol, str):
        mol = dm.to_mol(mol)

    if mol is None:
        raise ValueError("Molecule is None")

    n_rotatable_bonds = (
        n_rotatable_bonds if n_rotatable_bonds is not None else dm.descriptors.n_rotatable_bonds(mol)
    )
    tpsa = tpsa if tpsa is not None else dm.descriptors.tpsa(mol)

    return n_rotatable_bonds <= 10 and tpsa <= 140


def rule_of_reos(
    mol: Union[dm.Mol, str],
    mw: Optional[float] = None,
    clogp: Optional[float] = None,
    n_hba: Optional[float] = None,
    n_hbd: Optional[float] = None,
    charge: Optional[int] = None,
    n_rotatable_bonds: Optional[int] = None,
    n_heavy_atoms: Optional[int] = None,
    **kwargs: Any,
) -> bool:
    """
    Compute the REOS filter. The REOS filter is a filter designed to filter out unuseful compounds from HTS screening results.
    The filter is described in: [Waters & Namchuk (2003) Designing screens: how to make your hits a hit.](https://pubmed.ncbi.nlm.nih.gov/12669025/)

    It computes: `MW in [200, 500] & logP in [-5, 5] & HBA in [0, 10] & HBD in [0, 5] & charge in [-2, 2] & ROTBONDS in [0, 8] & NHeavyAtoms in [15, 50]`

    Args:
        mol: input molecule
        mw: precomputed molecular weight.
        clogp: precomputed cLogP.
        n_hba: precomputed number of HBA.
        n_hbd: precomputed number of HBD.
        charge: precomputed formal charge.
        n_rotatable_bonds: precomputed number of rotatable bonds in the molecule.
        n_heavy_atoms: precomputed number of heavy atoms in the molecule.
        **kwargs: Allow extra arguments for descriptors pre-computation.

    Returns:
        ror: True if molecule is compliant, False otherwise
    """
    if isinstance(mol, str):
        mol = dm.to_mol(mol)

    if mol is None:
        raise ValueError("Molecule is None")

    mw = mw if mw is not None else dm.descriptors.mw(mol)
    clogp = clogp if clogp is not None else dm.descriptors.clogp(mol)
    n_hba = n_hba if n_hba is not None else dm.descriptors.n_hba(mol)
    n_hbd = n_hbd if n_hbd is not None else dm.descriptors.n_hbd(mol)
    n_rotatable_bonds = (
        n_rotatable_bonds if n_rotatable_bonds is not None else dm.descriptors.n_rotatable_bonds(mol)
    )
    n_heavy_atoms = n_heavy_atoms if n_heavy_atoms is not None else dm.descriptors.n_heavy_atoms(mol)
    charge = charge if charge is not None else dm.descriptors.formal_charge(mol)

    return (
        in_range(mw, 200, 500)
        and in_range(clogp, -5, 5)
        and in_range(n_hba, 0, 10)
        and in_range(n_hbd, 0, 5)
        and in_range(n_rotatable_bonds, 0, 8)
        and in_range(n_heavy_atoms, 15, 50)
    )


def rule_of_chemaxon_druglikeness(
    mol: Union[dm.Mol, str],
    mw: Optional[float] = None,
    clogp: Optional[float] = None,
    n_hba: Optional[float] = None,
    n_hbd: Optional[float] = None,
    n_rotatable_bonds: Optional[int] = None,
    n_rings: Optional[int] = None,
    **kwargs: Any,
) -> bool:
    """
    Compute the druglikeness filter according to chemaxon:

    It computes: `MW < 400 & logP < 5 & HBA <= 10 & HBD <= 5 & ROTBONDS < 5 & RINGS > 0`

    Args:
        mol: input molecule
        mw: precomputed molecular weight.
        clogp: precomputed cLogP.
        n_hba: precomputed number of HBA.
        n_hbd: precomputed number of HBD.
        n_rotatable_bonds: precomputed number of rotatable bonds in the molecule.
        n_rings: precomputed number of rings in the molecule.
        **kwargs: Allow extra arguments for descriptors pre-computation.

    Returns:
        roc: True if molecule is compliant, False otherwise
    """

    if isinstance(mol, str):
        mol = dm.to_mol(mol)

    if mol is None:
        raise ValueError("Molecule is None")

    mw = mw if mw is not None else dm.descriptors.mw(mol)
    clogp = clogp if clogp is not None else dm.descriptors.clogp(mol)
    n_hba = n_hba if n_hba is not None else dm.descriptors.n_hba(mol)
    n_hbd = n_hbd if n_hbd is not None else dm.descriptors.n_hbd(mol)
    n_rotatable_bonds = (
        n_rotatable_bonds if n_rotatable_bonds is not None else dm.descriptors.n_rotatable_bonds(mol)
    )
    n_rings = n_rings if n_rings is not None else dm.descriptors.n_rings(mol)

    return mw < 400 and clogp < 5 and n_hbd <= 5 and n_hba <= 10 and n_rotatable_bonds < 5 and n_rings > 0


def rule_of_egan(
    mol: Union[dm.Mol, str],
    clogp: Optional[float] = None,
    tpsa: Optional[float] = None,
    **kwargs: Any,
) -> bool:
    """
    Compute passive intestinal absorption according to Egan Rules as described in:
    [Egan, William J., Kenneth M. Merz, and John J. Baldwin (2000) Prediction of drug absorption using multivariate statistics.](https://pubs.acs.org/doi/abs/10.1021/jm000292e)

    It computes: `TPSA in [0, 132] & logP in [-1, 6]`

    !!! abstract "Summary of the paper"
        The author built a multivariate statistics model of passive intestinal absorption with robust outlier detection.
        Outliers were identified as being actively transported. They chose PSA and AlogP98 (cLogP), based on consideration of the physical processes
        involved in membrane permeability and the interrelationships and redundancies between other available descriptors.

        **Compounds, which had been assayed for Caco-2 cell permeability, demonstrated a good rate of successful predictions (74−92%)**

    Args:
        mol: input molecule
        clogp: precomputed cLogP.
        tpsa: precomputed TPSA.
        **kwargs: Allow extra arguments for descriptors pre-computation.

    Returns:
        roe: True if molecule is compliant, False otherwise
    """
    if isinstance(mol, str):
        mol = dm.to_mol(mol)

    if mol is None:
        raise ValueError("Molecule is None")

    tpsa = tpsa if tpsa is not None else dm.descriptors.tpsa(mol)
    clogp = clogp if clogp is not None else dm.descriptors.clogp(mol)
    return in_range(tpsa, 0, 132) and in_range(clogp, -1, 6)


def rule_of_pfizer_3_75(
    mol: Union[dm.Mol, str],
    clogp: Optional[float] = None,
    tpsa: Optional[float] = None,
    **kwargs: Any,
) -> bool:
    """
    Compute Pfizer 3/75 Rule for invivo toxicity. It has been described in:
    [Hughes, et al. (2008) Physiochemical drug properties associated with in vivo toxicological outcomes](https://pubmed.ncbi.nlm.nih.gov/18691886/)
    [Price et al. (2009) Physicochemical drug properties associated with in vivo toxicological outcomes: a review](https://pubmed.ncbi.nlm.nih.gov/19519283/)

    It computes: `TPSA >= 75 & logP <= 3`

    !!! abstract "Summary of the paper"
        * In vivo toleration (IVT) studies on 245 preclinical Pfizer compounds found an increased likelihood of toxic events for less polar, more lipophilic compounds.
        * Compounds with `low clogP / high TPSA` are ∼ 2.5 times more likely **to not have any toxity issue** at a fixed concentration of 10 uM (total) or 1 uM (free);
        * Compounds with `high clogP / low TPSA` are ∼ 2.5 times more likely **to have a toxicity finding**; this represents an overall odds >= 6.

    Args:
        mol: input molecule
        clogp: precomputed cLogP.
        tpsa: precomputed TPSA.
        **kwargs: Allow extra arguments for descriptors pre-computation.

    Returns:
        rop: True if molecule is compliant, False otherwise

    """
    if isinstance(mol, str):
        mol = dm.to_mol(mol)

    if mol is None:
        raise ValueError("Molecule is None")

    tpsa = tpsa if tpsa is not None else dm.descriptors.tpsa(mol)
    clogp = clogp if clogp is not None else dm.descriptors.clogp(mol)

    return not (clogp > 3 and tpsa < 75)


def rule_of_gsk_4_400(
    mol: Union[dm.Mol, str],
    mw: Optional[float] = None,
    clogp: Optional[float] = None,
    **kwargs: Any,
) -> bool:
    """
    Compute GSK Rule (4/400) for druglikeness using interpretable ADMET rule of thumb. It has been described in:
    [Gleeson, M. Paul (2008). Generation of a set of simple, interpretable ADMET rules of thumb.](https://pubs.acs.org/doi/10.1021/jm701122q)

    It computes: `MW <= 400 & logP <= 4`.

    !!! abstract "Summary of the paper"
        - The rule are based on a set of consistent structure-property guides determined from an analysis of a number of key
            ADMET assays run within GSK: solubility, permeability, bioavailability, volume of distribution, plasma protein binding,
            CNS penetration, brain tissue binding, P-gp efflux, hERG inhibition, and cytochrome P450 1A2/2C9/2C19/2D6/3A4 inhibition.
        - Conclusion: _It is clear from the analyses reported herein that almost all ADMET parameters deteriorate with either increasing molecular weight,
            logP, or both, with ionization state playing either a beneficial or detrimental affect depending on the parameter in question._

    Args:
        mol: input molecule
        clogp: precomputed cLogP.
        **kwargs: Allow extra arguments for descriptors pre-computation.

    Returns:
        rog: True if molecule is compliant, False otherwise

    """
    if isinstance(mol, str):
        mol = dm.to_mol(mol)

    if mol is None:
        raise ValueError("Molecule is None")

    mw = mw if mw is not None else dm.descriptors.mw(mol)
    clogp = clogp if clogp is not None else dm.descriptors.clogp(mol)

    return mw <= 400 and clogp <= 4


def rule_of_oprea(
    mol: Union[dm.Mol, str],
    n_hba: Optional[float] = None,
    n_hbd: Optional[float] = None,
    n_rotatable_bonds: Optional[int] = None,
    n_rings: Optional[int] = None,
    **kwargs: Any,
) -> bool:
    """
    Computes Oprea's rule of drug likeness obtained by comparing drug vs non drug compounds across multiple datasets.
    The rules have been described in: [Oprea (2000) Property distribution of drug-related chemical databases](https://pubmed.ncbi.nlm.nih.gov/10756480/)

    It computes: `HBD in [0, 2] & HBA in [2, 9] & ROTBONDS in [2,8] and RINGS in [1, 4]`

    !!! abstract "Summary of the paper"
        Seventy percent of the `drug-like` compounds were found between the following limits: `0 ≤ HDO ≤ 2, 2 ≤ HAC ≤ 9, 2 ≤ RTB ≤ 8, and 1 ≤ RNG ≤ 4`

    Args:
        mol: input molecule
        n_hba: precomputed number of HBA.
        n_hbd: precomputed number of HBD.
        n_rotatable_bonds: precomputed number of rotatable bonds in the molecule.
        n_rings: precomputed number of rings in the molecule.
        **kwargs: Allow extra arguments for descriptors pre-computation.

    Returns
        roo: True if molecule is compliant, False otherwise
    """

    if isinstance(mol, str):
        mol = dm.to_mol(mol)

    if mol is None:
        raise ValueError("Molecule is None")

    n_hba = n_hba if n_hba is not None else dm.descriptors.n_hba(mol)
    n_hbd = n_hbd if n_hbd is not None else dm.descriptors.n_hbd(mol)
    n_rotatable_bonds = (
        n_rotatable_bonds if n_rotatable_bonds is not None else dm.descriptors.n_rotatable_bonds(mol)
    )
    n_rings = n_rings if n_rings is not None else dm.descriptors.n_rings(mol)

    return (
        in_range(n_hbd, 0, 2)
        and in_range(n_hba, 2, 9)
        and in_range(n_rotatable_bonds, 2, 8)
        and in_range(n_rings, 1, 4)
    )


def rule_of_xu(
    mol: Union[dm.Mol, str],
    n_hba: Optional[float] = None,
    n_hbd: Optional[float] = None,
    n_rotatable_bonds: Optional[int] = None,
    n_rings: Optional[int] = None,
    n_heavy_atoms: Optional[int] = None,
    **kwargs: Any,
) -> bool:
    """
    Computes Xu's rule of drug likeness as described in:
    [Xu & Stevenson (2000), Drug-like Index: A New Approach To Measure Drug-like Compounds and Their Diversity.](https://pubmed.ncbi.nlm.nih.gov/11045811/)

    It computes `HBD <= 5 & HBA <= 10 & ROTBONDS in [2, 35] & RINGS in [1, 7] & NHeavyAtoms in [10, 50]`.

    !!! note
        A compound's Drug Likeness Index is calculated based upon the knowledge derived from known drugs selected from Comprehensive Medicinal Chemistry (CMC) database.

    Args:
        mol: input molecule
        n_hba: precomputed number of HBA.
        n_hbd: precomputed number of HBD.
        n_rotatable_bonds: precomputed number of rotatable bonds in the molecule.
        n_rings: precomputed number of rings in the molecule.
        n_heavy_atoms: precomputed number of rings in the molecule.
        **kwargs: Allow extra arguments for descriptors pre-computation.

    Returns
        rox: True if molecule is compliant, False otherwise

    """

    if isinstance(mol, str):
        mol = dm.to_mol(mol)

    if mol is None:
        raise ValueError("Molecule is None")

    n_hba = n_hba if n_hba is not None else dm.descriptors.n_hba(mol)
    n_hbd = n_hbd if n_hbd is not None else dm.descriptors.n_hbd(mol)
    n_rotatable_bonds = (
        n_rotatable_bonds if n_rotatable_bonds is not None else dm.descriptors.n_rotatable_bonds(mol)
    )
    n_rings = n_rings if n_rings is not None else dm.descriptors.n_rings(mol)
    n_heavy_atoms = n_heavy_atoms if n_heavy_atoms is not None else dm.descriptors.n_heavy_atoms(mol)

    return (
        n_hba <= 10
        and n_hbd <= 5
        and in_range(n_rotatable_bonds, 2, 35)
        and in_range(n_rings, 1, 7)
        and in_range(n_heavy_atoms, 10, 50)
    )


def rule_of_cns(
    mol: Union[dm.Mol, str],
    mw: Optional[float] = None,
    clogp: Optional[float] = None,
    n_hba: Optional[float] = None,
    n_hbd: Optional[float] = None,
    tpsa: Optional[int] = None,
    **kwargs: Any,
) -> bool:
    """
    Computes drug likeness rule for CNS penetrant molecules as described in:
    [Jeffrey & Summerfield (2010) Assessment of the blood-brain barrier in CNS drug discovery.](https://pubmed.ncbi.nlm.nih.gov/19664709/)

    It computes: `MW in [135, 582]  & logP in [-0.2, 6.1] & TPSA in [3, 118] & HBD <= 3 & HBA <= 5`

    Args:
        mol: input molecule
        mw: precomputed molecular weight.
        clogp: precomputed logP.
        n_hba: precomputed number of HBA.
        n_hbd: precomputed number of HBD.
        tpsa: precomputed TPSA.
        **kwargs: Allow extra arguments for descriptors pre-computation.

    Returns:
        roc: True if molecule is compliant, False otherwise
    """

    if isinstance(mol, str):
        mol = dm.to_mol(mol)

    if mol is None:
        raise ValueError("Molecule is None")

    mw = mw if mw is not None else dm.descriptors.mw(mol)
    clogp = clogp if clogp is not None else dm.descriptors.clogp(mol)
    n_hbd = n_hbd if n_hbd is not None else dm.descriptors.n_hbd(mol)
    n_hba = n_hba if n_hba is not None else dm.descriptors.n_hba(mol)
    tpsa = tpsa if tpsa is not None else dm.descriptors.tpsa(mol)

    return (
        in_range(mw, 135, 582)
        and in_range(clogp, -0.2, 6.1)
        and in_range(tpsa, 3, 118)
        and n_hbd <= 3
        and n_hba <= 5
    )


def rule_of_respiratory(
    mol: Union[dm.Mol, str],
    mw: Optional[float] = None,
    clogp: Optional[float] = None,
    n_hba: Optional[float] = None,
    n_hbd: Optional[float] = None,
    tpsa: Optional[int] = None,
    n_rotatable_bonds: Optional[int] = None,
    n_rings: Optional[int] = None,
    **kwargs: Any,
) -> bool:
    """
    Computes drug likeness rule for Respiratory (nasal/inhalatory) molecules as described in:
    [Ritchie et al. (2009) Analysis of the Calculated Physicochemical Properties of Respiratory Drugs: Can We Design for Inhaled Drugs Yet?](https://pubmed.ncbi.nlm.nih.gov/19275169/)

    It computes: `MW in [240, 520]  & logP in [-2, 4.7] & HBONDS in [6, 12] & TPSA in [51, 135] & ROTBONDS in [3,8] & RINGS in [1,5]`

    Args:
        mol: input molecule
        mw: precomputed molecular weight.
        clogp: precomputed logP.
        n_hba: precomputed number of HBA.
        n_hbd: precomputed number of HBD.
        tpsa: precomputed TPSA.
        n_rotatable_bonds: precomputed number of rotatable bonds in the molecule.
        n_rings: precomputed number of rings.
        **kwargs: Allow extra arguments for descriptors pre-computation.

    Returns:
        roc: True if molecule is compliant, False otherwise
    """

    if isinstance(mol, str):
        mol = dm.to_mol(mol)

    if mol is None:
        raise ValueError("Molecule is None")

    mw = mw if mw is not None else dm.descriptors.mw(mol)
    clogp = clogp if clogp is not None else dm.descriptors.clogp(mol)
    n_hbd = n_hbd if n_hbd is not None else dm.descriptors.n_hbd(mol)
    n_hba = n_hba if n_hba is not None else dm.descriptors.n_hba(mol)
    n_hbonds = n_hbd + n_hba
    tpsa = tpsa if tpsa is not None else dm.descriptors.tpsa(mol)
    n_rotatable_bonds = (
        n_rotatable_bonds if n_rotatable_bonds is not None else dm.descriptors.n_rotatable_bonds(mol)
    )
    n_rings = n_rings if n_rings is not None else dm.descriptors.n_rings(mol)

    return (
        in_range(mw, 240, 520)
        and in_range(clogp, -2, 4.7)
        and in_range(n_hbonds, 6, 12)
        and in_range(tpsa, 51, 135)
        and in_range(n_rotatable_bonds, 3, 8)
        and in_range(n_rings, 1, 5)
    )


def rule_of_generative_design(
    mol: Union[dm.Mol, str],
    mw: Optional[float] = None,
    clogp: Optional[float] = None,
    n_lipinski_hba: Optional[float] = None,
    n_lipinski_hbd: Optional[float] = None,
    tpsa: Optional[float] = None,
    n_rotatable_bonds: Optional[int] = None,
    n_hetero_atoms: Optional[int] = None,
    charge: Optional[float] = None,
    **kwargs: Any,
) -> bool:
    """
    Compute druglikeness rule of generative design.

    This set of rules are proprietary to © Valence Discovery and have been curated to better prioritize molecules
    suggested by generative models for small molecules.

    It computes:

    ```
    MW in [200, 600] & logP < in [-3, 6] & HBD <= 7  & HBA <= 12 & TPSA in [40, 180] &
    ROTBONDS <= 15 & RIGID BONDS <= 30 & N_AROMATIC_RINGS <= 5 & N_FUSED_AROMATIC_RINGS_TOGETHER <= 2 &
    MAX_SIZE_RING_SYSTEM <= 18  & N_CARBONS in [3, 40] & N_HETEROATOMS in [1, 15] & CHARGE in [-2, 2] &
    N_ATOM_CHARGE <= 2 & N_TOTAL_ATOMS < 70 & N_HEAVY_METALS < 1
    ```

    Args:
        mol: input molecule
        mw: precomputed molecular weight.
        clogp: precomputed cLogP.
        n_lipinski_hba: precomputed number of HBA.
        n_lipinski_hbd: precomputed number of HBD.
        tpsa: precomputed TPSA.
        n_rotatable_bonds: precomputed number of rotatable bonds.
        n_hetero_atoms: precomputed number of heteroatoms.
        charge: precomputed charge.
        **kwargs: Allow extra arguments for descriptors pre-computation.

    """
    if isinstance(mol, str):
        mol = dm.to_mol(mol)
    mol = dm.sanitize_mol(mol)

    if mol is None:
        raise ValueError("Molecule is None")

    mw = mw if mw is not None else dm.descriptors.mw(mol)
    clogp = clogp if clogp is not None else dm.descriptors.clogp(mol)
    n_lipinski_hba = n_lipinski_hba if n_lipinski_hba is not None else dm.descriptors.n_lipinski_hba(mol)
    n_lipinski_hbd = n_lipinski_hbd if n_lipinski_hbd is not None else dm.descriptors.n_lipinski_hbd(mol)
    tpsa = tpsa if tpsa is not None else dm.descriptors.tpsa(mol)
    n_rotatable_bonds = (
        n_rotatable_bonds if n_rotatable_bonds is not None else dm.descriptors.n_rotatable_bonds(mol)
    )
    n_hetero_atoms = n_hetero_atoms if n_hetero_atoms is not None else dm.descriptors.n_hetero_atoms(mol)

    # reionize first before computing charge
    with dm.without_rdkit_log():
        standard_mol = dm.standardize_mol(mol, reionize=True, uncharge=False, stereo=False)

    charge = charge if charge is not None else dm.descriptors.formal_charge(standard_mol)
    num_charged_atom = dm.descriptors.n_charged_atoms(standard_mol)
    n_rigid_bonds = dm.descriptors.n_rigid_bonds(mol)
    ring_system = dm.compute_ring_system(mol, include_spiro=False)
    max_size_ring = 0 if len(ring_system) == 0 else max([len(x) for x in ring_system])
    n_carbons = len([at for at in mol.GetAtoms() if at.GetSymbol() == "C"])
    n_aromatic_rings = dm.descriptors.n_aromatic_rings(mol)
    n_fused_aro_rings = n_fused_aromatic_rings(mol, require_all_aromatic=True, pairwise=False)
    n_total_atoms = mol.GetNumAtoms()
    n_heavy_mets = n_heavy_metals(mol)

    return (
        in_range(mw, 200, 600)
        and in_range(clogp, -3, 6)
        and n_lipinski_hbd <= 7
        and n_lipinski_hba <= 12
        and in_range(tpsa, 40, 180)
        and n_rotatable_bonds <= 15
        and n_rigid_bonds <= 30
        and n_aromatic_rings <= 5
        and n_fused_aro_rings <= 2
        and max_size_ring <= 18
        and in_range(n_carbons, 3, 40)
        and in_range(n_hetero_atoms, 1, 15)
        and in_range(charge, -2, 2)
        and num_charged_atom <= 2
        and n_total_atoms < 70
        and n_heavy_mets < 1
    )


def rule_of_generative_design_strict(
    mol: Union[dm.Mol, str],
    mw: Optional[float] = None,
    clogp: Optional[float] = None,
    n_lipinski_hba: Optional[float] = None,
    n_lipinski_hbd: Optional[float] = None,
    tpsa: Optional[float] = None,
    n_rotatable_bonds: Optional[int] = None,
    n_hetero_atoms: Optional[int] = None,
    charge: Optional[float] = None,
    **kwargs: Any,
) -> bool:
    """
    Compute the strict version of the druglikeness rule of generative design, which includes long aliphatic chains.

    This set of rules are proprietary to © Valence Discovery and have been curated to better prioritize molecules
    suggested by generative models for small molecules.

    It computes:

    ```
    MW in [200, 600] & logP < in [-3, 6] & HBD <= 7  & HBA <= 12 & TPSA in [40, 180] &
    ROTBONDS <= 15 & RIGID BONDS <= 30 & N_AROMATIC_RINGS <= 5 & N_FUSED_AROMATIC_RINGS_TOGETHER <= 2 &
    MAX_SIZE_RING_SYSTEM <= 18  & N_CARBONS in [3, 40] & N_HETEROATOMS in [1, 15] & CHARGE in [-2, 2] &
    N_ATOM_CHARGE <= 2 & N_TOTAL_ATOMS < 70 & N_HEAVY_METALS < 1 & N_STEREO_CENTER <= 3 &
    HAS_NO_SPIDER_SIDE_CHAINS & FRACTION_RING_SYSTEM >= 0.25
    ```

    By default SPIDER_SIDE_CHAINS are defined as having at least 2 'chains' of  >=4 consecutif atoms in side chains (not part of any ring system)


    Args:
        mol: input molecule
        mw: precomputed molecular weight.
        clogp: precomputed cLogP.
        n_lipinski_hba: precomputed number of HBA.
        n_lipinski_hbd: precomputed number of HBD.
        tpsa: precomputed TPSA.
        n_rotatable_bonds: precomputed number of rotatable bonds.
        n_hetero_atoms: precomputed number of heteroatoms.
        charge: precomputed charge.
        **kwargs: Allow extra arguments for descriptors pre-computation.
    """

    mol = dm.to_mol(mol)
    mol = dm.sanitize_mol(mol)

    if mol is None:
        raise ValueError("Molecule is None")

    generative_rule = rule_of_generative_design(
        mol=mol,
        mw=mw,
        clogp=clogp,
        n_lipinski_hba=n_lipinski_hba,
        n_lipinski_hbd=n_lipinski_hbd,
        tpsa=tpsa,
        n_rotatable_bonds=n_rotatable_bonds,
        n_hetero_atoms=n_hetero_atoms,
        charge=charge,
        **kwargs,
    )

    n_stereo_center = dm.descriptors.n_stereo_centers(mol)
    has_spider_appendages = has_spider_chains(mol)
    good_fraction_ring_system = fraction_atom_in_scaff(mol) >= 0.25

    return (
        generative_rule and good_fraction_ring_system and not has_spider_appendages and n_stereo_center <= 3
    )
