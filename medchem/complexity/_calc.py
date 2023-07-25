import math

from rdkit.Chem.rdmolops import GetMolFrags
from rdkit.Chem.rdmolops import FindPotentialStereo
from rdkit.Chem import FindMolChiralCenters

import datamol as dm


def WhitlockCT(
    mol: dm.Mol,
    ringval: float = 4,
    unsatval: float = 2,
    heteroval: float = 1,
    chiralval: float = 2,
):
    """
    A chemically intuitive measure for molecular complexity. This complexity measure
    has been described in : H. W. Whitlock, J. Org. Chem., 1998, 63, 7982-7989.
    Benzyls, fenyls, etc. are not treated at all.

    On the zinc 15 commercially available dataset, the range of this score is [0, 172] with a median of 25

    Args:
        mol: The input molecule.
        ringval: The contribution of rings
        unsatval: The contribution of the unsaturated bond.
        heteroval: The contribution of the heteroatom.
        chiralval: The contribution of the chiral center.
    """

    # possibly TODO phenyls and others are considered as protecting group thus could be removed (keyword argument removePhenyl)
    mol = dm.to_mol(mol, add_hs=True)

    if mol is None:
        raise ValueError("Invalid molecule")

    n_bonds = 0
    n_unsat = 0
    n_hetero = 0
    n_chiral = len(FindMolChiralCenters(mol, includeUnassigned=True))
    n_arom = 0

    bonds = mol.GetBonds()
    atoms = mol.GetAtoms()

    for bond in bonds:
        n_bonds += 1
        btype = bond.GetBondType()
        if btype == dm.DOUBLE_BOND:
            n_unsat += 1
            n_bonds += 1
        elif btype == dm.TRIPLE_BOND:
            n_unsat += 2
            n_bonds += 2
        elif btype == dm.AROMATIC_BOND:
            n_arom += 1

    for atom in atoms:
        if atom.GetAtomicNum() != 1 and atom.GetAtomicNum() != 6:
            n_hetero += 1

    n_rings = len(bonds) - len(atoms) + len(GetMolFrags(mol))
    complexity = ringval * n_rings + unsatval * n_unsat + heteroval * n_hetero + chiralval * n_chiral

    return complexity


def BaroneCT(mol: dm.Mol, chiral: bool = False):
    """
    Compute the Barone complexity measure for a molecule as described in:
    [R. Barone and M. Chanon, J. Chem. Inf. Comput. Sci., 2001, 41 (2), pp 269–272](https://pubs.acs.org/doi/abs/10.1021/ci000145p)

    On zinc 15 commercially available dataset, the range of this score is [30, 4266] with a median of 538

    Args:
        mol: The input molecule.
        chiral: Whether to include chirality in the calculation.
    """
    mol = dm.to_mol(mol)
    cmpx = 0
    rinfo = mol.GetRingInfo()
    for aring in rinfo.AtomRings():
        cmpx += len(aring) * 6
    for atom in mol.GetAtoms():
        deg = atom.GetExplicitValence()
        if deg == 1:
            cmpx += 3
        elif deg == 2:
            cmpx += 6
        elif deg == 3:
            cmpx += 12
        else:
            cmpx += 24
        if atom.GetAtomicNum() != 6:
            cmpx += 6
        else:
            cmpx += 3
    if chiral:
        cmpx += 20 * len(FindMolChiralCenters(mol, includeUnassigned=True))
    return cmpx


_SMCM_ENs = {
    5: 0.851,
    6: 1,
    7: 1.149,
    8: 1.297,
    9: 1.446,
    15: 1.086,
    16: 1.235,
    17: 1.384,
    35: 1.244,
    53: 1.103,
}

_SMCM_BNs = {
    dm.SINGLE_BOND: {
        6: {
            6: 1.0,
            7: 0.857,
            8: 0.75,
            9: 0.667,
            15: 0.4,
            16: 0.375,
            17: 0.353,
            35: 0.171,
            53: 0.113,
        },
        7: {7: 0.735, 8: 0.643},
        8: {16: 0.281},
    },
    dm.DOUBLE_BOND: {
        6: {
            6: 0.5,
            7: 0.429,
            8: 0.375,
            16: 0.188,
        },
        7: {7: 0.367, 8: 0.321},
        8: {16: 0.141},
    },
    dm.TRIPLE_BOND: {6: {6: 0.333, 7: 0.286}},
    dm.AROMATIC_BOND: {
        6: {6: 0.667, 7: 0.571, 16: 0.25},
        7: {7: 0.49, 8: 0.423},
    },
}
"""
Table of bond
atom    atom    single    double    triple    aromatic
C        C      1.000     0.500     0.333     0.667
C        N      0.857     0.429     0.286     0.571
C        O      0.750     0.375
C        F      0.667
C        P      0.400
C        S      0.375     0.188               0.250
C        Cl     0.353
C        Br     0.171
C        I      0.113
N        N      0.735     0.367               0.490
N        O      0.643     0.321               0.423
O        S      0.281     0.141
"""


def _SMCM_GetBondScore(btype, a1, a2):
    ats = sorted([a1, a2])
    s = 0
    try:
        s = _SMCM_BNs[btype][ats[0]][ats[1]]
    except KeyError:
        pass
    return s


# SMARTS patterns for SMCM metric which lower structure complexity
_SMCM_SMARTSs = [
    dm.from_smarts(patt)
    for patt in (
        "[$([C!D1]-@C(=O)-@[N!a])]",  # peptides, cyclic
        "[$([C!D1]-!@C(=O)-!@[N!a])]",  # peptides, acyclic
        "[$([#6]=,:[*R!C!c]=,:[#6R]~@[*R!C!c])]",  # C attached to Hetero in ring
        "[$([A!C!H](-!@[#6])-!@[#6])]",  # hetero attached to 2 carbons not in ring
        "[$([#7]=[C!$(*a)]-!@[N!H0])]",  # amidine, guanidine
        "[$([#8]=[#6H0]-!@[N!H2])!$(NC[!#6])]",  # nonterminal amide
        "[#8]=[#6](-!@[NH1,NH0])-!@[N!H2]",  # nonterminal urea
        "[$(O=[CD3]([#6])[#6])!$([!#6]=CC=O)]",  # ketone, not diketo
        "[$([#16X2v2]([#6])[#6])]",  # thio ether
        "[$([#8X2H0][#6]=[#8D1])!$(O~C(~O)~O)]",  # carboxyl ester, not carbonate
        "[$([#8X2v2]([#6])[#6])]",  # ether oxygen
        "C1:c-@C-@N-@C-@C-@O-@C-@c:c-1",  # SMARTS from ICCB's Diversity-Oriented
        "[N,O,C]C(=O)C1=C[C@H](*)C[C@H](O)O1",  # Synthesis approach19
        "C[C@H]1O[C@H]~3O[C@H]~C~2~C~C[C@H]([C@@H]~1)[C@@H]~23",  #
        "C[C@@H]1O[C@@H]~2O[C@H][C@@H][C@@H]~3~C~C~C~1[C@@H]~23",  #
        "C[C@H]2C[C@]14~C~C~C~C[C@H]1Oc3cccc(CN2)c34",  #
        "[#6][C@H]1C[C@@H]([#6])O[C@@H](-a)O1",  #
        "C2=CC[C@@H]1C(=O)~*~*~C(=O)[C@@H]1[C@@H]2-a",  #
        "C2=CC[C@@H]1C(=O)~*~C(=O)[C@@H]1[C@@H]2-a",  #
        "a-[CH,CH2;R0;0*]",  # from ref 17
        "[R;0*]-[CH2R0,NHR0,OR0;0*]-[R]",  #
        "*-[CD3H,ND2;R0;0*](-a)-a",  #
        "[a]-&!@[a]",  #
        "[NR;0*]-[CD3R0;0*](=O)-[R]",  #
        "[NR;0*]-[CD2R0;0*]-[R]",  #
        "[NR;0*]-[CD2R0;0*]-[CD2,CD3,OD2,ND2,ND3,aD2,aD3]",  #
        "a-[NHR0]-[CR0;0*](=O)-[OR0,NR0,0*]",  #
        "[CR,NR]=[CR]-&!@[a]",  #
        "[$([#6](-!@[#7])-!@[#7!H0])]",  # NHCN not in ring
        "[$([#6](@[#7])@[#7!H0])]",  # NHCN in ring
        "[$([A!C!H](-@[#6])-@[#6])]",  # hetero attached to 2 carbons in ring
        "[$([#6](=[#8])-[#6]=[!#6])]",  # diketo, keto-Ryl
        "[$([#16X3v4,#16X4v6]([#6])[A])]",  # hypervalent sulfur
        "[$(C(~O)(~O)~O)]",  # carbonate
    )
]


def SMCM(mol: dm.Mol):
    """
    Compute synthetic and molecular complexity as described in:
    [TK Allu, TI Oprea, J. Chem. Inf. Model. 2005, 45(5), pp. 1237-1243](https://pubs.acs.org/doi/10.1021/ci0501387)

    On zinc 15 commercially available dataset, the range of this score is [1.93, 192.00] with a median of 42.23

    Args:
        mol: the input molecule
    """
    mol = dm.to_mol(mol)
    ats = mol.GetAtoms()
    bonds = mol.GetBonds()
    a_score = 0
    c = 0
    for a in ats:
        if a.GetAtomicNum() == 6:
            c += 1
        if a.GetAtomicNum() in _SMCM_ENs:
            a_score += _SMCM_ENs[a.GetAtomicNum()]
        else:
            a_score += 1
            # TODO log it or print this in the same way as RDKit
            print("Unsupported atom #{}".format(a.GetAtomicNum()))

    b_score = 0
    for b in bonds:
        b_score += _SMCM_GetBondScore(
            b.GetBondType(),
            b.GetBeginAtom().GetAtomicNum(),
            b.GetEndAtom().GetAtomicNum(),
        )

    chiral_a = len(FindPotentialStereo(mol)) * 2
    mot_score = sum([0] + [len(mol.GetSubstructMatches(patt)) for patt in _SMCM_SMARTSs])
    score = a_score + b_score + chiral_a - mot_score

    return round(score, 3)


def _AWC(k: int, atom: int, table: list, neighbors: list):
    """Compute walk count for atom"""
    if atom not in table[k]:
        table[k][atom] = 0
        for i in neighbors[atom]:
            table[k][atom] += _AWC(k - 1, i, table, neighbors)
    return table[k][atom]


def TWC(mol: dm.Mol, log10: bool = True):
    """Compute total walk count in a molecules as proxy for complexity. This score is described in:
    [Gerta Rucker and Christoph Rucker, J. Chem. Inf. Comput. Sci. 1993, 33, 683-695](https://pubs.acs.org/doi/pdf/10.1021/ci00015a005)

    The total walk count is defined as: $twc = \\frac{1}{2} \sum_{k=1}^{n-1} \sum_{i}^{Natoms} \\text{awc}(k,i)$
    where $\\text{awc}(k,i)$ is the number of walk of length `k` starting at atom `i`.

    On zinc 15 commercially available dataset, the range of this score is [1.20, 39.08] with a median of 10.65

    Args:
        mol: the input molecule
        log10: whether to return the log10 of the values
    """
    mol = dm.to_mol(mol)
    twc = 0
    table = []
    neighbors = []
    atoms = mol.GetAtoms()
    for k in range(len(atoms)):
        table.append({})
    for atom in atoms:
        neighbors.append([i.GetIdx() for i in atom.GetNeighbors()])
        table[0][atom.GetIdx()] = len(neighbors[atom.GetIdx()])

    for k in range(len(atoms) - 1):
        for i in atoms:
            twc += _AWC(k, i.GetIdx(), table, neighbors)
    twc /= 2
    if log10:
        try:
            return math.log10(twc)
        except ValueError:
            return float("nan")
    return twc
