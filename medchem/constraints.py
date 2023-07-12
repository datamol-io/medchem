from typing import Dict
from typing import Callable
from typing import List

from loguru import logger

import datamol as dm

from rdkit.Chem.rdchem import SubstructMatchParameters
from rdkit.Chem.rdmolfiles import MolFragmentToSmiles


class Constraints:
    """Complex query system for performing substructure matches with additional constraints

    !!! example

        ```python
        mol1 = dm.to_mol("CN(C)C(=O)c1cncc(C)c1")
        mol2 = dm.to_mol("c1ccc(cc1)-c1cccnc1")
        core = dm.from_smarts("c1cncc([*:1])c1")
        [atom.SetProp("query", "my_constraints") for atom in core.GetAtoms() if atom.GetAtomMapNum() == 1]
        constraint_fns = dict(my_constraints=lambda x: dm.descriptors.n_aromatic_atoms(x) > 0)
        constraint = Constraints(core, constraint_fns)
        matches = [constraint(mol1), constraint(mol2)] # False, True
        ```
    """

    def __init__(self, core: dm.Mol, constraint_fns: Dict[str, Callable], prop_name: str = "query"):
        """Initialize the constraint matcher

        Args:
            core: the scaffold/query molecule to match against. Needs to be a molecule
            constraint_fns: a dictionary of constraint functions that defines the required constraints after the substructure match
            prop_name: the property name to use in the match at each atom defined by the core
                for further matches against the constraints functions
        """
        self.core = core
        self.prop_name = prop_name
        self.constraint_fns = constraint_fns
        self.atom_to_query = self._initialize()

    def _initialize(self):
        """Initialize the constraint matcher"""
        atom_to_query = dict()
        for a in self.core.GetAtoms():
            if a.HasProp(self.prop_name) and str(a.GetProp(self.prop_name)) in self.constraint_fns:
                atom_to_query[a.GetIdx()] = a.GetProp(self.prop_name)
        return atom_to_query

    @staticmethod
    def validate(mol: dm.Mol, constraints: List["Constraints"]):
        """Validate a list of constraint object against a molecule

        Args:
            mol: the molecule object
            constraints: list of Contraints object to validate against the molecule
        """
        mol = dm.to_mol(mol)

        if mol is None:
            raise ValueError("Input molecule is None")

        for constraint in constraints:
            if not isinstance(constraint, Constraints):
                raise ValueError("Input constraint should be an instance of Constraints")

            if not constraint(mol):
                return False

        return True

    def get_matches(self, mol: dm.Mol, multiple: bool = True):
        """Get matches that respect the constraints in the molecules

        Args:
            mol: input molecule
            multiple: if True, return all the matches, if False, return the first match
        """
        if not isinstance(mol, dm.Mol):
            mol = dm.to_mol(mol)

        params = SubstructMatchParameters()
        params.setExtraFinalCheck(self._check_final)
        matches = mol.GetSubstructMatches(self.core, params)

        if multiple:
            return matches

        return matches[0]

    def has_match(self, mol: dm.Mol):
        """Check if input molecule respect the constraints

        Args:
            mol: input molecule
        """
        return len(self.get_matches(mol)) > 0

    def __call__(self, mol: dm.Mol):
        """Check if input molecule respect the constraints

        Args:
            mol: input molecule
        """
        return self.has_match(mol)

    def _check_final(self, mol: dm.Mol, vect: List[int]):
        """
        Perform a breadth first search over current matches.

        !!! note
            This function is not designed to be used directly but this setup is required to make things work
            Always use `validate` or `match


        Args:
            mol: input molecule
            vect: list of atom indexes to check

        """
        seen = [0] * mol.GetNumAtoms()
        for idx in vect:
            seen[idx] = 1

        # loop over the atoms we care about:
        for idx, qfn in self.atom_to_query.items():
            midx = vect[idx]
            atom = mol.GetAtomWithIdx(midx)
            # now do a breadth-first search from that atom, checking
            # all of its neighbors that aren't in the substructure
            # query:
            stack = [atom]
            connected_atoms = set()
            while stack:
                atom = stack.pop(0)
                connected_atoms.add(atom.GetIdx())
                for neigh in atom.GetNeighbors():
                    neigh_idx = neigh.GetIdx()
                    if not seen[neigh_idx] and neigh_idx not in connected_atoms:
                        stack.append(neigh)
            is_ok = False
            try:
                submol = dm.to_mol(MolFragmentToSmiles(mol, atomsToUse=connected_atoms))
                is_ok = self.constraint_fns[qfn](submol)
            except Exception as e:
                logger.error(e)
                # we can't raise here

            if not is_ok:
                return False
        return True
