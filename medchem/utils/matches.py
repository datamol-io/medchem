from __future__ import annotations
from collections.abc import Mapping
from typing import Dict
from typing import Callable
from typing import List
import datamol as dm
from rdkit import Chem


class Constraints:
    """Complex query system for matches with additional constraints

    !!! example
        >>> import datamol as dm
        >>> mol = dm.from_smiles('CC(=O)O')
        >>> core =
    """

    def __init__(
        self,
        core: dm.Mol,
        constraint_fns: Dict[Callable],
        prop_name: str = "query",
    ):
        """Initialize the constraint matcher
        Args:
            core: the scaffold/query molecule to match against. Needs to be a molecule
            constraint_fns: a dictionary of constraints functions
            prop_name: the property name to use in the match at each atom defined by the core
                for further matches against the constraints functions
        """
        self.core = core
        self.prop_name = prop_name
        self.constraint_fns = self._get_constraints(constraint_fns)
        self.atom_to_query = self._initialize()

    def _initialize(self):
        """Initialize the constraint matcher"""
        atom_to_query = dict()
        for a in self.core.GetAtoms():
            if (
                a.HasProp(self.prop_name)
                and a.GetProp(self.prop_name) in self.constraint_fns
            ):
                atom_to_query[a.GetIdx()] = a.GetProp(self.prop_name)
        return atom_to_query

    @classmethod
    def _get_constraints(cls, constraint_fns):
        """Get constraints list for current instance

        Args:
            constraint_fns: list of constraint functions (callable)
        """
        constraint_fns = dict()
        if constraint_fns and isinstance(constraint_fns, Mapping):
            constraint_fns.update(constraint_fns)
        return constraint_fns

    @staticmethod
    def validate(mol, constraints: List[Constraints]):
        """Validate a list of constraint object against a molecule
        Args:
            mol: the molecule object
            constraints: list of Contraints object to validate against the molecule
        """
        mol = dm.to_mol(mol)
        for cons in constraints:
            if not isinstance(cons, Constraints):
                raise ValueError(
                    "Input constraint should be an instance of Constraints"
                )
            valid = len(cons.match(mol)) > 0
            if not valid:
                return False
        return True

    def match(self, mol):
        """Check if input molecule respect the constraints

        Args:
            mol: input molecule
        """
        if not isinstance(mol, Chem.Mol):
            mol = dm.to_mol(mol)
        params = Chem.SubstructMatchParameters()
        params.setExtraFinalCheck(self)
        matches = mol.GetSubstructMatches(self.core, params)
        return matches

    def __call__(self, mol: dm.Mol, vect: List[int]):
        """
        Perform a breadth search over current matches.

        !!! note
            This function is not designed to be used directly but this setup is required to make things work
            Always use `validate` or `match


        Args:
            mol: input molecule
            vecte: list of atom indexes to check

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
                submol = dm.to_mol(
                    Chem.MolFragmentToSmiles(mol, atomsToUse=connected_atoms)
                )
                is_ok = all(self.constraint_fns[qfn].applyFunctions(submol))
            except:
                pass
            if not is_ok:
                return False
        return True
