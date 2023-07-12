import pytest


import datamol as dm

from medchem.constraints import Constraints


def test_constraints():
    def my_constraint(mol):
        # we want to either (have phenol) OR (have less than 7 atoms and not ring)
        return mol.HasSubstructMatch(dm.to_mol("Oc1ccccc1")) or (
            mol.GetNumAtoms() < 7 and dm.descriptors.n_rings(mol) < 1
        )

    smiles = [
        "CN(C)C(=O)c1cncc(C)c1",  # match, n_atoms < 7 and no ring
        "Cc1cncc(CC2CCCCC2)c1",  # not match, n_atoms < 7 but ring
        "Cc1cncc(c1)-c1ccc(O)cc1",  # match phenol
        "Cc1cncc(c1)-c1cccc2nc[nH]c12",  # no match n_atoms >= 7
    ]
    expected_results = []
    mols = [dm.to_mol(x) for x in smiles]
    core = dm.from_smarts("[C;H3]c1cncc([*:1])c1")

    assert core is not None

    # now let's set the constraints query
    for atom in core.GetAtoms():
        # we add a recursive query to check again on any match that starts with this atom position
        if atom.GetAtomMapNum() == 1:
            atom.SetProp("query", "my_constraint")

    constraint_fns = {"my_constraint": my_constraint}
    constraint = Constraints(core, constraint_fns)

    matches = [constraint(mol) for mol in mols]
    expected_results = [True, False, True, False]
    assert expected_results == matches

    one_match = constraint.get_matches(mols[0], multiple=False)
    assert isinstance(one_match[0], int)
    assert len(one_match) > 0


def test_constraints_validate():
    def my_constraint(mol):
        return mol.HasSubstructMatch(dm.to_mol("Oc1ccccc1")) or (
            mol.GetNumAtoms() < 7 and dm.descriptors.n_rings(mol) < 1
        )

    core = dm.from_smarts("[C;H3]c1cncc([*:1])c1")

    assert core is not None

    for atom in core.GetAtoms():
        if atom.GetAtomMapNum() == 1:
            atom.SetProp("query", "my_constraint")

    # Create the constraints object
    constraint_fns = {"my_constraint": my_constraint}
    constraint = Constraints(core, constraint_fns)

    mol = dm.to_mol("CN(C)C(=O)c1cncc(C)c1")

    assert Constraints.validate(mol, constraints=[constraint]) is True


def test_constraints_validate_fails():
    mol = dm.to_mol("CN(C)C(=O)c1cncc(C)c1")

    with pytest.raises(ValueError):
        Constraints.validate(mol, constraints=["xxxxxx"])  # type: ignore

    with pytest.raises(ValueError):
        Constraints.validate(None, constraints=[])  # type: ignore
