import datamol as dm
import pandas as pd

from rdkit.Chem import rdfiltercatalog  # type: ignore

from medchem.groups import ChemicalGroup
from medchem.groups import list_default_chemical_groups
from medchem.groups import list_functional_group_names
from medchem.groups import get_functional_group_map

from medchem.utils.loader import get_data_path


def test_list_default_chemical_groups():
    assert len(list_default_chemical_groups()) > 0
    assert all([isinstance(ll, str) for ll in list_default_chemical_groups()])

    with_hierarchy = list_default_chemical_groups(hierarchy=True)
    assert len(with_hierarchy) > 0
    assert all([isinstance(ll, str) for ll in with_hierarchy])
    assert all(["." in ll for ll in with_hierarchy])


def test_list_functional_group_names():
    group_names = list_functional_group_names()

    assert len(group_names) > 0
    assert all([isinstance(ll, str) for ll in group_names])


def test_get_functional_group_map():
    group_map = get_functional_group_map()
    assert len(group_map) > 0
    assert isinstance(group_map, dict)


def test_chemical_group_filter():
    group = ChemicalGroup(groups="rings_in_drugs")
    assert len(group.data) == 92
    group.filter(["azole", "3,4-dihydroquino-2(1H)-one"], fuzzy=False)
    assert len(group.data) == 1

    group = ChemicalGroup(groups="rings_in_drugs")

    assert len(group.data) == 92
    group.filter(["azole", "3,4-dihydroquino-2(1H)-one"], fuzzy=True)
    assert len(group.data) == 12

    group = ChemicalGroup(groups="rings_in_drugs")

    assert len(group.data) == 92
    group.filter([])
    assert len(group.data) == 92


def test_chemical_group_attributes():
    group = ChemicalGroup(groups="rings_in_drugs")
    assert isinstance(group.name, list)
    assert isinstance(group.smiles, list)
    assert isinstance(group.mols, list)
    assert isinstance(group.smarts, list)
    assert isinstance(group.mol_smarts, list)
    assert isinstance(group.dataframe, pd.DataFrame)


def test_chemical_group_list():
    group = ChemicalGroup(groups="rings_in_drugs")
    assert isinstance(group.list_groups(), list)
    assert isinstance(group.list_hierarchy_groups(), list)


def test_chemical_group_catalogs():
    group = ChemicalGroup(groups="rings_in_drugs")
    assert isinstance(group.get_catalog(), rdfiltercatalog.FilterCatalog)


def test_chemical_group():
    c_group = ChemicalGroup(groups="rings_in_drugs")
    rings_numbers = 92
    assert len(c_group) == rings_numbers
    assert len(c_group.smarts) == rings_numbers

    out_false = c_group.has_match("CCCCCCCCCC")
    assert out_false is False

    out_true = c_group.has_match("C1CCCCC1")
    assert out_true is True


def test_chemical_group_query():
    c_group = ChemicalGroup(groups="rings_in_drugs")
    mol = dm.to_mol("CCS(=O)(=O)N1CC(C1)(CC#N)N2C=C(C=N2)C3=C4C=CNC4=NC=N3")

    # we have two rings of interest for smarts
    out = c_group.get_matches(mol, use_smiles=False)

    assert out is not None
    assert out.shape[0] == 2
    assert set(out.group.unique()) == {"rings_in_drugs"}
    assert set(out.name.unique()) == {"diazine", "1H-pyrrole"}

    # however, if we use smiles, we would have 3
    # this is a bug that needs to be fixed
    out_smiles = c_group.get_matches(mol, use_smiles=True)
    assert out_smiles is not None
    assert out_smiles.shape[0] == 3
    assert set(out_smiles.group.unique()) == {"rings_in_drugs"}
    assert set(out_smiles.name.unique()) == {"diazine", "1H-pyrrole", "1H-pyrazole"}


def test_external_bank():
    c_group = ChemicalGroup(groups_db=get_data_path("smarts_bank.csv"))
    mol = dm.to_mol("CCS(=O)(=O)N1CC(C1)(CC#N)N2C=C(C=N2)C3=C4C=CNC4=NC=N3")
    out = c_group.get_matches(mol, use_smiles=False)
    expected_match = set(["HBA", "HBD", "Hydrogen", "SP3 Nitrogen", "SP2 Carbon"])

    assert out is not None
    assert expected_match.issubset(set(out["name"].values)) is True
