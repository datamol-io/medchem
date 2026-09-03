import datamol as dm
import pandas as pd

from medchem.utils.loader import get_data_path
from medchem.catalogs import merge_catalogs
from medchem.catalogs import NamedCatalogs
from medchem.catalogs import catalog_from_smarts


def test_catalog_merge():
    data = dm.data.freesolv()
    mols = data["smiles"].apply(dm.to_mol).tolist()

    tox = NamedCatalogs.tox(pains_a=False, pains_b=False, pains_c=False, nih=True, brenk=True)
    nih = NamedCatalogs.nih()
    brenk = NamedCatalogs.brenk()

    merged_tox = merge_catalogs(nih, brenk)

    toxic1 = [tox.HasMatch(m) for m in mols]
    toxic2 = [merged_tox.HasMatch(m) for m in mols]
    assert toxic1 == toxic2


def test_catalog_from_smarts():
    smarts_bank = pd.read_csv(get_data_path("smarts_bank.csv"))
    custom_catalog = catalog_from_smarts(
        smarts_bank["smarts"],
        smarts_bank["name"],
        entry_as_inds=False,
    )

    assert custom_catalog.GetNumEntries() == 30


def test_NamedCatalogs():
    for catalog_name in [
        "alarm_nmr",
        "alerts",
        "alphascreen",
        "bms",
        "bredt",
        "brenk",
        "carcinogen",
        "chelator",
        "chemical_groups",
        "dnabinder",
        "dundee",
        "electrophilic",
        "glaxo",
        "gst_hitters",
        "his_hitters",
        "hitters",
        "inpharmatica",
        "ld50_oral",
        "lint",
        "luciferase",
        "mlsmr",
        "nibr",
        "nih",
        "pains",
        "pains_a",
        "pains_b",
        "pains_c",
        "reactive_unstable_toxic",
        "schembl",
        "skin",
        "tox",
        "toxicophore",
        "unstable_graph",
        "zinc",
    ]:
        catalog = getattr(NamedCatalogs, catalog_name)()
        assert catalog.GetNumEntries() > 0


def test_toxicophore_michael_acceptor_requires_an_accepting_double_bond():
    catalog = NamedCatalogs.toxicophore()

    def matches(smiles):
        return {
            match.filterMatch.GetName().split("||")[1].split("_min(")[0]
            for match in catalog.GetFilterMatches(dm.to_mol(smiles))
        }

    assert "michael_acceptors_double" not in matches("C=CCN")
    assert "michael_acceptors_double" not in matches("C=CCS")
    assert "michael_acceptors_double" in matches("C=CC=N")
    assert "michael_acceptors_double" in matches("C=CC=S")
    assert "michael_acceptors_triple" not in matches("C#CCN")
    assert "michael_acceptors_triple" not in matches("C#CCS")
    assert "michael_acceptors_triple" in matches("C#CC=N")
    assert "michael_acceptors_triple" in matches("C#CC=S")
