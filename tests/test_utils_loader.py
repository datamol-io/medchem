import importlib.resources

import pandas as pd

from medchem.utils.loader import get_data_path
from medchem.utils.loader import get_grammar


def test_get_data_path():
    df = pd.read_csv(get_data_path("medchem_rule_list.csv"))
    assert df.columns.tolist() == ["name", "rules", "description"]


def test_get_grammar_default():
    with open(get_grammar()) as f:
        grammar_text = f.read()

    assert "bool_term (OR_OP bool_term)" in grammar_text


def test_get_grammar_as_string():
    grammar_text = get_grammar(as_string=True)
    assert "bool_term (OR_OP bool_term)" in grammar_text


def test_get_grammar_custom():
    path = importlib.resources.files("medchem.data").joinpath("grammar.lark")
    path = str(path)

    with open(get_grammar(path)) as f:
        grammar_text = f.read()

    assert "bool_term (OR_OP bool_term)" in grammar_text
