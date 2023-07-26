import pytest

import operator

import datamol as dm
import numpy as np

from lark import Lark, exceptions

from medchem.query import QueryParser
from medchem.query import QueryOperator
from medchem.query import QueryFilter
from medchem.query import EvaluableQuery

from medchem.utils.loader import get_grammar


def get_test_grammar():
    return Lark(get_grammar(as_string=True), parser="lalr", transformer=QueryParser())


def test_language():
    grammar = get_test_grammar()

    bad_queries = [
        "HASPROP(tpsa > 120 )",  # missing quotes, UnexpectedCharacters
        """HASPROP("tpsa" > 120 ) OR OR HASSUPERSTRUCTURE("Cc1cc2ccccc2cc1C")""",  # wrong syntax with 2 OR, UnexpectedCharacters
        """HASPROP("tpsa" > 120 ) OR TRUE("wrong")""",  # wrong syntax, UnexpectedCharacters
        """HASPROP("tpsa" > 120 ) OR HASSUPERSTRUCTURE("Cc1cc2ccccc2cc1C") OR""",  # wrong syntax, wrong terminal node (OR), KeyError
        """HASPROP("tpsa" > 120 ) OR  (HASSUPERSTRUCTURE("Cc1cc2ccccc2cc1C")""",  # wrong syntax, parenthesis mismatch, KeyError
        """HASPROP("tpsa" > 120 ) OR HASSUPERSTRUCTURE("Cc1cc2ccccc2cc1C") TRUE""",  # wrong syntax, wrong terminal node (TRUE) , KeyError
        """HASPROP("tpsa" > 120 ) OR () OR HASSUPERSTRUCTURE("Cc1cc2ccccc2cc1C")""",  # wrong syntax, empty statement, UnexpectedCharacters
        """HASPROP("tpsa" > 120) OR HASSUBSTRUCTURE("[OH]", True, >, 3)""",  # unexpected token '>' at wrong place, UnexpectedCharacters
        """HASPROP("tpsa" + 120) OR HASSUBSTRUCTURE("[OH]")"""  # unexpected token '+', UnexpectedCharacters
        """HASPROP("tpsa") OR HASSUBSTRUCTURE("[OH]")""",  # incomplete HASPROP statement UnexpectedCharacters
        """HASPROP("tpsa" > 120) OR HASSUBSTRUCTURE("[OH]" True max 3)""",  # mixing comma in substructure statement token, UnexpectedCharacters
        """HASPROP("tpsa" > 120) OR HASSUBSTRUCTURE("[OH]", min, True, 3)""",  # wrong order of arguments, UnexpectedCharacters
        """HASPROP("tpsa" > 120) OR BIND("herg")""",  # BIND is not a valid statement
        """HASALERT("tpsa" > 120)""",  # this is not a valid statement
    ]
    valid_queries = [
        """HASPROP("tpsa" > 120 ) OR HASSUPERSTRUCTURE("ls -lah /")""",  # This is valid query, but "ls -lah" will never be evaluated
        """HASPROP("tpsa" > 120 ) OR HASSUPERSTRUCTURE("Cc1cc2ccccc2cc1C")""",
        """HASPROP("tpsa" > 120 ) OR HASSUPERSTRUCTURE("Cc1cc2ccccc2cc1C") OR NOT TRUE""",  # this is valid, but should be simplified,
        """HASPROP("tpsa" > 120 ) OR HASSUPERSTRUCTURE("Cc1cc2ccccc2cc1C") OR TRUE""",
        """HASPROP("tpsa",  >, 120 ) OR HASSUPERSTRUCTURE("Cc1cc2ccccc2cc1C")""",
        """HASPROP("tpsa" > 120) OR HASSUBSTRUCTURE("[OH]")""",  # min match, not smarts, limit=1
        """HASPROP("tpsa" > 120) OR HASSUBSTRUCTURE("[OH]", True)""",  # min match, not smarts, limit=1
        """HASPROP("tpsa" > 120) OR HASSUBSTRUCTURE("[OH]", True, min)""",  # min match, not smarts, limit=1
        """HASPROP("tpsa" > 120) OR HASSUBSTRUCTURE("[OH]", True, "min")""",  # doesn't matter if min is quoted
        """HASPROP("tpsa" > 120) OR HASSUBSTRUCTURE("[OH]", True, min, 1)""",  # min match, not smarts, limit=1
        """HASPROP("tpsa" > 120) OR HASSUBSTRUCTURE("[OH]", True, 1)""",  # this is correct but bad practice
    ]
    for query in bad_queries:
        with pytest.raises((exceptions.UnexpectedToken, KeyError, exceptions.UnexpectedCharacters)):
            print(query)
            grammar.parse(query)

    valid_results = [grammar.parse(query) is not None for query in valid_queries]
    assert valid_results == [True] * len(valid_queries)


def test_query_equivalence():
    grammar = get_test_grammar()

    # the following should all be equal in the resulting query
    equivalent_queries = [
        """HASPROP("tpsa" = 120) OR ! HASSUBSTRUCTURE("CO")""",  # A OR ! B
        """HASPROP("tpsa", =, 120) OR ! HASSUBSTRUCTURE("CO")""",  # A OR ! B, comma does not matter
        """HASPROP("tpsa" =, 120) OR ! HASSUBSTRUCTURE("CO")""",  # A OR ! B, comma does not matter
        'HASPROP("tpsa" = 120) OR ! HASSUBSTRUCTURE("CO")',  # A OR ! B
        """(HASPROP("tpsa" = 120)) OR ! (HASSUBSTRUCTURE("CO"))""",  # (A) OR  ! (B)
        """(HASPROP("tpsa" = 120) OR ! HASSUBSTRUCTURE("CO"))""",  # (A OR ! B)
        """HASPROP("tpsa" = 120) OR NOT HASSUBSTRUCTURE("CO")""",  # A OR NOT B
        """HASPROP("tpsa" = 120) OR ~HASSUBSTRUCTURE("CO")""",  # A OR ~B
        """HASPROP("tpsa" == 120) OR ! HASSUBSTRUCTURE("CO")""",  # A OR ! B
        """HASPROP("tpsa" =   120)   OR ! HASSUBSTRUCTURE("CO")""",  # change in space
    ]
    different_query = """( ! HASPROP("tpsa" > 120) OR HASSUBSTRUCTURE("CO"))"""
    same_results = [grammar.parse(equivalent_queries[i]) for i in range(len(equivalent_queries))]

    expected_equivalent_results = [same_results[i] == same_results[0] for i in range(len(same_results))]
    expected_different_results = grammar.parse(different_query)

    assert expected_equivalent_results == [True] * len(same_results)
    assert expected_different_results != same_results[0]


def test_prop():
    mol = dm.to_mol("CC(=O)OC1=CC=CC=C1C(=O)O")
    tpsa = dm.descriptors.tpsa(mol)
    clogp = dm.descriptors.clogp(mol)

    assert QueryOperator.getprop(mol, "tpsa") == tpsa
    assert QueryOperator.getprop(mol, "clogp") == clogp

    with pytest.raises(Exception):
        QueryOperator.getprop(mol, "fake")

    with pytest.raises(Exception):
        QueryOperator.hasprop(mol, "tpsa", ">", 120)  # type: ignore

    assert QueryOperator.hasprop(mol, "tpsa", operator.le, 120) is True


def test_alert():
    mol = dm.to_mol("Oc1cscc1-c1ccc(O)cc1")  # contains tiophene hydroxy should match pains

    assert QueryOperator.hasalert(mol, "pains") is True
    assert QueryOperator.hasalert(mol, "brenk") is False

    with pytest.raises(Exception):
        QueryOperator.hasalert(mol, "fake")


def test_substructure():
    mol = "Oc1cscc1-c1ccc(OC)cc1"  # contains tiophene hydroxy should match pains
    substruct1 = "s1ccc(c1)-[#8;H1]"
    substruct2 = "CO"  # should match twice for SMILES and once for SMARTS

    assert QueryOperator.hassubstructure(mol, substruct1, is_smarts=True) is True
    assert QueryOperator.hassubstructure(mol, substruct2, is_smarts=True) is True
    assert QueryOperator.hassubstructure(mol, substruct2, is_smarts=False) is True
    assert QueryOperator.hassubstructure(mol, substruct2, True, "min", limit=1) is True
    assert QueryOperator.hassubstructure(mol, substruct2, True, "min", limit=2) is False
    assert QueryOperator.hassubstructure(mol, substruct2, True, None, limit=2) is False
    assert QueryOperator.hassubstructure(mol, substruct2, True, "max", limit=2) is True
    assert QueryOperator.hassubstructure(mol, substruct2, False, "min", limit=1) is True
    assert QueryOperator.hassubstructure(mol, substruct2, False, "min", limit=2) is True
    assert QueryOperator.hassubstructure(mol, substruct2, False, None, limit=2) is True
    assert QueryOperator.hassubstructure(mol, substruct2, False, "max", limit=2) is False


def test_rule():
    mol = "Oc1cscc1-c1ccc(OC)cc1"
    assert QueryOperator.matchrule(mol, "rule_of_five") is True
    assert QueryOperator.matchrule(mol, "rule_of_three") is False

    with pytest.raises(Exception):
        QueryOperator.matchrule(mol, "fake")


def test_group():
    mol = "Oc1cscc1-c1ccc(OC)cc1"
    assert QueryOperator.hasgroup(mol, "Ethers") is True
    assert QueryOperator.hasgroup(mol, "Ketones") is False

    # should normally match an alcohol, but how the group as defined, the specificity is high
    # so another more specific group should match and not alcohol for the tiophene hydroxy
    assert QueryOperator.hasgroup(mol, "Alcohols") is False
    assert QueryOperator.hasgroup(mol, "Hydroxy compounds: alcohols or phenols") is True

    with pytest.raises(Exception):
        QueryOperator.hasgroup(mol, "fake")


def test_similarity():
    mol1 = "Oc1cscc1-c1ccc(OC)cc1"
    mol2 = dm.to_mol("Oc1cscc1-c1ccc(O)cc1")
    dist = dm.cdist([mol1], [mol2])[0][0]

    assert QueryOperator.like(mol1, mol2, operator.le, 0.8) == True  # noqa: E712
    assert QueryOperator.similarity(mol1, mol2) == 1 - dist


def test_superstructure():
    query = "Oc1cscc1-c1ccc(OC)cc1"
    mol = "Oc1ccsc1"

    assert QueryOperator.hassuperstructure(mol, query) is True
    assert QueryOperator.hassuperstructure(query, mol) is False


def test_query_eval():
    mol = "Oc1cscc1-c1ccc(OC)cc1"
    grammar = Lark(get_grammar(as_string=True), parser="lalr", transformer=QueryParser())
    queries = [
        """HASPROP("tpsa" <= 120) AND ! HASALERT("pains") AND HASGROUP("Alcohols")""",  # False, does not match alchol
        """HASPROP("tpsa" <= 120) AND NOT HASALERT("pains") AND HASGROUP("Ethers")""",  # False has pains alerts
        """(HASPROP("tpsa" <= 120) AND NOT HASALERT("brenk")) AND HASGROUP("Ethers")""",  # True
        """HASPROP("tpsa" <= 120) AND MATCHRULE("rule_of_three")""",  # False
        """HASPROP("tpsa" <= 120) AND HASSUBSTRUCTURE("[OH]", True)""",  # True
        """HASPROP("tpsa" <= 120) AND HASSUBSTRUCTURE("[OH]", True, 3)""",  # False
    ]
    parsed_queries = [grammar.parse(q) for q in queries]
    expected_results = [False, False, True, False, True, False]

    observed_results = []
    observed_results_str = []
    for q in parsed_queries:
        evaluator = EvaluableQuery(q)
        observed_results.append(evaluator(mol))
        observed_results_str.append(evaluator(mol, exec=False))

    assert expected_results == observed_results


def test_query_filter():
    queries = [
        # complex query
        """(HASPROP("tpsa" < 100) AND HASPROP("clogp" < 3) AND ! HASALERT("pains")) OR (HASPROP("n_heavy_atoms" >= 10) AND (HASGROUP("Alcohols") OR HASSUBSTRUCTURE("[CX3](=[OX1])O", True, 1)))""",
        # is a rewriting of the above
        """(HASPROP("tpsa" < 100) AND HASPROP("clogp" < 3) AND ! HASALERT("pains")) OR (HASPROP("n_heavy_atoms" >= 10) AND (HASGROUP("Alcohols") OR HASSUBSTRUCTURE("[CX3](=[OX1])O", True, min, 1)))""",
        # is a rewriting of the above with spacing and a differentm yet equivalent bool expression
        """
        (
            HASPROP("tpsa" < 100)
            AND
            HASPROP("clogp" < 3)
            AND
            ! HASALERT("pains")
        )
        OR
        (
            HASPROP("n_heavy_atoms" >= 10)
            AND
            HASGROUP("Alcohols")
        )
        OR
        (
            HASPROP("n_heavy_atoms" >= 10)
            AND
            HASSUBSTRUCTURE("[CX3](=[OX1])O", True)
        )
        """,
        # always true
        """(HASPROP("tpsa" < 100) AND HASPROP("clogp" < 3)) OR TRUE""",
        # always false
        """(HASPROP("tpsa" < 100) AND HASPROP("clogp" < 3)) AND False""",
    ]

    query_filters = [QueryFilter(q) for q in queries]
    df = dm.cdk2()
    df["tpsa"] = df["mol"].apply(dm.descriptors.tpsa)
    df["clogp"] = df["mol"].apply(dm.descriptors.clogp)
    df["n_heavy_atoms"] = df["mol"].apply(dm.descriptors.n_heavy_atoms)
    df["has_carboxyl"] = df["mol"].apply(lambda x: QueryOperator.hassubstructure(x, "[CX3](=[OX1])O", True))
    df["has_pains"] = df["mol"].apply(lambda x: QueryOperator.hasalert(x, "pains"))
    df["has_alcohol"] = df["mol"].apply(lambda x: QueryOperator.hasgroup(x, "Alcohols"))

    tmp = df.query(
        "((tpsa < 100) & (clogp < 3) & ~has_pains) | (n_heavy_atoms >= 10 & (has_carboxyl | has_alcohol))"
    )
    df["expected_results"] = False
    df.loc[tmp.index, "expected_results"] = True

    computed_1 = query_filters[0](df["mol"].tolist(), progress=False)
    computed_2 = query_filters[1](df["mol"].tolist(), progress=False)
    computed_3 = query_filters[2](df["smiles"].tolist(), progress=False)
    all_true = query_filters[3](df["mol"].tolist(), progress=False)
    all_false = query_filters[4](df["mol"].tolist(), progress=False)

    np.testing.assert_array_equal(np.asarray(computed_1), np.asarray(computed_2))
    np.testing.assert_array_equal(np.asarray(computed_1), np.asarray(computed_3))
    np.testing.assert_array_equal(df.expected_results, np.asarray(computed_1))
    np.testing.assert_array_equal(np.asarray(all_true), [True] * len(df))
    np.testing.assert_array_equal(np.asarray(all_false), [False] * len(df))
