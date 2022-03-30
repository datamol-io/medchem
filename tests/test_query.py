import unittest as ut
import datamol as dm
import operator
from lark import Lark, Transformer, v_args, exceptions
from medchem.query.parser import QueryParser
from medchem.query.eval import QueryOperator
from medchem.utils.loader import get_grammar


class Test_QueryParser(ut.TestCase):
    grammar = Lark(
        get_grammar(as_string=True), parser="lalr", transformer=QueryParser()
    )

    def test_language(self):
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
            with self.assertRaises(
                (exceptions.UnexpectedToken, KeyError, exceptions.UnexpectedCharacters)
            ):
                print(query)
                self.grammar.parse(query)

        valid_results = [
            self.grammar.parse(query) is not None for query in valid_queries
        ]
        self.assertListEqual(valid_results, [True] * len(valid_queries))

    def test_query_equivalence(self):
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
        same_results = [
            self.grammar.parse(equivalent_queries[i])
            for i in range(len(equivalent_queries))
        ]

        expected_equivalent_results = [
            same_results[i] == same_results[0] for i in range(len(same_results))
        ]
        expected_different_results = self.grammar.parse(different_query)
        self.assertListEqual(expected_equivalent_results, [True] * len(same_results))
        self.assertNotEqual(expected_different_results, same_results[0])


class Test_QueryOperator(ut.TestCase):
    def test_prop(self):
        mol = dm.to_mol("CC(=O)OC1=CC=CC=C1C(=O)O")
        tpsa = dm.descriptors.tpsa(mol)
        clogp = dm.descriptors.clogp(mol)

        self.assertEqual(QueryOperator.getprop(mol, "tpsa"), tpsa)
        self.assertEqual(QueryOperator.getprop(mol, "clogp"), clogp)

        with self.assertRaises(Exception):
            QueryOperator.getprop(mol, "fake")

        with self.assertRaises(Exception):
            QueryOperator.hasprop(mol, "tpsa", ">", 120)

        self.assertTrue(QueryOperator.hasprop(mol, "tpsa", operator.le, 120))

    def test_alert(self):
        mol = dm.to_mol(
            "Oc1cscc1-c1ccc(O)cc1"
        )  # contains tiophene hydroxy should match pains
        self.assertTrue(QueryOperator.hasalert(mol, "pains"))
        self.assertFalse(QueryOperator.hasalert(mol, "brenk"))
        with self.assertRaises(Exception):
            QueryOperator.hasalert(mol, "fake")

    def test_substructure(self):
        mol = "Oc1cscc1-c1ccc(OC)cc1"  # contains tiophene hydroxy should match pains
        substruct1 = "s1ccc(c1)-[#8;H1]"
        substruct2 = "CO"  # should match twice for SMILES and once for SMARTS
        self.assertTrue(QueryOperator.hassubstructure(mol, substruct1, is_smarts=True))
        self.assertTrue(QueryOperator.hassubstructure(mol, substruct2, is_smarts=True))
        self.assertTrue(QueryOperator.hassubstructure(mol, substruct2, is_smarts=False))
        self.assertTrue(
            QueryOperator.hassubstructure(mol, substruct2, True, "min", limit=1)
        )
        self.assertFalse(
            QueryOperator.hassubstructure(mol, substruct2, True, "min", limit=2)
        )
        self.assertFalse(
            QueryOperator.hassubstructure(mol, substruct2, True, None, limit=2)
        )
        self.assertTrue(
            QueryOperator.hassubstructure(mol, substruct2, True, "max", limit=2)
        )

        self.assertTrue(
            QueryOperator.hassubstructure(mol, substruct2, False, "min", limit=1)
        )
        self.assertTrue(
            QueryOperator.hassubstructure(mol, substruct2, False, "min", limit=2)
        )
        self.assertTrue(
            QueryOperator.hassubstructure(mol, substruct2, False, None, limit=2)
        )
        self.assertFalse(
            QueryOperator.hassubstructure(mol, substruct2, False, "max", limit=2)
        )

    def test_rule(self):
        mol = "Oc1cscc1-c1ccc(OC)cc1"
        self.assertTrue(QueryOperator.matchrule(mol, "rule_of_five"))
        self.assertFalse(QueryOperator.matchrule(mol, "rule_of_three"))
        with self.assertRaises(Exception):
            QueryOperator.matchrule(mol, "fake")

    def test_group(self):
        mol = "Oc1cscc1-c1ccc(OC)cc1"
        self.assertTrue(QueryOperator.hasgroup(mol, "Ethers"))
        self.assertFalse(QueryOperator.hasgroup(mol, "Ketones"))
        # should normally match an alcohol, but how the group as defined, the specificity is high
        # so another more specific group should match and not alcohol for the tiophene hydroxy
        self.assertFalse(QueryOperator.hasgroup(mol, "Alcohols"))
        self.assertTrue(
            QueryOperator.hasgroup(mol, "Hydroxy compounds: alcohols or phenols")
        )
        with self.assertRaises(Exception):
            QueryOperator.hasgroup(mol, "fake")

    def test_similarity(self):
        mol1 = "Oc1cscc1-c1ccc(OC)cc1"
        mol2 = dm.to_mol("Oc1cscc1-c1ccc(O)cc1")
        dist = float(dm.cdist([mol1], [mol2]))
        self.assertTrue(QueryOperator.like(mol1, mol2, operator.le, 0.8))
        self.assertEqual(QueryOperator.similarity(mol1, mol2), 1 - dist)

    def test_superstructure(self):
        query = "Oc1cscc1-c1ccc(OC)cc1"
        mol = "Oc1ccsc1"
        self.assertTrue(QueryOperator.hassuperstructure(mol, query))
        self.assertFalse(QueryOperator.hassuperstructure(query, mol))


if __name__ == "__main__":
    ut.main()
