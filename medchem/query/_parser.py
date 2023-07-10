from lark import Transformer, v_args


class QueryParser(Transformer):
    """
    Query parser for the custom query language for molecule. This parses the input language, build a parseable and evaluable representation.
    The trick for lazy evaluation is to define custom guard with '`fn(*)`' around expression that needs to be evaluated.

    Note that you **SHOULD NOT HAVE TO INTERACT WITH THIS CLASS DIRECTLY**.

    Example:
        ```python
        import medchem
        import lark
        QUERY_GRAMMAR = medchem.utils.loader.get_grammar(as_string=True)
        QUERY_PARSER = Lark(QUERY_GRAMMAR, parser="lalr", transformer=QueryParser())
        # see how the string needs to be "quoted". This builds on the json quote requirements to avoid dealing with unwanted outcomes
        example = \"""(HASPROP("tpsa" > 120 ) | HASSUBSTRUCTURE("c1ccccc1")) AND NOT HASALERT("pains") OR HASSUBSTRUCTURE("[OH]", max)\"""
        t = QUERY_PARSER.parse(example)
        print(t)
        ((((`fn(getprop, prop='tpsa')` > 120.0) or `fn(hassubstructure, query='c1ccccc1', operator='None', limit=None, is_smarts=None)`) and not `fn(hasalert, alert='pains')`) or `fn(hassubstructure, query='[OH]', operator='max', limit=None, is_smarts=None)`)
        ```
    """

    @v_args(inline=True)
    def not_bool_factor(self, *args):
        """Define representation of a negation"""
        return " ".join([str(x) for x in args])

    @v_args(inline=True)
    def bool_expr(self, bool_term, *others):
        """Define how boolean expressions should be parsed"""
        or_clauses = " ".join(f"{other}" for other in others)
        return f"({bool_term} {or_clauses})"

    @v_args(inline=True)
    def bool_term(self, bool_factor, *others):
        """Define how boolean terms should be parsed"""
        and_clauses = " ".join(f"{other}" for other in others)
        return f"({bool_factor} {and_clauses})"

    @v_args(inline=True)
    def hassubstructure(self, value, is_smarts, operator, limit):
        """Format the substructure node in the query

        !!! note
            The parser does not enforce any validity on the argument and
            the underlying function is supposed to handle it.

        """
        return f"`fn(hassubstructure, query='{value}', is_smarts={is_smarts}, operator={operator}, limit={limit})`"

    @v_args(inline=True)
    def hassuperstructure(self, value):
        """Format the superstructure node in the query

        !!! note
            The parser does not enforce any validity on the argument and
            the underlying function is supposed to handle it.

        """
        return f"`fn(hassuperstructure, query='{value}')`"

    @v_args(inline=True)
    def hasalert(self, value):
        """Format the hasalert node in the query

        !!! note
            The parser does not enforce any validity on the argument and
            the underlying function is supposed to handle it.

        """
        return f"`fn(hasalert, alert='{value}')`"

    @v_args(inline=True)
    def hasgroup(self, value):
        """Format the hasgroup node in the query

        !!! note
            The parser does not enforce any validity on the argument and
            the underlying function is supposed to handle it.

        """
        return f"`fn(hasgroup, group='{value}')`"

    @v_args(inline=True)
    def matchrule(self, value):
        """Format the matchrule node in the query

        !!! note
            The parser does not enforce any validity on the argument and
            the underlying function is supposed to handle it.

        """
        return f"`fn(matchrule, rule='{value}')`"

    @v_args(inline=True)
    def hasprop(self, value, comparator, limit):
        """Format the hasprop node in the query

        !!! note
            The parser does not enforce any validity on the argument and
            the underlying function is supposed to handle it.

        """
        return f"(`fn(getprop, prop='{value}')` {comparator} {limit})"

    @v_args(inline=True)
    def like(self, value, comparator, limit):
        """Format the like node in the query

        !!! note
            The parser does not enforce any validity on the argument and
            the underlying function is supposed to handle it.

        """
        return f"(`fn(similarity, query='{value}')` {comparator} {limit})"

    # non function call
    @v_args(inline=True)
    def operator(self, value):
        if value is None:
            return "None"
        return f"'{value}'"

    @v_args(inline=True)
    def comparator(self, value):
        if value == "=":
            # equality fix
            value = "=="
        return f"{value}"

    def is_smarts(self, value):
        return bool(value)

    def TRUE(self, _):
        return True

    def FALSE(self, _):
        return False

    def MIN(self, _):
        return "min"

    def MAX(self, _):
        return "max"

    def NOT_OP(self, _):
        return "not"

    def AND_OP(self, _):
        return "and"

    def OR_OP(self, _):
        return "or"

    def ESCAPED_STRING(self, value):
        return value[1:-1]

    def INT(self, i):
        return int(i)

    def SIGNED_NUMBER(self, f):
        return float(f)
