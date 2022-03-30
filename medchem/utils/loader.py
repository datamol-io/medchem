from typing import Optional

import os
import fsspec
import medchem


def get_data(file=None):
    """Return the folder that contains the package specific data"""

    path = os.path.join(medchem.__path__[0], "data/")
    if file is not None:
        path = os.path.join(path, file)
    return path


def get_grammar(grammar: Optional[os.PathLike] = None, as_string: bool = False):
    """Return the default lark grammar file for queries

    Args:
        grammar: The path to the grammar file. If None, the default grammar
        as_string (bool, optional): If True, return the grammar as a string. Defaults to False.
    """
    if grammar is None:
        grammar = get_data("grammar.lark")
    if as_string:
        with fsspec.open(grammar, "r") as IN:
            grammar = IN.read()
    return grammar
