from typing import Optional

import os
import importlib.resources as importlib_resources

import fsspec


def get_data_path_path(filename: str, module: str = "medchem.data"):
    """Return the filepath of a data file."""

    path = importlib_resources.files(module).joinpath(filename)
    path = str(path)
    return path


def get_grammar(grammar: Optional[os.PathLike] = None, as_string: bool = False):
    """Return the default lark grammar file for queries

    Args:
        grammar: The path to the grammar file. If None, the default medchem grammar file is used.
        as_string: If True, return the grammar as a string. Defaults to False.
    """
    if grammar is None:
        grammar = importlib_resources.files("medchem.data").joinpath("grammar.lark")
        grammar = str(grammar)

    if as_string:
        with fsspec.open(grammar, "r") as f:
            grammar = f.read()

    return grammar
