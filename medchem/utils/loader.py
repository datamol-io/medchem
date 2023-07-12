from typing import Optional
from typing import cast
from typing import Union

import os
import io
import functools
import importlib.resources

import fsspec


@functools.lru_cache(maxsize=10)
def get_data_path(filename: str, module: str = "medchem.data"):
    """Return the filepath of an internal data file."""

    path = importlib.resources.files(module).joinpath(filename)
    return str(path)


def get_grammar(
    grammar: Optional[Union[os.PathLike, str]] = None,
    as_string: bool = False,
):
    """Return the default lark grammar file for the medchem query system

    Args:
        grammar: The path to the grammar file. If None, the default medchem grammar file is used.
        as_string: If True, return the grammar as a string. Defaults to False.
    """
    if grammar is None:
        _grammar = importlib.resources.files("medchem.data").joinpath("grammar.lark")
        _grammar = str(_grammar)
    else:
        _grammar = str(grammar)

    if as_string:
        with fsspec.open(_grammar, "r") as f:
            f = cast(io.TextIOBase, f)
            _grammar = f.read()

    return _grammar
