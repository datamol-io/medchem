from typing import Callable
from typing import Union

import types
import functools


def get_function_name(fn: Union[types.FunctionType, Callable, object]) -> str:
    """Get the name of a function"""

    if isinstance(fn, functools.partial):
        fn = fn.func
    if isinstance(fn, types.FunctionType):
        return fn.__name__
    elif isinstance(fn, object):
        return type(fn).__name__

    raise ValueError(f"Cannot find name of input function: {fn}")
