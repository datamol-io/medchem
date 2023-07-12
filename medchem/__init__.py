from ._version import __version__

from typing import TYPE_CHECKING

import os

import importlib


# The below lazy import logic is coming from openff-toolkit:
# https://github.com/openforcefield/openff-toolkit/blob/b52879569a0344878c40248ceb3bd0f90348076a/openff/toolkit/__init__.py#L44

# Dictionary of objects to lazily import; maps the object's name to its module path
_lazy_imports_obj = {
    # version
    "__version__": "medchem._version",
}

# Dictionary of modules to lazily import; maps the modules's name to its path
_lazy_imports_mod = {
    "utils": "medchem.utils",
    "groups": "medchem.groups",
    "catalogs": "medchem.catalogs",
    "constraints": "medchem.constraints",
    "complexity": "medchem.complexity",
    "rules": "medchem.rules",
    "structural": "medchem.structural",
    "functional": "medchem.functional",
    "query": "medchem.query",
}


def __getattr__(name):
    """Lazily import objects from _lazy_imports_obj or _lazy_imports_mod

    Note that this method is only called by Python if the name cannot be found
    in the current module."""
    obj_mod = _lazy_imports_obj.get(name)
    if obj_mod is not None:
        mod = importlib.import_module(obj_mod)
        return mod.__dict__[name]

    lazy_mod = _lazy_imports_mod.get(name)
    if lazy_mod is not None:
        return importlib.import_module(lazy_mod)

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    """Add _lazy_imports_obj and _lazy_imports_mod to dir(<module>)"""
    keys = (*globals().keys(), *_lazy_imports_obj.keys(), *_lazy_imports_mod.keys())
    return sorted(keys)


if TYPE_CHECKING or os.environ.get("MEDCHEM_DISABLE_LAZY_LOADING", "0") == "1":
    # These types are imported lazily at runtime, but we need to tell type
    # checkers what they are.

    from . import utils
    from . import groups
    from . import catalogs
    from . import constraints
    from . import complexity
    from . import rules
    from . import structural
    from . import functional
    from . import query
