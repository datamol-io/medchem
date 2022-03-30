# `medchem.query`

This module helps build a filter based on a query language that can be parsed. 
By default, the default query parser will be used, which contains the following instructions that can be orchestrated using boolean operation (`or`, `and`, `not` and parenthesis)

## Example

```python
import datamol as dm
from medchem.query.eval import QueryFilter

query = """HASPROP("tpsa" < 120) AND HASSUBSTRUCTURE("[OH]", True)"""
chemical_filter = QueryFilter(query, parser="lalr")
mols = dm.data.cdk2().mol[:10]
chemical_filter(mols, n_jobs=-1) # [False, False, False, False, False, True, True, True, False, False]
```

## Syntax

Any string provided as `query` argument needs to be quoted (similar to json) to avoid ambiguity in parsing. 
* An example of valid query is `"""(HASPROP("tpsa" > 120 ) | HASSUBSTRUCTURE("c1ccccc1")) AND NOT HASALERT("pains") OR HASSUBSTRUCTURE("[OH]", max, 2)"""`.
* Examples of invalid queries are 
  * `"""HASPROP("tpsa" > 120) OR HASSUBSTRUCTURE("[OH]", True, >, 3)"""` : unexpected wrong operator `>`
  * `"""HASPROP(tpsa > 120)"""` : tpsa is not quoted
  * `"""HASPROP("tpsa") > 120"""` : this is not part of the language specification
  * `"""(HASPROP("tpsa" > 120) AND HASSUBSTRUCTURE("[OH]", True, max, 3 )"""`: mismatching parenthesis `(`

* `"""HASPROP("tpsa" > 120) OR HASSUBSTRUCTURE("CO")"""`, `"""(HASPROP("tpsa" > 120)) OR (HASSUBSTRUCTURE("CO"))"""` and `"""(HASPROP("tpsa" > 120) OR HASSUBSTRUCTURE("CO"))"""` are equivalent

  
### HASALERT
check whether a molecule has an `alert` from a catalog
```python 
# alert is one supported alert catalog by `medchem`. For example `pains`
HASALERT(alert:str) 
```

### HASGROUP
check whether a molecule has a specific functional group from a catalog

```python 
# group is one supported functional group provided by `medchem`
HASGROUP(group:str) 
```


### MATCHRULE
check whether a molecule match a predefined druglikeness `rule` from a catalog
```python 
# rule is one supported rule provided by `medchem`. For example `rule_of_five`
MATCHRULE(rule:str) 
```

### HASSUPERSTRUCTURE
check whether a molecule has `query` as superstructure
```python 
# query is a SMILES
HASSUPERSTRUCTURE(query:str) 
```

### HASSUBSTRUCTURE
Check whether a molecule has `query` as substructure. 
**Note that providing the comma separator `,` is _mandatory_ here as each variable is an argument.**

```python
# query is a SMILES or a SMARTS, operator is defined below, is_smarts is a boolean

HASSUBSTRUCTURE(query:str, is_smarts:Optional[bool], operator:Optional[str], limit:Optional[int])

# which correspond to setting this default values
HASSUBSTRUCTURE(query:str, is_smarts=False, operator="min", limit=1)
# same as
HASSUBSTRUCTURE(query:str, is_smarts=None, operator=None, limit=None)
```

Not providing optional arguments is allowed, but they need to be provided in the exact same order shown above. Thus:

* `HASSUBSTRUCTURE("CO")`
* `HASSUBSTRUCTURE("CO", False)`
* `HASSUBSTRUCTURE("CO", False, min)`
* `HASSUBSTRUCTURE("CO", False, min, 1)`
  
are all `valid` and `equivalent` (given their default values)

Furthermore, since the correct argument map can be inferred when no ambiguity arises, the following `are valid but discouraged`

* `HASSUBSTRUCTURE("CO", False, 1)`
* `HASSUBSTRUCTURE("CO", min, 1)`

Whereas, this is invalid:
* `HASSUBSTRUCTURE("CO", min, False, 1)`
  

### HASPROP
Check whether a molecule has `prop` as property within a defined limit.
**Any comma `,` provided between arguments will be ignored**

```python
# prop is a valid datamol.descriptors property, comparator is a required comparator operator and defined below
HASPROP(prop:str comparator:str limit:float)
```

### LIKE
Check whether a molecule is similar enough to another molecule.
**Any comma `,` provided between arguments will be ignored**

```python
# query is a SMILES
LIKE(query:str  comparator:str limit:float)
```

### Basic operators:

* comparator: one of `=` `==`, `!=`, `<`, `>`, `<=`,  `>=`
* misc: the following misc values are accepted and parsed `true`, `false`, `True`, `False`, `TRUE`, `FALSE`
* operator (can be quoted or unquoted):
  * MIN: `min`, `MIN`
  * MAX: `max`, `MAX`
* boolean operator: 
  * AND operator : `AND` or `&` or `&&` or `and`
  * OR operator : `OR`  or `|` or `||` or `or`
  * NOT operator : `NOT` or  `!` or  `~` or `not`



## API

::: medchem.query.parser

---

::: medchem.query.eval