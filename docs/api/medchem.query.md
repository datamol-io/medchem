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

Any string provided in the query needs to be quoted (similar to json) to avoid ambiguity in parsing. 
* An example of valid queries is `"""(HASPROP("tpsa" > 120 ) | HASSUBSTRUCTURE("c1ccccc1")) AND NOT HASALERT("pains") OR HASSUBSTRUCTURE("[OH]", max, 2)"""`.
* An example of invalid queries are 
  * `"""HASPROP("tpsa" > 120) OR HASSUBSTRUCTURE("[OH]", True, >, 3)"""` : unexpected wrong operator `>`
  * `"""(HASPROP("tpsa" > 120) AND HASSUBSTRUCTURE("[OH]", True, max, 3 )"""`: mismatching parenthesis `(`

* `"""HASPROP("tpsa" > 120) OR HASSUBSTRUCTURE("CO")"""`, `"""(HASPROP("tpsa" > 120)) OR (HASSUBSTRUCTURE("CO"))"""` and `"""(HASPROP("tpsa" > 120) OR HASSUBSTRUCTURE("CO"))"""` are equivalent

  
### HASALERT
check whether a molecule has an `alert` from a catalog
```python 
# alert is one supported alert catalog by `medchem`. For example `pains`
HASALERT(alert:str) 
```

### HASSUPERSTRUCTURE
check whether a molecule has `query` as superstructure
```python 
# query is a SMILES
HASSUPERSTRUCTURE(query:str) 
```

### HASSUBSTRUCTURE
Check whether a molecule has `query` as substructure

```python
# query is a SMILES or a SMARTS, operator is defined below, is_smarts is a boolean
HASSUBSTRUCTURE(query:str, is_smarts:Optional[bool], operator:Optional[str], limit:Optional[int])
```

### HASPROP
Check whether a molecule has `prop` as property within a defined limit.

```python
# prop is a valid datamol.descriptors property, comparator is a required comparator operator and defined below
HASPROP(prop:str, comparator:str, limit:float)
```

### LIKE
Check whether a molecule is similar enough to another molecule.

```python
# query is a SMILES
LIKE(query:str, comparator:str, limit:float)
```

### Basic operators:

* comparator: one of `=` `==`, `!=`, `<`, `>`, `<=`,  `>=`
* misc: the following misc values are accepted and parsed `true`, `false`, `True`, `False`, `TRUE`, `FALSE`
* operator: one of `min`, `max`
* boolean operator: 
  * AND operator : `AND` or `&` or `&&` or `and`
  * OR operator : `OR`  or `|` or `||` or `or`
  * NOT operator : `NOT` or  `!` or  `~` or `not`



::: medchem.query.parser

---

::: medchem.query.eval