# Medchem

Package for applying common medchem filters to a dataset of molecules.

## Summary

This package contains various implementation of medchem rules collected from various sources that may be applied as filters on generated or screened molecules. It centralizes all common filters used at Valence Discovery.

Although the list is as exhaustive as possible, filtering rules mainly depends on the drug discovery programs.

It should be noted that **systematically applying all filters is to be avoided**. For example, "PAINS C" filters are usually not very relevant, another example is the filtering are very strict and could flag important substructure for a project (example some ZBGs).

## Available Filters

The following filters are available:

### **Eli Lilly Medchem Rules**

These are python binding of the implementation of Eli Lilly Medchem Rules published under "Rules for Identifying Potentially Reactive or Promiscuous Compounds" by Robert F. Bruns and Ian W. Watson, J. Med. Chem. 2012, 55, 9763--9772 as ACS Author choice, i.e. open access at [doi 10.1021/jm301008n](https://doi.org/10.1021/jm301008n).

These rules are used in `medchem.filter.demerit_filter` function.

### NIBR filters

Rules used by Novartis to build their new screening deck. The rules are published under "Evolution of Novartis' small molecule screening deck design" by Schuffenhauer, A. et al. J. Med. Chem. (2020), https://dx.doi.org/10.1021/acs.jmedchem.0c01332.

These rules are used in lead filtering as `medchem.filter.lead.screening_filter`

### Common tox and assay interference rules

These are pure rdkit filtering rules based on PAINS, BRENK, NIH and ZINC filters. There are used in lead filtering as `medchem.filter.lead.common_filter`

### ChEMBL filters

These are alerts rules from the ChEMBL database that have been collected from various Pharma groups and commons assays. The rule set are:

| Rule Set                                                | Number of Alerts |
| ------------------------------------------------------- | ---------------: |
| BMS                                                     |              180 |
| Dundee                                                  |              105 |
| Glaxo                                                   |               55 |
| Inpharmatica                                            |               91 |
| LINT                                                    |               57 |
| MLSMR                                                   |              116 |
| [PAINS](https://pubs.acs.org/doi/abs/10.1021/jm901137j) |              479 |
| SureChEMBL                                              |              166 |

There are used in lead filtering as `medchem.filter.lead.alert_filter`

#### Generic filters

These are generic filters based on specific molecular property such as number of atoms, size of macrocycles, etc. They are available at `medchem.filter.generic`

## Installation

### conda

```bash
conda install -c invivoai medchem
```

### Source

This package requires : `gcc` and `g++` for compilation. Use your OS package manager or conda:

```bash
conda install -c conda-forge c-compiler cxx-compiler makez zlib
```

Clone the repo and install it locally

```bash
git clone https://github.com/valence-platform/medchem.git
cd medchem
pip install -e .
```

This will build (compile the C source) and install the package. If you are having trouble with pip, use setup.py:

```bash
python setup.py install # "python setup.py build" should not be necessary

```

### pip

You can install directly from git using pip too:

```bash
pip install git+https://github.com/valence-platform/medchem.git
```

### Troubleshooting

In the rare case where you `$CPATH` is not well configured, you will get some `error: stray \number in program`.
Please set your CPATH to empy with `export CPATH=` before compiling/installing.

### Module

You can import the package and run the filters of interest. For more information see the [Getting Started](docs/tutorials/getting-started.ipynb) tutorial.

```python
from medchem.filter.lead import common_filter
from medchem.filter.lead import demerit_filter

test_config = {
    'output': 'test',
    'min_atoms': 7,
    'soft_max_atoms': 30,
    'hard_max_atoms': 50,
    'smarts': [],
    'nodemerit': False,
    'dthresh': 160,
    'odm': [],
    'okiso': False,
    'noapdm': False
} # optional config dict

smiles_list = [
    'Cc1cnc(CNc(cccc2-c3cn(CC(C4)CC4O)c4ncnc(N)c34)c2F)s1'
    'Cc1cnc(CNc(cccc2-c3cn(C[C@H](C4)C[C@@H]4O)c4ncnc(N)c34)c2F)s1',
    'Nc1c(c(-c2cccc(NCc3nccs3)c2F)cn2Cc3cn(C4CC5(COCC5)OCC4)nn3)c2ncn1',
    'C[C@@H](c1cn(CC2CCCCC2)nn1)N(Cc1nccs1)c(cccc1-c2c[nH]c3ncnc(N)c23)c1F',
    'Nc1c(c(-c2cccc(NCc3nccs3)c2F)cn2Cc3cn(CCN4CCOCC4)nn3)c2ncn1',
    'Cc1cnc(C(Nc(cccc2-c3cn(C[C@H](C4)C[C@@H]4O)c4ncnc(N)c34)c2F)=O)s1',
    'Nc1c(c(-c(cccc2NCc3nccs3)c2F)cn2Cc3cn(-c(cc4)cc5c4scc5)nn3)c2ncn1',
    'Nc1c(c(-c2cccc(NCc3nccs3)c2F)cn2Cc3cn(CCC4CC4)nn3)c2ncn1',
    'CC(C)Cn1ncnc1Cn1nnc(CNc2c(c(-c3cccc(NCc4nccs4)c3F)c[nH]3)c3ncn2)c1',
    'Nc1c(c(-c2cccc(NCc3nc(CO)cs3)c2F)cn2CC3CCC3)c2ncn1',
    'Nc1c(c(-c2cccc(NCc3nccs3)c2F)cn2Cc3cn(C4CCCC4)nn3)c2ncn1'
]
res_demerits = demerit_filter(smiles_list, **test_config)
res_common = common_filter(smiles_list)
```

### Command line

You can also use the provided binary: `chemfilter --help`. This will only apply the demerits (Eli Lilly) filters.

```bash
Usage: chemfilter [OPTIONS]

Options:
  -n, --min-atoms INTEGER        Use lower atom count cutoff
  -Cs, --soft-max-atoms INTEGER  Apply soft upper atom count cuttof
  -Ch, --hard-max-atoms INTEGER  Apply hard upper atom count cuttof
  --min-num-rings INTEGER        Minimum number of rings accepted
  --max-num-rings INTEGER        Maximum number of rings accepted
  --max-size-rings INTEGER       Maximum ring size (number of atom per ring)
  --max-size-chain INTEGER       Threshold for long carbon chain (7)
  --smarts TEXT                  Optional smarts to reject. File or smarts
                                 string

  --nodemerit                    Use hard rejections only, do not apply any
                                 demerits

  --dthresh INTEGER              Demerit threshold. For relaxed rules, use 160
                                 demerit cutoff

  -o, --output TEXT              output file where to write result  [required]
  --odm TEXT                     Optional demerits to omit to apply. File or
                                 demerits

  --okiso                        Allow isotopic atoms to pass through
  --noapdm                       Do not append demerit reasons
  --smcol TEXT                   Name of the smiles columns
  --alert-filter                 Whether to run lead filtering (alerts) on the molecules
  --alerts TEXT                  Alerts to use for lead filtering (multiple allowed)
  --allow-non-interesting        Allow molecules with non interesting atoms
                                 only to pass

  -i, --input-file TEXT          Input csv or smi files. Header expected and
                                 first column should always be the smiles if
                                 smiles column name is not provided.
                                 [required]

  --help                         Show this message and exit.
```

## Maintainers

@maclandrol - Emmanuel Noutahi
