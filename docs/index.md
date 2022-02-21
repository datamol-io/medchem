Medchem
===================

Package for applying common medchem filters to a dataset of molecules.

## Summary

This package contains various implementation of medchem rules collected from various sources that may be applied as filters on generated or screened molecules. It centralizes all common filters used at Valence Discovery.

Although the list is as exhaustive as possible, filtering rules mainly depends on the drug discovery programs. 

It should be noted that **systematically applying all filters is to be avoided**. For example, "PAINS C" filters are usually not very relevant, another example is the filtering are very strict and could flag important substructure for a project (example some ZBGs).


## Available Filters

The following filters are available:

#### **Eli Lilly Medchem Rules**

These are python binding of the implementation of Eli Lilly Medchem Rules published under "Rules for Identifying Potentially Reactive or Promiscuous Compounds" by Robert F. Bruns and Ian W. Watson, J. Med. Chem. 2012, 55, 9763--9772 as ACS Author choice, i.e. open access at [doi 10.1021/jm301008n](https://doi.org/10.1021/jm301008n).

These rules are used in `medchem.filter.lilly_demerit_filter` function and are the main offering of this package.

#### NIBR filters

Rules used by Novartis to build their new screening deck. The rules are published under "Evolution of Novartis' small molecule screening deck design" by Schuffenhauer, A. et al. J. Med. Chem. (2020), https://dx.doi.org/10.1021/acs.jmedchem.0c01332. 

These rules are used in lead filtering as `medchem.filter.lead.screening_filter`

#### Common tox and assay interference rules

These are filtering rules based on PAINS, BRENK, NIH and ZINC and any other catalog provided by medchem. There are used in lead filtering as `medchem.filter.lead.alert_filter` and you need to provide the list of catalog you want to use.


#### Bredt filters

These are filters based on the Bredt's rules for unstable chemistry.There are used in lead filtering as `medchem.filter.lead.bredt_filter`.


#### ChEMBL filters

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

There are used in lead filtering as `medchem.filter.lead.chembl_filter`

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
conda install -c conda-forge c-compiler cxx-compiler
# conda install gcc_linux-64 
# conda install gxx_linux-64
```

Clone the repo and install it locally
```bash
git clone https://github.com/valence-platform/medchem.git
cd medchem 
pip install . # Alternatively you can install a develop version
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


## Maintainers

@maclandrol - Emmanuel Noutahi
