# Medchem

Package for applying common medchem filters to a dataset of molecules.

## Summary

This package contains various implementation of medchem rules collected from various sources that may be applied as filters on generated or screened molecules. It centralizes all common filters used at Valence Discovery.

Although the list is as exhaustive as possible, filtering rules mainly depends on the drug discovery programs.

It should be noted that **systematically applying all filters is to be avoided**. For example, "PAINS C" filters are usually not very relevant, another example is the filtering are very strict and could flag important substructure for a project (example some ZBGs).

## Installation

```bash
micromamba install -c conda-forge medchem
```

## Available Filters

The following filters are available:

### **Eli Lilly Medchem Rules**

These are python binding of the implementation of Eli Lilly Medchem Rules published under "Rules for Identifying Potentially Reactive or Promiscuous Compounds" by Robert F. Bruns and Ian W. Watson, J. Med. Chem. 2012, 55, 9763--9772 as ACS Author choice, i.e. open access at [doi 10.1021/jm301008n](https://doi.org/10.1021/jm301008n).

These rules are used in `medchem.filter.lilly_demerit_filter` function and are the main offering of this package.

### NIBR filters

Rules used by Novartis to build their new screening deck. The rules are published under "Evolution of Novartis' small molecule screening deck design" by Schuffenhauer, A. et al. J. Med. Chem. (2020), <https://dx.doi.org/10.1021/acs.jmedchem.0c01332>.

These rules are used in lead filtering as `medchem.filter.lead.screening_filter`

### Bredt filters

These are filters based on the Bredt's rules for unstable chemistry.There are used in lead filtering as `medchem.filter.lead.bredt_filter`.

### Alerts filters

These are alerts rules from the ChEMBL database curation scheme and public litterature on promiscuous compounds on commons assays. The rule set are:

| name                          | # alerts | source      |
| :---------------------------- | -------: | :---------- |
| Glaxo                         |       55 | ChEMBL      |
| Dundee                        |      105 | ChEMBL      |
| BMS                           |      180 | ChEMBL      |
| PAINS                         |      481 | ChEMBL      |
| SureChEMBL                    |      166 | ChEMBL      |
| MLSMR                         |      116 | ChEMBL      |
| Inpharmatica                  |       91 | ChEMBL      |
| LINT                          |       57 | ChEMBL      |
| Alarm-NMR                     |       75 | Litterature |
| AlphaScreen-Hitters           |        6 | Litterature |
| GST-Hitters                   |       34 | Litterature |
| HIS-Hitters                   |       19 | Litterature |
| LuciferaseInhibitor           |        3 | Litterature |
| DNABinder                     |       78 | Litterature |
| Chelator                      |       55 | Litterature |
| Frequent-Hitter               |       15 | Litterature |
| Electrophilic                 |      119 | Litterature |
| Genotoxic-Carcinogenicity     |      117 | Litterature |
| LD50-Oral                     |       20 | Litterature |
| Non-Genotoxic-Carcinogenicity |       22 | Litterature |
| Reactive-Unstable-Toxic       |      335 | Litterature |
| Skin                          |      155 | Litterature |
| Toxicophore                   |      154 | Litterature |

There are used in lead filtering through `medchem.filter.lead.alert_filter`

### Generic filters

These are generic filters based on specific molecular property such as number of atoms, size of macrocycles, etc. They are available at `medchem.filter.generic`
