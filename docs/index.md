# Medchem

_Medchem - Molecular filtering for drug discovery._

Medchem is a Python library that proposes multiple molecular medchem filters to a wide range of use cases relevant in a drug discovery context.

## Installation

```bash
micromamba install -c conda-forge medchem
```

## Getting Started

The best way to get started is by going through the tutorials.

## Usage Notice

While medchem filters, alerts and rules are powerfull ways to triage a list of drug-like compounds, it's **important** to always keep in mind to **never apply those filters blindly**. It's impossible for such filters to take into account the full diversity of the drug-like chemical space and by applying them blindly you might discard interesting compounds or let pass toxic or unwanted compounds.

> _Don't blindly apply Medchem filters; you may miss gems or allow toxins._

## Acknowledgement

Medchem proposes a large list of mechem filters, alerts and rules. All of those has been built over the years by the scientific community and we would like to thanks everyone who has contributed to those collections and filtering methods.

### Eli Lilly Medchem Rules

Originally proposed in ["Rules for Identifying Potentially Reactive or Promiscuous Compounds" in 2012](https://doi.org/10.1021/jm301008n) by Robert F. Bruns and Ian A. Watson. Medchem is re-using the implementation from <https://github.com/IanAWatson/Lilly-Medchem-Rules>.

The Medchem implementation can be used from `medchem.structural.lilly_demerits`.

### NIBR Filters

Rules used by Novartis to build their new screening deck. First proposed in ["Evolution of Novartis' small molecule screening deck design" by Schuffenhauer, A. et al. J. Med. Chem. (2020)](https://dx.doi.org/10.1021/acs.jmedchem.0c01332).

The Medchem implementation ca be found at `medchem.structural.NIBRFilters()`.

### Structural Alerts Filters

These are alerts rules from the ChEMBL database curation scheme and public litterature on promiscuous compounds on commons assays. We thanks [Patrick Walters](https://twitter.com/wpwalters) for putting together this list of alerts.

The original implementation by Patrick Walters can be found at <https://github.com/PatWalters/rd_filters> with [its companion blog post](https://practicalcheminformatics.blogspot.com/2018/08/filtering-chemical-libraries.html).

The Medchem implementation can be found at `medchem.structural.CommonAlertsFilters()`.

### RDKit Catalogs

RDKit contains also contains a large list of filtering catalogs on which Medchem relies for the filtering logic.

Those catalogs can be found at `medchem.catalogs.NamedCatalogs`.
