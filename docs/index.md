# Medchem

_Medchem - Molecular filtering for drug discovery._

Medchem is a Python library that proposes molecular filters and prioritization rules for a wide range of use cases relevant in drug discovery.

## Installation

```bash
micromamba install -c conda-forge medchem

# or using pip
pip install medchem
```

## Getting Started

The best way to get started is by going through [**the tutorials**](./tutorials/Basic_Concepts.ipynb).

## Usage Notice

When using medchem filters, alerts, and rules to triage a list of drug-like compounds, it is **essential to avoid blindly applying these filters**. It's impossible for such filters to take into account the full diversity of the drug-like chemical space and by blindly applying them you might discard interesting compounds or allow toxic or unwanted compounds for your particular context.

> _Don't blindly apply Medchem filters; you may miss gems or allow toxins._

## Acknowledgement

Medchem incorporates a comprehensive collection of medchem filters, alerts, and rules that have been developed by the scientific community over the years. We extend our gratitude to all the contributors to these collections and filtering methods.

### Eli Lilly Medchem Rules

Originally proposed in ["Rules for Identifying Potentially Reactive or Promiscuous Compounds" in 2012](https://doi.org/10.1021/jm301008n) by Robert F. Bruns and Ian A. Watson. Medchem is re-using the implementation from <https://github.com/IanAWatson/Lilly-Medchem-Rules>.

The Medchem implementation is accessible through `medchem.structural.lilly_demerits`.

### NIBR Filters

Rules used by the Novartis Institutes for BioMedical Research to build their screening deck. The initial proposal can be found in ["Evolution of Novartis' small molecule screening deck design" by Schuffenhauer, A. et al. J. Med. Chem. (2020)](https://dx.doi.org/10.1021/acs.jmedchem.0c01332).

The Medchem implementation can be found at `medchem.structural.NIBRFilters()`.

### Common Structural Alert Filters

These alert rules are derived from the ChEMBL database curation scheme and public literature on promiscuous compounds in common assays. We extend our thanks [Patrick Walters](https://twitter.com/wpwalters) for putting together this list of alerts.

The original implementation by Patrick Walters can be found at <https://github.com/PatWalters/rd_filters> along with its accompanying blog post [Filtering Chemical Libraries](https://practicalcheminformatics.blogspot.com/2018/08/filtering-chemical-libraries.html).

The Medchem implementation can be found at `medchem.structural.CommonAlertsFilters()`.

### RDKit Catalogs

RDKit also contains a large list of catalogs on which Medchem relies for our filtering logic.

Those catalogs can be found at `medchem.catalogs.NamedCatalogs`.
