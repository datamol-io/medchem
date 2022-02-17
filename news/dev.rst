**Added:**

* Added  `medchem.catalog.list_named_catalogs`
* Added  `medchem.catalog.list_chemical_groups`
* Added  `medchem.group.ChemicalGroup`
* Added  `medchem.filter.lead.chemical_group_filters`

**Changed:**

* Moved `medchem.novartis.NovartisFilters` to `medchem.alerts`
* Renamed `medchem.alerts.AlertFilters` to `medchem.alerts.ChEMBLFilters`
* Renamed `medchem.filter.lead.common_filter` to `medchem.filter.lead.chembl_filter`
* Changed  `medchem.filter.lead.alert_filter` to take either a list of catalog or named catalogs
* Renamed  `medchem.catalog.merge` to `medchem.catalog.merge_catalogs`

**Deprecated:**

* <news item>

**Removed:**

* Removed `medchem.novartis`

**Fixed:**

* Bug in catalog merging that does not support `FilterCatalog` but only `FilterCatalogParams.FilterCatalogs`

**Security:**

* <news item>
