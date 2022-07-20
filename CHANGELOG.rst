==================
medchem Change Log
==================

.. current developments

v1.3.6
====================

**Added:**

* Alternative rule of generative design

**Changed:**

* Behaviour of default rule of generative design

**Fixed:**

* Typo in docs



v1.3.5
====================

**Added:**

* Rule of generative design + test
* Various new descriptors
* Update smarts utils with both attachment standardization and new set of utils



v1.3.4
====================

**Fixed:**

* Non-kekulized molecules used for bredt filters



v1.3.3
====================

**Changed:**

* Make `lead.catalog_filter` much faster using batching.



v1.3.2
====================

**Added:**

* Add complexity filters



v1.3.1
====================

**Added:**

* Query system for complex filtering

**Changed:**

* Use datamol property functions
* update env to require `datamol>=0.7.1`

**Removed:**

* Local property function

**Fixed:**

* Error in documentation path



v1.3.0
====================

**Added:**

* SMARTS utils for efficient construction of complex smarts query
* `medchem.rules` module for an extensive list of physchem rules for molecular properties
* New molecular alerts and rules for toxicity and target/assay specific filtering
* Start building a collection for smarts bank to be reused across projects

**Changed:**

* medchem.utils into a module
* Rename ChemblFilter to AlertFilter (again), due to name inconsistency. 
* Rename `lead.alert_filter` to `lead.catalog_filter`
* Rename `lead.chembl_filter` to `lead.alert_filter`
* Chemical Group now has a name property instead of uipac and can support custom data bank



v1.2.0
====================

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

**Removed:**

* Removed `medchem.novartis`

**Fixed:**

* Docs
* Bug in catalog merging that does not support `FilterCatalog` but only `FilterCatalogParams.FilterCatalogs`



v1.1.7
====================

**Added:**

* Support for bredt's rule
* Support for illegal/nasty molecular graphs
* Parallelization for demerits score
* Move to threads as main scheduler because of rdkit



v1.1.6
====================

**Fixed:**

* Explicit conversion to int type for df index



v1.1.5
====================



v1.1.4
====================

**Added:**

* add `.pains_a`, `.pains_b` and `.pains_c` to `NamedCatalogs`.

**Fixed:**

* Remove an annoying `RuntimeWarnings` when important `medchem.catalogs`



v1.1.3
====================

**Fixed:**

* Fix bugs in `medchem.catalogs` module preventing loading some catalogs (the alert ones).



v1.1.2
====================



v1.1.2
====================



v1.1.2
====================



v1.1.2
====================

**Changed:**

* Catalog class more general now



v1.1.1
====================



v1.1.1
====================



v1.1.0
====================


