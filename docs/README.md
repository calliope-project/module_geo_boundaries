# Home

Welcome to the documentation of the `module_geo_boundaries` data module!
This module combines country regions from three datasources (GADM, NUTS, or Overture Maps), and Exclusive Economic Zones from the Marine regions database.

Please consult the [`clio` documentation](https://clio.readthedocs.io/) for more information.


## Overview

![rulegraph](./rulegraph.png)

The analysis of the module is structured as follows.

1. The configuration file is read to identify the countries and regional aggregation (`subtype` in the configuration) to process.
1. Individual country files are downloaded and harmonised to fit a standardised schema.
Contested regions are removed at this stage.
1. Land is clipped using maritime Exclusive Economic Zones (EEZ).
1. Each polygon is clipped using its neighbours to minimise overlapping polygons.

> [!TIP]
> The subtype naming matches that of the source database. For example, NUTS uses 0, 1, 2 and 3 (NUTS0, NUTS1, NUTS2, etc.).
> Please consult data sources below for more details.

## Configuration

See [the configuration README](./../config/README.md).

## Outputs

See the [interface file](./../INTERFACE.yaml).

## Data sources

We encourage users to cite both the datasets requested and our workflow.

- GADM 4.1. (2018). Global Administrative Areas (GADM).
    - License: GADM data is freely available for academic and non-commercial use. <https://gadm.org/license.html>.
- NUTS (various years). Nomenclature of territorial units for statistics (NUTS).
    - License: reuse is authorised provided the source is acknowledged. <https://ec.europa.eu/eurostat/statistics-explained/index.php?title=Copyright/licence_policy>
- Overture Maps Divisions database (most recent version). Overture Maps Foundation.
    - License: ODbL. See <https://docs.overturemaps.org/attribution/> and <https://opendatacommons.org/licenses/odbl/summary/> for details.
- Marine Regions World EEZ v12 (2023). Flanders Marine Institute (MarineRegions.org).
    - License: CC-By. See <https://www.marineregions.org/disclaimer.php>.
