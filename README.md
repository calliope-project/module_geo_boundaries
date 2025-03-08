# clio - Geo-boundaries module

A data module to create arbitrary regional boundaries for energy systems modelling.

A modular `snakemake` workflow built for [`clio`](https://clio.readthedocs.io/) data modules.

## Development

We use [`pixi`](https://pixi.sh/) for as our package manager for development.
Once installed, run the following to clone this repo and install all dependencies.

```shell
git clone git@github.com:calliope-project/geo_boundaries_module.git
cd geo_boundaries_module
pixi install --all
```

For testing, simply run:

```shell
pixi run test
```

To view the documentation locally, use:

```shell
pixi run serve-docs
```

To test a minimal example of a workflow using this module:

```shell
pixi shell    # activate this project's environment
cd tests/integration/  # navigate to the integration example
snakemake --use-conda  # run the workflow!
```
