"""Generate documentation for all the existing modules."""

from pathlib import Path
from textwrap import indent

from cffconvert import Citation
from clio_tools.data_module import ModuleInterface
from yaml import safe_dump, safe_load

MODULE_PATH = Path("./")
DOCS_FILE_PATH = MODULE_PATH / "docs" / "specification.md"

TEMPLATE = """
# Specification

This specification is based on `clio` standards.
For more information, please consult the [`clio` documentation](https://clio.readthedocs.io/).

## Input - Output structure

Here is a general summary of the required **resources**, **results** and **wildcards** needed to use this module.

- `resources/` are files needed by the module to operate.
    - `user/` resource files must be provided before executing the module.
    - `automatic/` resources will be downloaded at the start of the module's execution.
- `results/` are files produced by the module. We only list the most relevant ones here, others might be generated in intermediate steps.
- `wildcards` are filename patterns that can alter the module's behaviour without changing its configuration.

???+ info "Interfacing summary"

    ```yaml
{interface}
    ```

???+ info "Visual summary"

    ```mermaid
{flowchart}
    ```


## Configuration

We recommend to start with the following configuration.
All configuration options can be found in the module's schema.

???+ example "Example configuration"

    ```yaml
    module_geo_boundaries_module:
        --8<-- "{example_config}"
    ```

??? info "Configuration schema"

    ```yaml
    --8<-- "{schema}"
    ```

## How to use

1. Add the necessary module [configuration](#configuration) to your `snakemake` workflow.
2. Include this module in your workflow specifying a version `tag:` (to protect against future changes) and a `prefix:` (to redirect where the module will place files).

    ??? example "Importing this module"

        ```python
        module geo_boundaries_module:
          snakefile: github(
            "calliope-project/geo_boundaries_module",
            path="workflow/Snakefile",
            tag="vX.Y.Z"
          )
          config: config["module_geo_boundaries_module"]
          prefix: "results/geo_boundaries_module"

        use rule * from geo_boundaries_module as module_geo_boundaries_module_*
        ```

For more information, please consult the [`snakemake` documentation](https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html).

## Attribution

!!! quote "Citation"

    {citation}

??? info "Contributors"

    --8<-- "{authors}"

??? info "License"

    ```txt
    --8<-- "{license}"
    ```

"""


def on_pre_build(config):
    """Generate standard documentation page per module.

    Automatically called by mkdocs if the hook is configured.
    """
    create_module_docfile()


def on_post_build(config):
    """Remove automatically generated files."""
    DOCS_FILE_PATH.unlink()


def create_module_docfile():
    """Save a fully documented page for the requested module."""
    name = "geo_boundaries_module"

    authors = MODULE_PATH / "AUTHORS"
    example_config = MODULE_PATH / "config/example.yaml"
    schema = MODULE_PATH / "workflow/internal/config.schema.yaml"
    citation = MODULE_PATH / "CITATION.cff"
    license = MODULE_PATH / "LICENSE"
    interface = MODULE_PATH / "INTERFACE.yaml"
    assert all(
        [
            file.exists()
            for file in [authors, example_config, schema, citation, license, interface]
        ]
    )

    flowchart = ModuleInterface.from_yaml(interface).to_mermaid_flowchart(name)

    text = TEMPLATE.format(
        interface=indent(safe_dump(safe_load(interface.read_text())), "\t"),
        flowchart=indent(flowchart, "\t"),
        example_config=example_config,
        schema=schema,
        citation=get_apa_citation(citation),
        authors=authors,
        license=license,
    )

    with open(DOCS_FILE_PATH, "w") as file:
        file.write(text)


def get_apa_citation(citation_file: Path) -> str:
    """Return APA citation if the .cff file is correct."""
    cff = Citation(citation_file.read_text())
    return cff.as_apalike()
