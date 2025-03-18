"""Download division area data from the Overture Maps foundation."""

from typing import TYPE_CHECKING, Any

import duckdb
import pycountry

if TYPE_CHECKING:
    snakemake: Any

OVERTURE_LINK = (
    "s3://overturemaps-us-west-2/release/{version}/theme=divisions/type=division_area/*"
)


def download_country_overture(country: str, subtype: str, version: str, path: str):
    """Download country division areas from Overture maps.

    Uses duckdb for remote interfacing and 'larger than memory' file generation.
    """
    # Prepare variables for the request
    country_a2 = pycountry.countries.get(alpha_3=country).alpha_2

    # Setup SQL connection to the remote dataset
    connection = duckdb.connect()
    for extension in ["spatial", "httpfs"]:
        connection.install_extension(extension)
        connection.load_extension(extension)
    connection.sql("SET s3_region='us-west-2'")

    # Request country dataset with added metadata
    connection.sql(
        f"""
        COPY (
            SELECT
                '{country}' || '_' || 'overture' || '_' || id AS shape_id,
                '{country}' AS country_id,
                class AS shape_class,
                geometry,
                'overture' AS parent,
                subtype AS parent_subtype,
                id AS parent_id,
                names.primary AS parent_name
            FROM
                read_parquet(
                    '{OVERTURE_LINK.format(version=version)}',
                    filename=true,
                    hive_partitioning=true
                )
            WHERE
                country == '{country_a2}'
                AND subtype == '{subtype}'
        )
        TO '{path}'
        WITH (
            FORMAT parquet,
            COMPRESSION zstd
        );
        """
    )


if __name__ == "__main__":
    download_country_overture(
        country=snakemake.wildcards.country,
        subtype=snakemake.wildcards.subtype,
        version=snakemake.params.version,
        path=snakemake.output.path,
    )
