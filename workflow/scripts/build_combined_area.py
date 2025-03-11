"""Combine country shapes and marine regions into one harmonized dataset."""

from typing import TYPE_CHECKING, Any

import geopandas as gpd
import pandas as pd

if TYPE_CHECKING:
    snakemake: Any


STANDARD_COLS = [
    "shape_id",
    "country_id",
    "class",
    "geometry",
    "parent",
    "parent_subtype",
    "parent_id",
    "parent_name",
]


def build_combined_area(
    country_files: list[str], marine_file: str, crs: str, combined_file: str
):
    """Create a single file with the requested geographical scope."""
    combined = gpd.GeoDataFrame(columns=STANDARD_COLS)
    combined = combined.set_crs(crs)
    marine = gpd.read_parquet(marine_file)
    marine = marine.to_crs(crs)
    # Combine land and marine boundary for each country
    for file in country_files:
        # Fetch the country file and ensure crs is compatible
        country_land = gpd.read_parquet(file)
        country_id = country_land["country_id"].unique()
        if len(country_id) != 1:
            raise ValueError(
                f"Country file {file} should be a single country. Found {country_id}."
            )
        country_id = country_id[0]
        country_land = country_land.to_crs(crs)

        # if the country has maritime boundaries, clip them
        if country_id in marine["country_id"].unique():
            country_marine = marine[marine["country_id"] == country_id]
            country_land["geometry"] = country_land.difference(
                country_marine.union_all(), align=False
            )
            combined = pd.concat(
                [combined, country_land, country_marine], ignore_index=True
            )
        else:
            combined = pd.concat([combined, country_land], ignore_index=True)

    combined = combined.reset_index(drop=True)
    # TODO: improve standardisation of the datasets
    combined["parent_id"] = combined["parent_id"].apply(str)
    combined.to_parquet(combined_file)


if __name__ == "__main__":
    # build_combined_area()
    build_combined_area(
        country_files=snakemake.input.countries,
        marine_file=snakemake.input.marine,
        crs=snakemake.params.crs,
        combined_file=snakemake.output.combined,
    )
