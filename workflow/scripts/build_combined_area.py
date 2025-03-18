"""Combine country shapes and marine regions into one harmonized dataset."""

from typing import TYPE_CHECKING, Any

import geopandas as gpd
import pandas as pd
from _schema import schema
from pyproj import CRS

if TYPE_CHECKING:
    snakemake: Any


def _remove_overlaps(
    gdf: gpd.GeoDataFrame, buffer: float, projected_crs: str
) -> gpd.GeoDataFrame:
    """Remove overlaps between regional shapes.

    Applies a buffer to all shapes, and then clips each shape using its neighbors.

    Args:
        gdf (gpd.GeoDataFrame): dataframe with regional shapes.
        buffer (float): buffer radius (unit depends on the projected CRS used).
        projected_crs (str): CRS to use. Must be projected.

    Returns:
        gpd.GeoDataFrame: buffered dataframe.
    """
    # Buffering requires a projected CRS
    assert CRS(projected_crs).is_projected
    original_crs = gdf.crs
    projected = gdf.to_crs(projected_crs)

    # Buffer size is divided by two since it is applied at both sides of a border
    buffered = projected.buffer(buffer / 2)

    for index, row in projected.iterrows():
        minx, miny, maxx, maxy = row.geometry.bounds

        # Find neighbouring regions and only use those for the calculation
        neighbours = buffered.cx[minx:maxx, miny:maxy]
        neighbours = neighbours[neighbours.index != index]
        new_geometry = row["geometry"].difference(neighbours.union_all())

        # Ensure the resulting shape has no bowties
        if not new_geometry.is_valid:
            new_geometry = new_geometry.buffer(0)
            assert new_geometry.is_valid
        projected.loc[index, "geometry"] = new_geometry

    return projected.to_crs(original_crs)


def _combine_shapes(
    country_files: list[str], marine_file: str, geographic_crs: str
) -> gpd.GeoDataFrame:
    """Merge each all countries and maritime boundaries into one file.

    Args:
        country_files (list[str]): List of standardised country files to combine.
        marine_file (str): Standardised file with maritime boundaries.
        geographic_crs (str): CRS to use. Must be geographic.

    Raises:
        ValueError: Country file is not a unique country.

    Returns:
        gpd.GeoDataFrame: Combined dataframe using the given CRS.
    """
    assert CRS(geographic_crs).is_geographic
    combined = gpd.GeoDataFrame(columns=schema.columns)
    combined = combined.set_crs(geographic_crs)

    marine = gpd.read_parquet(marine_file)
    marine = marine.to_crs(geographic_crs)
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
        country_land = country_land.to_crs(geographic_crs)

        if country_id in marine["country_id"].unique():
            # clip maritime boundaries
            country_marine = marine[marine["country_id"] == country_id]
            country_land["geometry"] = country_land.difference(
                country_marine.union_all(), align=False
            )
            # only add uncontested maritime boundaries
            uncontested = country_marine[~country_marine["contested"]].drop(
                "contested", axis="columns"
            )
            combined = pd.concat(
                [combined, country_land, uncontested], ignore_index=True
            )
        else:
            combined = pd.concat([combined, country_land], ignore_index=True)

    combined = combined.reset_index(drop=True)
    return combined


def build_combined_area(
    country_files: list[str],
    marine_file: str,
    crs: dict[str, str],
    buffer: float,
    combined_file: str,
):
    """Create a single file with the requested geographical scope."""
    combined = _combine_shapes(country_files, marine_file, crs["geographic"])
    if buffer > 0:
        combined = _remove_overlaps(combined, buffer, crs["projected"])

    combined = schema(combined)
    combined.to_parquet(combined_file)


if __name__ == "__main__":
    # build_combined_area()
    build_combined_area(
        country_files=snakemake.input.countries,
        marine_file=snakemake.input.marine,
        crs=snakemake.params.crs,
        buffer=snakemake.params.buffer,
        combined_file=snakemake.output.combined,
    )
