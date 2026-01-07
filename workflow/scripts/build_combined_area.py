"""Combine country shapes and marine regions into one harmonized dataset."""

import sys
from typing import TYPE_CHECKING, Any

import _schemas
import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd
from pyproj import CRS

if TYPE_CHECKING:
    snakemake: Any
sys.stderr = open(snakemake.log[0], "w")


def plot_combined_area(combined_file: str, path: str):
    """Generate a nice figure of the resulting file."""
    gdf = gpd.read_parquet(combined_file)
    ax = gdf.plot(figsize=(10, 10), column="shape_class")
    ax.set_xlabel("longitude")
    ax.set_ylabel("latitude")
    ax.set_title("Combined regions")
    plt.savefig(path)


def _remove_overlaps(gdf: gpd.GeoDataFrame, projected_crs: str) -> gpd.GeoDataFrame:
    """Remove overlaps between regional shapes and clip shapes using neighbors.

    Args:
        gdf (gpd.GeoDataFrame): dataframe with regional shapes.
        projected_crs (str): CRS to use. Must be projected.

    Returns:
        gpd.GeoDataFrame: dataframe in the projected CRS.
    """
    # Buffering requires a projected CRS
    assert CRS(projected_crs).is_projected
    projected = gdf.to_crs(projected_crs)

    # A buffer of 0 resolves floating point mismatches that occur during geospatial operations
    buffered = projected.buffer(0)

    for index, row in projected.iterrows():
        minx, miny, maxx, maxy = row.geometry.bounds

        # Find neighbouring regions and only use those for the calculation
        neighbours = buffered.cx[minx:maxx, miny:maxy]
        neighbours = neighbours[neighbours.index != index]
        new_geometry = row["geometry"].difference(neighbours.union_all())

        if not new_geometry.is_valid:
            new_geometry = new_geometry.buffer(0)
            assert new_geometry.is_valid, "Invalid bowties could not be corrected."
        assert new_geometry is not None
        projected.loc[index, "geometry"] = new_geometry
    return projected


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

    frames = []
    marine = gpd.read_parquet(marine_file)
    marine = marine.to_crs(geographic_crs)
    # Combine land and marine boundary for each country
    for file in country_files:
        # Fetch the country file and ensure crs is compatible
        country_land = gpd.read_parquet(file).to_crs(geographic_crs)
        country_id = country_land["country_id"].unique()
        if len(country_id) != 1:
            raise ValueError(
                f"Country file {file} should be a single country. Found {country_id}."
            )
        country_id = country_id[0]

        if country_id in marine["country_id"].unique():
            # clip maritime boundaries
            country_marine = marine[marine["country_id"] == country_id]
            marine_geom = country_marine.geometry.union_all()
            country_land = country_land.copy()
            country_land.geometry = country_land.geometry.difference(marine_geom)

            # only add uncontested maritime boundaries
            uncontested = country_marine[~country_marine["contested"]].drop(
                "contested", axis="columns"
            )
            frames.extend([country_land, uncontested])
        else:
            frames.append(country_land)

    combined = gpd.GeoDataFrame(
        pd.concat(frames, ignore_index=True), crs=geographic_crs
    )
    return combined


def build_combined_area(
    country_files: list[str], marine_file: str, crs: dict[str, str], combined_file: str
):
    """Create a single file with the requested geographical scope."""
    combined = _combine_shapes(country_files, marine_file, crs["geographic"])
    combined = _remove_overlaps(combined, crs["projected"])

    combined = combined.to_crs(crs["geographic"])
    # A buffer of 0 resolves floating point mismatches that occur during CRS conversion
    combined["geometry"] = combined.buffer(0)
    combined = _schemas.ShapesSchema.validate(combined)
    combined.reset_index(drop=True).to_parquet(combined_file)


if __name__ == "__main__":
    build_combined_area(
        country_files=snakemake.input.countries,
        marine_file=snakemake.input.marine,
        crs=snakemake.params.crs,
        combined_file=snakemake.output.combined,
    )
    plot_combined_area(
        combined_file=snakemake.output.combined, path=snakemake.output.plot
    )
