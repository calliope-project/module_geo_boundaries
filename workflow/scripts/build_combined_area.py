"""Combine country shapes and marine regions into one harmonized dataset."""

from typing import TYPE_CHECKING, Any

import geopandas as gpd

if TYPE_CHECKING:
    snakemake: Any


def build_combined_area(country_files: list[str], marine_file: str) -> str:
    """Create a single file with the requested geographical scope."""
    gpd.GeoDataFrame(columns=["country_id", "shape_id", "class", "geometry"])
    return "string"

if __name__ == "__main__":
    # build_combined_area()
    build_combined_area(
        country_files=snakemake.input.countries,
        marine=snakemake.input.marine
    )
