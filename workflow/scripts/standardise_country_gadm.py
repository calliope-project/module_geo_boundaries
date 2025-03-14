"""Transform a GADM country dataset into a cross-compatible shape."""

from typing import TYPE_CHECKING, Any

import geopandas as gpd

if TYPE_CHECKING:
    snakemake: Any


def standardise_country_gadm(
    input_path: str, country_id: str, subtype: str, output_path: str
):
    """Transformation of GADM dataset to clio.

    Args:
        input_path (str): path to input file.
        country_id (str): ISO alpha 3 code.
        subtype (str): regional definition in the parent dataset (0, 1, 2).
        output_path (str): path to output file.
    """
    gdf = gpd.read_parquet(input_path)
    standardised = gpd.GeoDataFrame(
        {
            "shape_id": gdf[f"GID_{subtype}"].apply(lambda x: f"{country_id}_gadm_{x}"),
            "country_id": gdf["GID_0"],
            "class": "land",
            "geometry": gdf["geometry"],
            "parent": "gadm",
            "parent_subtype": subtype,
            "parent_id": gdf[f"GID_{subtype}"],
            "parent_name": gdf[f"NAME_{subtype}"],
        }
    )
    standardised.to_parquet(output_path)


if __name__ == "__main__":
    standardise_country_gadm(
        input_path=snakemake.input.raw,
        country_id=snakemake.params.country_id,
        subtype=snakemake.params.subtype,
        output_path=snakemake.output.standardised,
    )
