"""Download data from the Marineregions database.

https://www.marineregions.org/
"""

import sys
from pathlib import Path
from shutil import unpack_archive
from tempfile import TemporaryDirectory
from typing import TYPE_CHECKING, Any

import geopandas as gpd
import requests
from _schema import shape_schema

if TYPE_CHECKING:
    snakemake: Any
sys.stderr = open(snakemake.log[0], "w")

MARINE_URL = "https://www.marineregions.org/download_file.php"
FILE = "World_EEZ_v12_20231025_gpkg"


def transform_to_clio(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """Transform the MarineRegions dataset for better clio compatibility.

    - Removes geopolitically contested areas
    - Adds common naming conventions.

    Args:
        gdf (gpd.GeoDataFrame): A marine regions geo-dataframe.

    Returns:
        gpd.GeoDataFrame: standardised dataframe.
    """
    standardised = gpd.GeoDataFrame(
        {
            "shape_id": gdf.apply(
                lambda x: "_".join(
                    [str(x["ISO_TER1"]), "marineregions", str(x["MRGID"])]
                ),
                axis="columns",
            ),
            "country_id": gdf["ISO_TER1"],
            "shape_class": "maritime",
            "geometry": gdf["geometry"],
            "parent": "marineregions",
            "parent_subtype": "eez",
            "parent_id": gdf["MRGID"],
            "parent_name": gdf["GEONAME"],
        }
    )
    # Remove cases without territorial ISO code
    standardised = standardised[~standardised["country_id"].isna()]
    # Check that the base columns fit the schema
    standardised = shape_schema.validate(standardised)
    # Extra: identify contested areas and potential attribution conflicts
    standardised["contested"] = gdf["POL_TYPE"].apply(
        lambda x: True if x in ["Joint regime", "Overlapping claim"] else False
    )
    return standardised


def download_marine_eez_area(path: str):
    """Download data from the marineregions database.

    Adapted from PyPSA-Eur code (MIT licensed).
    https://github.com/PyPSA/pypsa-eur/blob/v2025.01.0/rules/retrieve.smk
    """
    name = "geo-boundaries-module"
    org = "calliope-project"
    response = requests.post(
        MARINE_URL,
        params={"name": FILE + ".zip"},
        data={
            "name": name,
            "organisation": org,
            "email": f"{name}@{org}.org",
            "country": "Netherlands",
            "user_category": "academia",
            "purpose_category": "Research",
            "agree": "1",
        },
    )

    with TemporaryDirectory() as tmpdirname:
        tmp_path = Path(tmpdirname)
        with open(tmp_path / "tmp.zip", "w+b") as file:
            file.write(response.content)
        unpack_archive(tmp_path / "tmp.zip", tmp_path)
        gdf = gpd.read_file(tmp_path / FILE / "eez_v12.gpkg")

    gdf = transform_to_clio(gdf)
    gdf.to_parquet(path)


if __name__ == "__main__":
    download_marine_eez_area(snakemake.output.path)
