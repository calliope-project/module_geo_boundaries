"""Download data from the Marineregions database.

https://www.marineregions.org/
"""

from pathlib import Path
from shutil import unpack_archive
from tempfile import TemporaryDirectory
from typing import TYPE_CHECKING, Any

import geopandas as gpd
import requests

if TYPE_CHECKING:
    snakemake: Any

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
    # remove contested areas and potential attribution conflicts
    std_gdf = gdf[~gdf["POL_TYPE"].isin(["Joint regime", "Overlapping claim"])]

    # Standardise
    std_gdf = std_gdf.rename(
        columns={
            "ISO_TER1": "country_id",
            "ISO_SOV1": "sovereign_country_id",
            "MRGID": "parent_id",
            "GEONAME": "parent_name",
        }
    )
    std_gdf["shape_id"] = std_gdf.apply(
        lambda x: "_".join([str(x["country_id"]), "marineregions", str(x["parent_id"])]), axis="columns"
    )
    std_gdf["class"] = "maritime"
    std_gdf["parent"] = "marineregions"
    std_gdf["parent_subtype"] = "eez"
    # slim the data
    std_gdf = std_gdf[["shape_id", "country_id", "class", "geometry", "parent", "parent_subtype", "parent_id", "parent_name"]]
    return std_gdf


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
