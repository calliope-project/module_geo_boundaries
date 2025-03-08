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


def download_marine_shape(path: str):
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

    gdf.to_parquet(path)


if __name__ == "__main__":
    download_marine_shape(snakemake.output.path)
