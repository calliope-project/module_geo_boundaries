"""Download data from the GADM database.

Built for version 4.1 of the dataset.

https://gadm.org/index.html
"""

import sys
from typing import TYPE_CHECKING, Any

import geopandas as gpd
import requests

if TYPE_CHECKING:
    snakemake: Any


GADM_URL = (
    "https://geodata.ucdavis.edu/gadm/gadm4.1/json/gadm41_{country}_{subtype}.json{zip}"
)
GADM_CRS = "EPSG:4326"


def download_country_gadm(country: str, subtype: str) -> gpd.GeoDataFrame:
    """Download country and save it to parquet."""
    session = requests.Session()
    gdf: gpd.GeoDataFrame | None = None

    for zip_ext in (".zip", ""):
        uri = GADM_URL.format(country=country, subtype=subtype, zip=zip_ext)
        try:
            r = session.get(uri, stream=True, timeout=30)
            r.raise_for_status()
            # URL is valid, read and format
            gdf = gpd.read_file(uri).to_crs(GADM_CRS)
            break
        except Exception:
            continue

    if gdf is None:
        raise RuntimeError(f"Could not fetch GADM request for {country!r}:{subtype!r}.")

    return gdf


def main():
    """Main snakemake process."""
    country = download_country_gadm(
        snakemake.wildcards.country, snakemake.wildcards.subtype
    )
    country.to_parquet(snakemake.output.path)


if __name__ == "__main__":
    sys.stderr = open(snakemake.log[0], "w")
    main()
