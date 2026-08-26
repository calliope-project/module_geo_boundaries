"""Download data from the GADM database.

https://gadm.org/index.html
"""

import sys
import tempfile
from pathlib import Path
from typing import TYPE_CHECKING, Any

import geopandas as gpd
from _utils import DownloadTimeouts, download_file, read_geojson_file

if TYPE_CHECKING:
    snakemake: Any


URL = "https://geodata.ucdavis.edu/gadm/gadm{release}/json/gadm{nodot}_{country}_{subtype}.json{zip}"
CRS = "EPSG:4326"
SUPPORTED = ("4.1",)


def download_country_gadm(
    release: str,
    country: str,
    subtype: str,
    timeouts: DownloadTimeouts,
    geojson_max_obj_size_mb: int,
) -> gpd.GeoDataFrame:
    """Attempts to download country GADM data in .json or zipped json."""
    last_error: Exception | None = None

    if release not in SUPPORTED:
        raise ValueError(f"GADM {release=} is not supported.")

    for zip_ext in (".zip", ""):
        url = URL.format(
            release=release,
            nodot=release.replace(".", ""),
            country=country,
            subtype=subtype,
            zip=zip_ext,
        )
        try:
            with tempfile.TemporaryDirectory() as tmp_dir:
                tmp_path = Path(tmp_dir) / f"download.json{zip_ext}"

                download_file(url, tmp_path, timeouts)
                gdf = read_geojson_file(tmp_path, geojson_max_obj_size_mb)
                if gdf.empty:
                    raise RuntimeError(f"Downloaded empty GADM file from {url!r}.")
                return gdf.to_crs(CRS)

        except Exception as exc:
            last_error = exc
    raise RuntimeError(
        f"Could not fetch GADM request for {country!r}:{subtype!r}."
    ) from last_error


def main():
    """Main snakemake process."""
    timeouts = DownloadTimeouts(**snakemake.params.timeouts)
    country = download_country_gadm(
        snakemake.wildcards.release,
        snakemake.wildcards.country,
        snakemake.wildcards.subtype,
        timeouts,
        snakemake.params.geojson_max_obj_size_mb,
    )
    country.to_parquet(snakemake.output.path)


if __name__ == "__main__":
    sys.stderr = open(snakemake.log[0], "w", buffering=1)
    main()
