"""Download data from the geoBoundaries repository.

https://github.com/wmgeolab/geoBoundaries
"""

import sys
import tempfile
from pathlib import Path
from typing import TYPE_CHECKING, Any

import geopandas as gpd
from _utils import DownloadTimeouts, download_file, read_geojson_file

if TYPE_CHECKING:
    snakemake: Any


URL = (
    "https://github.com/wmgeolab/geoBoundaries/raw/refs/tags/v{release}/releaseData/"
    "{release_type}/{country}/ADM{subtype}/geoBoundaries-{country}-ADM{subtype}.geojson"
)
CRS = "EPSG:4326"


def download_country_geoboundaries(
    country: str,
    subtype: str,
    release: str,
    release_type: str,
    timeouts: DownloadTimeouts,
    geojson_max_obj_size_mb: int,
) -> gpd.GeoDataFrame:
    """Download country data from geoBoundaries."""
    geojson_url = URL.format(
        release=release, release_type=release_type, country=country, subtype=subtype
    )

    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)

        geojson_path = tmp_path / "download.geojson"
        download_file(geojson_url, geojson_path, timeouts)

        gdf = read_geojson_file(geojson_path, geojson_max_obj_size_mb)
        if gdf.empty:
            raise RuntimeError(
                f"Downloaded empty geoBoundaries file from {geojson_url!r}."
            )

    return gdf.to_crs(CRS)


def main() -> None:
    """Main snakemake process."""
    timeouts = DownloadTimeouts(**snakemake.params.timeouts)

    country = download_country_geoboundaries(
        snakemake.wildcards.country,
        snakemake.wildcards.subtype,
        snakemake.wildcards.release,
        snakemake.wildcards.release_type,
        timeouts,
        snakemake.params.geojson_max_obj_size_mb,
    )
    country.to_parquet(snakemake.output.path)


if __name__ == "__main__":
    sys.stderr = open(snakemake.log[0], "w", buffering=1)
    main()
