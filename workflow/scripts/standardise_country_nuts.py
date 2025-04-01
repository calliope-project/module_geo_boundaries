"""Download a NUTS map from the GISCO website."""

import sys
from typing import TYPE_CHECKING, Any

import geopandas as gpd
import pycountry
import requests

if TYPE_CHECKING:
    snakemake: Any
sys.stderr = open(snakemake.log[0], "w")

NUTS_TO_ISO2 = {
    "EL": "GR",  # Greece
    "UK": "GB",  # United Kingdom
}
URL = "https://gisco-services.ec.europa.eu/distribution/v2/nuts/geojson/NUTS_RG_{resolution}_{year}_{espg}_LEVL_{level}.geojson"

def _nuts_to_iso_3(code):
    iso2 = NUTS_TO_ISO2.get(code, None)
    if not iso2:
        iso2 = code
    return pycountry.countries.get(alpha_2=iso2).alpha_3
