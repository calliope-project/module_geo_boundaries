"""Rules to used to download resource files.

Small transformations might be performed to make the data easier to work with.
"""


rule download_country_overture:
    message:
        "Download '{wildcards.country}_{wildcards.subtype}' dataset from Overture Maps."
    params:
        version=internal["resources"]["overture_release"],
    output:
        path="resources/automatic/countries/overture_{country}_{subtype}.parquet",
    conda:
        "../envs/shape.yaml"
    script: "../scripts/download_country_overture.py"


rule download_country_gadm:
    message:
        "Download '{wildcards.country}_{wildcards.subtype}' dataset from GADM."
    params:
        admin="{country}",
        content_level=lambda wc: int(wc.subtype)
    output:
        path=temp("resources/automatic/countries/raw_gadm_{country}_{subtype}.parquet")
    wrapper: "master/geo/pygadm/item"


rule standardise_country_gadm:
    message:
        "Standardise '{wildcards.country}_{wildcards.subtype}' dataset."
    params:
        country_id= lambda wc: str(wc.country),
        subtype= lambda wc: str(wc.subtype)
    input:
        raw="resources/automatic/countries/raw_gadm_{country}_{subtype}.parquet"
    output:
        standardised= "resources/automatic/countries/gadm_{country}_{subtype}.parquet"
    conda:
        "../envs/shape.yaml"
    script: "../scripts/standardise_country_gadm.py"


rule download_marine_eez_area:
    message:
        "Download global exclusive economic zone (eez) polygons."
    output:
        path="resources/automatic/marineregions/eez.parquet",
    conda:
        "../envs/shape.yaml"
    script: "../scripts/download_marine_eez_area.py"
