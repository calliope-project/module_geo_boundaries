"""Rules to used to download automatic resource files."""


rule download_country_area:
    message:
        "Download country area polygons: {wildcards.country}_{wildcards.subtype}."
    params:
        version=internal["resources"]["overture_release"],
    output:
        path="resources/automatic/overture/{country}_{subtype}.parquet",
    conda:
        "../envs/shape.yaml"
    script: "../scripts/download_country_area.py"


rule download_marine_eez_area:
    message:
        "Download global exclusive economic zone (eez) polygons."
    output:
        path="resources/automatic/marineregions/eez.parquet",
    conda:
        "../envs/shape.yaml"
    script: "../scripts/download_marine_eez_area.py"
