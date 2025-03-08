"""Rules to used to download automatic resource files."""


rule download_country_shapes:
    message:
        "Download country shape: {wildcards.country}-{wildcards.subtype}."
    params:
        version=internal["resources"]["overture_release"],
    output:
        path="resources/automatic/overture/{country}_{subtype}.parquet",
    conda:
        "../envs/shape.yaml"
    script: "../scripts/download_country_shapes.py"


rule download_marine_eez_shapes:
    message:
        "Download global exclusive economic zone (eez) data."
    output:
        path="resources/automatic/marineregions/eez.parquet",
    conda:
        "../envs/shape.yaml"
    script: "../scripts/download_marine_shapes.py"
