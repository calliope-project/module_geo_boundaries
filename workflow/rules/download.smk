"""Rules to used to download resource files."""


rule download_duckdb_extensions:
    output:
        path="<resources>/automatic/overture/duckdb_extensions.txt",
    log:
        "<logs>/download_duckdb_extensions.log",
    localrule: True
    conda:
        "../envs/module.yaml"
    threads: 1
    message:
        "Downloading DuckDB extensions."
    script:
        "../scripts/download_duckdb_extensions.py"


rule download_geoboundaries:
    output:
        path="<resources>/automatic/geoboundaries/download/{release}/{country}_{subtype}_{release_type}.parquet",
    log:
        "<logs>/geoboundaries/download/{release}/{country}_{subtype}_{release_type}.log",
    localrule: True
    conda:
        "../envs/module.yaml"
    params:
        timeouts=internal["timeouts"],
        geojson_max_obj_size_mb=get_gdal_config()["geojson_max_obj_size_mb"],
    message:
        "Downloading geoBoundaries {wildcards.release}: {wildcards.country}_{wildcards.subtype}_{wildcards.release_type}."
    script:
        "../scripts/download_geoboundaries.py"


rule download_gadm:
    output:
        path="<resources>/automatic/gadm/download/{release}/{country}_{subtype}.parquet",
    log:
        "<logs>/gadm/download/{release}/{country}_{subtype}.log",
    localrule: True
    conda:
        "../envs/module.yaml"
    params:
        timeouts=internal["timeouts"],
        geojson_max_obj_size_mb=get_gdal_config()["geojson_max_obj_size_mb"],
    message:
        "Downloading GADM {wildcards.release}: {wildcards.country}_{wildcards.subtype}."
    script:
        "../scripts/download_gadm.py"


rule download_nuts:
    output:
        path="<resources>/automatic/nuts/download/{release}/{subtype}_{resolution}.parquet",
    log:
        "<logs>/nuts/download/{release}/{subtype}_{resolution}.log",
    localrule: True
    conda:
        "../envs/module.yaml"
    params:
        timeouts=internal["timeouts"],
    message:
        "Downloading NUTS {wildcards.release}: {wildcards.subtype}_{wildcards.resolution}."
    script:
        "../scripts/download_nuts.py"
