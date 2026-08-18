"""Rules to used to harmonise country files."""


rule harmonise_geoboundaries:
    input:
        raw=rules.download_geoboundaries.output.path,
    output:
        path="<resources>/automatic/geoboundaries/harmonise/{release}/{country}_{subtype}_{release_type}.parquet",
    log:
        "<logs>/geoboundaries/harmonise/{release}/{country}_{subtype}_{release_type}.log",
    conda:
        "../envs/module.yaml"
    message:
        "Harmonising geoBoundaries {wildcards.release}: {wildcards.country}_{wildcards.subtype}_{wildcards.release_type}."
    script:
        "../scripts/harmonise_geoboundaries.py"


rule download_harmonised_overture:
    input:
        duckdb_extensions=rules.download_duckdb_extensions.output.path,
    output:
        path="<resources>/automatic/overture/harmonise/{release}/{country}_{subtype}.parquet",
    log:
        "<logs>/overture/download_and_harmonise/{release}/{country}_{subtype}.log",
    localrule: True
    conda:
        "../envs/module.yaml"
    message:
        "Downloading|harmonising Overture {wildcards.release}: {wildcards.country}_{wildcards.subtype}."
    script:
        "../scripts/download_harmonised_overture.py"


rule harmonise_gadm:
    input:
        raw=rules.download_gadm.output.path,
    output:
        standardised="<resources>/automatic/gadm/harmonise/{release}/{country}_{subtype}.parquet",
    log:
        "<logs>/gadm/harmonise/{release}/{country}_{subtype}.log",
    conda:
        "../envs/module.yaml"
    message:
        "Harmonising GADM {wildcards.release}: {wildcards.country}_{wildcards.subtype}."
    script:
        "../scripts/harmonise_gadm.py"


rule harmonise_nuts:
    input:
        raw=rules.download_nuts.output.path,
    output:
        path="<resources>/automatic/nuts/harmonise/{release}/{country}_{subtype}_{resolution}.parquet",
    log:
        "<logs>/nuts/harmonise/{release}/{country}_{subtype}_{resolution}.log",
    conda:
        "../envs/module.yaml"
    message:
        "Harmonising NUTS {wildcards.release}: {wildcards.subtype}_{wildcards.resolution}."
    script:
        "../scripts/harmonise_nuts.py"


rule download_harmonised_eez:
    output:
        path="<resources>/automatic/eez/single/{release}/{eez}.parquet",
        plot=report(
            "<resources>/automatic/eez/single/{release}/{eez}.png",
            caption="../report/download_harmonised_eez.rst",
            category="Module Geo-Boundaries",
            subcategory="EEZ area",
        ),
    log:
        "<logs>/eez/harmonise/{release}/{eez}.log",
    localrule: True
    conda:
        "../envs/module.yaml"
    params:
        timeouts=internal["timeouts"],
    message:
        "Downloading|harmonising MarineRegions {wildcards.release}: {wildcards.eez}."
    script:
        "../scripts/download_harmonised_eez.py"
