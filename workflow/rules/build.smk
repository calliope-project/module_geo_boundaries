"""Rules used to construct the final dataset."""


rule build_combined_area:
    message:
        "Combine land and marine polygons."
    params:
        crs=internal["standardisation"]["crs"],
    input:
        countries=[
            f"resources/automatic/countries/{c['source']}_{c['country_id']}_{c['subtype']}.parquet"
            for c in config["countries"]
        ],
        marine="resources/automatic/marineregions/eez.parquet",
    output:
        combined="results/shapes.parquet",
    conda:
        "../envs/shape.yaml"
    script:
        "../scripts/build_combined_area.py"
