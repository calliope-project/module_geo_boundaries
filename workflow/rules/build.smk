"""Rules used to construct the final dataset."""


rule build_combined_area:
    message:
        "Combine land and marine polygons."
    params:
        crs=internal["standardisation"]["crs"],
    input:
        countries=[
            f"resources/automatic/countries/{country['source']}_{country['country_id']}_{country['subtype']}.parquet"
            for country in config["countries"]
        ],
        marine="resources/automatic/marineregions/eez.parquet",
    output:
        combined="results/shapes.parquet",
    conda:
        "../envs/shape.yaml"
    script:
        "../scripts/build_combined_area.py"
