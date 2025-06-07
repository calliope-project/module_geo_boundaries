"""Rules used to construct the final dataset."""


rule build_combined_area:
    message:
        "Combine land and marine polygons."
    params:
        buffer=config["buffer"],
        crs=config["crs"],
    input:
        countries=[
            f"resources/automatic/countries/{data['source']}_{country}_{data['subtype']}.parquet"
            for country, data in config["countries"].items()
        ],
        marine="resources/automatic/marineregions/eez.parquet",
        schema=workflow.source_path("../internal/shape.schema.yaml"),
    output:
        combined="results/shapes.parquet",
        plot=report(
            "results/shapes.png",
            caption="../report/results.rst",
            category="Combined area",
        ),
    log:
        "logs/build_combined_area.log",
    conda:
        "../envs/shape.yaml"
    script:
        "../scripts/build_combined_area.py"
