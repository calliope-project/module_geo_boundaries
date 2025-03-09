"""Rules used to construct the final dataset."""

rule build_combined_area:
    message:
        "Combine land and marine polygons."
    input:
        countries= [f"resources/automatic/overture/{c['name']}_{c['subtype']}.parquet" for c in config["countries"]],
        marine= "resources/automatic/marineregions/eez.parquet"
    output:
        shapes= "results/shapes.parquet"
    conda:
        "../envs/shape.yaml"
    shell:
        "touch results/shapes.parquet"
