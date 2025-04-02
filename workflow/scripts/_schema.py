"""General schema that the produced dataset should follow."""

from pandera import Check, Column, DataFrameSchema

shape_schema = DataFrameSchema(
    {
        "shape_id": Column(
            str, unique=True, description="A unique identifier for this shape."
        ),
        "country_id": Column(str, description="Country ISO alpha-3 code."),
        "shape_class": Column(
            str,
            checks=Check.isin(["land", "maritime"]),
            description="Identifier of the shape's context.",
        ),
        "geometry": Column("geometry", description="Shape (multi)polygon."),
        "parent": Column(
            str,
            checks=Check.isin(["gadm", "overture", "marineregions", "nuts"]),
            description="Parent dataset.",
        ),
        "parent_subtype": Column(
            str, description="Region disaggregation level in the parent dataset."
        ),
        "parent_id": Column(str, description="Unique id in the parent dataset."),
        "parent_name": Column(
            str, description="Human-readable name in the parent dataset."
        ),
    },
    coerce=True,
    strict=True,
)
