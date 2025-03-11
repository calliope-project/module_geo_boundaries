"""General schema that the produced dataset should follow."""

import pandera as pa

schema = pa.DataFrameSchema(
    {
        "shape_id": pa.Column(str),
        "country_id": pa.Column(str),
        "class": pa.Column(str, pa.Check.isin(["land", "maritime"])),
        "geometry": pa.Column("geometry"),
        "parent": pa.Column(str, pa.Check.isin(["gadm", "overture", "marineregions"])),
        "parent_subtype": pa.Column(str),
        "parent_id": pa.Column(str),
        "parent_name": pa.Column(str),
    },
    coerce=True,
    strict=True
)
