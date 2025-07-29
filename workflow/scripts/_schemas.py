import pandera.pandas as pa
from pandera.typing.geopandas import GeoSeries
from pandera.typing.pandas import Series


class ShapesSchema(pa.DataFrameModel):
    """Schema for geographic shapes."""

    shape_id: Series[str] = pa.Field(unique=True)
    "A unique identifier for this shape."
    country_id: Series[str]
    "Country ISO alpha-3 code."
    shape_class: Series[str] = pa.Field(isin=["land", "maritime"])
    "Identifier of the shape's context."
    geometry: GeoSeries
    "Shape (multi)polygon."
    parent: Series[str] = pa.Field(isin=["gadm", "overture", "marineregions", "nuts"])
    "Parent dataset."
    parent_subtype: Series[str]
    "Region disaggregation level in the parent dataset."
    parent_id: Series[str]
    "Unique id in the parent dataset."
    parent_name: Series[str]
    "Human-readable name in the parent dataset."

    class Config:
        # top-level schema options from your YAML
        coerce = True
        strict = True
