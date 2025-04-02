"""General schema that the produced dataset should follow."""

from pandera import DataFrameModel, Field
from pandera.typing import Series
from pandera.typing.geopandas import GeoDataFrame, GeoSeries  # noqa: F401


class ShapeSchema(DataFrameModel):
    shape_id: Series[str] = Field(unique=True)
    """A unique identifier for this shape."""
    country_id: Series[str]
    """Country ISO alpha-3 code."""
    shape_class: Series[str] = Field(isin=["land", "maritime"])
    """Identifier of the shape's context."""
    geometry: GeoSeries
    """Shape (multi)polygon."""
    parent: Series[str] = Field(isin=["gadm", "overture", "marineregions", "nuts"])
    """Parent dataset."""
    parent_subtype: Series[str]
    """Region disaggregation level in the parent dataset."""
    parent_id: Series[str]
    """Unique id in the parent dataset."""
    parent_name: Series[str]
    """Human-readable name in the parent dataset."""

    class Config:
        coerce = True
        strict = True
