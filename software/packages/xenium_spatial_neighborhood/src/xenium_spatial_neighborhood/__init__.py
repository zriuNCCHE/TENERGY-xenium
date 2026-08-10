"""Spatial cellular-neighborhood analysis for Xenium-style AnnData objects."""

from .spatial_neighborhood_analyzer import (
    SpatialNeighborhoodAnalyzer,
    SpatialNeighborhoodCalculator,
    SpatialNeighborhoodConfig,
)

__all__ = [
    "SpatialNeighborhoodAnalyzer",
    "SpatialNeighborhoodCalculator",
    "SpatialNeighborhoodConfig",
]

__version__ = "0.2.4"
