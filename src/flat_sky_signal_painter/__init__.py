"""Flat-sky astrophysical signal painters."""

from .catalog import MovingLensCatalog, PositionCatalog, read_csv
from .geometry import MapGeometry
from .ksz import paint_ksz_gaussian
from .masks import paint_circular_masks
from .moving_lens import paint_moving_lens, single_moving_lens_template
from .profiles import (
    gaussian_profile,
    nfw_deflection_angle,
    truncated_nfw_profile,
    tsz_spectral_factor,
)
from .sources import (
    paint_dust_test_sources,
    paint_localization_disks,
    paint_point_sources,
)
from .tsz import paint_tsz_gaussian

__all__ = [
    "MapGeometry",
    "PositionCatalog",
    "MovingLensCatalog",
    "read_csv",
    "paint_moving_lens",
    "paint_ksz_gaussian",
    "single_moving_lens_template",
    "paint_tsz_gaussian",
    "paint_point_sources",
    "paint_dust_test_sources",
    "paint_localization_disks",
    "paint_circular_masks",
    "truncated_nfw_profile",
    "nfw_deflection_angle",
    "tsz_spectral_factor",
    "gaussian_profile",
]
