"""Square flat-sky geometry."""

from __future__ import annotations
from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class MapGeometry:
    """Square tangent-plane geometry.

    Coordinates are angular offsets in radians around a configurable center.
    """

    npix: int = 256
    size_deg: float = 4.0
    center_ra_rad: float = 0.0
    center_dec_rad: float = 0.0

    def __post_init__(self):
        if self.npix < 2:
            raise ValueError("npix must be >= 2")
        if self.size_deg <= 0:
            raise ValueError("size_deg must be positive")

    @property
    def size_rad(self) -> float:
        return np.deg2rad(self.size_deg)

    @property
    def pixel_size_rad(self) -> float:
        return self.size_rad / (self.npix - 1)

    @property
    def resolution_arcmin(self) -> float:
        return np.rad2deg(self.pixel_size_rad) * 60.0

    @property
    def extent_deg(self):
        h = self.size_deg / 2.0
        return (-h, h, -h, h)

    def coordinate_vectors(self):
        half = 0.5 * self.size_rad
        ra = np.linspace(self.center_ra_rad - half, self.center_ra_rad + half, self.npix)
        dec = np.linspace(self.center_dec_rad - half, self.center_dec_rad + half, self.npix)
        return ra, dec

    def coordinate_grids(self):
        ra, dec = self.coordinate_vectors()
        ra_grid, dec_grid = np.meshgrid(ra, dec, indexing="xy")
        return ra_grid, dec_grid

    def nearest_pixel(self, ra_rad: float, dec_rad: float):
        """Return (row, column) of the nearest in-bounds map pixel."""
        ra0 = self.center_ra_rad - 0.5 * self.size_rad
        dec0 = self.center_dec_rad - 0.5 * self.size_rad
        col = int(np.rint((ra_rad - ra0) / self.pixel_size_rad))
        row = int(np.rint((dec_rad - dec0) / self.pixel_size_rad))
        if not (0 <= row < self.npix and 0 <= col < self.npix):
            return None
        return row, col

    def radial_distance(self, ra_rad: float, dec_rad: float):
        ra_grid, dec_grid = self.coordinate_grids()
        return np.hypot(ra_grid - ra_rad, dec_grid - dec_rad)
