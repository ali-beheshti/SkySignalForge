"""Circular halo-mask painter."""

from __future__ import annotations
import numpy as np
from tqdm.auto import tqdm

from .geometry import MapGeometry


def paint_circular_masks(
    ra_rad,
    dec_rad,
    radius_rad,
    geometry: MapGeometry,
    *,
    combine: str = "sum",
    show_progress: bool = True,
):
    """Paint per-object circular masks.

    ``combine='sum'`` counts overlapping disks; ``combine='union'`` returns a binary union mask.
    """
    ra = np.atleast_1d(np.asarray(ra_rad, dtype=float))
    dec = np.atleast_1d(np.asarray(dec_rad, dtype=float))
    radii = np.atleast_1d(np.asarray(radius_rad, dtype=float))
    if not (len(ra) == len(dec) == len(radii)):
        raise ValueError("positions and radii must have equal length")
    if np.any(radii < 0):
        raise ValueError("mask radii must be non-negative")
    if combine not in {"sum", "union"}:
        raise ValueError("combine must be 'sum' or 'union'")

    out = np.zeros((geometry.npix, geometry.npix), dtype=float)
    iterator = range(len(ra))
    if show_progress:
        iterator = tqdm(iterator, desc="Masks")

    for i in iterator:
        inside = geometry.radial_distance(ra[i], dec[i]) < radii[i]
        if combine == "sum":
            out[inside] += 1.0
        else:
            out[inside] = 1.0
    return out
