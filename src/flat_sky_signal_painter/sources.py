"""Point-like and catalog-impulse source painters."""

from __future__ import annotations
import numpy as np
from tqdm.auto import tqdm

from .geometry import MapGeometry


def paint_point_sources(
    ra_rad,
    dec_rad,
    geometry: MapGeometry,
    *,
    amplitudes=None,
    show_progress: bool = True,
):
    """Paint catalog objects into their nearest map pixels.

    Unit amplitudes produce an occupancy map; supplying amplitudes produces a weighted catalog layer.
    """
    ra = np.atleast_1d(np.asarray(ra_rad, dtype=float))
    dec = np.atleast_1d(np.asarray(dec_rad, dtype=float))
    if len(ra) != len(dec):
        raise ValueError("ra_rad and dec_rad must have equal length")

    if amplitudes is None:
        amp = np.ones(len(ra), dtype=float)
    else:
        amp = np.atleast_1d(np.asarray(amplitudes, dtype=float))
        if len(amp) != len(ra):
            raise ValueError("amplitudes must match source count")

    out = np.zeros((geometry.npix, geometry.npix), dtype=float)
    iterator = range(len(ra))
    if show_progress:
        iterator = tqdm(iterator, desc="Point sources")

    for i in iterator:
        index = geometry.nearest_pixel(ra[i], dec[i])
        if index is not None:
            out[index] += amp[i]
    return out


def paint_dust_test_sources(
    ra_rad,
    dec_rad,
    geometry: MapGeometry,
    *,
    amplitudes_uk=None,
    default_amplitude_uk: float = 4.25,
    show_progress: bool = True,
):
    """Paint the catalog-position dust/CIB test layer in microkelvin.

    This is a compact impulse-source layer for dust/CIB contamination tests. It is intentionally not a full CIB spectral model.
    """
    ra = np.atleast_1d(np.asarray(ra_rad, dtype=float))
    if amplitudes_uk is None:
        amplitudes_uk = np.full(len(ra), default_amplitude_uk)
    return paint_point_sources(
        ra,
        dec_rad,
        geometry,
        amplitudes=amplitudes_uk,
        show_progress=show_progress,
    )


def paint_localization_disks(
    ra_rad,
    dec_rad,
    geometry: MapGeometry,
    *,
    radius_arcmin: float = 0.25,
    show_progress: bool = True,
):
    """Paint unit disks around catalog locations for localization diagnostics."""
    ra = np.atleast_1d(np.asarray(ra_rad, dtype=float))
    dec = np.atleast_1d(np.asarray(dec_rad, dtype=float))
    radius_rad = np.deg2rad(radius_arcmin / 60.0)
    out = np.zeros((geometry.npix, geometry.npix), dtype=float)

    iterator = range(len(ra))
    if show_progress:
        iterator = tqdm(iterator, desc="Localization disks")

    for i in iterator:
        out[geometry.radial_distance(ra[i], dec[i]) < radius_rad] += 1.0
    return out
