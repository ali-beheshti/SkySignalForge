"""Gaussian thermal-SZ painter."""

from __future__ import annotations
import numpy as np
from tqdm.auto import tqdm

from .constants import TCMB_UK
from .geometry import MapGeometry
from .profiles import gaussian_profile, tsz_spectral_factor


def paint_tsz_gaussian(
    ra_rad,
    dec_rad,
    profile_amplitude,
    geometry: MapGeometry,
    *,
    frequency_ghz: float = 150.0,
    intrinsic_fwhm_arcmin: float = 5.83,
    amplitude_scale: float = 1.328,
    show_progress: bool = True,
):
    """Paint Gaussian tSZ temperature profiles.

    Parameters
    ----------
    profile_amplitude
        Dimensionless per-object profile amplitudes. This corresponds to the
        Per-object ``tszGaussAmp`` profile amplitude.
    intrinsic_fwhm_arcmin
        Gaussian FWHM before any external beam treatment. Default: 5.83 arcmin.
    amplitude_scale
        Multiplicative calibration factor. Default: 1.328.

    Returns
    -------
    numpy.ndarray
        Thermodynamic temperature map in microkelvin.
    """
    ra = np.atleast_1d(np.asarray(ra_rad, dtype=float))
    dec = np.atleast_1d(np.asarray(dec_rad, dtype=float))
    amp = np.atleast_1d(np.asarray(profile_amplitude, dtype=float))
    if not (len(ra) == len(dec) == len(amp)):
        raise ValueError("ra_rad, dec_rad, and profile_amplitude must have equal length")

    sigma_rad = np.deg2rad(intrinsic_fwhm_arcmin / 60.0) / np.sqrt(8.0 * np.log(2.0))
    freq_factor = float(tsz_spectral_factor(frequency_ghz))
    out = np.zeros((geometry.npix, geometry.npix), dtype=float)

    iterator = range(len(ra))
    if show_progress:
        iterator = tqdm(iterator, desc="tSZ")

    for i in iterator:
        radius = geometry.radial_distance(ra[i], dec[i])
        out += (
            amplitude_scale
            * amp[i]
            * gaussian_profile(radius, sigma_rad)
            * TCMB_UK
            * freq_factor
        )
    return out
