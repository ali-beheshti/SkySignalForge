"""Gaussian kinetic-SZ painter."""

from __future__ import annotations

import numpy as np
from tqdm.auto import tqdm

from .constants import C_KMS, TCMB_UK
from .geometry import MapGeometry
from .profiles import gaussian_profile


def paint_ksz_gaussian(
    ra_rad,
    dec_rad,
    tau_amplitude,
    v_los_kms,
    geometry: MapGeometry,
    *,
    intrinsic_fwhm_arcmin: float = 5.83,
    show_progress: bool = True,
):
    """Paint Gaussian kinetic-SZ temperature profiles.

    The thermodynamic-temperature shift is

        Delta T / T_CMB = -tau * v_los / c,

    with positive ``v_los_kms`` defined as motion away from the observer.

    Parameters
    ----------
    ra_rad, dec_rad
        Angular source positions in radians.
    tau_amplitude
        Dimensionless peak optical-depth amplitude for each object.
    v_los_kms
        Line-of-sight peculiar velocity in km/s.
    geometry
        Output flat-sky map geometry.
    intrinsic_fwhm_arcmin
        Gaussian FWHM used for the optical-depth profile.

    Returns
    -------
    numpy.ndarray
        kSZ temperature map in microkelvin.
    """
    ra = np.atleast_1d(np.asarray(ra_rad, dtype=float))
    dec = np.atleast_1d(np.asarray(dec_rad, dtype=float))
    tau0 = np.atleast_1d(np.asarray(tau_amplitude, dtype=float))
    vlos = np.atleast_1d(np.asarray(v_los_kms, dtype=float))

    if not (len(ra) == len(dec) == len(tau0) == len(vlos)):
        raise ValueError("positions, tau_amplitude, and v_los_kms must have equal length")
    if np.any(tau0 < 0):
        raise ValueError("tau_amplitude must be non-negative")

    sigma_rad = np.deg2rad(intrinsic_fwhm_arcmin / 60.0) / np.sqrt(8.0 * np.log(2.0))
    out = np.zeros((geometry.npix, geometry.npix), dtype=float)

    iterator = range(len(ra))
    if show_progress:
        iterator = tqdm(iterator, desc="kSZ")

    for i in iterator:
        radius = geometry.radial_distance(ra[i], dec[i])
        tau = tau0[i] * gaussian_profile(radius, sigma_rad)
        out += -TCMB_UK * tau * (vlos[i] / C_KMS)

    return out
