"""Analytic profiles used by the sky painters."""

from __future__ import annotations
import numpy as np

from .constants import (
    BOLTZMANN_J_K,
    G_MPC_KMS2_MSUN,
    PLANCK_J_HZ,
    TCMB_K,
)


def truncated_nfw_profile(x_over_rs, concentration: float, *, inner_cutoff=1.0e-2):
    """Dimensionless projected deflection profile of a truncated NFW halo."""
    x = np.asarray(x_over_rs, dtype=np.complex128)
    c = float(concentration)
    if c <= 0:
        raise ValueError("concentration must be positive")

    xr = x.real
    # Avoid the removable numerical singularity at exactly x=1.
    exact_one = np.isclose(xr, 1.0, rtol=0.0, atol=1e-12)
    if np.any(exact_one):
        x = x.copy()
        x[exact_one] = 1.0 + 1e-8 + 0j
        xr = x.real

    inside = xr < c
    resolved = xr > inner_cutoff
    unresolved = xr < inner_cutoff
    out = np.zeros(x.shape, dtype=np.complex128)

    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        m = inside & resolved
        if np.any(m):
            xm = x[m]
            sqc = np.sqrt(c**2 - xm**2)
            sq1 = np.sqrt(1 - xm**2)
            out[m] = (
                -xm / (1 + c) / (c + sqc)
                + (np.log(1 + c) + np.log(xm) - np.log(c + sqc)) / xm
                - np.log(-c - xm**2 + sq1 * sqc) / xm / sq1
                + (
                    np.log(1 + c)
                    + np.log(-xm)
                    - 2j * np.log(-xm).imag
                ) / xm / sq1
            )

        out[inside & unresolved] = 0.0 + 0.0j

        m = (~inside) & resolved
        if np.any(m):
            xm = x[m]
            out[m] = ((-c / (1 + c)) + np.log(1 + c)) / xm

    return out.real


def nfw_deflection_angle(
    x_over_rs,
    *,
    concentration: float,
    rho_s_msun_mpc3: float,
    rs_mpc: float,
    scale_factor: float,
    inner_cutoff=1.0e-2,
):
    """Dimensionless deflection angle beta for the truncated NFW profile."""
    profile = truncated_nfw_profile(x_over_rs, concentration, inner_cutoff=inner_cutoff)
    c_kms = 299_792.458
    prefactor = (
        16.0 * np.pi * G_MPC_KMS2_MSUN * rho_s_msun_mpc3 * rs_mpc**2
        / (scale_factor * c_kms**2)
    )
    return prefactor * profile


def tsz_spectral_factor(frequency_ghz, *, tcmb_k: float = TCMB_K):
    """Non-relativistic tSZ thermodynamic-temperature factor g(x).

    g(x) = x coth(x/2) - 4
    """
    nu_hz = np.asarray(frequency_ghz, dtype=float) * 1e9
    x = PLANCK_J_HZ * nu_hz / (BOLTZMANN_J_K * tcmb_k)
    return x / np.tanh(x / 2.0) - 4.0


def gaussian_profile(radius_rad, sigma_rad):
    """Unit-amplitude circular Gaussian profile."""
    if sigma_rad <= 0:
        raise ValueError("sigma_rad must be positive")
    r = np.asarray(radius_rad, dtype=float)
    return np.exp(-0.5 * (r / sigma_rad) ** 2)
