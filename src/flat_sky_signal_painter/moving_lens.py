"""Moving-lens temperature painter."""

from __future__ import annotations
import numpy as np
from tqdm.auto import tqdm

from .catalog import MovingLensCatalog
from .constants import C_KMS, TCMB_UK
from .geometry import MapGeometry
from .profiles import nfw_deflection_angle


def single_moving_lens_template(halo: dict[str, float], geometry: MapGeometry):
    """Paint one halo's moving-lens temperature dipole in microkelvin."""
    ra_grid, dec_grid = geometry.coordinate_grids()
    dx_dec = dec_grid - halo["dec_rad"]
    dx_ra = ra_grid - halo["ra_rad"]

    radius_rad = np.hypot(dx_ra, dx_dec)
    polar_angle = np.arctan2(dx_dec, dx_ra)
    r_comoving_mpc = radius_rad * halo["comoving_distance_mpc"]
    x_over_rs = r_comoving_mpc / halo["rs_mpc"]

    beta = nfw_deflection_angle(
        x_over_rs,
        concentration=halo["concentration"],
        rho_s_msun_mpc3=halo["rho_s_msun_mpc3"],
        rs_mpc=halo["rs_mpc"],
        scale_factor=1.0 / (1.0 + halo["redshift"]),
    )

    v_theta = halo["v_theta_kms"]
    v_phi = halo["v_phi_kms"]
    speed = np.hypot(v_theta, v_phi)
    if speed == 0:
        return np.zeros_like(beta)

    velocity_angle = np.arctan2(v_theta, v_phi)
    return beta * (speed / C_KMS) * np.cos(velocity_angle - polar_angle) * TCMB_UK


def paint_moving_lens(catalog: MovingLensCatalog, geometry: MapGeometry, *, show_progress=True):
    """Paint the summed moving-lens signal from a halo catalog."""
    out = np.zeros((geometry.npix, geometry.npix), dtype=float)
    iterator = range(len(catalog))
    if show_progress:
        iterator = tqdm(iterator, desc="Moving lens")
    for i in iterator:
        out += single_moving_lens_template(catalog.row(i), geometry)
    return out
