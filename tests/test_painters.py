import numpy as np

from flat_sky_signal_painter import (
    MapGeometry,
    MovingLensCatalog,
    paint_circular_masks,
    paint_moving_lens,
    paint_point_sources,
    paint_tsz_gaussian,
)


def ml_catalog(vtheta=400.0, vphi=200.0):
    return MovingLensCatalog(
        ra_rad=[0.0],
        dec_rad=[0.0],
        redshift=[0.35],
        concentration=[4.5],
        rs_mpc=[0.32],
        rho_s_msun_mpc3=[1.1e15],
        comoving_distance_mpc=[1200.0],
        v_theta_kms=[vtheta],
        v_phi_kms=[vphi],
    )


def test_zero_velocity_zero_moving_lens():
    geo = MapGeometry(npix=48, size_deg=0.6)
    m = paint_moving_lens(ml_catalog(0.0, 0.0), geo, show_progress=False)
    assert np.allclose(m, 0.0)


def test_moving_lens_velocity_linearity():
    geo = MapGeometry(npix=48, size_deg=0.6)
    m1 = paint_moving_lens(ml_catalog(200.0, 100.0), geo, show_progress=False)
    m2 = paint_moving_lens(ml_catalog(400.0, 200.0), geo, show_progress=False)
    assert np.allclose(m2, 2.0 * m1, rtol=1e-10, atol=1e-14)


def test_point_source_total_amplitude():
    geo = MapGeometry(npix=64, size_deg=1.0)
    amps = np.array([1.0, 2.0, 3.0])
    m = paint_point_sources(
        np.deg2rad([-0.2, 0.0, 0.2]),
        np.deg2rad([0.0, 0.2, -0.2]),
        geo,
        amplitudes=amps,
        show_progress=False,
    )
    assert np.isclose(m.sum(), amps.sum())


def test_union_mask_is_binary():
    geo = MapGeometry(npix=64, size_deg=1.0)
    m = paint_circular_masks(
        np.deg2rad([0.0, 0.05]),
        np.deg2rad([0.0, 0.0]),
        np.deg2rad([0.12, 0.12]),
        geo,
        combine="union",
        show_progress=False,
    )
    assert set(np.unique(m)).issubset({0.0, 1.0})


def test_tsz_shape_and_sign_at_150():
    geo = MapGeometry(npix=64, size_deg=1.0)
    m = paint_tsz_gaussian(
        [0.0], [0.0], [1.0e-6], geo,
        frequency_ghz=150.0,
        show_progress=False,
    )
    assert m.shape == (64, 64)
    assert np.min(m) < 0.0
