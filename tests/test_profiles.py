import numpy as np

from flat_sky_signal_painter.profiles import (
    truncated_nfw_profile,
    tsz_spectral_factor,
)


def test_nfw_profile_is_finite():
    x = np.geomspace(0.02, 20.0, 400)
    y = truncated_nfw_profile(x, 4.5)
    assert y.shape == x.shape
    assert np.all(np.isfinite(y))


def test_nfw_inner_cutoff():
    y = truncated_nfw_profile(np.array([0.0, 0.005]), 4.5)
    assert np.all(y == 0.0)


def test_tsz_is_negative_at_150_ghz():
    assert tsz_spectral_factor(150.0) < 0.0


def test_tsz_null_is_near_217_ghz():
    assert abs(float(tsz_spectral_factor(217.0))) < 0.05
