import numpy as np
import pytest

from opticspy.fourier_optics import transfer_functions as tf, spatial_filtering as sf


def test_pupil_function_is_binary_and_centered():
    pupil = tf.pupil_function_circular(128, 1.0, 0.3)
    assert set(np.unique(np.real(pupil))) <= {0.0, 1.0}
    center = pupil[64, 64]
    assert np.real(center) == 1.0


def test_incoherent_psf_peak_at_center():
    pupil = tf.pupil_function_circular(128, 1.0, 0.3)
    psf = tf.incoherent_psf(pupil)
    peak_idx = np.unravel_index(np.argmax(psf), psf.shape)
    assert peak_idx == (64, 64)


def test_otf_at_zero_frequency_is_one():
    pupil = tf.pupil_function_circular(64, 1.0, 0.3)
    otf = tf.optical_transfer_function(pupil)
    center = otf[32, 32]
    assert np.abs(center) == pytest.approx(1.0, rel=1e-6)


def test_incoherent_cutoff_twice_coherent_cutoff():
    NA, wl = 0.5, 550e-9
    fc_coh = tf.diffraction_limited_cutoff_coherent(NA, wl)
    fc_incoh = tf.diffraction_limited_cutoff_incoherent(NA, wl)
    assert fc_incoh == pytest.approx(2 * fc_coh)


def test_rayleigh_and_abbe_resolution_same_order():
    NA, wl = 0.9, 500e-9
    r_rayleigh = tf.rayleigh_resolution(wl, NA)
    r_abbe = tf.abbe_resolution(wl, NA)
    assert 0.5 * r_abbe < r_rayleigh < 2 * r_abbe


def test_four_f_system_identity_with_no_filter():
    rng = np.random.default_rng(0)
    field = rng.normal(size=(32, 32)) + 1j * rng.normal(size=(32, 32))
    out, _ = sf.four_f_system(field, filter_mask=None)
    assert np.allclose(out, field, atol=1e-9)


def test_four_f_low_pass_reduces_high_frequency_content():
    N = 64
    x = np.linspace(-1, 1, N)
    X, Y = np.meshgrid(x, x)
    # a high-frequency checkerboard-like pattern
    field = np.sign(np.sin(40 * np.pi * X)).astype(complex)
    mask = sf.low_pass_filter_mask(N, cutoff_fraction=0.05)
    out, _ = sf.four_f_system(field, filter_mask=mask)
    # low-pass filtering a fast oscillation should strongly reduce its variance
    assert np.var(np.real(out)) < 0.5 * np.var(np.real(field))


def test_high_pass_mask_is_complement_of_low_pass():
    N = 32
    lp = sf.low_pass_filter_mask(N, 0.3)
    hp = sf.high_pass_filter_mask(N, 0.3)
    assert np.allclose(lp + hp, np.ones((N, N)))
