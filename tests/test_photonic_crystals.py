import numpy as np
import pytest

from opticspy.photonics import photonic_crystals as pc


def test_no_layers_reduces_to_single_interface_fresnel():
    from opticspy.electromagnetic import fresnel
    n_inc, n_sub = 1.0, 1.5
    result = pc.multilayer_reflectance([], 500e-9, n_incident=n_inc, n_substrate=n_sub)
    expected_R = fresnel.normal_incidence_reflectance(n_inc, n_sub)
    assert result['R'] == pytest.approx(expected_R, rel=1e-6)


def test_energy_conservation_lossless_stack():
    layers = [(2.3, 100e-9), (1.45, 150e-9), (2.3, 100e-9)]
    result = pc.multilayer_reflectance(layers, 550e-9, n_incident=1.0, n_substrate=1.5)
    assert result['R'] + result['T'] == pytest.approx(1.0, rel=1e-6)


def test_quarter_wave_stack_high_reflectance_at_design_wavelength():
    n_high, n_low = 2.3, 1.45
    wl0 = 550e-9
    layers = pc.quarter_wave_stack(n_high, n_low, wl0, n_periods=15)
    result = pc.multilayer_reflectance(layers, wl0, n_incident=1.0, n_substrate=1.5)
    assert result['R'] > 0.99


def test_quarter_wave_stack_reflectance_increases_with_periods():
    n_high, n_low, wl0 = 2.3, 1.45, 550e-9
    Rs = []
    for N in [1, 5, 15]:
        layers = pc.quarter_wave_stack(n_high, n_low, wl0, n_periods=N)
        Rs.append(pc.multilayer_reflectance(layers, wl0, n_substrate=1.5)['R'])
    assert Rs[0] < Rs[1] < Rs[2]


def test_reflectance_spectrum_peaks_near_design_wavelength():
    n_high, n_low, wl0 = 2.3, 1.45, 550e-9
    layers = pc.quarter_wave_stack(n_high, n_low, wl0, n_periods=10)
    wavelengths = np.linspace(400e-9, 750e-9, 300)
    R = pc.reflectance_spectrum(layers, wavelengths, n_substrate=1.5)
    peak_wl = wavelengths[np.argmax(R)]
    assert peak_wl == pytest.approx(wl0, rel=0.05)


def test_stopband_edges_bracket_high_reflectance_region():
    n_high, n_low, wl0 = 2.3, 1.45, 550e-9
    lo, hi = pc.bragg_stopband_edges(n_high, n_low, wl0)
    assert lo < wl0 < hi
    layers = pc.quarter_wave_stack(n_high, n_low, wl0, n_periods=20)
    R_mid = pc.multilayer_reflectance(layers, wl0, n_substrate=1.5)['R']
    R_outside = pc.multilayer_reflectance(layers, hi * 1.15, n_substrate=1.5)['R']
    assert R_mid > R_outside
