import numpy as np
import pytest

from opticspy.modulation import electrooptics as eo, acoustooptics as ao


# --------------------------- electro-optics ---------------------------

def test_half_wave_voltage_gives_pi_phase_shift():
    n0, r, L, d, wl = 2.2, 30.8e-12, 0.02, 5e-6, 1550e-9  # LiNbO3-like
    V_pi = eo.half_wave_voltage(n0, r, L, d, wl)
    phase = eo.phase_shift_pockels(n0, r, V_pi, d, L, wl)
    assert abs(phase) == pytest.approx(np.pi, rel=1e-9)


def test_mach_zehnder_transmission_at_zero_and_vpi():
    V_pi = 5.0
    T0 = eo.mach_zehnder_eo_transmission(0.0, V_pi, bias_phase=0.0)
    T_pi = eo.mach_zehnder_eo_transmission(V_pi, V_pi, bias_phase=0.0)
    assert T0 == pytest.approx(1.0, rel=1e-9)
    assert T_pi == pytest.approx(0.0, abs=1e-9)


def test_index_change_sign_convention_consistent():
    n0, r = 2.2, 30.8e-12
    dn_pos = eo.index_change_pockels(n0, r, 1e6)
    dn_neg = eo.index_change_pockels(n0, r, -1e6)
    assert dn_pos == pytest.approx(-dn_neg)


def test_eo_bandwidth_infinite_when_velocity_matched():
    bw = eo.eo_modulator_bandwidth_limit(0.01, n_microwave=2.2, n_optical=2.2)
    assert np.isinf(bw)


# --------------------------- acousto-optics ---------------------------

def test_klein_cook_regime_classification():
    Q_thin = ao.klein_cook_parameter(wavelength=633e-9, interaction_length=1e-3,
                                      acoustic_wavelength=200e-6, n=1.5)
    Q_thick = ao.klein_cook_parameter(wavelength=633e-9, interaction_length=5e-2,
                                       acoustic_wavelength=20e-6, n=1.5)
    assert ao.regime(Q_thin) in ("Raman-Nath", "intermediate")
    assert ao.regime(Q_thick) == "Bragg"


def test_bragg_angle_small_for_typical_aom():
    wl, lam_ac = 633e-9, 200e-6
    theta_B = ao.bragg_angle(wl, lam_ac, n=1.5)
    assert 0 < theta_B < np.radians(5)


def test_acoustic_wavelength_relation():
    v, f = 4000.0, 100e6
    lam = ao.acoustic_wavelength(v, f)
    assert lam == pytest.approx(v / f)


def test_diffracted_frequency_shift_up_and_down():
    f0, fa = 4.74e14, 100e6
    f_up = ao.diffracted_order_frequency_shift(f0, fa, order=1)
    f_down = ao.diffracted_order_frequency_shift(f0, fa, order=-1)
    assert f_up == pytest.approx(f0 + fa)
    assert f_down == pytest.approx(f0 - fa)


def test_normalized_bragg_efficiency_full_transfer_at_pi_over_2():
    eta = ao.normalized_bragg_efficiency(np.pi / 2)
    assert eta == pytest.approx(1.0, rel=1e-9)


def test_raman_nath_efficiency_zero_order_at_zero_parameter():
    eta0 = ao.raman_nath_diffraction_efficiency(0, 0.0)
    eta1 = ao.raman_nath_diffraction_efficiency(1, 0.0)
    assert eta0 == pytest.approx(1.0)
    assert eta1 == pytest.approx(0.0)
