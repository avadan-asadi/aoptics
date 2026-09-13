import numpy as np
import pytest

from opticspy.nonlinear_optics_boyd import two_level_atom as tla
from opticspy.nonlinear_optics_boyd import self_action as sa
from opticspy.nonlinear_optics_boyd import stimulated_scattering as ss
from opticspy.nonlinear_optics_boyd import multiphoton_absorption as tpa


# --------------------------- two-level atom ---------------------------

def test_absorption_unsaturated_limit():
    alpha0 = 10.0
    alpha = tla.saturated_absorption_coefficient(alpha0, intensity=0.0, saturation_intensity=1.0)
    assert alpha == pytest.approx(alpha0)


def test_absorption_halved_at_saturation_intensity_on_resonance():
    alpha0 = 10.0
    alpha = tla.saturated_absorption_coefficient(alpha0, intensity=1.0, saturation_intensity=1.0,
                                                  detuning=0.0, linewidth=1.0)
    assert alpha == pytest.approx(alpha0 / 2)


def test_power_broadening_increases_linewidth():
    g0 = 1.0
    g = tla.power_broadened_linewidth(g0, intensity=3.0, saturation_intensity=1.0)
    assert g == pytest.approx(g0 * 2.0)


def test_inversion_limits():
    w_low = tla.steady_state_population_inversion(rabi_frequency=1e-6, T1=1.0, T2=1.0, detuning=0.0)
    w_high = tla.steady_state_population_inversion(rabi_frequency=1e6, T1=1.0, T2=1.0, detuning=0.0)
    assert w_low == pytest.approx(-1.0, abs=1e-4)
    assert w_high == pytest.approx(0.0, abs=1e-4)


def test_susceptibility_shape_zero_detuning_pure_absorption():
    disp, absn = tla.two_level_susceptibility_shape(detuning=0.0, linewidth=1.0)
    assert disp == pytest.approx(0.0)
    assert absn == pytest.approx(1.0)


# --------------------------- self-action ---------------------------

def test_critical_power_matches_marburger_formula():
    wl, n0, n2 = 800e-9, 1.33, 4.1e-16 * 1e-4  # water n2 in m^2/W (from cm^2/W)
    Pcr = sa.critical_power_self_focusing(wl, n0, n2)
    expected = 3.77 * wl ** 2 / (8 * np.pi * n0 * n2)
    assert Pcr == pytest.approx(expected)
    assert Pcr > 0


def test_is_self_focusing_threshold():
    Pcr = 1e6
    assert sa.is_self_focusing(2e6, Pcr)
    assert not sa.is_self_focusing(0.5e6, Pcr)


def test_collapse_distance_infinite_below_threshold():
    z = sa.self_focusing_collapse_distance(power=0.5e6, critical_power=1e6, rayleigh_range=1.0)
    assert np.isinf(z)


def test_collapse_distance_finite_and_positive_above_threshold():
    z = sa.self_focusing_collapse_distance(power=5e6, critical_power=1e6, rayleigh_range=0.69)
    assert 0 < z < np.inf


def test_collapse_distance_known_example():
    # From the literature example: air at 1030 nm, Pcr=4.9 GW, zR=690 m
    Pcr = 4.9e9
    zR = 690.0
    P = 10 * Pcr
    z = sa.self_focusing_collapse_distance(P, Pcr, zR)
    # sanity: should be a positive distance much less than zR for P >> Pcr
    assert 0 < z < zR


def test_nonlinear_phase_scales_linearly_with_length():
    phi1 = sa.nonlinear_phase_kerr(n2=1e-20, intensity=1e13, length=1e-3, wavelength=800e-9)
    phi2 = sa.nonlinear_phase_kerr(n2=1e-20, intensity=1e13, length=2e-3, wavelength=800e-9)
    assert phi2 == pytest.approx(2 * phi1)


# --------------------------- stimulated scattering ---------------------------

def test_effective_length_low_loss_limit():
    L = 1.0
    Leff = ss.effective_length(L, alpha=1e-9)
    assert Leff == pytest.approx(L, rel=1e-5)


def test_stokes_below_pump_anti_stokes_above():
    wp, wv = 2.0e15, 4.0e13
    wS = ss.stokes_frequency_shift(wp, wv)
    wAS = ss.anti_stokes_frequency_shift(wp, wv)
    assert wS < wp < wAS


def test_srs_threshold_scales_inversely_with_gain():
    Pth1 = ss.srs_threshold_power(gain_coefficient=1e-13, effective_area=50e-12, effective_length_value=1000)
    Pth2 = ss.srs_threshold_power(gain_coefficient=2e-13, effective_area=50e-12, effective_length_value=1000)
    assert Pth2 == pytest.approx(Pth1 / 2)


def test_sbs_threshold_reasonable_order_of_magnitude_for_fiber():
    # matches the classic ~mW-level SBS threshold for a long fiber span
    Pth = ss.sbs_threshold_power(gain_coefficient=5e-11, effective_area=46.5e-12,
                                  effective_length_value=7846)
    assert 1e-3 < Pth < 1.0


def test_sbs_threshold_increases_with_pump_linewidth():
    base = ss.sbs_threshold_power(1e-11, 50e-12, 1000)
    broadened = ss.sbs_threshold_power(1e-11, 50e-12, 1000, pump_linewidth=50e6, brillouin_linewidth=25e6)
    assert broadened > base


def test_brillouin_backscatter_shift_positive():
    nu_B = ss.brillouin_frequency_shift(refractive_index=1.45, acoustic_velocity=5960.0,
                                         wavelength=1550e-9, theta=np.pi)
    assert nu_B > 0
    # typical silica-fiber Brillouin shift is ~10-11 GHz at 1550 nm
    assert 8e9 < nu_B < 13e9


# --------------------------- multiphoton absorption ---------------------------

def test_tpa_reduces_transmission_below_linear():
    I_out_linear = 1.0 * np.exp(-0.0)  # no linear loss
    I_out_tpa = tpa.intensity_after_tpa(intensity_in=1.0, beta_tpa=0.5, length=1.0, alpha_linear=0.0)
    assert I_out_tpa < I_out_linear


def test_tpa_zero_beta_reduces_to_beer_lambert():
    alpha = 0.3
    I_out = tpa.intensity_after_tpa(intensity_in=2.0, beta_tpa=0.0, length=1.5, alpha_linear=alpha)
    expected = 2.0 * np.exp(-alpha * 1.5)
    assert I_out == pytest.approx(expected, rel=1e-6)


def test_tpa_transmission_decreases_with_intensity():
    T_low, _ = tpa.nonlinear_transmission_tpa(intensity_in=0.1, beta_tpa=1.0, length=1.0)
    T_high, _ = tpa.nonlinear_transmission_tpa(intensity_in=10.0, beta_tpa=1.0, length=1.0)
    assert T_high < T_low


def test_saturable_absorber_transmission_increases_with_intensity():
    T_low = tpa.saturable_absorber_transmission(0.01, alpha0=0.5, saturation_intensity=1.0)
    T_high = tpa.saturable_absorber_transmission(100.0, alpha0=0.5, saturation_intensity=1.0)
    assert T_high > T_low
    assert T_high == pytest.approx(1.0, abs=0.01)
