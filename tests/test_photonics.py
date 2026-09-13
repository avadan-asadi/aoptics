import numpy as np
import pytest

from opticspy.photonics import waveguides as wg
from opticspy.photonics import fibers
from opticspy.photonics import resonators as res
from opticspy.photonics import nonlinear_optics as nl


# --------------------------- planar waveguides ---------------------------

def test_slab_waveguide_fundamental_mode_index_between_core_and_clad():
    wavelength = 1.0e-6
    n1, n2, d = 1.5, 1.45, 2.0e-6
    modes = wg.solve_te_modes(wavelength, d, n1, n2)
    assert len(modes) >= 1
    assert n2 < modes[0] < n1


def test_slab_waveguide_more_modes_for_thicker_core():
    wavelength = 1.0e-6
    n1, n2 = 1.5, 1.45
    n_thin = wg.number_of_te_modes(wavelength, 1.0e-6, n1, n2)
    n_thick = wg.number_of_te_modes(wavelength, 10.0e-6, n1, n2)
    assert n_thick > n_thin


def test_slab_waveguide_single_mode_thin_core():
    wavelength = 1.55e-6
    n1, n2, d = 1.46, 1.45, 1.0e-6
    modes = wg.solve_te_modes(wavelength, d, n1, n2)
    assert len(modes) == 1


def test_mode_confinement_factor_increases_with_v_number():
    wavelength = 1.0e-6
    n1, n2 = 1.5, 1.45
    k0 = 2 * np.pi / wavelength
    d_thin, d_thick = 0.3e-6, 3.0e-6
    modes_thin = wg.solve_te_modes(wavelength, d_thin, n1, n2)
    modes_thick = wg.solve_te_modes(wavelength, d_thick, n1, n2)
    gamma_thin = wg.mode_confinement_factor(modes_thin[0], k0, n1, n2, d_thin)
    gamma_thick = wg.mode_confinement_factor(modes_thick[0], k0, n1, n2, d_thick)
    assert 0 < gamma_thin < gamma_thick <= 1.0


# --------------------------- fibers ---------------------------

def test_v_number_single_mode_cutoff_consistency():
    wavelength = 1.55e-6
    n_core, n_clad = 1.4504, 1.4447  # typical SMF-28-like
    # find core radius giving V just below and above the 2.405 cutoff
    a_smf = 4.1e-6
    V = fibers.v_number(wavelength, a_smf, n_core, n_clad)
    assert fibers.is_single_mode(wavelength, a_smf, n_core, n_clad) == (V < 2.405)


def test_numerical_aperture_positive_and_consistent():
    NA = fibers.numerical_aperture(1.45, 1.44)
    assert NA > 0
    assert NA == pytest.approx(np.sqrt(1.45**2 - 1.44**2))


def test_mode_field_radius_reasonable_range():
    wavelength = 1.55e-6
    a = 4.1e-6
    n_core, n_clad = 1.4504, 1.4447
    w = fibers.mode_field_radius(wavelength, a, n_core, n_clad)
    assert a < w < 3 * a  # MFD is somewhat larger than the core for typical SMF


def test_sellmeier_silica_matches_known_index_near_1550nm():
    n = fibers.sellmeier_index_silica(1.55e-6)
    assert n == pytest.approx(1.444, abs=0.002)  # well-known fused-silica index near 1.55 um


def test_sellmeier_silica_matches_known_index_near_589nm():
    n = fibers.sellmeier_index_silica(589e-9)
    assert n == pytest.approx(1.458, abs=0.003)  # standard "fused silica ~1.458" figure at visible wavelengths


def test_zero_dispersion_wavelength_of_silica_near_1270nm():
    zeros = fibers.zero_dispersion_wavelength(fibers.sellmeier_index_silica, 1.0e-6, 1.4e-6)
    assert len(zeros) >= 1
    assert any(1.2e-6 < z < 1.35e-6 for z in zeros)  # silica ZDW is ~1.27-1.32 um


def test_group_index_greater_than_phase_index_normal_dispersion():
    wavelength = 700e-9  # below the ZDW -> normal dispersion region -> n_g > n
    n = fibers.sellmeier_index_silica(wavelength)
    ng = fibers.group_index(fibers.sellmeier_index_silica, wavelength)
    assert ng > n


# --------------------------- resonators ---------------------------

def test_finesse_high_reflectivity_mirrors():
    F = res.finesse(0.99, 0.99)
    assert F > 100


def test_free_spectral_range_known_value():
    L = 0.15  # m -> FSR = c/(2L) = 1 GHz
    fsr = res.free_spectral_range(L)
    assert fsr == pytest.approx(299_792_458.0 / 0.3, rel=1e-9)


def test_airy_transmission_peak_is_unity_at_resonance():
    T = res.airy_transmission(0.0, 0.9, 0.9)
    assert T == pytest.approx(1.0, rel=1e-9)


def test_airy_transmission_bounded():
    delta = np.linspace(0, 4 * np.pi, 500)
    T = res.airy_transmission(delta, 0.8, 0.85)
    assert np.all(T <= 1.0 + 1e-9) and np.all(T >= 0)


def test_symmetric_confocal_stability():
    L = 0.5
    R = L  # confocal
    assert res.is_stable(L, R, R)


def test_plane_parallel_cavity_is_marginally_stable_boundary():
    L = 0.5
    g1, g2 = res.g_parameters(L, np.inf, np.inf)
    assert g1 == pytest.approx(1.0) and g2 == pytest.approx(1.0)
    assert res.is_stable(L, np.inf, np.inf)  # g1*g2 = 1, boundary of stability


def test_hemispherical_unstable_cavity_detected():
    L = 2.0
    R1, R2 = 0.5, np.inf  # L > R1 with a flat mirror -> unstable
    assert not res.is_stable(L, R1, R2)


def test_confocal_resonator_waist_matches_closed_form():
    # NOTE: an EXACTLY confocal cavity (R1=R2=L) has a degenerate round-trip
    # matrix (M = -Identity, B=0 at every reference plane), which is a
    # genuine, well-known singularity of the single-round-trip ABCD
    # self-consistency law for this special geometry. We test the
    # near-confocal limit instead, where the general machinery is regular
    # and should converge to the confocal closed-form result.
    L = 0.5
    wavelength = 1.0e-6
    R = L * (1 + 1e-4)
    w_mirror, q, stable = res.resonator_waist_at_mirror1(L, R, R, wavelength)
    assert stable
    expected_mirror = res.symmetric_confocal_waist_at_mirror(L, wavelength)
    assert w_mirror == pytest.approx(expected_mirror, rel=2e-3)


def test_confocal_resonator_mirror_spot_matches_closed_form():
    L = 0.5
    wavelength = 1.0e-6
    R = L * (1 + 1e-4)
    w_mirror, q, stable = res.resonator_waist_at_mirror1(L, R, R, wavelength)
    assert stable
    expected = res.symmetric_confocal_waist_at_mirror(L, wavelength)
    assert w_mirror == pytest.approx(expected, rel=2e-3)


def test_near_planar_symmetric_resonator_is_stable_and_has_finite_waist():
    L = 0.3
    R = 10.0  # nearly planar but still curved (g close to but below 1)
    wavelength = 1.0e-6
    assert res.is_stable(L, R, R)
    w_center, stable = res.resonator_waist_at_center(L, R, R, wavelength)
    assert stable
    assert np.isfinite(w_center) and w_center > 0


# --------------------------- nonlinear optics ---------------------------

def test_shg_phase_matched_efficiency_is_unity():
    eta = nl.shg_phase_matching_efficiency(0.0, 0.01)
    assert eta == pytest.approx(1.0)


def test_shg_efficiency_drops_at_coherence_length():
    delta_k = 1000.0
    Lcoh = nl.coherence_length_shg(delta_k)
    eta = nl.shg_phase_matching_efficiency(delta_k, 2 * Lcoh)
    assert eta == pytest.approx(0.0, abs=1e-9)  # first null of the sinc^2 at L = Lcoh... check 2x


def test_shg_efficiency_array_matches_scalar():
    dk = np.linspace(-500, 500, 21)
    L = 0.001
    arr = nl.shg_phase_matching_efficiency_array(dk, L)
    scalars = np.array([nl.shg_phase_matching_efficiency(d, L) for d in dk])
    assert np.allclose(arr, scalars)


def test_nonlinear_length_and_phase_consistency():
    gamma, P = 1.5, 0.2  # 1/(W*m), W
    L_NL = nl.nonlinear_length(gamma, P)
    phi = nl.nonlinear_phase_spm(gamma, P, L_NL)
    assert phi == pytest.approx(1.0, rel=1e-9)  # by definition, phase = 1 rad at L = L_NL


def test_soliton_order_fundamental():
    # choose parameters so N=1 exactly: gamma*P*T0^2 = |beta2|
    beta2 = -20e-27  # s^2/m (typical anomalous dispersion)
    T0 = 1e-12
    gamma = 1.5e-3  # 1/(W*m)
    P = abs(beta2) / (gamma * T0 ** 2)
    N = nl.soliton_order(gamma, P, T0, beta2)
    assert N == pytest.approx(1.0, rel=1e-6)


def test_fwm_phase_mismatch_zero_for_degenerate_case():
    k = 1.2345
    assert nl.four_wave_mixing_phase_mismatch(k, k, k, k) == pytest.approx(0.0)
