import numpy as np
import pytest

from aoptics.nonlinear_optics_boyd import susceptibility as sus
from aoptics.nonlinear_optics_boyd import coupled_wave_mixing as cwm


# --------------------------- susceptibility ---------------------------

def test_d_coefficient_round_trip():
    chi2 = 4.5e-12
    d = sus.d_coefficient_from_chi2(chi2)
    assert sus.chi2_from_d_coefficient(d) == pytest.approx(chi2)


def test_chi1_from_index_matches_definition():
    n = 1.5
    assert sus.chi1_from_refractive_index(n) == pytest.approx(n ** 2 - 1)


def test_millers_rule_round_trip():
    chi1_1, chi1_2, chi1_3 = 1.2, 1.3, 1.25
    chi2 = 5e-12
    delta = sus.millers_delta(chi2, chi1_1, chi1_2, chi1_3)
    chi2_est = sus.estimate_chi2_miller_rule(delta, chi1_1, chi1_2, chi1_3)
    assert chi2_est == pytest.approx(chi2, rel=1e-9)


def test_symmetric_tensor_detection():
    # a fully (Kleinman) symmetric synthetic tensor
    base = np.random.default_rng(0).normal(size=(3, 3, 3))
    sym = base + base.transpose(1, 2, 0) + base.transpose(2, 0, 1) \
        + base.transpose(0, 2, 1) + base.transpose(1, 0, 2) + base.transpose(2, 1, 0)
    assert sus.is_symmetric_chi2_tensor(sym)


def test_asymmetric_tensor_detected_as_not_symmetric():
    rng = np.random.default_rng(1)
    t = rng.normal(size=(3, 3, 3))
    assert not sus.is_symmetric_chi2_tensor(t)


def test_total_polarization_matches_manual_expansion():
    E = 1e6
    chi1, chi2, chi3 = 2.0, 1e-11, 1e-20
    P = sus.total_polarization(E, chi1, chi2, chi3)
    expected = sus.EPSILON_0 * (chi1 * E + chi2 * E ** 2 + chi3 * E ** 3)
    assert P == pytest.approx(expected)


# --------------------------- coupled-wave mixing (SFG/OPA) ---------------------------

def test_sfg_energy_manley_rowe_conserved():
    omega1, omega2 = 1.5e15, 1.8e15
    omega3 = omega1 + omega2
    n1 = n2 = n3 = 1.5
    deff = 1e-12
    z, A1, A2, A3 = cwm.simulate_three_wave_mixing(
        A1_0=1.0, A2_0=1.0, A3_0=0.0, z_max=0.1,
        omega1=omega1, omega2=omega2, omega3=omega3,
        n1=n1, n2=n2, n3=n3, deff=deff, delta_k=0.0, n_steps=1000)
    f1, f2, f3 = cwm.manley_rowe_photon_fluxes(A1, A2, A3, omega1, omega2, omega3)
    # total photon number-like invariant f1+f3 and f2+f3 should be conserved
    inv13 = f1 + f3
    inv23 = f2 + f3
    assert np.allclose(inv13, inv13[0], rtol=1e-6)
    assert np.allclose(inv23, inv23[0], rtol=1e-6)


def test_sfg_generates_sum_frequency_from_two_inputs():
    omega1, omega2 = 1.5e15, 1.8e15
    omega3 = omega1 + omega2
    z, A1, A2, A3 = cwm.simulate_three_wave_mixing(
        A1_0=1.0, A2_0=1.0, A3_0=0.0, z_max=0.1,
        omega1=omega1, omega2=omega2, omega3=omega3,
        n1=1.5, n2=1.5, n3=1.5, deff=2e-12, delta_k=0.0, n_steps=2000)
    assert np.abs(A3[-1]) > 1e-7  # sum-frequency wave grows measurably from zero


def test_opa_amplifies_seed_with_strong_pump():
    omega1 = omega2 = 1.5e15  # near-infrared, degenerate signal/idler
    omega3 = omega1 + omega2
    z, A1, A2, A3 = cwm.simulate_three_wave_mixing(
        A1_0=0.01, A2_0=0.01, A3_0=50.0, z_max=3000.0,
        omega1=omega1, omega2=omega2, omega3=omega3,
        n1=1.5, n2=1.5, n3=1.5, deff=2e-12, delta_k=0.0, n_steps=2000)
    assert np.abs(A1[-1]) > np.abs(A1[0])  # signal amplified by the strong pump
    assert np.abs(A3[-1]) == pytest.approx(np.abs(A3[0]), rel=1e-3)  # pump ~undepleted


def test_phase_mismatch_reduces_conversion():
    omega1, omega2 = 1.5e15, 1.8e15
    omega3 = omega1 + omega2
    kwargs = dict(A1_0=1.0, A2_0=1.0, A3_0=0.0, z_max=0.1,
                  omega1=omega1, omega2=omega2, omega3=omega3,
                  n1=1.5, n2=1.5, n3=1.5, deff=2e-12, n_steps=1000)
    _, _, _, A3_matched = cwm.simulate_three_wave_mixing(delta_k=0.0, **kwargs)
    _, _, _, A3_mismatched = cwm.simulate_three_wave_mixing(delta_k=200.0, **kwargs)
    assert np.abs(A3_matched[-1]) > np.abs(A3_mismatched[-1])

