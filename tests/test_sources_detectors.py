import numpy as np
import pytest

from opticspy.sources_detectors import lasers, photodetectors as pd


# --------------------------- lasers ---------------------------

def test_below_threshold_photon_number_stays_small():
    tau_2, tau_photon, gain_coeff = 1e-9, 1e-11, 1e-7
    R_th = lasers.threshold_pump_rate(tau_2, tau_photon, gain_coeff)
    t, N, phi = lasers.simulate_laser(0.3 * R_th, tau_2, tau_photon, gain_coeff,
                                       phi0=1.0, t_max=50 * tau_2, n_steps=4000)
    assert phi[-1] < 10 * phi[0]  # no macroscopic buildup below threshold


def test_above_threshold_reaches_predicted_steady_state():
    tau_2, tau_photon, gain_coeff = 1e-9, 1e-11, 1e-7
    R_th = lasers.threshold_pump_rate(tau_2, tau_photon, gain_coeff)
    R_pump = 3.0 * R_th
    t, N, phi = lasers.simulate_laser(R_pump, tau_2, tau_photon, gain_coeff,
                                       phi0=1.0, t_max=80 * tau_2, n_steps=6000)
    expected_phi_ss = lasers.steady_state_photon_number(R_pump, tau_2, tau_photon, gain_coeff)
    assert phi[-1] == pytest.approx(expected_phi_ss, rel=0.05)
    N_th = lasers.threshold_inversion(tau_photon, gain_coeff)
    assert N[-1] == pytest.approx(N_th, rel=0.05)  # gain clamping


def test_steady_state_zero_below_threshold():
    tau_2, tau_photon, gain_coeff = 1e-9, 1e-11, 1e-7
    R_th = lasers.threshold_pump_rate(tau_2, tau_photon, gain_coeff)
    assert lasers.steady_state_photon_number(0.5 * R_th, tau_2, tau_photon, gain_coeff) == 0.0


def test_output_power_scales_linearly_above_threshold():
    tau_2, tau_photon, gain_coeff = 1e-9, 1e-11, 1e-7
    R_th = lasers.threshold_pump_rate(tau_2, tau_photon, gain_coeff)
    hv = 3e-19
    phi1 = lasers.steady_state_photon_number(2 * R_th, tau_2, tau_photon, gain_coeff)
    phi2 = lasers.steady_state_photon_number(4 * R_th, tau_2, tau_photon, gain_coeff)
    P1 = lasers.output_power(phi1, tau_photon, hv)
    P2 = lasers.output_power(phi2, tau_photon, hv)
    assert P2 > P1


# --------------------------- photodetectors ---------------------------

def test_responsivity_known_value_at_1550nm():
    # eta=1 ideal detector at 1550 nm: R = q*lambda/(h*c) ~ 1.25 A/W
    R = pd.responsivity(1.0, 1550e-9)
    assert R == pytest.approx(1.25, rel=0.02)


def test_responsivity_scales_with_quantum_efficiency():
    R_full = pd.responsivity(1.0, 900e-9)
    R_half = pd.responsivity(0.5, 900e-9)
    assert R_half == pytest.approx(0.5 * R_full)


def test_quantum_efficiency_from_absorption_limits():
    assert pd.quantum_efficiency_from_absorption(1e6, 0.0) == pytest.approx(0.0)
    assert pd.quantum_efficiency_from_absorption(1e6, 1e-3) == pytest.approx(1.0, abs=1e-6)


def test_shot_noise_increases_with_current_and_bandwidth():
    i1 = pd.shot_noise_current(1e-6, 1e6)
    i2 = pd.shot_noise_current(1e-6, 4e6)
    assert i2 == pytest.approx(2 * i1, rel=1e-6)  # sqrt(bandwidth) scaling


def test_thermal_noise_decreases_with_load_resistance():
    i1 = pd.thermal_noise_current(300, 1e6, 1e3)
    i2 = pd.thermal_noise_current(300, 1e6, 1e5)
    assert i2 < i1


def test_snr_increases_with_optical_power():
    snr1 = pd.signal_to_noise_ratio(1e-6, 0.8, 1550e-9, 1e6, 300, 1e3)
    snr2 = pd.signal_to_noise_ratio(1e-3, 0.8, 1550e-9, 1e6, 300, 1e3)
    assert snr2 > snr1
