"""
opticspy.sources_detectors.lasers
====================================
Four-level laser rate equations: population inversion N(t) and cavity
photon number phi(t) coupled through stimulated emission, integrated
numerically (RK4). Includes the threshold pump rate and above-threshold
steady state derived directly from the same rate equations (gain
clamping), so the closed-form results are guaranteed consistent with
the numerical model rather than being independently-recalled formulas.

Reference: Saleh & Teich, "Fundamentals of Photonics", Ch. 15-16
(Laser Amplifiers / Lasers).
"""

import numpy as np


def rate_equations_rhs(N, phi, R_pump, tau_2, tau_photon, gain_coeff):
    """
    dN/dt = R_pump - N/tau_2 - gain_coeff*N*phi
    dphi/dt = gain_coeff*N*phi - phi/tau_photon
    `gain_coeff` absorbs sigma*c/V (stimulated-emission rate coefficient).
    """
    dN = R_pump - N / tau_2 - gain_coeff * N * phi
    dphi = gain_coeff * N * phi - phi / tau_photon
    return dN, dphi


def simulate_laser(R_pump, tau_2, tau_photon, gain_coeff, N0=0.0, phi0=1.0,
                    t_max=None, n_steps=2000):
    """
    RK4 integration of the four-level laser rate equations. `phi0`
    should be small but nonzero (seed photon, representing spontaneous
    emission) so the photon number can build up from below threshold.
    Returns (t, N(t), phi(t)).
    """
    if t_max is None:
        t_max = 20 * tau_2
    t = np.linspace(0, t_max, n_steps)
    dt = t[1] - t[0]
    N, phi = N0, phi0
    N_arr, phi_arr = np.zeros(n_steps), np.zeros(n_steps)
    N_arr[0], phi_arr[0] = N, phi
    for i in range(1, n_steps):
        k1N, k1p = rate_equations_rhs(N, phi, R_pump, tau_2, tau_photon, gain_coeff)
        k2N, k2p = rate_equations_rhs(N + 0.5 * dt * k1N, phi + 0.5 * dt * k1p,
                                       R_pump, tau_2, tau_photon, gain_coeff)
        k3N, k3p = rate_equations_rhs(N + 0.5 * dt * k2N, phi + 0.5 * dt * k2p,
                                       R_pump, tau_2, tau_photon, gain_coeff)
        k4N, k4p = rate_equations_rhs(N + dt * k3N, phi + dt * k3p,
                                       R_pump, tau_2, tau_photon, gain_coeff)
        N = N + (dt / 6) * (k1N + 2 * k2N + 2 * k3N + k4N)
        phi = phi + (dt / 6) * (k1p + 2 * k2p + 2 * k3p + k4p)
        phi = max(phi, 0.0)
        N_arr[i], phi_arr[i] = N, phi
    return t, N_arr, phi_arr


def threshold_inversion(tau_photon, gain_coeff):
    """N_th such that gain_coeff*N_th*tau_photon = 1 (round-trip gain = loss)."""
    return 1.0 / (gain_coeff * tau_photon)


def threshold_pump_rate(tau_2, tau_photon, gain_coeff):
    """R_th = N_th / tau_2 (steady-state pump rate needed to sustain N_th below threshold)."""
    return threshold_inversion(tau_photon, gain_coeff) / tau_2


def steady_state_photon_number(R_pump, tau_2, tau_photon, gain_coeff):
    """
    Above-threshold steady-state photon number, derived from the rate
    equations themselves (gain clamping: N stays pinned at N_th above
    threshold): phi_ss = (R_pump - R_th) * tau_photon.
    Returns 0 for R_pump below threshold.
    """
    R_th = threshold_pump_rate(tau_2, tau_photon, gain_coeff)
    return max(0.0, (R_pump - R_th) * tau_photon)


def output_power(phi_ss, tau_photon, photon_energy, output_coupling_fraction=1.0):
    """
    Power leaving the cavity through the loss channel with photon
    lifetime `tau_photon`: P = phi_ss * photon_energy / tau_photon,
    scaled by the fraction of total loss attributable to the output
    coupler (1.0 = idealized, all loss is useful output).
    """
    return output_coupling_fraction * phi_ss * photon_energy / tau_photon


def slope_efficiency(tau_photon, gain_coeff, tau_2, photon_energy, pump_photon_energy,
                      output_coupling_fraction=1.0):
    """
    d(P_out)/d(P_pump) above threshold (differential/slope efficiency),
    obtained analytically by differentiating `output_power` and
    `steady_state_photon_number` with respect to R_pump (R_pump = P_pump/pump_photon_energy):
        slope = output_coupling_fraction * (photon_energy / pump_photon_energy)
    (Every above-threshold pump photon added increases phi_ss by exactly
    tau_photon, per the linear relation above -- the ratio of output to
    pump photon energy sets the quantum-limited slope efficiency.)
    """
    return output_coupling_fraction * photon_energy / pump_photon_energy
