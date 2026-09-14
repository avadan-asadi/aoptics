"""
aoptics.photonics.nonlinear_optics
=====================================
Core quantities of chi(2) (second-harmonic generation, phase matching)
and chi(3) (Kerr effect / self-phase modulation, four-wave mixing)
nonlinear optics.

Reference: Saleh & Teich, "Fundamentals of Photonics", Ch. 21 (Nonlinear
Optics); Agrawal, "Nonlinear Fiber Optics", Ch. 2 (SPM parameters).
"""

import numpy as np


# ---------------------------------------------------------------------------
# chi(2): second-harmonic generation & phase matching
# ---------------------------------------------------------------------------

def phase_mismatch(k_2omega, k_omega):
    """Delta k = k(2*omega) - 2*k(omega), for collinear SHG."""
    return k_2omega - 2 * k_omega


def coherence_length_shg(delta_k):
    """SHG coherence length L_coh = pi / |Delta k| (crystal length over
    which the second-harmonic field stays in phase with the driving
    nonlinear polarization)."""
    return np.pi / np.abs(delta_k)


def shg_phase_matching_efficiency(delta_k, length):
    """
    Relative SHG conversion efficiency vs. the ideal phase-matched case,
    for undepleted-pump, plane-wave SHG:
        eta_rel(Delta k, L) = sinc^2(Delta k * L / 2)
    using sinc(u) = sin(pi u)/(pi u) (np.sinc convention with the
    argument pre-divided by pi to match the physics sinc(x)=sin(x)/x form).
    """
    x = delta_k * length / 2.0
    return (np.sin(x) / x) ** 2 if x != 0 else 1.0


def shg_phase_matching_efficiency_array(delta_k, length):
    x = np.asarray(delta_k) * length / 2.0
    out = np.ones_like(x, dtype=float)
    nz = x != 0
    out[nz] = (np.sin(x[nz]) / x[nz]) ** 2
    return out


def shg_conversion_efficiency_undepleted(kappa, length, delta_k=0.0):
    """
    Normalized SHG power-conversion efficiency for the undepleted-pump
    regime with a nonlinear coupling coefficient `kappa` (units 1/(m*sqrt(W)),
    absorbing all the material/geometric prefactors -- deff, omega, n, c,
    epsilon0, effective area -- that are system-specific):
        P_2omega/P_omega = (kappa * sqrt(P_omega) * L)^2 * sinc^2(Delta k L/2)
    Provided as a normalized building block; supply `kappa` from your
    material/beam parameters for an absolute efficiency estimate.
    """
    x = delta_k * length / 2.0
    sinc2 = 1.0 if x == 0 else (np.sin(x) / x) ** 2
    return (kappa * length) ** 2 * sinc2


# ---------------------------------------------------------------------------
# chi(3): Kerr effect, self-phase modulation, four-wave mixing
# ---------------------------------------------------------------------------

def nonlinear_phase_spm(gamma, power, length):
    """
    Self-phase-modulation nonlinear phase shift accumulated over fiber
    `length`, with the standard fiber nonlinear parameter
    gamma = 2*pi*n2/(wavelength*A_eff)  (units 1/(W*m)):
        phi_NL = gamma * P * L
    """
    return gamma * power * length


def nonlinear_length(gamma, power):
    """Fiber nonlinear length L_NL = 1/(gamma * P)."""
    return 1.0 / (gamma * power)


def kerr_index_change(n2, intensity):
    """Kerr-effect refractive index change Delta n = n2 * I (n2 in m^2/W, I in W/m^2)."""
    return n2 * intensity


def fiber_gamma(n2, wavelength, A_eff):
    """Fiber nonlinear parameter gamma = 2*pi*n2 / (wavelength * A_eff)."""
    return 2 * np.pi * n2 / (wavelength * A_eff)


def spm_induced_chirp(t, intensity_envelope):
    """
    Instantaneous SPM-induced angular-frequency chirp
        delta_omega(t) = -d(phi_NL)/dt, proportional to -d|A(t)|^2/dt,
    computed by numerical differentiation of the supplied intensity
    envelope |A(t)|^2 (already time-normalized by gamma*L so its
    derivative directly gives the chirp) -- avoids assuming a specific
    pulse shape.
    """
    t = np.asarray(t, dtype=float)
    I = np.asarray(intensity_envelope, dtype=float)
    return -np.gradient(I, t)


def four_wave_mixing_phase_mismatch(k1, k2, k3, k4):
    """Phase mismatch for a general four-wave-mixing process
    (2 pumps + signal -> new wave, or non-degenerate FWM): Delta k = k3+k4-k1-k2."""
    return k3 + k4 - k1 - k2


def soliton_order(gamma, power, T0, beta2):
    """
    Soliton order N for pulse propagation in an anomalous-dispersion
    fiber (beta2 < 0), N^2 = gamma * P0 * T0^2 / |beta2|
    (N=1 is the fundamental soliton).
    """
    return np.sqrt(gamma * power * T0 ** 2 / abs(beta2))


def dispersion_length(T0, beta2):
    """Fiber dispersion length L_D = T0^2 / |beta2|."""
    return T0 ** 2 / abs(beta2)

