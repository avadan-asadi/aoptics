"""
aoptics.nonlinear_optics_boyd.coupled_wave_mixing
======================================================
The coupled-amplitude equations for chi(2) three-wave mixing
(sum-frequency/difference-frequency generation, optical parametric
amplification), integrated numerically (RK4), plus the Manley-Rowe
relations -- the photon-flux conservation laws obeyed by any lossless
three-wave-mixing interaction, used here as a numerical self-consistency
check on the integrated coupled-wave equations rather than as an
independently asserted formula.

Coupled-amplitude equations (undepleted-pump-free, general form; omega3
= omega1 + omega2, Delta_k = k3 - k1 - k2):
    dA1/dz = i*(omega1*deff)/(n1*c) * A3 * conj(A2) * exp(-i*Delta_k*z)
    dA2/dz = i*(omega2*deff)/(n2*c) * A3 * conj(A1) * exp(-i*Delta_k*z)
    dA3/dz = i*(omega3*deff)/(n3*c) * A1 * A2        * exp(+i*Delta_k*z)

Reference: Boyd, "Nonlinear Optics", Ch. 2 (Wave-Equation Description
of Nonlinear Optical Interactions).
"""

import numpy as np

C_LIGHT = 299_792_458.0


def coupled_wave_rhs(A1, A2, A3, z, omega1, omega2, omega3, n1, n2, n3, deff, delta_k):
    dA1 = 1j * (omega1 * deff) / (n1 * C_LIGHT) * A3 * np.conj(A2) * np.exp(-1j * delta_k * z)
    dA2 = 1j * (omega2 * deff) / (n2 * C_LIGHT) * A3 * np.conj(A1) * np.exp(-1j * delta_k * z)
    dA3 = 1j * (omega3 * deff) / (n3 * C_LIGHT) * A1 * A2 * np.exp(1j * delta_k * z)
    return dA1, dA2, dA3


def simulate_three_wave_mixing(A1_0, A2_0, A3_0, z_max, omega1, omega2, omega3,
                                n1, n2, n3, deff, delta_k=0.0, n_steps=2000):
    """
    RK4 integration of the three-wave-mixing coupled-amplitude
    equations from z=0 to z_max. Returns (z, A1(z), A2(z), A3(z)).
    With A3_0 large and A1_0/A2_0 small (seeded or from vacuum noise, a
    small nonzero seed) this models optical parametric amplification;
    with A1_0, A2_0 large and A3_0=0, it models sum-frequency generation.
    """
    z = np.linspace(0, z_max, n_steps)
    dz = z[1] - z[0]
    A1, A2, A3 = complex(A1_0), complex(A2_0), complex(A3_0)
    out1, out2, out3 = np.zeros(n_steps, complex), np.zeros(n_steps, complex), np.zeros(n_steps, complex)
    out1[0], out2[0], out3[0] = A1, A2, A3
    args = (omega1, omega2, omega3, n1, n2, n3, deff, delta_k)
    for i in range(1, n_steps):
        zi = z[i - 1]
        k1 = coupled_wave_rhs(A1, A2, A3, zi, *args)
        k2 = coupled_wave_rhs(A1 + 0.5 * dz * k1[0], A2 + 0.5 * dz * k1[1], A3 + 0.5 * dz * k1[2],
                               zi + 0.5 * dz, *args)
        k3 = coupled_wave_rhs(A1 + 0.5 * dz * k2[0], A2 + 0.5 * dz * k2[1], A3 + 0.5 * dz * k2[2],
                               zi + 0.5 * dz, *args)
        k4 = coupled_wave_rhs(A1 + dz * k3[0], A2 + dz * k3[1], A3 + dz * k3[2],
                               zi + dz, *args)
        A1 = A1 + (dz / 6) * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
        A2 = A2 + (dz / 6) * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        A3 = A3 + (dz / 6) * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2])
        out1[i], out2[i], out3[i] = A1, A2, A3
    return z, out1, out2, out3


def manley_rowe_photon_fluxes(A1, A2, A3, omega1, omega2, omega3):
    """
    Photon-flux-like quantities |A_j|^2/omega_j; the Manley-Rowe
    relations state these satisfy
        d/dz(|A1|^2/omega1) = d/dz(|A2|^2/omega2) = -d/dz(|A3|^2/omega3)
    i.e. for every photon created at omega3, one photon each is
    destroyed at omega1 and omega2 (and vice versa for DFG/OPA).
    """
    return np.abs(A1) ** 2 / omega1, np.abs(A2) ** 2 / omega2, np.abs(A3) ** 2 / omega3


def conversion_efficiency(A3, A3_input_equivalent):
    """Fractional power conversion into the sum/difference-frequency wave."""
    return np.abs(A3) ** 2 / np.abs(A3_input_equivalent) ** 2

