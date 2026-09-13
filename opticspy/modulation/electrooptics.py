"""
opticspy.modulation.electrooptics
====================================
The linear electro-optic (Pockels) effect: field-induced refractive-
index change, phase retardation in a Pockels cell, the half-wave
voltage, and the transmission of a Mach-Zehnder electro-optic amplitude
modulator built from two phase modulators.

Reference: Saleh & Teich, "Fundamentals of Photonics", Ch. 20
(Electro-Optics).
"""

import numpy as np


def index_change_pockels(n0, r_coefficient, E_field):
    """Pockels-effect index change: Delta n = -0.5 * n0^3 * r * E."""
    return -0.5 * n0 ** 3 * r_coefficient * E_field


def phase_shift_pockels(n0, r_coefficient, voltage, electrode_gap, length, wavelength):
    """
    Phase shift accumulated over interaction length `length` in a
    Pockels cell with applied voltage `voltage` across gap
    `electrode_gap`: Delta_phi = (2*pi/wavelength) * Delta_n * length,
    with Delta_n from `index_change_pockels` and E = voltage/electrode_gap.
    """
    E = voltage / electrode_gap
    dn = index_change_pockels(n0, r_coefficient, E)
    return 2 * np.pi * dn * length / wavelength


def half_wave_voltage(n0, r_coefficient, length, electrode_gap, wavelength):
    """
    V_pi: the voltage giving a pi phase shift, found by solving
    phase_shift_pockels(V_pi) = pi:
        V_pi = wavelength * electrode_gap / (n0^3 * r_coefficient * length)
    """
    return wavelength * electrode_gap / (n0 ** 3 * r_coefficient * length)


def mach_zehnder_eo_transmission(voltage, V_pi, bias_phase=0.0):
    """
    Intensity transmission of a push-pull Mach-Zehnder electro-optic
    amplitude modulator: T(V) = cos^2( (pi/2)*(V/V_pi) + bias_phase/2 ).
    """
    phase = np.pi * voltage / V_pi + bias_phase
    return 0.5 * (1 + np.cos(phase))


def eo_modulator_bandwidth_limit(electrode_length, n_microwave, n_optical, c=299_792_458.0):
    """
    Velocity-mismatch (walk-off) 3-dB bandwidth estimate for a
    traveling-wave electro-optic modulator of electrode length L, from
    the microwave/optical group-index mismatch:
        f_3dB ~ 1.4 * c / (pi * L * |n_microwave - n_optical|)
    (standard traveling-wave EOM bandwidth estimate).
    """
    mismatch = abs(n_microwave - n_optical)
    if mismatch == 0:
        return np.inf
    return 1.4 * c / (np.pi * electrode_length * mismatch)
