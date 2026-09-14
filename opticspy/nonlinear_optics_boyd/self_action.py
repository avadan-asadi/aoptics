"""
aoptics.nonlinear_optics_boyd.self_action
==============================================
Processes resulting from the intensity-dependent refractive index
n = n0 + n2*I: self-focusing (critical power, and the empirical
Marburger/Dawes-Marburger collapse-distance formula fit to numerical
simulations of the nonlinear Schrodinger equation), and self-phase
modulation.

Reference: Boyd, "Nonlinear Optics", Ch. 4 (The Intensity-Dependent
Refractive Index) and Ch. 7 (Processes Resulting from the Intensity-
Dependent Refractive Index).
"""

import numpy as np


def intensity_dependent_index(n0, n2, intensity):
    """n(I) = n0 + n2*I (Kerr effect, third-order intensity-dependent index)."""
    return n0 + n2 * intensity


def critical_power_self_focusing(wavelength, n0, n2):
    """
    Marburger critical power for self-focusing of a Gaussian beam:
        Pcr = 3.77 * wavelength^2 / (8*pi*n0*n2)
    (numerically determined coefficient 3.77, from Marburger's analysis
    of the nonlinear Schrodinger equation for a Gaussian beam -- widely
    cited and reproduced, e.g. in laser-damage and filamentation literature).
    """
    return 3.77 * wavelength ** 2 / (8 * np.pi * n0 * n2)


def is_self_focusing(power, critical_power):
    return power > critical_power


def self_focusing_collapse_distance(power, critical_power, rayleigh_range):
    """
    Marburger (Dawes-Marburger) empirical self-focusing collapse
    distance, fit to numerical solutions of the nonlinear Schrodinger
    equation for a collimated Gaussian beam:
        z_f = 0.367 * z_R / sqrt[ (sqrt(P/Pcr) - 0.852)^2 - 0.0219 ]
    Only defined (real) for P somewhat above Pcr; returns np.inf
    otherwise (no collapse within the paraxial/CW model).
    """
    ratio = power / critical_power
    arg = (np.sqrt(ratio) - 0.852) ** 2 - 0.0219
    if arg <= 0:
        return np.inf
    return 0.367 * rayleigh_range / np.sqrt(arg)


def nonlinear_phase_kerr(n2, intensity, length, wavelength):
    """Self-phase-modulation nonlinear phase: phi_NL = (2*pi/wavelength) * n2 * I * L."""
    return 2 * np.pi * n2 * intensity * length / wavelength


def self_focusing_length_estimate(n0, n2, intensity, beam_radius):
    """
    Simple ray-optics (thin-lens) estimate of the self-focusing length,
    treating the Kerr-induced index gradient as a lens of focal length
    f ~ beam_radius / (2*sqrt(2*n2*intensity/n0)) -- an order-of-magnitude
    estimate distinct from (and generally shorter/less accurate than)
    the Marburger formula above, included for comparison.
    """
    return beam_radius / (2 * np.sqrt(2 * n2 * intensity / n0))

