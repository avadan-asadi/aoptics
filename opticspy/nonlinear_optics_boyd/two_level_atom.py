"""
opticspy.nonlinear_optics_boyd.two_level_atom
=================================================
Nonlinear optics in the two-level approximation: the saturated,
power-broadened Lorentzian absorption/dispersion lineshape that results
from the steady-state solution of the optical Bloch equations for a
two-level atom driven by a monochromatic field -- the central result
used throughout saturable-absorber, laser-gain-medium, and
electromagnetically-induced-transparency-adjacent discussions.

Reference: Boyd, "Nonlinear Optics", Ch. 6 (Nonlinear Optics in the
Two-Level Approximation), and Ch. 3 (semiclassical theory).
"""

import numpy as np


def saturation_parameter(intensity, saturation_intensity, detuning=0.0, linewidth=1.0):
    """s = (I/Isat) / (1 + (2*Delta/Gamma)^2), the effective saturation
    parameter including detuning-dependent power broadening."""
    return (intensity / saturation_intensity) / (1 + (2 * detuning / linewidth) ** 2)


def saturated_absorption_coefficient(alpha0, intensity, saturation_intensity,
                                      detuning=0.0, linewidth=1.0):
    """
    alpha(I, Delta) = alpha0 / (1 + (2*Delta/Gamma)^2 + I/Isat)
    the standard saturated/power-broadened Lorentzian absorption
    coefficient of a two-level medium.
    """
    return alpha0 / (1 + (2 * detuning / linewidth) ** 2 + intensity / saturation_intensity)


def power_broadened_linewidth(linewidth0, intensity, saturation_intensity):
    """Gamma' = Gamma0 * sqrt(1 + I/Isat) (power broadening of the resonance)."""
    return linewidth0 * np.sqrt(1 + intensity / saturation_intensity)


def two_level_susceptibility_shape(detuning, linewidth, saturation_param=0.0):
    """
    Normalized (dimensionless-shape) real and imaginary parts of the
    two-level susceptibility (dispersion and absorption), both sharing
    the same saturated-Lorentzian denominator:
        chi'' (absorption)  ~  1 / D
        chi'  (dispersion)  ~  (2*Delta/Gamma) / D
    with D = 1 + (2*Delta/Gamma)^2 + saturation_param.
    Returns (chi_real_shape, chi_imag_shape).
    """
    x = 2 * detuning / linewidth
    D = 1 + x ** 2 + saturation_param
    return x / D, 1.0 / D


def steady_state_population_inversion(rabi_frequency, T1, T2, detuning=0.0):
    """
    Steady-state population inversion w_ss = -1/(1+s) of a driven
    two-level atom (w=-1: fully in the ground state, unsaturated;
    w -> 0: equal populations, strongly saturated), with the Bloch-
    equation saturation parameter s = Omega^2*T1*T2 / (1+(Delta*T2)^2).
    """
    s = (rabi_frequency ** 2 * T1 * T2) / (1 + (detuning * T2) ** 2)
    return -1.0 / (1 + s)


def saturation_intensity_estimate(photon_energy, cross_section, upper_state_lifetime):
    """
    Order-of-magnitude two-level saturation intensity estimate:
        Isat = h*nu / (sigma * tau)
    (photon energy divided by the absorption cross section and the
    upper-state lifetime) -- a standard, dimensionally-transparent
    estimate used throughout laser physics.
    """
    return photon_energy / (cross_section * upper_state_lifetime)
