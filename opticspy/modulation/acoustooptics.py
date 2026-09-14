"""
aoptics.modulation.acoustooptics
====================================
Acousto-optic diffraction of light by a traveling sound wave: the
Klein-Cook parameter distinguishing the thin-grating (Raman-Nath) and
thick-grating (Bragg) regimes, the Bragg diffraction angle/condition,
and the acousto-optic frequency shift (Doppler shift of the diffracted
order by the acoustic frequency).

Reference: Saleh & Teich, "Fundamentals of Photonics", Ch. 19
(Acousto-Optics).
"""

import numpy as np


def klein_cook_parameter(wavelength, interaction_length, acoustic_wavelength, n=1.0):
    """
    Klein-Cook parameter Q = 2*pi*wavelength*L / (n * acoustic_wavelength^2).
    Q << 1: Raman-Nath (thin-grating) regime, multiple diffraction orders.
    Q >> 1: Bragg (thick-grating) regime, single dominant diffraction order.
    """
    return 2 * np.pi * wavelength * interaction_length / (n * acoustic_wavelength ** 2)


def regime(Q):
    if Q < 0.3:
        return "Raman-Nath"
    elif Q > 10:
        return "Bragg"
    return "intermediate"


def bragg_angle(wavelength, acoustic_wavelength, n=1.0):
    """
    Bragg angle (angle of incidence, inside the medium, measured from
    the acoustic wavefronts) satisfying the Bragg condition:
        sin(theta_B) = wavelength / (2 * n * acoustic_wavelength)
    """
    return np.arcsin(wavelength / (2 * n * acoustic_wavelength))


def acoustic_wavelength(acoustic_velocity, acoustic_frequency):
    return acoustic_velocity / acoustic_frequency


def diffracted_order_frequency_shift(optical_frequency, acoustic_frequency, order=1):
    """
    Frequency of the diffracted light in order `order` (+1 = upshifted,
    -1 = downshifted), from momentum/energy conservation with the
    traveling acoustic phonon: f_diffracted = f_optical + order*f_acoustic.
    """
    return optical_frequency + order * acoustic_frequency


def raman_nath_diffraction_efficiency(order, raman_nath_parameter):
    """
    Raman-Nath regime diffraction efficiency into order m:
        eta_m = J_m(raman_nath_parameter)^2
    where the Raman-Nath parameter is proportional to sqrt(acoustic
    intensity) (supplied by the caller) and J_m is the m-th order
    Bessel function of the first kind (the standard thin-grating result).
    """
    from scipy.special import jv
    return jv(order, raman_nath_parameter) ** 2


def normalized_bragg_efficiency(coupling_parameter):
    """
    Bragg-regime diffraction efficiency vs. a normalized coupling
    parameter xi (proportional to sqrt(acoustic power) x interaction
    length; xi=pi/2 gives 100% diffraction): eta = sin^2(xi)
    (the standard two-wave coupled-mode Bragg-diffraction result).
    """
    return np.sin(coupling_parameter) ** 2

