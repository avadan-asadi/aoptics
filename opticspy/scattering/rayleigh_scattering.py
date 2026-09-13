"""
opticspy.scattering.rayleigh_scattering
==========================================
Rayleigh scattering by particles much smaller than the wavelength
(dipole scattering): the scattering cross section and its characteristic
1/wavelength^4 law, the dipole angular scattering (phase) function, and
the degree of linear polarization of singly-scattered light (the
textbook explanation of blue-sky polarization).

Reference: Born & Wolf, "Principles of Optics", Ch. 14 (Scattering from
inhomogeneous media, in the small-particle limit); Hecht, "Optics",
sky-polarization discussion; Bohren & Huffman, "Absorption and
Scattering of Light by Small Particles".
"""

import numpy as np


def rayleigh_cross_section(wavelength, particle_diameter, n_particle, n_medium=1.0):
    """
    Rayleigh scattering cross section for a small dielectric sphere:
        sigma = (2*pi^5/3) * (d^6/wavelength^4) * ((m^2-1)/(m^2+2))^2
    with m = n_particle/n_medium the relative refractive index, valid
    for particle_diameter << wavelength.
    """
    m = n_particle / n_medium
    factor = ((m ** 2 - 1) / (m ** 2 + 2)) ** 2
    return (2 * np.pi ** 5 / 3) * (particle_diameter ** 6 / wavelength ** 4) * factor


def rayleigh_phase_function(theta, unpolarized=True):
    """
    Angular scattering (phase) function for a dipole scatterer, from
    unpolarized incident light: f(theta) = (3/4)*(1+cos^2(theta)) / (4*pi)
    normalized so that its integral over 4*pi steradians is 1.
    """
    return (3.0 / 4.0) * (1 + np.cos(theta) ** 2) / (4 * np.pi)


def rayleigh_polarization_degree(theta_scatter):
    """
    Degree of linear polarization of Rayleigh-scattered light from
    initially unpolarized illumination, at scattering angle theta:
        DoP(theta) = sin^2(theta) / (1 + cos^2(theta))
    Maximum (DoP=1) at theta = 90 deg -- the classic blue-sky-
    polarization result.
    """
    return np.sin(theta_scatter) ** 2 / (1 + np.cos(theta_scatter) ** 2)


def relative_scattering_vs_wavelength(wavelengths, reference_wavelength=550e-9):
    """Relative Rayleigh scattering strength vs wavelength, normalized to 1 at
    `reference_wavelength`: scales as 1/wavelength^4 (blue skies, red sunsets)."""
    return (reference_wavelength / np.asarray(wavelengths)) ** 4


def scattering_mean_free_path(number_density, cross_section):
    """Mean free path between scattering events: 1 / (N * sigma)."""
    return 1.0 / (number_density * cross_section)
