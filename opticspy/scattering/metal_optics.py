"""
opticspy.scattering.metal_optics
===================================
Reflection from absorbing media (metals), described by a complex
refractive index n_tilde = n + i*k. Normal-incidence reflectance,
optical penetration (skin) depth, and oblique-incidence reflectance via
the existing Fresnel-equation machinery, which already supports complex
indices (so metal reflection at any angle simply reuses
`opticspy.electromagnetic.fresnel` with a complex n2).

Reference: Born & Wolf, "Principles of Optics", Ch. 14 (Optics of
metals).
"""

import numpy as np

from opticspy.electromagnetic import fresnel


def complex_refractive_index(n, k):
    return n + 1j * k


def normal_incidence_reflectance(n, k, n_incident=1.0):
    """R = [(n-n_i)^2 + k^2] / [(n+n_i)^2 + k^2] at normal incidence."""
    return ((n - n_incident) ** 2 + k ** 2) / ((n + n_incident) ** 2 + k ** 2)


def skin_depth(wavelength, k):
    """Optical penetration (skin) depth: delta = wavelength / (4*pi*k)."""
    return wavelength / (4 * np.pi * k)


def absorption_coefficient(wavelength, k):
    """Power absorption coefficient alpha = 4*pi*k/wavelength (Beer-Lambert alpha)."""
    return 4 * np.pi * k / wavelength


def oblique_reflectance_metal(theta_i, n, k, n_incident=1.0):
    """
    Reflectance at oblique incidence on a metal (complex index n+ik),
    reusing the Fresnel amplitude-coefficient machinery directly (it
    already handles complex n2 via a complex sqrt for cos(theta_t)).
    Only the R_s/R_p/R_unpolarized entries of the returned dict are
    physically meaningful here: for a semi-infinite absorbing medium
    there is no transmitted beam in the usual sense (the field decays
    exponentially with depth), so the T_s/T_p entries -- derived
    assuming a real Poynting-vector ratio -- should be ignored; use
    `absorptance` (= 1 - R) for the energy going into the metal instead.
    """
    n2 = complex_refractive_index(n, k)
    return fresnel.reflectance_transmittance(theta_i, n_incident, n2)


def absorptance(theta_i, n, k, n_incident=1.0, polarization='unpolarized'):
    """
    Fraction of incident power absorbed by a semi-infinite metal,
    A = 1 - R (valid since there is no transmitted channel): the
    correct energy balance for an absorbing substrate, as opposed to
    the R+T=1 balance used for lossless dielectrics.
    """
    r = oblique_reflectance_metal(theta_i, n, k, n_incident)
    key = {'s': 'R_s', 'p': 'R_p', 'unpolarized': 'R_unpolarized'}[polarization]
    return 1.0 - r[key]
