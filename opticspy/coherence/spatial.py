"""
opticspy.coherence.spatial
============================
Spatial coherence of light from an incoherent extended source, governed
by the Van Cittert-Zernike theorem: the complex degree of spatial
coherence between two points in an observation plane a distance z from
an incoherent quasi-monochromatic source equals the (normalized) Fourier
transform of the source's intensity distribution, evaluated at spatial
frequency (separation)/(wavelength * z).

Reference: Born & Wolf, "Principles of Optics", Ch. 10.4 (Van
Cittert-Zernike theorem); Goodman, "Statistical Optics", Ch. 5.
"""

import numpy as np
from scipy.special import j1


def jinc(x):
    """jinc(x) = 2*J1(x)/x, with jinc(0) = 1 (the Fourier transform of a
    uniform circular disk -- same function that gives the Airy pattern)."""
    x = np.asarray(x, dtype=float)
    out = np.ones_like(x)
    nz = x != 0
    out[nz] = 2.0 * j1(x[nz]) / x[nz]
    return out


def van_cittert_zernike_circular(source_diameter, wavelength, distance, r):
    """
    Complex degree of spatial coherence gamma(r) between two observation
    points separated by r, due to a uniform incoherent circular source of
    diameter `source_diameter` at distance `distance`:
        gamma(r) = jinc( pi * source_diameter * r / (wavelength * distance) )
    (real and non-negative until the first zero, matching the classic
    Michelson stellar-interferometer result).
    """
    v = np.pi * source_diameter * np.asarray(r) / (wavelength * distance)
    return jinc(v)


def van_cittert_zernike_slit(source_width, wavelength, distance, r):
    """
    Complex degree of spatial coherence for a uniform incoherent slit
    (1-D) source of full width `source_width`:
        gamma(r) = sinc( source_width * r / (wavelength * distance) )
    using the normalized sinc(u) = sin(pi u)/(pi u) convention
    (np.sinc), consistent with Young's double-slit visibility results.
    """
    u = source_width * np.asarray(r) / (wavelength * distance)
    return np.sinc(u)


def van_cittert_zernike_gaussian(source_sigma, wavelength, distance, r):
    """
    Complex degree of spatial coherence for a Gaussian incoherent source
    of intensity standard deviation `source_sigma`: the Fourier transform
    of a Gaussian is a Gaussian, gamma(r) = exp(-2 pi^2 source_sigma^2 r^2
    / (wavelength*distance)^2).
    """
    r = np.asarray(r, dtype=float)
    arg = (2 * np.pi ** 2 * source_sigma ** 2 * r ** 2) / (wavelength * distance) ** 2
    return np.exp(-arg)


def coherence_radius_circular_source(wavelength, distance, source_diameter):
    """Transverse coherence radius (first zero of the VCZ jinc function):
    rho_c = 1.22 * wavelength * distance / source_diameter."""
    return 1.22 * wavelength * distance / source_diameter


def coherence_width_slit_source(wavelength, distance, source_width):
    """First zero of the VCZ sinc function for a slit source: wavelength*distance/source_width."""
    return wavelength * distance / source_width


def young_double_slit_visibility(slit_separation, wavelength, distance,
                                  source_diameter, source_shape='circular'):
    """
    Fringe visibility observed in a Young's double-slit experiment
    illuminated by an incoherent extended source, per the Van
    Cittert-Zernike theorem: V = |gamma12(slit_separation)|.
    """
    if source_shape == 'circular':
        gamma = van_cittert_zernike_circular(source_diameter, wavelength, distance, slit_separation)
    elif source_shape == 'slit':
        gamma = van_cittert_zernike_slit(source_diameter, wavelength, distance, slit_separation)
    else:
        raise ValueError("source_shape must be 'circular' or 'slit'")
    return np.abs(gamma)
