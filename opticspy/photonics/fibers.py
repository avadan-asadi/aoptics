"""
aoptics.photonics.fibers
===========================
Step-index optical fiber parameters: V-number, numerical aperture,
mode count, single-mode cutoff, mode-field diameter, and chromatic
dispersion (material dispersion of fused silica via the Sellmeier
equation, computed numerically -- no hard-coded dispersion-parameter
formula -- plus the standard fiber V-number and single-mode criterion).

Reference: Saleh & Teich, "Fundamentals of Photonics", Ch. 9 (Fiber
Optics).
"""

import numpy as np

C_LIGHT = 299_792_458.0  # m/s

# Sellmeier coefficients for fused silica (Malitson, J. Opt. Soc. Am. 55,
# 1205 (1965)), wavelength in micrometers.
_SIO2_B = [0.6961663, 0.4079426, 0.8974794]
_SIO2_LAM = [0.0684043, 0.1162414, 9.896161]  # micrometers


def numerical_aperture(n_core, n_clad):
    return np.sqrt(n_core ** 2 - n_clad ** 2)


def v_number(wavelength, core_radius, n_core, n_clad):
    """Fiber normalized frequency V = (2*pi*a/wavelength) * NA."""
    return (2 * np.pi * core_radius / wavelength) * numerical_aperture(n_core, n_clad)


def is_single_mode(wavelength, core_radius, n_core, n_clad):
    """Single-mode (only LP01 guided) below the V = 2.405 cutoff (first zero of J0)."""
    return v_number(wavelength, core_radius, n_core, n_clad) < 2.405


def number_of_modes_multimode(wavelength, core_radius, n_core, n_clad):
    """Weakly-guiding step-index multimode approximation: M ~ V^2 / 2."""
    V = v_number(wavelength, core_radius, n_core, n_clad)
    return V ** 2 / 2.0


def mode_field_radius(wavelength, core_radius, n_core, n_clad):
    """
    Marcuse's empirical formula for the (near-Gaussian) mode-field radius
    w of the LP01 mode of a step-index fiber, valid for 1.2 < V < 2.4:
        w/a = 0.65 + 1.619/V^1.5 + 2.879/V^6
    """
    V = v_number(wavelength, core_radius, n_core, n_clad)
    ratio = 0.65 + 1.619 / V ** 1.5 + 2.879 / V ** 6
    return ratio * core_radius


def sellmeier_index_silica(wavelength):
    """Refractive index of fused silica via the 3-term Sellmeier equation.
    `wavelength` in meters (scalar or array)."""
    lam_um = np.asarray(wavelength, dtype=float) * 1e6
    n2 = np.ones_like(lam_um)
    for B, L in zip(_SIO2_B, _SIO2_LAM):
        n2 = n2 + B * lam_um ** 2 / (lam_um ** 2 - L ** 2)
    return np.sqrt(n2)


def group_index(index_fn, wavelength, dlam=1e-11):
    """
    Group index n_g = n - lambda * dn/dlambda, computed by numerical
    differentiation of any supplied index_fn(wavelength) (e.g.
    `sellmeier_index_silica`) -- avoids hard-coding a dispersion formula.
    """
    wavelength = np.asarray(wavelength, dtype=float)
    n = index_fn(wavelength)
    dn = (index_fn(wavelength + dlam) - index_fn(wavelength - dlam)) / (2 * dlam)
    return n - wavelength * dn


def material_dispersion_parameter(index_fn, wavelength, dlam=1e-11):
    """
    Material dispersion parameter D_mat = -(wavelength/c) * d^2n/dlambda^2,
    in s/m^2 (multiply by 1e6 to get the conventional ps/(nm.km)),
    computed by numerical second-differentiation of index_fn.
    """
    wavelength = np.asarray(wavelength, dtype=float)
    n_plus = index_fn(wavelength + dlam)
    n_0 = index_fn(wavelength)
    n_minus = index_fn(wavelength - dlam)
    d2n = (n_plus - 2 * n_0 + n_minus) / dlam ** 2
    return -(wavelength / C_LIGHT) * d2n


def zero_dispersion_wavelength(index_fn, lam_min, lam_max, n_points=2000):
    """
    Find the wavelength(s) where the material dispersion parameter
    crosses zero, by scanning `material_dispersion_parameter` for sign
    changes between lam_min and lam_max (meters) and refining with
    linear interpolation.
    """
    lams = np.linspace(lam_min, lam_max, n_points)
    D = np.array([material_dispersion_parameter(index_fn, l) for l in lams])
    zeros = []
    for i in range(len(lams) - 1):
        if D[i] == 0:
            zeros.append(lams[i])
        elif D[i] * D[i + 1] < 0:
            frac = -D[i] / (D[i + 1] - D[i])
            zeros.append(lams[i] + frac * (lams[i + 1] - lams[i]))
    return zeros


def pulse_broadening_material(length, spectral_width, wavelength, index_fn=sellmeier_index_silica):
    """
    Pulse broadening (s) from material dispersion alone over fiber
    `length` (m), for a source of rms spectral width `spectral_width`
    (m), at `wavelength`: Delta_t = D_mat(wavelength) * length * spectral_width.
    """
    D = material_dispersion_parameter(index_fn, wavelength)
    return abs(D) * length * spectral_width

