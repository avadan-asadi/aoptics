"""
aoptics.photonics.waveguides
===============================
Guided modes of a symmetric planar (slab) dielectric waveguide: core of
thickness d and index n1, cladding index n2 < n1, computed by solving
the exact TE/TM eigenvalue (dispersion) equations numerically -- this
avoids relying on memorized closed-form mode counts and instead finds
every physically allowed mode directly from the transcendental
equation, for any d, wavelength, n1, n2.

Reference: Saleh & Teich, "Fundamentals of Photonics", Ch. 8 (Guided-
Wave Optics).
"""

import numpy as np
from scipy.optimize import brentq


def _te_dispersion_residual(beta, k0, n1, n2, d, parity):
    """
    Residual of the symmetric-slab TE transcendental equation, zero at a
    guided mode. kappa = transverse wavenumber in the core, gamma = decay
    constant in the cladding.
        even modes: kappa*tan(kappa*d/2) = gamma
        odd modes:  kappa*cot(kappa*d/2) = -gamma   (i.e. -kappa*cot(...) = gamma)
    """
    kappa2 = (n1 * k0) ** 2 - beta ** 2
    gamma2 = beta ** 2 - (n2 * k0) ** 2
    if kappa2 <= 0 or gamma2 <= 0:
        return np.nan
    kappa = np.sqrt(kappa2)
    gamma = np.sqrt(gamma2)
    if parity == 'even':
        return kappa * np.tan(kappa * d / 2) - gamma
    else:
        return -kappa / np.tan(kappa * d / 2) - gamma


def solve_te_modes(wavelength, d, n1, n2, n_search=4000):
    """
    Find all guided TE mode effective indices n_eff = beta/k0 of a
    symmetric slab waveguide by scanning beta in (n2*k0, n1*k0) for sign
    changes of the dispersion residual (both even and odd families) and
    refining each with Brent's method.

    Returns a sorted list of n_eff values (highest first = fundamental mode).
    """
    k0 = 2 * np.pi / wavelength
    beta_lo, beta_hi = n2 * k0 * (1 + 1e-9), n1 * k0 * (1 - 1e-9)
    betas = np.linspace(beta_lo, beta_hi, n_search)

    roots = []
    for parity in ('even', 'odd'):
        vals = np.array([_te_dispersion_residual(b, k0, n1, n2, d, parity) for b in betas])
        finite = np.isfinite(vals)
        for i in range(len(betas) - 1):
            if not (finite[i] and finite[i + 1]):
                continue
            if vals[i] == 0:
                roots.append(betas[i])
            elif vals[i] * vals[i + 1] < 0:
                try:
                    root = brentq(_te_dispersion_residual, betas[i], betas[i + 1],
                                   args=(k0, n1, n2, d, parity))
                    roots.append(root)
                except ValueError:
                    continue
    n_eff = sorted({round(r / k0, 12) for r in roots}, reverse=True)
    return n_eff


def number_of_te_modes(wavelength, d, n1, n2):
    return len(solve_te_modes(wavelength, d, n1, n2))


def waveguide_v_number(wavelength, d, n1, n2):
    """
    Slab-waveguide normalized frequency (analogous to the fiber V-number):
        V = (pi * d / wavelength) * sqrt(n1^2 - n2^2)
    Guided-mode count for a symmetric slab scales as roughly V/(pi/2)
    rounded up (each additional half-pi of V adds one more even/odd mode);
    `number_of_te_modes` gives the exact count from the dispersion relation.
    """
    NA = np.sqrt(n1 ** 2 - n2 ** 2)
    return (np.pi * d / wavelength) * NA


def mode_confinement_factor(n_eff, k0, n1, n2, d):
    """
    Fraction of the modal power confined to the core, for a symmetric
    slab TE mode of the given effective index (computed by numerically
    integrating |E(x)|^2 over the core vs. the full transverse profile).
    """
    beta = n_eff * k0
    kappa = np.sqrt(max((n1 * k0) ** 2 - beta ** 2, 0.0))
    gamma = np.sqrt(max(beta ** 2 - (n2 * k0) ** 2, 0.0))

    x_core = np.linspace(-d / 2, d / 2, 2000)
    x_clad = np.linspace(d / 2, d / 2 + 8 / max(gamma, 1e-9), 2000)

    # even-mode field shape cos(kappa x) in core, matched exponential decay outside
    E_core = np.cos(kappa * x_core)
    E_clad_edge = np.cos(kappa * d / 2)
    E_clad = E_clad_edge * np.exp(-gamma * (x_clad - d / 2))

    P_core = np.trapezoid(E_core ** 2, x_core)
    P_clad = 2 * np.trapezoid(E_clad ** 2, x_clad)  # two claddings by symmetry
    return P_core / (P_core + P_clad)

