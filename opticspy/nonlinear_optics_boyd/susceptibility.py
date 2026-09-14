"""
aoptics.nonlinear_optics_boyd.susceptibility
=================================================
The nonlinear optical susceptibility: the perturbative expansion of the
induced polarization in powers of the applied field, the conventional
d-coefficient used for second-order (chi(2)) processes, and Miller's
rule, an empirical scaling relation used to estimate chi(2) of a new
material from its linear refractive indices.

Convention: SI units throughout (P in C/m^2, E in V/m, chi(n) in
(m/V)^(n-1)). Boyd's text uses both SI and Gaussian (esu) units
side-by-side; to avoid an unnecessary source of unit-conversion error
this module works in SI only.

Reference: Boyd, "Nonlinear Optics", Ch. 1 (The Nonlinear Optical
Susceptibility).
"""

import numpy as np

EPSILON_0 = 8.8541878128e-12  # F/m


def linear_polarization(chi1, E):
    return EPSILON_0 * chi1 * E


def second_order_polarization(chi2, E1, E2):
    """P(2) = eps0 * chi2 * E1 * E2 (degenerate case E1=E2=E gives P(2)=eps0*chi2*E^2)."""
    return EPSILON_0 * chi2 * E1 * E2


def third_order_polarization(chi3, E1, E2, E3):
    return EPSILON_0 * chi3 * E1 * E2 * E3


def total_polarization(E, chi1, chi2=0.0, chi3=0.0):
    """P = eps0*(chi1*E + chi2*E^2 + chi3*E^3), the standard perturbative expansion."""
    return EPSILON_0 * (chi1 * E + chi2 * E ** 2 + chi3 * E ** 3)


def d_coefficient_from_chi2(chi2):
    """
    Conventional nonlinear-optical d-coefficient: chi(2) = 2*d
    (Boyd's convention, matched to how d is quoted for common crystals
    such as KDP, BBO, LiNbO3).
    """
    return chi2 / 2.0


def chi2_from_d_coefficient(d):
    return 2.0 * d


def millers_delta(chi2, chi1_omega1, chi1_omega2, chi1_omega3):
    """
    Miller's delta: Delta = chi(2)(omega3=omega1+omega2) /
                             [chi1(omega1) * chi1(omega2) * chi1(omega3)]
    An empirical quantity found to be roughly constant (same order of
    magnitude) across many different nonlinear crystals -- used to
    estimate chi(2) for a new material from its linear susceptibilities.
    """
    return chi2 / (chi1_omega1 * chi1_omega2 * chi1_omega3)


def estimate_chi2_miller_rule(millers_delta_value, chi1_omega1, chi1_omega2, chi1_omega3):
    """Inverse of `millers_delta`: estimate chi(2) given an assumed (typical) Miller's delta."""
    return millers_delta_value * chi1_omega1 * chi1_omega2 * chi1_omega3


def chi1_from_refractive_index(n):
    """Linear susceptibility from refractive index: chi1 = n^2 - 1 (non-magnetic medium)."""
    return n ** 2 - 1


def is_symmetric_chi2_tensor(chi2_tensor, atol=1e-8):
    """
    Check full permutation (Kleinman) symmetry of a 3x3x3 chi(2) tensor
    component array chi2_ijk: valid in the limit that all frequencies
    involved are far from any resonance (dispersion of chi(2) negligible),
    in which case the Cartesian indices i,j,k may be freely permuted.
    """
    t = np.asarray(chi2_tensor)
    perms = [(0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)]
    for p in perms:
        if not np.allclose(t, np.transpose(t, p), atol=atol):
            return False
    return True

