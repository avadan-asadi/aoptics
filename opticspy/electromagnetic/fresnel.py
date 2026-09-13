"""
opticspy.electromagnetic.fresnel
=================================
Fresnel reflection/transmission coefficients for a plane interface
between two homogeneous isotropic media, including total internal
reflection (handled with complex cos(theta_t)).

Sign/phase convention: fields ~ exp(i(kz - wt)); s-polarization (TE,
E perpendicular to the plane of incidence) and p-polarization (TM, E
parallel to the plane of incidence).

Reference: Born & Wolf, "Principles of Optics", Ch. 1 (Basic properties
of the electromagnetic field); Hecht, "Optics"; Saleh & Teich, Ch. 6.
"""

import numpy as np


def _cos_theta_t(theta_i, n1, n2):
    """Complex cosine of the transmission angle (handles TIR)."""
    sin_t = (n1 / n2) * np.sin(theta_i)
    cos_t = np.sqrt(1.0 - sin_t ** 2 + 0j)
    return cos_t


def snell_angle(theta_i, n1, n2):
    """Real transmission angle (radians), NaN beyond the critical angle."""
    sin_t = (n1 / n2) * np.sin(theta_i)
    with np.errstate(invalid='ignore'):
        theta_t = np.where(np.abs(sin_t) <= 1.0, np.arcsin(np.clip(sin_t, -1, 1)), np.nan)
    return theta_t


def fresnel_coefficients(theta_i, n1, n2):
    """
    Complex amplitude reflection/transmission coefficients.

    Returns
    -------
    r_s, r_p, t_s, t_p : complex (or complex ndarray, if theta_i is an array)
    """
    theta_i = np.asarray(theta_i, dtype=float)
    cos_i = np.cos(theta_i)
    cos_t = _cos_theta_t(theta_i, n1, n2)

    r_s = (n1 * cos_i - n2 * cos_t) / (n1 * cos_i + n2 * cos_t)
    r_p = (n2 * cos_i - n1 * cos_t) / (n2 * cos_i + n1 * cos_t)
    t_s = (2 * n1 * cos_i) / (n1 * cos_i + n2 * cos_t)
    t_p = (2 * n1 * cos_i) / (n2 * cos_i + n1 * cos_t)
    return r_s, r_p, t_s, t_p


def reflectance_transmittance(theta_i, n1, n2):
    """
    Power reflectance/transmittance for s and p polarization, plus the
    unpolarized (average) values. T accounts for the beam cross-section
    and index change: T = (n2 cos(theta_t) / n1 cos(theta_i)) * |t|^2.
    """
    theta_i = np.asarray(theta_i, dtype=float)
    r_s, r_p, t_s, t_p = fresnel_coefficients(theta_i, n1, n2)
    cos_i = np.cos(theta_i)
    cos_t = _cos_theta_t(theta_i, n1, n2)

    R_s = np.abs(r_s) ** 2
    R_p = np.abs(r_p) ** 2

    factor = np.real((n2 * cos_t) / (n1 * cos_i))
    T_s = factor * np.abs(t_s) ** 2
    T_p = factor * np.abs(t_p) ** 2
    # Beyond critical angle, cos_t is purely imaginary -> factor real part is 0 -> T=0 automatically.

    R_unpolarized = 0.5 * (R_s + R_p)
    T_unpolarized = 0.5 * (T_s + T_p)
    return dict(R_s=R_s, R_p=R_p, T_s=T_s, T_p=T_p,
                R_unpolarized=R_unpolarized, T_unpolarized=T_unpolarized)


def brewster_angle(n1, n2):
    """Angle of incidence at which R_p = 0 (light polarized in-plane fully transmitted)."""
    return np.arctan2(n2, n1)


def critical_angle(n1, n2):
    """Total-internal-reflection onset angle; returns NaN if n1 <= n2 (no TIR possible)."""
    if n1 <= n2:
        return np.nan
    return np.arcsin(n2 / n1)


def normal_incidence_reflectance(n1, n2):
    """R at theta_i = 0 (polarization-independent)."""
    return ((n1 - n2) / (n1 + n2)) ** 2


def tir_phase_shift(theta_i, n1, n2, polarization='s'):
    """
    Phase shift (radians) acquired by the reflected wave under total
    internal reflection (theta_i > critical angle). `polarization` is
    's' or 'p'.
    """
    r_s, r_p, _, _ = fresnel_coefficients(theta_i, n1, n2)
    r = r_s if polarization == 's' else r_p
    return np.angle(r)
