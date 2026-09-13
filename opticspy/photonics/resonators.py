"""
opticspy.photonics.resonators
===============================
Fabry-Perot resonators (Airy function transmission, finesse, free
spectral range, linewidth, Q factor) and general two-mirror laser
resonators: stability via the g-parameters, and the self-consistent
Gaussian-beam mode found from the cavity's round-trip ABCD matrix
(reusing `opticspy.geometrical.ray_tracing`).

Reference: Saleh & Teich, "Fundamentals of Photonics", Ch. 7 (Beam
Optics / resonators) and Ch. 9-10 discussions of the Fabry-Perot etalon;
Yariv, "Quantum Electronics", for the ABCD self-consistency law.
"""

import numpy as np

from opticspy.geometrical import ray_tracing as rt

C_LIGHT = 299_792_458.0  # m/s


# ---------------------------------------------------------------------------
# Fabry-Perot etalon / cavity
# ---------------------------------------------------------------------------

def free_spectral_range(cavity_length, n=1.0):
    """FSR (Hz) of a linear two-mirror cavity of optical length n*L: c/(2 n L)."""
    return C_LIGHT / (2 * n * cavity_length)


def finesse(R1, R2):
    """Fabry-Perot finesse from power reflectivities R1, R2: pi*(R1 R2)^(1/4) / (1 - sqrt(R1 R2))."""
    r = np.sqrt(R1 * R2)
    return np.pi * r ** 0.5 / (1 - r)


def linewidth(cavity_length, R1, R2, n=1.0):
    """FWHM linewidth (Hz) of a cavity resonance: FSR / Finesse."""
    return free_spectral_range(cavity_length, n) / finesse(R1, R2)


def quality_factor(frequency, cavity_length, R1, R2, n=1.0):
    """Q = frequency / linewidth."""
    return frequency / linewidth(cavity_length, R1, R2, n)


def airy_transmission(delta, R1, R2):
    """
    Fabry-Perot (Airy function) power transmission as a function of the
    single-pass round-trip phase `delta` (radians):
        T(delta) = (1-R1)(1-R2) / [ (1-sqrt(R1 R2))^2 + 4 sqrt(R1 R2) sin^2(delta/2) ]
    """
    r = np.sqrt(R1 * R2)
    return (1 - R1) * (1 - R2) / ((1 - r) ** 2 + 4 * r * np.sin(delta / 2) ** 2)


def round_trip_phase(frequency, cavity_length, n=1.0):
    """Round-trip phase delta = 2 * (2*pi*frequency*n*L/c) = 4*pi*n*L*frequency/c."""
    return 4 * np.pi * n * cavity_length * frequency / C_LIGHT


# ---------------------------------------------------------------------------
# Two-mirror resonator stability and Gaussian mode (ABCD self-consistency)
# ---------------------------------------------------------------------------

def g_parameters(L, R1, R2):
    """Resonator g-parameters g_i = 1 - L/R_i (mirror sign convention: R>0
    if the mirror's center of curvature faces the cavity interior)."""
    g1 = 1 - L / R1 if np.isfinite(R1) else 1.0
    g2 = 1 - L / R2 if np.isfinite(R2) else 1.0
    return g1, g2


def is_stable(L, R1, R2):
    """Resonator stability criterion: 0 <= g1*g2 <= 1."""
    g1, g2 = g_parameters(L, R1, R2)
    return 0.0 <= g1 * g2 <= 1.0


def round_trip_matrix(L, R1, R2):
    """
    Build the cavity round-trip ABCD matrix (starting just after mirror
    1, propagating to mirror 2, reflecting, back to mirror 1, reflecting)
    by composing `opticspy.geometrical.ray_tracing` building blocks.
    """
    M_prop = rt.free_space_matrix(L)
    M1 = rt.mirror_matrix(R1)
    M2 = rt.mirror_matrix(R2)
    # 1 -> 2 -> reflect at 2 -> 2 -> 1 -> reflect at 1
    return rt.combine_matrices(M_prop, M2, M_prop, M1)


def gaussian_mode_from_round_trip(M, wavelength):
    """
    Self-consistent Gaussian-beam complex parameter q at the reference
    plane of round-trip ABCD matrix M = [[A,B],[C,D]], via Kogelnik's
    ABCD self-consistency law:
        1/q = (D-A)/(2B) + i * sqrt(4 - (A+D)^2) / (2*|B|)      (stable cavity)
    Returns (waist w0 at the reference plane, complex q, stable: bool).
    """
    A, B = M[0, 0], M[0, 1]
    D = M[1, 1]
    stability_param = (A + D) / 2.0
    if abs(stability_param) > 1.0:
        return None, None, False
    inv_q_imag = np.sqrt(1.0 - stability_param ** 2) / abs(B)
    inv_q_real = (D - A) / (2 * B)
    inv_q = inv_q_real + 1j * inv_q_imag
    q = 1.0 / inv_q
    # Convention used throughout this module: 1/q = 1/R + i*wavelength/(pi w^2)
    # (sign fixed by inv_q_imag >= 0 above), so w^2 = wavelength / (pi * Im(1/q)).
    w0 = np.sqrt(wavelength / (np.pi * inv_q.imag))
    return w0, q, True


def resonator_waist_at_mirror1(L, R1, R2, wavelength):
    """
    Fundamental-mode spot size at the plane just after mirror 1, for a
    two-mirror cavity of length L and mirror radii R1, R2 (np.inf for a
    flat mirror), found via the round-trip ABCD matrix and Kogelnik's
    self-consistency law.
    """
    M = round_trip_matrix(L, R1, R2)
    w0, q, stable = gaussian_mode_from_round_trip(M, wavelength)
    return w0, q, stable


def propagate_q(q0, distance):
    """Free-space propagation of the complex Gaussian beam parameter: q(z) = q0 + z."""
    return q0 + distance


def waist_from_q(q, wavelength):
    """
    Beam radius w from the complex parameter q, using the sign
    convention adopted throughout this module (matching
    `gaussian_mode_from_round_trip`): 1/q = 1/R + i*wavelength/(pi w^2),
    so w = sqrt(wavelength / (pi * Im(1/q))).
    """
    inv_q = 1.0 / q
    return np.sqrt(wavelength / (np.pi * inv_q.imag))


def resonator_waist_at_center(L, R1, R2, wavelength):
    """Fundamental-mode spot size at the longitudinal center of a symmetric
    (R1 = R2) two-mirror cavity, propagating the mirror-plane q by L/2."""
    w_mirror, q_mirror, stable = resonator_waist_at_mirror1(L, R1, R2, wavelength)
    if not stable:
        return None, False
    q_center = propagate_q(q_mirror, L / 2)
    return waist_from_q(q_center, wavelength), True


def symmetric_confocal_waist_center(L, wavelength):
    """
    Closed-form fundamental-mode waist at the CENTER of a symmetric
    confocal resonator (R1 = R2 = L): w0^2 = wavelength*L / (2*pi)
    (classic textbook result, used as a cross-check of the general ABCD
    self-consistency machinery).
    """
    return np.sqrt(wavelength * L / (2 * np.pi))


def symmetric_confocal_waist_at_mirror(L, wavelength):
    """Closed-form spot size AT THE MIRRORS of a symmetric confocal
    resonator: w(mirror)^2 = wavelength*L/pi = 2 * w0_center^2."""
    return np.sqrt(wavelength * L / np.pi)
