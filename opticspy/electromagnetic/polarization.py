"""
opticspy.electromagnetic.polarization
======================================
Jones calculus (fully polarized light) and Stokes/Mueller calculus
(partially polarized light), plus the Poincare-sphere representation.

Convention: fields ~ exp(i(kz - wt)). A Jones vector (Ex, Ey) represents
the complex transverse field amplitudes at z=0 as a function of time,
E(t) = Re[(Ex, Ey) exp(-i w t)]. With this convention, (1, -i)/sqrt(2)
traces out a field vector that rotates counter-clockwise when viewed
looking into the oncoming beam (the "right circular" convention used
throughout this module) -- see `circular_right`/`circular_left` and the
matching Stokes-S3 sign in `jones_to_stokes`.

Reference: Saleh & Teich, "Fundamentals of Photonics", Ch. 6 (Polarization
Optics); Born & Wolf, Ch. 1; Goldstein, "Polarized Light".
"""

import numpy as np

# ---------------------------------------------------------------------------
# Jones vectors
# ---------------------------------------------------------------------------

def linear_horizontal():
    return np.array([1.0, 0.0], dtype=complex)


def linear_vertical():
    return np.array([0.0, 1.0], dtype=complex)


def linear_at_angle(theta):
    """Linearly polarized light, transmission axis at angle theta from horizontal."""
    return np.array([np.cos(theta), np.sin(theta)], dtype=complex)


def circular_right():
    """Right-circular, per this module's convention (see module docstring)."""
    return np.array([1.0, -1j], dtype=complex) / np.sqrt(2)


def circular_left():
    return np.array([1.0, 1j], dtype=complex) / np.sqrt(2)


def elliptical(a, b, phase=0.0, angle=0.0):
    """
    Elliptical polarization with semi-axis amplitudes a (x') and b (y')
    in a frame rotated by `angle`, and relative phase `phase` of the
    y' component (b) with respect to x' (a).
    """
    v = np.array([a, b * np.exp(1j * phase)], dtype=complex)
    v = v / np.linalg.norm(v)
    c, s = np.cos(angle), np.sin(angle)
    R = np.array([[c, -s], [s, c]])
    return R @ v


def normalize_jones(v):
    return np.asarray(v, dtype=complex) / np.linalg.norm(v)


# ---------------------------------------------------------------------------
# Jones matrices
# ---------------------------------------------------------------------------

def rotation_matrix(theta):
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, s], [-s, c]], dtype=complex)


def linear_polarizer(theta=0.0):
    """Ideal linear polarizer, transmission axis at angle theta."""
    J0 = np.array([[1.0, 0.0], [0.0, 0.0]], dtype=complex)
    return _rotate_jones(J0, theta)


def _rotate_jones(J0, theta):
    """Rotate a Jones matrix defined in its natural frame to angle theta."""
    R = rotation_matrix(theta)
    Rinv = rotation_matrix(-theta)
    return Rinv @ J0 @ R


def wave_plate(retardance, fast_axis_angle=0.0):
    """
    General retarder: `retardance` (radians) is the phase delay of the
    slow axis relative to the fast axis. fast_axis_angle measured from
    horizontal.
    """
    J0 = np.array([[1.0, 0.0], [0.0, np.exp(1j * retardance)]], dtype=complex)
    return _rotate_jones(J0, fast_axis_angle)


def quarter_wave_plate(fast_axis_angle=0.0):
    return wave_plate(np.pi / 2, fast_axis_angle)


def half_wave_plate(fast_axis_angle=0.0):
    return wave_plate(np.pi, fast_axis_angle)


def rotator(theta):
    """Optically-active rotator (rotates the polarization ellipse by theta)."""
    return rotation_matrix(-theta)


def apply_jones(matrix, vector):
    return matrix @ np.asarray(vector, dtype=complex)


def jones_intensity(vector):
    return float(np.sum(np.abs(vector) ** 2))


# ---------------------------------------------------------------------------
# Stokes parameters / Mueller calculus
# ---------------------------------------------------------------------------

def jones_to_stokes(vector):
    """Stokes vector [S0,S1,S2,S3] of a fully-polarized Jones vector."""
    Ex, Ey = vector
    S0 = np.real(Ex * np.conj(Ex) + Ey * np.conj(Ey))
    S1 = np.real(Ex * np.conj(Ex) - Ey * np.conj(Ey))
    S2 = np.real(Ex * np.conj(Ey) + Ey * np.conj(Ex))
    S3 = np.real(1j * (Ex * np.conj(Ey) - Ey * np.conj(Ex)))
    return np.array([S0, S1, S2, S3])


def degree_of_polarization(stokes):
    S0, S1, S2, S3 = stokes
    return np.sqrt(S1 ** 2 + S2 ** 2 + S3 ** 2) / S0 if S0 > 0 else 0.0


def poincare_coordinates(stokes):
    """Normalized (s1,s2,s3) coordinates on/inside the Poincare sphere."""
    S0, S1, S2, S3 = stokes
    if S0 == 0:
        return np.array([0.0, 0.0, 0.0])
    return np.array([S1, S2, S3]) / S0


def mueller_rotation(theta):
    c2, s2 = np.cos(2 * theta), np.sin(2 * theta)
    return np.array([
        [1, 0, 0, 0],
        [0, c2, s2, 0],
        [0, -s2, c2, 0],
        [0, 0, 0, 1],
    ], dtype=float)


def mueller_polarizer(theta=0.0):
    """Ideal linear polarizer Mueller matrix, transmission axis at theta."""
    c2, s2 = np.cos(2 * theta), np.sin(2 * theta)
    M = 0.5 * np.array([
        [1, c2, s2, 0],
        [c2, c2 ** 2, c2 * s2, 0],
        [s2, c2 * s2, s2 ** 2, 0],
        [0, 0, 0, 0],
    ], dtype=float)
    return M


def mueller_retarder(retardance, fast_axis_angle=0.0):
    """General retarder Mueller matrix."""
    d = retardance
    M0 = np.array([
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, np.cos(d), -np.sin(d)],
        [0, 0, np.sin(d), np.cos(d)],
    ], dtype=float)
    R = mueller_rotation(fast_axis_angle)
    Rinv = mueller_rotation(-fast_axis_angle)
    return Rinv @ M0 @ R


def mueller_rotator(theta):
    return mueller_rotation(-theta)


def mueller_depolarizer(factor):
    """Ideal diagonal depolarizer: scales S1,S2,S3 by `factor` in [0,1]."""
    return np.diag([1.0, factor, factor, factor]).astype(float)


def apply_mueller(matrix, stokes):
    return matrix @ np.asarray(stokes, dtype=float)


def jones_to_mueller(J):
    """
    Convert a (non-depolarizing) 2x2 complex Jones matrix to its
    equivalent 4x4 Mueller matrix, via M = A (J kron J*) A^{-1}
    with the standard Jones-Stokes transform matrix A.
    """
    A = np.array([
        [1, 0, 0, 1],
        [1, 0, 0, -1],
        [0, 1, 1, 0],
        [0, 1j, -1j, 0],
    ], dtype=complex)
    Ainv = np.linalg.inv(A)
    K = np.kron(J, np.conj(J))
    M = A @ K @ Ainv
    return np.real(M)
