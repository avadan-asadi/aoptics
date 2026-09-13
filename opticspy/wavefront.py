"""
opticspy.wavefront
==================
Wavefront analysis using Zernike polynomials.

Features:
    - Zernike polynomial generation (OSA/ANSI indexing)
    - Wavefront decomposition and reconstruction
    - Strehl ratio and RMS wavefront error
    - Shack-Hartmann simulation
    - Adaptive optics correction
"""

import numpy as np
from math import factorial


# ---------------------------------------------------------------------------
# Zernike polynomials (OSA/ANSI standard)
# ---------------------------------------------------------------------------

def _zernike_radial(n, m, r):
    """Radial Zernike polynomial R_n^|m|(r)."""
    m = abs(m)
    R = np.zeros_like(r, dtype=float)
    for s in range((n - m) // 2 + 1):
        coeff = ((-1)**s * factorial(n - s) /
                 (factorial(s) * factorial((n + m) // 2 - s) *
                  factorial((n - m) // 2 - s)))
        R += coeff * r**(n - 2 * s)
    return R


def zernike(n, m, N=256, normalize=True):
    """
    Compute a single Zernike polynomial Z_n^m on a unit disk.

    Parameters
    ----------
    n, m       : int – radial and azimuthal orders  (|m| <= n, n-|m| even)
    N          : int – grid size
    normalize  : bool – apply OSA normalization

    Returns
    -------
    Z : ndarray (N, N) – Zernike polynomial (NaN outside unit disk)
    """
    x1 = np.linspace(-1, 1, N)
    X, Y = np.meshgrid(x1, x1)
    r = np.sqrt(X**2 + Y**2)
    theta = np.arctan2(Y, X)

    Z = np.full((N, N), np.nan)
    mask = r <= 1.0

    R = _zernike_radial(n, m, r[mask])

    if m > 0:
        Z_val = R * np.cos(m * theta[mask])
    elif m < 0:
        Z_val = R * np.sin(-m * theta[mask])
    else:
        Z_val = R

    if normalize:
        norm = np.sqrt((2 * (n + 1)) / (1 + (m == 0)))
        Z_val *= norm

    Z[mask] = Z_val
    return Z


def zernike_basis(n_max=4, N=256):
    """
    Generate all Zernike polynomials up to radial order n_max.

    Returns
    -------
    basis : list of (n, m, Z) tuples
    """
    basis = []
    for n in range(n_max + 1):
        for m in range(-n, n + 1, 2 if n > 0 else 1):
            if (n - abs(m)) % 2 == 0:
                Z = zernike(n, m, N)
                basis.append((n, m, Z))
    return basis


def wavefront_decompose(wavefront, n_max=4):
    """
    Decompose a wavefront into Zernike coefficients.

    Parameters
    ----------
    wavefront : ndarray (N, N) – phase map [rad], may contain NaN outside pupil
    n_max     : int – maximum radial order

    Returns
    -------
    coefficients : dict {(n,m): coeff}
    reconstructed: ndarray (N, N)
    """
    N = wavefront.shape[0]
    basis = zernike_basis(n_max, N)
    mask = ~np.isnan(wavefront)
    wf_vals = wavefront[mask]

    Z_matrix = np.column_stack([Z[mask] for _, _, Z in basis])
    coeffs, _, _, _ = np.linalg.lstsq(Z_matrix, wf_vals, rcond=None)

    coeff_dict = {(n, m): c for (n, m, _), c in zip(basis, coeffs)}
    reconstructed = np.zeros((N, N))
    for (n, m, Z), c in zip(basis, coeffs):
        reconstructed += c * np.where(np.isnan(Z), 0, Z)

    return coeff_dict, reconstructed


def strehl_ratio(wavefront_rms_rad):
    """
    Maréchal approximation: S ≈ exp(-(2π σ/λ)²) with σ in radians.

    Parameters
    ----------
    wavefront_rms_rad : float – RMS wavefront error in radians

    Returns
    -------
    strehl : float  (1.0 = perfect, 0.0 = fully aberrated)
    """
    return float(np.exp(-wavefront_rms_rad**2))


def wavefront_rms(wavefront):
    """RMS wavefront error [rad] ignoring NaN."""
    return float(np.sqrt(np.nanmean(wavefront**2)))


# ---------------------------------------------------------------------------
# Shack-Hartmann simulator
# ---------------------------------------------------------------------------

def shack_hartmann_simulate(wavefront, lenslet_pitch=500e-6, f_lenslet=10e-3,
                             wavelength=633e-9, pupil_diameter=10e-3):
    """
    Simulate a Shack-Hartmann wavefront sensor.

    Parameters
    ----------
    wavefront      : ndarray (N, N) – phase [rad]
    lenslet_pitch  : float – lenslet array pitch [m]
    f_lenslet      : float – lenslet focal length [m]

    Returns
    -------
    centroids_x, centroids_y : ndarray – spot displacements [m]
    slopes_x, slopes_y       : ndarray – local wavefront slopes [rad/m]
    """
    N = wavefront.shape[0]
    dx = pupil_diameter / N
    n_lenslets = int(pupil_diameter / lenslet_pitch)

    centroids_x = []
    centroids_y = []
    slopes_x = []
    slopes_y = []

    step = N // n_lenslets
    for i in range(n_lenslets):
        for j in range(n_lenslets):
            row_s, row_e = i * step, (i + 1) * step
            col_s, col_e = j * step, (j + 1) * step
            sub = wavefront[row_s:row_e, col_s:col_e]
            valid = ~np.isnan(sub)
            if valid.sum() < step * step // 4:
                continue

            # Local slope from gradient
            gy, gx = np.gradient(np.where(valid, sub, 0))
            sx = np.nanmean(gx[valid]) / dx
            sy = np.nanmean(gy[valid]) / dx

            # Spot displacement at focal plane
            cx = sx * f_lenslet / (2 * np.pi / wavelength)
            cy = sy * f_lenslet / (2 * np.pi / wavelength)
            centroids_x.append(cx)
            centroids_y.append(cy)
            slopes_x.append(sx)
            slopes_y.append(sy)

    return (np.array(centroids_x), np.array(centroids_y),
            np.array(slopes_x), np.array(slopes_y))


# ---------------------------------------------------------------------------
# Named aberrations
# ---------------------------------------------------------------------------

ABERRATION_NAMES = {
    (0, 0): 'Piston',
    (1, -1): 'Tilt Y', (1, 1): 'Tilt X',
    (2, -2): 'Astigmatism 45°', (2, 0): 'Defocus', (2, 2): 'Astigmatism 0°',
    (3, -3): 'Trefoil', (3, -1): 'Coma Y', (3, 1): 'Coma X', (3, 3): 'Trefoil',
    (4, 0): 'Spherical aberration',
}


def generate_aberrated_wavefront(coefficients, N=256):
    """
    Build a wavefront from a dict of Zernike coefficients.

    Parameters
    ----------
    coefficients : dict {(n, m): amplitude_rad}

    Returns
    -------
    wavefront : ndarray (N, N)
    """
    wavefront = np.zeros((N, N))
    for (n, m), amp in coefficients.items():
        Z = zernike(n, m, N)
        wavefront += amp * np.where(np.isnan(Z), 0, Z)
    return wavefront
