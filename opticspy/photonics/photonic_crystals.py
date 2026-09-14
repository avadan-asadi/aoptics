"""
aoptics.photonics.photonic_crystals
=======================================
One-dimensional photonic crystals (multilayer dielectric stacks /
distributed Bragg reflectors), via the standard characteristic-matrix
(transfer-matrix) method of thin-film optics: each layer contributes a
2x2 matrix relating tangential E and H fields at its two interfaces,
and the stack's overall matrix gives the reflectance/transmittance.

Reference: Born & Wolf, "Principles of Optics", Ch. 1.6 (Reflection and
refraction at the boundary of a stratified medium / multilayer films);
Saleh & Teich, "Fundamentals of Photonics", Ch. 7 discussion of
periodic layered media (also the basis of distributed Bragg reflectors
used in VCSEL / fiber-Bragg-grating-adjacent devices).
"""

import numpy as np


def _layer_admittance(n, theta, polarization):
    """Optical admittance eta = n*cos(theta) for TE (s), n/cos(theta) for TM (p)."""
    if polarization == 's':
        return n * np.cos(theta)
    return n / np.cos(theta)


def _snell_angles(n_list, theta0):
    """Propagate the angle of incidence through each layer via Snell's law
    (n0 sin(theta0) = n_j sin(theta_j), complex-safe)."""
    n0 = n_list[0]
    sin0 = n0 * np.sin(theta0)
    thetas = [np.arcsin(np.clip(sin0 / n, -1, 1) + 0j) for n in n_list]
    return thetas


def layer_matrix(n, d, wavelength, theta, polarization='s'):
    """Characteristic matrix of a single homogeneous layer of index n,
    thickness d, at angle theta (inside the layer) and the given polarization."""
    delta = 2 * np.pi * n * d * np.cos(theta) / wavelength
    eta = _layer_admittance(n, theta, polarization)
    return np.array([
        [np.cos(delta), 1j * np.sin(delta) / eta],
        [1j * eta * np.sin(delta), np.cos(delta)],
    ], dtype=complex)


def multilayer_reflectance(layers, wavelength, n_incident=1.0, n_substrate=1.5,
                            theta0=0.0, polarization='s'):
    """
    Reflectance/transmittance of a stack of layers [(n1,d1), (n2,d2), ...]
    sandwiched between a semi-infinite incident medium (n_incident) and
    substrate (n_substrate), via the characteristic-matrix method.

    Returns dict(R=..., T=...).
    """
    n_list = [n_incident] + [n for n, d in layers] + [n_substrate]
    thetas = _snell_angles(n_list, theta0)
    theta_inc, theta_sub = thetas[0], thetas[-1]

    M = np.eye(2, dtype=complex)
    for (n, d), theta in zip(layers, thetas[1:-1]):
        M = M @ layer_matrix(n, d, wavelength, theta, polarization)

    eta_inc = _layer_admittance(n_incident, theta_inc, polarization)
    eta_sub = _layer_admittance(n_substrate, theta_sub, polarization)

    B, C = M[0, 0], M[0, 1]
    Cc, D = M[1, 0], M[1, 1]
    # standard thin-film formula: Y = (C_row2)/(row1) combination
    denom = eta_inc * (B + C * eta_sub) + (Cc + D * eta_sub)
    r = (eta_inc * (B + C * eta_sub) - (Cc + D * eta_sub)) / denom
    t = 2 * eta_inc / denom

    R = np.abs(r) ** 2
    T = np.real(eta_sub) / np.real(eta_inc) * np.abs(t) ** 2
    return dict(R=float(R), T=float(T))


def quarter_wave_stack(n_high, n_low, center_wavelength, n_periods,
                        n_incident=1.0, n_substrate=None):
    """
    Build a quarter-wave Bragg mirror layer list: alternating high/low
    index layers each of optical thickness wavelength/4 at the design
    (center) wavelength -- the standard distributed Bragg reflector (DBR).
    """
    d_high = center_wavelength / (4 * n_high)
    d_low = center_wavelength / (4 * n_low)
    layers = []
    for _ in range(n_periods):
        layers.append((n_high, d_high))
        layers.append((n_low, d_low))
    return layers


def reflectance_spectrum(layers, wavelengths, n_incident=1.0, n_substrate=1.5,
                          theta0=0.0, polarization='s'):
    R = np.array([multilayer_reflectance(layers, wl, n_incident, n_substrate,
                                          theta0, polarization)['R'] for wl in wavelengths])
    return R


def bragg_stopband_edges(n_high, n_low, center_wavelength):
    """
    Approximate fractional stopband width of an idealized infinite
    quarter-wave stack (normal incidence):
        Delta_wavelength / wavelength0 = (4/pi) * arcsin[(n_high-n_low)/(n_high+n_low)]
    (standard photonic-bandgap-width result for a 1D Bragg stack).
    Returns (wavelength_low, wavelength_high) of the stopband.
    """
    frac = (4.0 / np.pi) * np.arcsin((n_high - n_low) / (n_high + n_low))
    half = 0.5 * frac * center_wavelength
    return center_wavelength - half, center_wavelength + half

