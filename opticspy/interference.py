"""
opticspy.interference
=====================
Two-beam and multi-beam interference patterns, fringe analysis,
and phase reconstruction.

Features:
    - Two-beam interference (any pair of beams)
    - Multi-beam interference (Fabry-Perot)
    - Off-axis holography
    - Phase-shifting interferometry (PSI)
    - Speckle pattern simulation
    - Digital holographic reconstruction
"""

import numpy as np
from .beams import make_grid, gaussian_beam, plane_wave, laguerre_gaussian


# ---------------------------------------------------------------------------
# Core interference engine
# ---------------------------------------------------------------------------

def two_beam_interference(field1, field2):
    """
    Compute intensity and phase of interference between two fields.

    I = |E1 + E2|^2

    Returns
    -------
    intensity    : ndarray – interference intensity
    phase_map    : ndarray – phase of combined field [rad]
    visibility   : float   – fringe visibility (Michelson contrast)
    """
    combined = field1 + field2
    intensity = np.abs(combined)**2
    phase_map = np.angle(combined)

    I1 = np.abs(field1)**2
    I2 = np.abs(field2)**2
    Imax = np.max(intensity)
    Imin = np.min(intensity)
    visibility = (Imax - Imin) / (Imax + Imin + 1e-30)

    return intensity, phase_map, visibility


def multi_beam_interference(fields):
    """
    Interference of an arbitrary list of complex fields.

    Parameters
    ----------
    fields : list of ndarray

    Returns
    -------
    intensity : ndarray
    phase_map : ndarray
    """
    combined = sum(fields)
    return np.abs(combined)**2, np.angle(combined)


def fabry_perot_interference(field, R1=0.95, R2=0.95, n_roundtrips=20,
                              wavelength=633e-9, cavity_length=0.01, L=1e-3):
    """
    Fabry-Perot etalon multi-beam interference.

    Parameters
    ----------
    R1, R2      : float – mirror reflectivities
    n_roundtrips: int   – number of round trips to sum
    cavity_length: float – [m]

    Returns
    -------
    intensity : ndarray
    """
    N = field.shape[0]
    k = 2 * np.pi / wavelength
    T1 = np.sqrt(1 - R1)
    T2 = np.sqrt(1 - R2)
    r1 = np.sqrt(R1)
    r2 = np.sqrt(R2)

    combined = np.zeros_like(field)
    for m in range(n_roundtrips):
        phase = np.exp(1j * 2 * k * cavity_length * m)
        amp = T1 * T2 * (r1 * r2)**m
        combined += amp * phase * field

    return np.abs(combined)**2


# ---------------------------------------------------------------------------
# Specific interference configurations
# ---------------------------------------------------------------------------

def double_slit_interference(N=512, L=5e-3, wavelength=633e-9,
                               slit_sep=0.5e-3, slit_width=0.05e-3,
                               screen_dist=1.0):
    """
    Young's double-slit interference pattern (analytical).

    Returns
    -------
    x_screen : ndarray – screen coordinates [m]
    intensity : ndarray – 1-D intensity pattern
    """
    x_screen = np.linspace(-L, L, N)
    k = 2 * np.pi / wavelength
    d = slit_sep
    a = slit_width

    # Single-slit envelope * double-slit fringes
    beta = k * a * x_screen / (2 * screen_dist)
    delta = k * d * x_screen / screen_dist
    sinc_env = np.where(beta == 0, 1.0, np.sin(beta) / beta)**2
    fringe = np.cos(delta / 2)**2
    intensity = sinc_env * fringe
    return x_screen, intensity


def mach_zehnder(field, phase_obj=None, N=512, L=1e-3, wavelength=633e-9):
    """
    Mach-Zehnder interferometer simulation.

    Parameters
    ----------
    field     : complex ndarray – input beam
    phase_obj : ndarray or None – phase introduced in one arm [rad]

    Returns
    -------
    intensity_out : ndarray – output port intensity
    phase_map     : ndarray – differential phase [rad]
    """
    arm1 = field.copy()
    arm2 = field.copy()
    if phase_obj is not None:
        arm2 = arm2 * np.exp(1j * phase_obj)

    combined = (arm1 + arm2) / np.sqrt(2)
    return np.abs(combined)**2, np.angle(combined)


def vortex_interference(p=0, l=1, N=512, L=1e-3, w0=100e-6,
                         wavelength=633e-9):
    """
    Interference of LG beam with a co-axial plane wave –
    reveals the spiral phase structure (fork fringes).

    Returns
    -------
    intensity : ndarray
    phase_map : ndarray
    """
    from .beams import laguerre_gaussian, gaussian_beam
    lg = laguerre_gaussian(p, l, N, L, w0, wavelength)
    ref = gaussian_beam(N, L, w0 * 3, wavelength)   # broad reference
    return two_beam_interference(lg, ref)


# ---------------------------------------------------------------------------
# Phase-Shifting Interferometry
# ---------------------------------------------------------------------------

def phase_shifting_interferometry(field_obj, field_ref, N_steps=4):
    """
    Phase-shifting interferometry (PSI) – extract wrapped phase from
    N_steps interferograms with equal phase shifts.

    Parameters
    ----------
    field_obj : complex ndarray – object beam
    field_ref : complex ndarray – reference beam
    N_steps   : int – number of phase steps (3 or 4 recommended)

    Returns
    -------
    wrapped_phase : ndarray [rad]
    interferograms: list of ndarray – individual intensity patterns
    """
    shifts = np.linspace(0, 2 * np.pi, N_steps, endpoint=False)
    interferograms = []
    for shift in shifts:
        ref_shifted = field_ref * np.exp(1j * shift)
        I = np.abs(field_obj + ref_shifted)**2
        interferograms.append(I)

    # 4-step algorithm: phi = arctan2(I4-I2, I1-I3)
    if N_steps == 4:
        wrapped_phase = np.arctan2(
            interferograms[3] - interferograms[1],
            interferograms[0] - interferograms[2]
        )
    else:
        # Generalized least-squares
        A = np.array([[np.cos(s), np.sin(s), 1] for s in shifts])
        Iarr = np.array([I.ravel() for I in interferograms])
        coeffs, _, _, _ = np.linalg.lstsq(A, Iarr, rcond=None)
        wrapped_phase = np.arctan2(-coeffs[1], coeffs[0]).reshape(
            field_obj.shape)

    return wrapped_phase, interferograms


# ---------------------------------------------------------------------------
# Speckle
# ---------------------------------------------------------------------------

def speckle_pattern(N=512, L=1e-3, wavelength=633e-9,
                     roughness_scale=10e-6, seed=None):
    """
    Simulate a fully developed speckle pattern from a rough surface.

    Parameters
    ----------
    roughness_scale : float – spatial scale of surface roughness [m]
    seed            : int   – random seed for reproducibility

    Returns
    -------
    intensity : ndarray – speckle intensity pattern
    field     : complex ndarray – speckle field
    """
    rng = np.random.default_rng(seed)
    # Random phase on rough surface
    surface_phase = rng.uniform(0, 2 * np.pi, (N, N))
    x, y, r, _ = make_grid(N, L)
    # Gaussian envelope
    w = roughness_scale * 10
    envelope = np.exp(-(r**2) / (2 * w**2))
    field = envelope * np.exp(1j * surface_phase)

    # Fraunhofer diffraction → speckle at far field
    speckle_field = np.fft.fftshift(np.fft.fft2(field))
    intensity = np.abs(speckle_field)**2
    return intensity, speckle_field


# ---------------------------------------------------------------------------
# Off-axis holography
# ---------------------------------------------------------------------------

def record_hologram(object_field, reference_field, tilt_angle_deg=5.0,
                    wavelength=633e-9, L=1e-3):
    """
    Simulate off-axis hologram recording.

    Parameters
    ----------
    tilt_angle_deg : float – reference beam tilt [degrees]

    Returns
    -------
    hologram : ndarray – recorded intensity (hologram)
    """
    N = object_field.shape[0]
    x, y, _, _ = make_grid(N, L)
    k = 2 * np.pi / wavelength
    theta = np.deg2rad(tilt_angle_deg)
    tilted_ref = reference_field * np.exp(1j * k * np.sin(theta) * x)
    hologram = np.abs(object_field + tilted_ref)**2
    return hologram


def reconstruct_hologram(hologram, wavelength=633e-9, z_rec=0.1, L=1e-3):
    """
    Reconstruct object field from off-axis hologram by back-propagation.

    Returns
    -------
    recon : complex ndarray
    """
    from .propagation import fresnel_propagate
    recon = fresnel_propagate(hologram.astype(complex), z_rec, wavelength, L)
    return recon
