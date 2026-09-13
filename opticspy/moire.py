"""
moire.py - شبیه سازی الگوهای ماره
Moiré pattern simulation with various grating types.
"""
import numpy as np


def linear_grating(x, y, pitch=1e-3, angle=0.0, phase=0.0, grating_type='cosine'):
    """Single linear grating transmission function.
    
    Parameters
    ----------
    pitch        : grating period [m]
    angle        : orientation angle [rad]
    phase        : phase offset [rad]
    grating_type : 'cosine', 'binary', 'sinusoidal', 'blazed'
    
    Returns
    -------
    T : transmission (0 to 1)
    """
    # Rotated coordinate
    u = x * np.cos(angle) + y * np.sin(angle)
    xi = 2 * np.pi * u / pitch + phase
    
    if grating_type == 'cosine' or grating_type == 'sinusoidal':
        T = 0.5 * (1 + np.cos(xi))
    elif grating_type == 'binary':
        T = (np.mod(u, pitch) < pitch/2).astype(float)
    elif grating_type == 'blazed':
        T = np.mod(u / pitch, 1.0)
    else:
        T = 0.5 * (1 + np.cos(xi))
    return T


def circular_grating(x, y, pitch=1e-3, center=(0,0), phase=0.0):
    """Circular (zone-plate style) grating."""
    cx, cy = center
    r = np.sqrt((x - cx)**2 + (y - cy)**2)
    xi = 2 * np.pi * r / pitch + phase
    return 0.5 * (1 + np.cos(xi))


def radial_grating(x, y, N_spokes=20, phase=0.0):
    """Radial (spoke) grating."""
    phi = np.arctan2(y, x)
    xi = N_spokes * phi + phase
    return 0.5 * (1 + np.cos(xi))


def moire_pattern(x, y, pitch1=1e-3, angle1=0.0, pitch2=None, angle2=None,
                   phase1=0.0, phase2=0.0, grating_type='cosine'):
    """Two-grating Moiré pattern.
    
    Parameters
    ----------
    pitch1/2 : grating periods [m] (pitch2 defaults to pitch1)
    angle1/2 : orientation angles [rad]
    
    Returns
    -------
    I_moire    : Moiré intensity pattern
    T1, T2     : individual grating transmissions
    """
    if pitch2 is None:
        pitch2 = pitch1
    if angle2 is None:
        angle2 = angle1 + np.radians(5)  # default 5-degree rotation
    
    T1 = linear_grating(x, y, pitch1, angle1, phase1, grating_type)
    T2 = linear_grating(x, y, pitch2, angle2, phase2, grating_type)
    I_moire = T1 * T2
    return I_moire, T1, T2


def moire_from_circular(x, y, pitch1=1e-3, pitch2=None, phase_diff=0.0):
    """Moiré from two concentric circular gratings."""
    if pitch2 is None:
        pitch2 = pitch1 * 1.1
    T1 = circular_grating(x, y, pitch1)
    T2 = circular_grating(x, y, pitch2, phase=phase_diff)
    return T1 * T2, T1, T2


def moire_period(pitch1, pitch2, angle_diff):
    """Analytical Moiré period for two gratings.
    
    Parameters
    ----------
    pitch1, pitch2 : grating periods [m]
    angle_diff     : angle between gratings [rad]
    
    Returns
    -------
    Lambda_moire : Moiré period [m]
    """
    # Using wave vector approach
    k1 = 1 / pitch1
    k2 = 1 / pitch2
    kx1 = k1
    ky1 = 0
    kx2 = k2 * np.cos(angle_diff)
    ky2 = k2 * np.sin(angle_diff)
    dkx = kx1 - kx2
    dky = ky1 - ky2
    k_moire = np.sqrt(dkx**2 + dky**2)
    return 1 / k_moire if k_moire != 0 else np.inf
