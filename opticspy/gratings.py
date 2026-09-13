"""
gratings.py - انواع توری های اپتیکی
Grating types, efficiencies, and spectral analysis.
"""
import numpy as np


def grating_equation(wavelength, pitch, m, n_in=1.0, angle_in=0.0):
    """Grating equation: n*sin(theta_m) = sin(theta_in) + m*lambda/d.
    
    Parameters
    ----------
    wavelength : [m]
    pitch      : grating period [m]
    m          : diffraction order
    n_in       : refractive index of incident medium
    angle_in   : angle of incidence [rad]
    
    Returns
    -------
    angle_m : diffraction angle [rad] or None if evanescent
    """
    sin_m = n_in * np.sin(angle_in) + m * wavelength / pitch
    if abs(sin_m) > 1.0:
        return None   # evanescent order
    return np.arcsin(sin_m)


def blaze_condition(m, blaze_angle, wavelength_blaze, pitch):
    """Check blaze condition for given order."""
    angle_m = grating_equation(wavelength_blaze, pitch, m)
    if angle_m is None:
        return False
    # Blaze when reflection angle = blaze angle
    return abs(angle_m - 2*blaze_angle) < 1e-10


def sinusoidal_grating_field(x, y, pitch=10e-6, angle=0.0, depth=1.0,
                              phase=0.0, grating_type='amplitude'):
    """Generate sinusoidal grating complex field.
    
    Parameters
    ----------
    depth        : modulation depth (0-1 for amplitude, 0-2pi for phase)
    grating_type : 'amplitude', 'phase', 'complex'
    """
    u = x * np.cos(angle) + y * np.sin(angle)
    xi = 2 * np.pi * u / pitch + phase
    
    if grating_type == 'amplitude':
        T = 0.5 * (1 + depth * np.cos(xi))
        return T.astype(complex)
    elif grating_type == 'phase':
        return np.exp(1j * depth * np.cos(xi))
    else:  # complex
        A = 0.5 * (1 + depth * np.cos(xi))
        phi = depth * np.sin(xi)
        return A * np.exp(1j * phi)


def zone_plate(x, y, f=0.5, wavelength=632.8e-9, N_zones=50,
               zone_type='binary'):
    """Fresnel zone plate / diffractive lens.
    
    Parameters
    ----------
    f          : focal length [m]
    N_zones    : number of Fresnel zones
    zone_type  : 'binary', 'phase', 'continuous'
    """
    r2 = x**2 + y**2
    # Zone plate phase
    phi = np.pi * r2 / (wavelength * f)
    
    if zone_type == 'binary':
        T = (np.mod(phi, 2*np.pi) < np.pi).astype(float)
        return T
    elif zone_type == 'phase':
        T = np.exp(-1j * phi)
        return T
    else:  # continuous (lens-like)
        return np.exp(-1j * np.pi * r2 / (wavelength * f))


def echelle_grating_efficiency(m, blaze_angle, wavelength, pitch):
    """Approximate scalar efficiency of an echelle grating order."""
    # Littrow mount efficiency
    theta_B = blaze_angle
    lambda_B = 2 * pitch * np.sin(theta_B) / m  # blaze wavelength
    delta_lambda = wavelength - lambda_B
    # Sinc-squared approximation
    u = np.pi * m * delta_lambda / lambda_B
    eta = np.sinc(u / np.pi)**2
    return eta


def transmission_grating(x, y, pitch=10e-6, duty_cycle=0.5, angle=0.0):
    """Binary amplitude transmission grating."""
    u = x * np.cos(angle) + y * np.sin(angle)
    return (np.mod(u, pitch) < duty_cycle * pitch).astype(complex)


def volume_grating_coupled_wave(Kg, K, n0, n1, wavelength, d, angle_in):
    """Coupled-wave analysis (Kogelnik) for volume holographic grating.
    
    Parameters
    ----------
    Kg  : grating vector magnitude [m^-1]
    K   : incident wave vector [m^-1]
    n0  : average refractive index
    n1  : index modulation amplitude
    d   : grating thickness [m]
    angle_in : Bragg angle [rad]
    
    Returns
    -------
    eta : diffraction efficiency
    """
    nu = np.pi * n1 * d / (wavelength * np.cos(angle_in))
    xi = 0.0  # assume perfect Bragg condition
    eta = np.sin(np.sqrt(nu**2 + xi**2))**2 / (1 + (xi/nu)**2)
    return eta
