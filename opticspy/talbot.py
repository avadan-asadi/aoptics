"""
talbot.py - پدیده تالبوت و خود-تصویرسازی
Talbot effect and fractional Talbot carpets.
"""
import numpy as np


def talbot_distance(pitch, wavelength):
    """Compute Talbot distance z_T = 2*d^2/lambda.
    
    Parameters
    ----------
    pitch      : grating period [m]
    wavelength : [m]
    
    Returns
    -------
    z_T : Talbot distance [m]
    """
    return 2 * pitch**2 / wavelength


def talbot_carpet(pitch=10e-6, wavelength=532e-9, N=512, z_periods=2,
                  duty_cycle=0.5, N_periods=10):
    """Simulate Talbot carpet (longitudinal cross-section).
    
    Parameters
    ----------
    pitch       : grating period [m]
    wavelength  : [m]
    N           : grid points
    z_periods   : number of Talbot distances to simulate
    duty_cycle  : grating fill factor (0-1)
    N_periods   : number of grating periods
    
    Returns
    -------
    carpet : 2D intensity array (z vs x)
    x      : transverse coordinate [m]
    z      : longitudinal coordinate [m]
    z_T    : Talbot distance [m]
    """
    z_T = talbot_distance(pitch, wavelength)
    
    L = N_periods * pitch
    x = np.linspace(-L/2, L/2, N)
    z_max = z_periods * z_T
    z = np.linspace(1e-6, z_max, N)
    
    # Create grating field (binary amplitude)
    grating = (np.mod(x, pitch) < duty_cycle * pitch).astype(complex)
    
    carpet = np.zeros((N, N), dtype=float)
    dx = x[1] - x[0]
    k = 2 * np.pi / wavelength
    
    # FFT-based Fresnel propagation for each z slice
    G = np.fft.fft(grating)
    fx = np.fft.fftfreq(N, d=dx)
    
    for i, zi in enumerate(z):
        # Transfer function (paraxial)
        H = np.exp(-1j * np.pi * wavelength * zi * (fx**2))
        E = np.fft.ifft(G * H)
        carpet[i, :] = np.abs(E)**2
    
    return carpet, x, z, z_T


def fractional_talbot(x, z, pitch, wavelength, order_p=1, order_q=2):
    """Analytical fractional Talbot field at z = (p/q)*z_T."""
    z_T = talbot_distance(pitch, wavelength)
    z_frac = (order_p / order_q) * z_T
    
    k = 2 * np.pi / wavelength
    # Sum over Fourier components
    E = np.zeros_like(x, dtype=complex)
    N_terms = 20
    for m in range(-N_terms, N_terms+1):
        phase_grating = 2 * np.pi * m * x / pitch
        phase_talbot  = -np.pi * wavelength * z_frac * m**2 / pitch**2
        E += np.exp(1j * (phase_grating + phase_talbot))
    return E
