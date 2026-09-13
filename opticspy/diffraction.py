"""
opticspy.diffraction
====================
Diffraction gratings, Talbot effect, Moire patterns.
"""

import numpy as np
from .beams import make_grid
from .propagation import fresnel_propagate


def grating_1d(N=512, L=5e-3, period=500e-6, gtype='binary_amp',
               orientation=0.0, duty_cycle=0.5, phase_depth=np.pi):
    x, y, _, _ = make_grid(N, L)
    xi = x * np.cos(orientation) + y * np.sin(orientation)
    u = (xi % period) / period

    if gtype == 'binary_amp':
        return (u < duty_cycle).astype(complex)
    elif gtype == 'binary_phase':
        return np.where(u < duty_cycle, np.exp(1j * phase_depth), 1.0+0j)
    elif gtype == 'sinusoidal':
        return 0.5 * (1 + np.cos(2 * np.pi * xi / period)) + 0j
    elif gtype == 'blazed':
        return np.exp(1j * phase_depth * u)
    elif gtype == 'triangular':
        return np.exp(1j * phase_depth * 2 * np.abs(u - 0.5))
    else:
        raise ValueError(f"Unknown grating type: {gtype}")


def grating_2d(N=512, L=5e-3, period_x=500e-6, period_y=500e-6, gtype='binary_amp'):
    g1 = grating_1d(N, L, period_x, gtype, 0)
    g2 = grating_1d(N, L, period_y, gtype, np.pi/2)
    return g1 * g2


def talbot_carpet(N=512, L=2e-3, period=200e-6, wavelength=633e-9,
                  gtype='binary_amp', num_z=200):
    grating = grating_1d(N, L, period, gtype)
    z_T = 2 * period**2 / wavelength
    z_vals = np.linspace(0, 2 * z_T, num_z)
    z_norm = z_vals / z_T
    carpet = np.zeros((num_z, N))
    for i, z in enumerate(z_vals):
        if z == 0:
            prop = grating
        else:
            prop = fresnel_propagate(grating, z, wavelength, L)
        carpet[i, :] = np.abs(prop[N // 2, :])**2
    return carpet, z_norm, z_vals


def moire_pattern(N=512, L=5e-3, period1=300e-6, period2=310e-6,
                  angle1=0.0, angle2=0.05, gtype='sinusoidal'):
    g1 = np.abs(grating_1d(N, L, period1, gtype, angle1))**2
    g2 = np.abs(grating_1d(N, L, period2, gtype, angle2))**2
    return g1 * g2, g1, g2


def fresnel_zone_plate(N=512, L=5e-3, wavelength=633e-9, focal_length=0.1,
                       zone_type='amplitude'):
    _, _, r, _ = make_grid(N, L)
    zone_index = np.floor(r**2 / (wavelength * focal_length)).astype(int)
    is_odd = (zone_index % 2 == 1)
    if zone_type == 'amplitude':
        return np.where(is_odd, 0.0+0j, 1.0+0j)
    else:
        return np.where(is_odd, np.exp(1j * np.pi), 1.0+0j)


def circular_aperture_diffraction(radius=0.5e-3, N=512, L=5e-3, wavelength=633e-9, z=0.5):
    from .beams import plane_wave
    pw = plane_wave(N, L, wavelength)
    _, _, r, _ = make_grid(N, L)
    aperture = (r <= radius).astype(complex)
    field_out = fresnel_propagate(pw * aperture, z, wavelength, L)
    return np.abs(field_out)**2, field_out


def grating_efficiency(grating, N_orders=5):
    g1d = grating[grating.shape[0]//2, :] if grating.ndim == 2 else grating
    spectrum = np.fft.fft(g1d)
    power_total = np.sum(np.abs(spectrum)**2) + 1e-30
    N = len(g1d)
    orders = list(range(-N_orders, N_orders+1))
    efficiencies = [np.abs(spectrum[m % N])**2 / power_total for m in orders]
    return orders, efficiencies
