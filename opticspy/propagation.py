"""
opticspy.propagation - Beam propagation methods
"""
import numpy as np
from .beams import make_grid


def angular_spectrum(field, z, wavelength=633e-9, L=1e-3):
    N = field.shape[0]; dx = 2*L/N; k = 2*np.pi/wavelength
    fx = np.fft.fftfreq(N, d=dx); FX, FY = np.meshgrid(fx, fx)
    kz_sq = k**2 - (2*np.pi*FX)**2 - (2*np.pi*FY)**2
    mask = kz_sq >= 0
    kz = np.where(mask, np.sqrt(np.maximum(kz_sq, 0)), 0)
    H = np.where(mask, np.exp(1j*kz*z), 0)
    return np.fft.ifft2(np.fft.fft2(field)*H)


def fresnel_propagate(field, z, wavelength=633e-9, L=1e-3):
    N = field.shape[0]; dx = 2*L/N; k = 2*np.pi/wavelength
    fx = np.fft.fftfreq(N, d=dx); FX, FY = np.meshgrid(fx, fx)
    H = np.exp(1j*k*z)*np.exp(-1j*np.pi*wavelength*z*(FX**2+FY**2))
    return np.fft.ifft2(np.fft.fft2(field)*H)


def fraunhofer(field, wavelength=633e-9, f=0.1, L=1e-3):
    N = field.shape[0]; dx = 2*L/N
    field_far = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(field))) * dx**2/(1j*wavelength*f)
    du = wavelength*f/(N*dx)
    u = (np.arange(N)-N//2)*du; x_far, y_far = np.meshgrid(u, u)
    return field_far, x_far, y_far


def abcd_gaussian(w0, z0, M, wavelength=633e-9):
    A,B,C,D = M[0][0], M[0][1], M[1][0], M[1][1]
    zR = np.pi*w0**2/wavelength
    q_in = z0 + 1j*zR
    q_out = (A*q_in + B)/(C*q_in + D)
    w_out = np.sqrt(-wavelength/(np.pi*np.imag(1/q_out)))
    return w_out, q_out


def talbot_carpet(grating_field, wavelength=633e-9, period=100e-6, num_planes=50, L=1e-3):
    z_T = 2*period**2/wavelength
    z_vals = np.linspace(0, 2*z_T, num_planes)
    N = grating_field.shape[0]
    carpet = np.zeros((num_planes, N))
    for i, z in enumerate(z_vals):
        prop = grating_field if z==0 else fresnel_propagate(grating_field, z, wavelength, L)
        carpet[i, :] = np.abs(prop[N//2, :])**2
    return carpet, z_vals


def bpm_propagate(field, n_profile, dz, steps, wavelength=633e-9, L=1e-3):
    N = field.shape[0]; dx = 2*L/N; k0 = 2*np.pi/wavelength
    fx = np.fft.fftfreq(N, d=dx); FX, FY = np.meshgrid(fx, fx)
    H_free = np.exp(-1j*np.pi*wavelength*dz*(FX**2+FY**2))
    fields = [field.copy()]; current = field.copy()
    for s in range(steps):
        n = n_profile(s*dz) if callable(n_profile) else n_profile
        current = current*np.exp(1j*k0*(n-1.0)*dz)
        current = np.fft.ifft2(np.fft.fft2(current)*H_free)
        fields.append(current.copy())
    return fields
