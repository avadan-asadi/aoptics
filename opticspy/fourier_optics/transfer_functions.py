"""
aoptics.fourier_optics.transfer_functions
=============================================
Diffraction-limited imaging in Fourier-optics language: the pupil
function, coherent/incoherent point-spread functions, the coherent
transfer function (CTF, simply the scaled pupil function) and the
optical transfer function (OTF, the normalized autocorrelation of the
pupil function) with its modulus the modulation transfer function
(MTF), plus classic resolution criteria.

Reference: Goodman, "Introduction to Fourier Optics", Ch. 6 (Frequency
Analysis of Optical Imaging Systems); Saleh & Teich, "Fundamentals of
Photonics", Ch. 4 (Fourier Optics).
"""

import numpy as np


def pupil_function_circular(N, L, aperture_radius):
    """Binary circular pupil function on an N x N grid spanning [-L, L]."""
    x = np.linspace(-L, L, N)
    X, Y = np.meshgrid(x, x)
    R = np.sqrt(X ** 2 + Y ** 2)
    return (R <= aperture_radius).astype(complex)


def coherent_psf(pupil):
    """
    Coherent amplitude point-spread function h(x,y): the Fourier
    transform of the pupil function (Fraunhofer/focal-plane relation),
    normalized to unit peak.
    """
    H = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(pupil)))
    return H / np.max(np.abs(H))


def incoherent_psf(pupil):
    """Incoherent PSF |h(x,y)|^2, normalized to unit peak (unit volume left to the caller)."""
    h = coherent_psf(pupil)
    psf = np.abs(h) ** 2
    return psf / np.max(psf)


def coherent_transfer_function(pupil):
    """
    The coherent transfer function H(fx,fy) IS the (suitably rescaled)
    pupil function itself -- a foundational Fourier-optics result. This
    function simply returns the pupil array as the CTF for clarity of use.
    """
    return pupil


def optical_transfer_function(pupil):
    """
    OTF(fx,fy) = autocorrelation of the pupil function, normalized to 1
    at zero spatial frequency -- computed via the FFT-based correlation
    theorem (avoids assuming a specific closed-form pupil shape).
    """
    P = pupil
    F = np.fft.fft2(P)
    autocorr = np.fft.ifft2(F * np.conj(F))
    autocorr = np.fft.fftshift(autocorr)
    return autocorr / np.max(np.abs(autocorr))


def modulation_transfer_function(pupil):
    return np.abs(optical_transfer_function(pupil))


def diffraction_limited_cutoff_coherent(NA, wavelength):
    """Coherent-imaging cutoff spatial frequency: f_c = NA / wavelength."""
    return NA / wavelength


def diffraction_limited_cutoff_incoherent(NA, wavelength):
    """Incoherent-imaging (intensity) cutoff spatial frequency: f_c = 2*NA / wavelength
    (twice the coherent cutoff, since the OTF is the pupil autocorrelation)."""
    return 2 * NA / wavelength


def rayleigh_resolution(wavelength, NA):
    """Rayleigh two-point resolution criterion: 0.61 * wavelength / NA."""
    return 0.61 * wavelength / NA


def abbe_resolution(wavelength, NA):
    """Abbe diffraction limit: wavelength / (2 * NA)."""
    return wavelength / (2 * NA)

