"""
aoptics.fourier_optics.spatial_filtering
============================================
The coherent 4f imaging/spatial-filtering system: an input transparency
illuminated coherently, Fourier-transformed by a first lens onto a
Fourier plane where a mask can filter spatial frequencies, then
inverse-Fourier-transformed by a second lens back into an image.

Reference: Goodman, "Introduction to Fourier Optics", Ch. 7 (Frequency-
Domain Description of Optical Systems: spatial filtering); Saleh &
Teich, "Fundamentals of Photonics", Ch. 4.
"""

import numpy as np


def four_f_system(input_field, filter_mask=None):
    """
    Coherent 4f system: FT of the input (Fourier plane), optionally
    multiplied by `filter_mask` (same shape, complex- or real-valued),
    then inverse FT back to the image plane. Uses centered FFTs so the
    Fourier-plane array directly represents spatial frequency, with DC
    at the center -- the standard textbook idealization of a 4f system
    (unit magnification, scaling constants dropped).
    """
    F = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(input_field)))
    if filter_mask is not None:
        F = F * filter_mask
    image = np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(F)))
    return image, F


def low_pass_filter_mask(N, cutoff_fraction):
    """Circular low-pass mask, `cutoff_fraction` in (0,1] of the Nyquist/array half-width."""
    x = np.linspace(-1, 1, N)
    X, Y = np.meshgrid(x, x)
    R = np.sqrt(X ** 2 + Y ** 2)
    return (R <= cutoff_fraction).astype(complex)


def high_pass_filter_mask(N, cutoff_fraction):
    return 1.0 - low_pass_filter_mask(N, cutoff_fraction)


def vertical_slit_filter_mask(N, half_width_fraction):
    """Pass only a vertical strip of spatial frequencies (blocks horizontal
    frequency content) -- a classic spatial-filtering demonstration
    (e.g. removing horizontal scan-line noise)."""
    x = np.linspace(-1, 1, N)
    X, Y = np.meshgrid(x, x)
    return (np.abs(X) <= half_width_fraction).astype(complex)


def phase_contrast_mask(N, phase_shift=np.pi / 2, dc_radius_fraction=0.02):
    """
    Zernike phase-contrast filter: applies a phase shift `phase_shift`
    to only the low-frequency (DC) region of radius `dc_radius_fraction`,
    converting a pure-phase object into visible intensity contrast --
    the classic phase-contrast microscopy trick.
    """
    x = np.linspace(-1, 1, N)
    X, Y = np.meshgrid(x, x)
    R = np.sqrt(X ** 2 + Y ** 2)
    mask = np.ones((N, N), dtype=complex)
    mask[R <= dc_radius_fraction] = np.exp(1j * phase_shift)
    return mask

