"""
aoptics.fourier_optics
==========================
Fourier-optics description of imaging systems: pupil function,
coherent/incoherent PSF, coherent transfer function (CTF), optical/
modulation transfer function (OTF/MTF), classic resolution criteria,
and the coherent 4f spatial-filtering system.

Reference: Goodman, "Introduction to Fourier Optics"; Saleh & Teich,
"Fundamentals of Photonics", Ch. 4.
"""

from . import transfer_functions, spatial_filtering

__all__ = ["transfer_functions", "spatial_filtering"]

