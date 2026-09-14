"""
aoptics.photonics
=====================
Guided-wave and nonlinear photonics: planar waveguide modes, optical
fiber parameters and dispersion, laser resonators (Fabry-Perot and
Gaussian-mode ABCD analysis), and chi(2)/chi(3) nonlinear optics.

Reference: Saleh & Teich, "Fundamentals of Photonics", Chapters 7-9
(Beam/Guided-Wave/Fiber Optics) and Ch. 21 (Nonlinear Optics).
"""

from . import waveguides, fibers, resonators, nonlinear_optics, photonic_crystals

__all__ = ["waveguides", "fibers", "resonators", "nonlinear_optics", "photonic_crystals"]

