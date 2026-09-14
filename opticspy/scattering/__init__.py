"""
aoptics.scattering
======================
Scattering from small particles (Rayleigh regime) and reflection from
absorbing media (optics of metals, via a complex refractive index).

Reference: Born & Wolf, "Principles of Optics", Ch. 14.
"""

from . import rayleigh_scattering, metal_optics

__all__ = ["rayleigh_scattering", "metal_optics"]

