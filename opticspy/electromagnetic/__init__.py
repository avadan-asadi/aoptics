"""
opticspy.electromagnetic
=========================
Electromagnetic/vector aspects of light: Fresnel reflection/transmission
at interfaces, polarization optics (Jones and Stokes/Mueller calculus),
and uniaxial crystal optics (birefringence, walk-off).

Reference: Born & Wolf, "Principles of Optics", Chapters 1 and 14;
Saleh & Teich, "Fundamentals of Photonics", Chapter 6.
"""

from . import fresnel, polarization, crystal_optics

__all__ = ["fresnel", "polarization", "crystal_optics"]
