"""
aoptics.modulation
======================
Electro-optic (Pockels effect, phase/amplitude modulators) and
acousto-optic (Raman-Nath/Bragg diffraction) modulation of light.

Reference: Saleh & Teich, "Fundamentals of Photonics", Ch. 19-20.
"""

from . import electrooptics, acoustooptics

__all__ = ["electrooptics", "acoustooptics"]

