"""
opticspy.geometrical
=====================
Geometrical (ray) optics: paraxial ray tracing via ABCD matrices, cardinal
points of optical systems, real (skew-free, meridional) ray tracing through
spherical/conic surfaces via Snell's law, and third-order (Seidel)
aberration theory.

References
----------
Born & Wolf, "Principles of Optics", 7th ed., Chapters 3-5.
Saleh & Teich, "Fundamentals of Photonics", 3rd ed., Chapter 1 (Ray Optics)
and Chapter 6 (Polarization Optics, geometry parts).
"""

from . import ray_tracing, aberrations

__all__ = ["ray_tracing", "aberrations"]
