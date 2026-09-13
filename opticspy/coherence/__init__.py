"""
opticspy.coherence
====================
Coherence theory and statistical optics: temporal coherence
(Wiener-Khinchin theorem, coherence time/length, interferograms),
spatial coherence (Van Cittert-Zernike theorem, coherence radius,
extended-source visibility), and speckle statistics.

Reference: Goodman, "Statistical Optics"; Born & Wolf, "Principles of
Optics", Ch. 10.
"""

from . import temporal, spatial, speckle

__all__ = ["temporal", "spatial", "speckle"]
