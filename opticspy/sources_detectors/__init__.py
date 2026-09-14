"""
aoptics.sources_detectors
=============================
Laser rate equations (threshold, gain clamping, output power) and
photodetector figures of merit (responsivity, shot/thermal noise, NEP).

Reference: Saleh & Teich, "Fundamentals of Photonics", Ch. 15-16
(Laser amplifiers / Lasers) and Ch. 18 (Photodetectors).
"""

from . import lasers, photodetectors

__all__ = ["lasers", "photodetectors"]

