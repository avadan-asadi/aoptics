"""
aoptics – کتابخانه جامع پایتون برای اپتیک و کوانتوم اپتیک
=============================================================
"""
__version__ = "2.1.0"
__author__ = "aoptics Team"
__license__ = "MIT"

from . import beams, propagation, interference, diffraction, wavefront, quantum, utils
from . import geometrical, electromagnetic, coherence, photonics, qoptics
from . import fourier_optics, sources_detectors, modulation, scattering
from . import nonlinear_optics_boyd

__all__ = ["beams","propagation","interference","diffraction","wavefront","quantum","utils",
           "geometrical","electromagnetic","coherence","photonics","qoptics",
           "fourier_optics","sources_detectors","modulation","scattering",
           "nonlinear_optics_boyd"]

