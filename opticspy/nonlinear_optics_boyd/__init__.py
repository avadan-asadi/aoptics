"""
opticspy.nonlinear_optics_boyd
=================================
Nonlinear optics following the chapter structure of Robert W. Boyd's
"Nonlinear Optics": the nonlinear susceptibility and Miller's rule,
coupled three-wave-mixing amplitude equations and the Manley-Rowe
relations, the saturated two-level-atom response, self-action effects
(intensity-dependent refractive index, self-focusing, self-phase
modulation), stimulated Raman/Brillouin scattering, and two-photon
absorption / saturable absorption.

This complements (and, for chi(2)/chi(3) basics, partially overlaps
with) `opticspy.photonics.nonlinear_optics`, which follows Saleh &
Teich's treatment instead; the two are kept separate since they follow
different textbooks' notation and scope.

Reference: Boyd, "Nonlinear Optics", 3rd/4th ed. (full text).
"""

from . import (susceptibility, coupled_wave_mixing, two_level_atom,
               self_action, stimulated_scattering, multiphoton_absorption)

__all__ = ["susceptibility", "coupled_wave_mixing", "two_level_atom",
           "self_action", "stimulated_scattering", "multiphoton_absorption"]
