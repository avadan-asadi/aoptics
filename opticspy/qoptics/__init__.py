"""
opticspy.qoptics
===================
Advanced quantum optics: field operators in a truncated Fock basis,
the Jaynes-Cummings model (vacuum Rabi oscillations, collapse &
revival), the Lindblad master equation (cavity decay, spontaneous
emission), the quantum beamsplitter and Hong-Ou-Mandel interference,
quantum coherence functions (g^(2)), and entanglement/quantum-
information tools -- extending `opticspy.quantum`.

Reference: Gerry & Knight, "Introductory Quantum Optics" (full text).
"""

from . import operators, jaynes_cummings, master_equation, beamsplitter
from . import coherence_functions, entanglement, quantum_information

__all__ = ["operators", "jaynes_cummings", "master_equation", "beamsplitter",
           "coherence_functions", "entanglement", "quantum_information"]
