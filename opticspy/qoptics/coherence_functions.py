"""
aoptics.qoptics.coherence_functions
=======================================
Quantum degree-of-second-order-coherence g^(2)(0) and related photon-
statistics measures for a single-mode field state, computed directly
from the photon-number distribution (works for any state -- Fock,
coherent, thermal, squeezed, cat -- given as a state vector or density
matrix in the Fock basis).

Reference: Gerry & Knight, "Introductory Quantum Optics", Ch. 3 (Photon
statistics: g^(2), photon antibunching) and Ch. 5 (correlation
functions).
"""

import numpy as np


def _photon_number_distribution(state_or_rho):
    if state_or_rho.ndim == 1:
        return np.abs(state_or_rho) ** 2
    return np.real(np.diag(state_or_rho))


def g2_zero_delay(state_or_rho):
    """
    g^(2)(0) = <a^dagger a^dagger a a> / <a^dagger a>^2
             = sum_n n(n-1) P_n / (sum_n n P_n)^2

    g2 < 1: photon antibunching (nonclassical, e.g. a Fock state);
    g2 = 1: coherent state (Poissonian, classical-like);
    g2 > 1: photon bunching (e.g. thermal/chaotic light).
    """
    P_n = _photon_number_distribution(state_or_rho)
    n = np.arange(len(P_n))
    mean_n = np.sum(n * P_n)
    mean_n_n1 = np.sum(n * (n - 1) * P_n)
    if mean_n == 0:
        return np.nan
    return float(mean_n_n1 / mean_n ** 2)


def is_antibunched(state_or_rho):
    return g2_zero_delay(state_or_rho) < 1.0


def fano_factor(state_or_rho):
    """Fano factor F = Var(n)/<n> (F=1 Poissonian/coherent, F<1 sub-Poissonian, F>1 super-Poissonian)."""
    P_n = _photon_number_distribution(state_or_rho)
    n = np.arange(len(P_n))
    mean_n = np.sum(n * P_n)
    var_n = np.sum(n ** 2 * P_n) - mean_n ** 2
    return float(var_n / mean_n) if mean_n > 0 else np.nan

