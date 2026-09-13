"""
opticspy.qoptics.operators
=============================
Field operators in a truncated (dimension-N) single-mode Fock basis:
annihilation/creation operators, the number operator, quadratures, and
the displacement/squeeze operators built by exponentiating their
generators (matching how coherent/squeezed states are defined in
opticspy.quantum).

Reference: Gerry & Knight, "Introductory Quantum Optics", Ch. 2-3
(quantization of the field, quadrature operators, displaced and
squeezed states).
"""

import numpy as np
from scipy.linalg import expm


def annihilation_operator(dim):
    """a |n> = sqrt(n) |n-1>, in a truncated N-dimensional Fock basis."""
    n = np.arange(1, dim)
    return np.diag(np.sqrt(n), k=1)


def creation_operator(dim):
    return annihilation_operator(dim).conj().T


def number_operator(dim):
    return np.diag(np.arange(dim, dtype=complex))


def identity_operator(dim):
    return np.eye(dim, dtype=complex)


def quadrature_x(dim):
    """X = (a + a^dagger)/sqrt(2) (dimensionless quadrature)."""
    a = annihilation_operator(dim)
    return (a + a.conj().T) / np.sqrt(2)


def quadrature_p(dim):
    """P = (a - a^dagger)/(i*sqrt(2))."""
    a = annihilation_operator(dim)
    return (a - a.conj().T) / (1j * np.sqrt(2))


def commutator(A, B):
    return A @ B - B @ A


def displacement_operator(alpha, dim):
    """D(alpha) = exp(alpha a^dagger - alpha* a), the displacement operator
    (D(alpha)|0> is the coherent state |alpha>)."""
    a = annihilation_operator(dim)
    adag = a.conj().T
    return expm(alpha * adag - np.conj(alpha) * a)


def squeeze_operator(xi, dim):
    """S(xi) = exp( (xi* a^2 - xi a^dagger^2)/2 ), xi = r*exp(i*phi)."""
    a = annihilation_operator(dim)
    adag = a.conj().T
    return expm(0.5 * (np.conj(xi) * a @ a - xi * adag @ adag))


def expectation(op, state):
    """<state|op|state> for a state vector, or Tr(op*rho) for a density matrix."""
    state = np.asarray(state)
    if state.ndim == 1:
        return complex(np.conj(state) @ op @ state)
    return complex(np.trace(op @ state))


def variance(op, state):
    return expectation(op @ op, state) - expectation(op, state) ** 2


def heisenberg_uncertainty_product(state, dim):
    """<Delta X^2><Delta P^2>, should be >= 1/4 for any physical state
    (with the X,P normalization used here, the vacuum saturates 1/4)."""
    X, P = quadrature_x(dim), quadrature_p(dim)
    varX = np.real(variance(X, state))
    varP = np.real(variance(P, state))
    return varX * varP
