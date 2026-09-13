"""
opticspy.qoptics.jaynes_cummings
===================================
The Jaynes-Cummings model: a single two-level atom coupled to a single
quantized cavity mode under the rotating-wave approximation. Builds the
full Hamiltonian on the (field dim) x (2-level atom) Hilbert space and
propagates the state numerically, exhibiting vacuum Rabi oscillations
and, for a coherent-state field, the collapse-and-revival phenomenon.

Basis convention: atom states ordered as (|g>, |e>) (ground, excited);
combined basis index = field_index * 2 + atom_index (atom_index=0 -> g,
1 -> e), i.e. field \\otimes atom via np.kron(field_op, atom_op).

Reference: Gerry & Knight, "Introductory Quantum Optics", Ch. 4
(The Jaynes-Cummings model) and Ch. 8 (vacuum Rabi, collapse & revival).
"""

import numpy as np
from scipy.linalg import expm

from .operators import annihilation_operator, creation_operator, identity_operator

SIGMA_MINUS = np.array([[0, 1], [0, 0]], dtype=complex)  # |g><e|  (e -> g)
SIGMA_PLUS = SIGMA_MINUS.conj().T                          # |e><g|  (g -> e)
SIGMA_Z = np.array([[-1, 0], [0, 1]], dtype=complex)       # +1 for |e>, -1 for |g>
ATOM_GROUND = np.array([1, 0], dtype=complex)
ATOM_EXCITED = np.array([0, 1], dtype=complex)


def hamiltonian(field_dim, omega_cavity, omega_atom, g, hbar=1.0):
    """
    Full (non-interaction-picture) Jaynes-Cummings Hamiltonian:
        H = hbar*omega_c * a^dagger a  (x) I_atom
          + hbar*omega_a/2 * I_field  (x) sigma_z
          + hbar*g * ( a (x) sigma_+ + a^dagger (x) sigma_- )
    """
    a = annihilation_operator(field_dim)
    adag = a.conj().T
    I_field = identity_operator(field_dim)
    I_atom = np.eye(2, dtype=complex)

    H_field = hbar * omega_cavity * np.kron(adag @ a, I_atom)
    H_atom = hbar * omega_atom / 2.0 * np.kron(I_field, SIGMA_Z)
    H_int = hbar * g * (np.kron(a, SIGMA_PLUS) + np.kron(adag, SIGMA_MINUS))
    return H_field + H_atom + H_int


def initial_state(field_state, atom_state):
    """Combined field (x) atom state vector (field_state normalized, atom_state in {g,e} basis)."""
    return np.kron(field_state, atom_state)


def evolve(H, psi0, tlist):
    """
    Propagate psi0 under the time-independent Hamiltonian H, returning
    an array of state vectors at each time in tlist (via matrix
    exponentials -- fine for the small Hilbert space dimensions used here).
    """
    dim = H.shape[0]
    psis = np.zeros((len(tlist), dim), dtype=complex)
    for i, t in enumerate(tlist):
        U = expm(-1j * H * t)
        psis[i] = U @ psi0
    return psis


def atomic_inversion(psis, field_dim):
    """
    <sigma_z>(t) = P_e(t) - P_g(t), reshaping each combined state vector
    into (field_dim, 2) to trace out the field.
    """
    W = np.zeros(len(psis))
    for i, psi in enumerate(psis):
        psi_r = psi.reshape(field_dim, 2)
        P_g = np.sum(np.abs(psi_r[:, 0]) ** 2)
        P_e = np.sum(np.abs(psi_r[:, 1]) ** 2)
        W[i] = P_e - P_g
    return W


def mean_photon_number(psis, field_dim):
    n = np.arange(field_dim)
    N = np.zeros(len(psis))
    for i, psi in enumerate(psis):
        psi_r = psi.reshape(field_dim, 2)
        P_n = np.abs(psi_r[:, 0]) ** 2 + np.abs(psi_r[:, 1]) ** 2
        N[i] = np.sum(n * P_n)
    return N


def vacuum_rabi_frequency(g, n=0):
    """Rabi frequency for the n-photon manifold on resonance: 2*g*sqrt(n+1)."""
    return 2 * g * np.sqrt(n + 1)


def dressed_state_energies(n, omega, g, hbar=1.0):
    """
    On-resonance (omega_cavity = omega_atom = omega) dressed-state
    (JC ladder) energies for manifold n>=0:
        E_n_+/- = hbar*omega*(n + 1/2) +/- hbar*g*sqrt(n+1)
    """
    E0 = hbar * omega * (n + 0.5)
    splitting = hbar * g * np.sqrt(n + 1)
    return E0 - splitting, E0 + splitting


def collapse_and_revival_time(g, nbar):
    """
    Approximate revival time for a coherent-state field with mean photon
    number nbar, on resonance: t_revival ~ 2*pi*sqrt(nbar) / g
    (standard estimate, Gerry & Knight Ch. 8).
    """
    return 2 * np.pi * np.sqrt(nbar) / g
