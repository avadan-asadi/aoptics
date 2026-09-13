"""
opticspy.qoptics.beamsplitter
================================
The quantum beamsplitter acting on two truncated-Fock-basis modes, built
as the unitary exp[theta*(a^dagger b - a b^dagger)] (a Schwinger/SU(2)
beamsplitter generator) on the combined (dim x dim) two-mode Hilbert
space, with theta=pi/4 the symmetric (50/50) beamsplitter. Used here to
demonstrate the Hong-Ou-Mandel effect: two indistinguishable single
photons entering the two input ports of a balanced beamsplitter always
leave together (photon bunching), giving zero coincidence probability
at the outputs -- a genuinely quantum (two-photon interference) effect
with no classical analogue.

Reference: Gerry & Knight, "Introductory Quantum Optics", Ch. 6
(Beam splitters and interferometers) and the Hong-Ou-Mandel discussion
therein.
"""

import numpy as np
from scipy.linalg import expm

from .operators import annihilation_operator, identity_operator


def _two_mode_ops(dim):
    a = annihilation_operator(dim)
    I = identity_operator(dim)
    A = np.kron(a, I)
    B = np.kron(I, a)
    return A, B


def beamsplitter_unitary(theta, dim, phi=0.0):
    """
    Beamsplitter unitary on the combined two-mode Fock space (dimension
    dim*dim): U = exp[theta*(a^dagger b e^{-i phi} - a b^dagger e^{i phi})].
    theta=pi/4 -> 50/50 (symmetric) beamsplitter; theta=0 -> identity.
    """
    A, B = _two_mode_ops(dim)
    Adag, Bdag = A.conj().T, B.conj().T
    G = theta * (Adag @ B * np.exp(-1j * phi) - A @ Bdag * np.exp(1j * phi))
    return expm(G)


def two_mode_fock_state(n_a, n_b, dim):
    """|n_a, n_b> as a vector in the dim*dim combined Fock basis."""
    from .operators import identity_operator as _I  # local import, avoid unused warning
    va = np.zeros(dim, dtype=complex); va[n_a] = 1.0
    vb = np.zeros(dim, dtype=complex); vb[n_b] = 1.0
    return np.kron(va, vb)


def output_probabilities(U, psi_in, dim):
    """Return a (dim,dim) array of |<n_a,n_b|U|psi_in>|^2 output Fock probabilities."""
    psi_out = U @ psi_in
    probs = np.abs(psi_out.reshape(dim, dim)) ** 2
    return probs


def splitting_ratio(theta, dim=8, alpha=1.0):
    """
    Empirically verify the beamsplitter's power transmission/reflection
    by sending a coherent state into port a (vacuum into port b) and
    measuring the output mean photon numbers -- self-validates the
    T=cos^2(theta), R=sin^2(theta) convention numerically rather than
    assuming it.
    """
    from opticspy.quantum import coherent_state
    U = beamsplitter_unitary(theta, dim)
    psi_a = coherent_state(alpha, dim)
    vac = np.zeros(dim, dtype=complex); vac[0] = 1.0
    psi_in = np.kron(psi_a, vac)
    probs = output_probabilities(U, psi_in, dim)
    n = np.arange(dim)
    n_a_out = np.sum(probs.sum(axis=1) * n)
    n_b_out = np.sum(probs.sum(axis=0) * n)
    return n_a_out, n_b_out


def hong_ou_mandel_coincidence_probability(dim=6, theta=np.pi / 4):
    """
    Probability of a |1,1> coincidence at the OUTPUT of a beamsplitter
    of angle `theta`, given a |1,1> input (one indistinguishable photon
    in each input port). For the balanced case theta=pi/4 this should
    vanish (the Hong-Ou-Mandel dip): the two photons always bunch into
    |2,0> or |0,2>.
    """
    U = beamsplitter_unitary(theta, dim)
    psi_in = two_mode_fock_state(1, 1, dim)
    probs = output_probabilities(U, psi_in, dim)
    return probs[1, 1]


def hong_ou_mandel_visibility_vs_distinguishability(indistinguishability, dim=6):
    """
    Simple two-path (partially distinguishable photons) HOM model:
    mixes the fully-indistinguishable balanced-beamsplitter result
    (zero |1,1> coincidence) with the fully-distinguishable classical
    result (25% |1,1> coincidence from independent binomial splitting)
    in proportion to `indistinguishability` in [0, 1].
    """
    p_indist = hong_ou_mandel_coincidence_probability(dim, theta=np.pi / 4)
    p_dist = 0.25  # classical, distinguishable-photon coincidence rate at a 50/50 BS
    return indistinguishability * p_indist + (1 - indistinguishability) * p_dist
