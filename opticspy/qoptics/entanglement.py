"""
opticspy.qoptics.entanglement
================================
Quantum-information tools for finite-dimensional systems: partial
trace, von Neumann entropy, purity, fidelity, and Schmidt decomposition
of a bipartite pure state -- plus re-exports of the Bell-state and
concurrence helpers already in `opticspy.quantum` so this subpackage is
a single entry point for entanglement/quantum-information work.

Reference: Gerry & Knight, "Introductory Quantum Optics", Ch. 7
(Entanglement, mixed states, and the density operator).
"""

import numpy as np

from opticspy.quantum import bell_state, concurrence  # re-exported

__all__ = ["bell_state", "concurrence", "partial_trace", "von_neumann_entropy",
           "purity", "fidelity", "schmidt_decomposition", "schmidt_number",
           "density_matrix_from_state"]


def density_matrix_from_state(psi):
    psi = np.asarray(psi, dtype=complex)
    return np.outer(psi, psi.conj())


def partial_trace(rho, dims, keep):
    """
    Partial trace of a bipartite density matrix `rho` with subsystem
    dimensions `dims` (e.g. [2,2] for two qubits), keeping only the
    subsystem index in `keep` (0 or 1).
    """
    d0, d1 = dims
    rho4 = rho.reshape(d0, d1, d0, d1)
    if keep == 0:
        return np.einsum('ijkj->ik', rho4)
    elif keep == 1:
        return np.einsum('ijil->jl', rho4)
    raise ValueError("keep must be 0 or 1 for a bipartite system")


def von_neumann_entropy(rho, base=2):
    """S(rho) = -Tr(rho log rho), via the eigenvalues of rho (base-2 -> bits)."""
    eigvals = np.linalg.eigvalsh(rho)
    eigvals = eigvals[eigvals > 1e-14]
    if base == 2:
        return float(-np.sum(eigvals * np.log2(eigvals)))
    return float(-np.sum(eigvals * np.log(eigvals)))


def purity(rho):
    """Tr(rho^2); 1 for a pure state, 1/d for the maximally mixed state of dimension d."""
    return float(np.real(np.trace(rho @ rho)))


def fidelity(rho1, rho2):
    """
    Uhlmann fidelity F(rho1,rho2) = [Tr sqrt( sqrt(rho1) rho2 sqrt(rho1) )]^2.
    For a pure state rho1 = |psi><psi|, this reduces to <psi|rho2|psi>.
    A small regularization is added before the matrix square root to
    avoid numerical warnings/instability for (near-)singular rho1.
    """
    from scipy.linalg import sqrtm
    d = rho1.shape[0]
    eps = 1e-12
    rho1_reg = (1 - eps) * rho1 + eps * np.eye(d) / d
    sqrt_rho1 = sqrtm(rho1_reg)
    inner_arg = sqrt_rho1 @ rho2 @ sqrt_rho1
    inner_arg = (1 - eps) * inner_arg + eps * np.eye(d) / d
    inner = sqrtm(inner_arg)
    return float(np.real(np.trace(inner)) ** 2)


def schmidt_decomposition(psi, dim_a, dim_b):
    """
    Schmidt decomposition of a bipartite pure state vector `psi`
    (length dim_a*dim_b): returns (schmidt_coefficients, U, Vh) from the
    SVD of psi reshaped into a dim_a x dim_b matrix, with
    psi = sum_k coeff_k * |u_k>|v_k>  (u_k = columns of U, v_k = rows of Vh).
    """
    M = np.asarray(psi, dtype=complex).reshape(dim_a, dim_b)
    U, S, Vh = np.linalg.svd(M)
    return S, U, Vh


def schmidt_number(psi, dim_a, dim_b):
    """Inverse participation ratio 1/sum(p_k^2) of the Schmidt weights p_k=coeff_k^2
    (a continuous measure of entanglement: 1 for a product state)."""
    S, _, _ = schmidt_decomposition(psi, dim_a, dim_b)
    S = S / np.linalg.norm(S)
    p = S ** 2
    return 1.0 / np.sum(p ** 2)
