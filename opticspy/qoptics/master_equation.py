"""
aoptics.qoptics.master_equation
===================================
A general-purpose Lindblad master equation integrator,

    drho/dt = -i/hbar [H, rho] + sum_k ( L_k rho L_k^dagger
                                          - 1/2 {L_k^dagger L_k, rho} )

integrated with fixed-step 4th-order Runge-Kutta, plus convenience
collapse operators for the two archetypal open-system problems of
quantum optics: cavity photon loss and spontaneous emission of a
two-level atom.

Reference: Gerry & Knight, "Introductory Quantum Optics", Ch. 8
(Interaction of the atom with a damped cavity field / master equation).
"""

import numpy as np

from .operators import annihilation_operator


def lindblad_rhs(rho, H, collapse_ops, hbar=1.0):
    drho = -1j / hbar * (H @ rho - rho @ H)
    for L in collapse_ops:
        Ld = L.conj().T
        drho += L @ rho @ Ld - 0.5 * (Ld @ L @ rho + rho @ Ld @ L)
    return drho


def solve_lindblad(H, rho0, tlist, collapse_ops=(), hbar=1.0):
    """
    Fixed-step RK4 integration of the Lindblad master equation.
    Returns an array of density matrices, one per entry of `tlist`
    (which should be uniformly spaced).
    """
    tlist = np.asarray(tlist, dtype=float)
    dt = tlist[1] - tlist[0]
    rho = rho0.astype(complex).copy()
    out = [rho.copy()]
    for _ in range(len(tlist) - 1):
        k1 = lindblad_rhs(rho, H, collapse_ops, hbar)
        k2 = lindblad_rhs(rho + 0.5 * dt * k1, H, collapse_ops, hbar)
        k3 = lindblad_rhs(rho + 0.5 * dt * k2, H, collapse_ops, hbar)
        k4 = lindblad_rhs(rho + dt * k3, H, collapse_ops, hbar)
        rho = rho + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        # re-hermitize / re-normalize to control accumulated numerical drift
        rho = 0.5 * (rho + rho.conj().T)
        rho = rho / np.real(np.trace(rho))
        out.append(rho.copy())
    return np.array(out)


def cavity_decay_collapse_operator(kappa, field_dim):
    """Photon-loss (amplitude damping) collapse operator L = sqrt(kappa) * a."""
    return np.sqrt(kappa) * annihilation_operator(field_dim)


def spontaneous_emission_collapse_operator(gamma, sigma_minus):
    """Two-level-atom spontaneous emission collapse operator L = sqrt(gamma) * sigma_-."""
    return np.sqrt(gamma) * sigma_minus


def photon_number_decay(kappa, n0, tlist):
    """
    Closed-form mean photon number under pure cavity loss (no driving),
    <n>(t) = n0 * exp(-kappa*t) -- used as an analytic cross-check of
    `solve_lindblad` with a single amplitude-damping collapse operator.
    """
    return n0 * np.exp(-kappa * np.asarray(tlist))

