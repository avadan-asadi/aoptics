"""
aoptics.quantum
================
Quantum optics simulation: Fock/coherent/squeezed/cat states,
Wigner/Husimi functions, photon statistics, Jaynes-Cummings, HOM, Bell states.
"""

import numpy as np
from scipy.special import factorial, genlaguerre


def fock_state(n, dim=30):
    s = np.zeros(dim, dtype=complex); s[min(n, dim-1)] = 1.0; return s


def coherent_state(alpha, dim=60):
    n = np.arange(dim)
    return (np.exp(-0.5*np.abs(alpha)**2) * alpha**n / np.sqrt(factorial(n, exact=False))).astype(complex)


def squeezed_vacuum(r, phi=0.0, dim=60):
    state = np.zeros(dim, dtype=complex)
    for n in range(0, dim-1, 2):
        c = (np.sqrt(factorial(n, exact=False)) /
             (2**(n//2) * factorial(n//2, exact=False)))
        state[n] = np.sqrt(1/np.cosh(r)) * c * (-np.exp(1j*phi)*np.tanh(r))**(n//2)
    return state / np.linalg.norm(state)


def thermal_state(nbar, dim=60):
    n = np.arange(dim)
    return np.diag((nbar**n / (1+nbar)**(n+1)).astype(complex))


def cat_state(alpha, phi=0.0, dim=80):
    s = coherent_state(alpha, dim) + np.exp(1j*phi)*coherent_state(-alpha, dim)
    return s / np.linalg.norm(s)


def wigner_function(state_or_rho, xvec=None, pvec=None):
    if xvec is None: xvec = np.linspace(-6, 6, 100)
    if pvec is None: pvec = np.linspace(-6, 6, 100)
    rho = np.outer(state_or_rho, state_or_rho.conj()) if state_or_rho.ndim == 1 else state_or_rho
    dim = rho.shape[0]
    X, P = np.meshgrid(xvec, pvec)
    alpha = (X + 1j*P) / np.sqrt(2)
    W = np.zeros_like(X, dtype=float)
    for n in range(dim):
        for m in range(dim):
            if np.abs(rho[n,m]) < 1e-12: continue
            r2 = np.abs(alpha)**2
            if m >= n:
                pref = (-1)**n * np.sqrt(factorial(n,exact=False)/factorial(m,exact=False))
                L = genlaguerre(n, m-n)(4*r2)
                elem = pref * (2*alpha)**(m-n) * np.exp(-2*r2) * L
            else:
                pref = (-1)**m * np.sqrt(factorial(m,exact=False)/factorial(n,exact=False))
                L = genlaguerre(m, n-m)(4*r2)
                elem = np.conj(pref * (2*alpha)**(n-m) * np.exp(-2*r2) * L)
            W += np.real(rho[n,m] * elem)
    return W * 2/np.pi


def photon_statistics(state_or_rho, max_n=None):
    diag = np.abs(state_or_rho)**2 if state_or_rho.ndim == 1 else np.real(np.diag(state_or_rho))
    if max_n: diag = diag[:max_n]
    n_vals = np.arange(len(diag))
    P_n = diag / (diag.sum() + 1e-30)
    mean_n = np.sum(n_vals * P_n)
    var_n = np.sum(n_vals**2 * P_n) - mean_n**2
    mandel_Q = (var_n - mean_n) / (mean_n + 1e-30)
    return n_vals, P_n, mean_n, var_n, mandel_Q


def hong_ou_mandel_visibility(tau, tau_c=1e-12):
    return 1 - 0.5 * np.exp(-(tau/tau_c)**2)


def bell_state(state='phi+'):
    s = 1/np.sqrt(2)
    return {'phi+': np.array([s,0,0,s]), 'phi-': np.array([s,0,0,-s]),
            'psi+': np.array([0,s,s,0]), 'psi-': np.array([0,s,-s,0])}[state].astype(complex)


def concurrence(rho_2qubit):
    sy = np.array([[0,-1j],[1j,0]])
    YY = np.kron(sy, sy)
    M = rho_2qubit @ (YY @ rho_2qubit.conj() @ YY)
    ev = np.sqrt(np.maximum(np.sort(np.real(np.linalg.eigvals(M)))[::-1], 0))
    return float(max(0, ev[0]-ev[1]-ev[2]-ev[3]))

