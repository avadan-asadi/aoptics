import numpy as np
import pytest

from aoptics.qoptics import operators as ops
from aoptics.qoptics import jaynes_cummings as jc
from aoptics.qoptics import master_equation as me
from aoptics.qoptics import beamsplitter as bs
from aoptics.qoptics import coherence_functions as cf
from aoptics.qoptics import entanglement as ent
from aoptics.quantum import fock_state, coherent_state, thermal_state, bell_state


# --------------------------- operators ---------------------------

def test_commutator_a_adag_equals_identity():
    dim = 15
    a = ops.annihilation_operator(dim)
    adag = ops.creation_operator(dim)
    comm = ops.commutator(a, adag)
    # [a, a^dagger] = 1 exactly, except at the truncation boundary
    assert np.allclose(np.diag(comm)[:-1], 1.0, atol=1e-10)


def test_number_operator_eigenvalues():
    dim = 10
    n_op = ops.number_operator(dim)
    assert np.allclose(np.diag(n_op), np.arange(dim))


def test_vacuum_saturates_heisenberg_uncertainty():
    dim = 20
    vac = fock_state(0, dim)
    product = ops.heisenberg_uncertainty_product(vac, dim)
    assert product == pytest.approx(0.25, rel=1e-6)


def test_displacement_operator_creates_coherent_state():
    dim = 40
    alpha = 1.3 + 0.4j
    vac = fock_state(0, dim)
    D = ops.displacement_operator(alpha, dim)
    disp_vac = D @ vac
    target = coherent_state(alpha, dim)
    overlap = np.abs(np.vdot(disp_vac, target)) ** 2
    assert overlap == pytest.approx(1.0, rel=1e-3)


def test_coherent_state_has_symmetric_quadrature_variances():
    dim = 40
    alpha = 1.0 + 0.5j
    state = coherent_state(alpha, dim)
    X, P = ops.quadrature_x(dim), ops.quadrature_p(dim)
    varX = np.real(ops.variance(X, state))
    varP = np.real(ops.variance(P, state))
    assert varX == pytest.approx(0.5, rel=1e-3)
    assert varP == pytest.approx(0.5, rel=1e-3)


# --------------------------- Jaynes-Cummings ---------------------------

def test_vacuum_rabi_oscillation_matches_analytic_cosine():
    # Atom excited, field in vacuum, on resonance: P_e(t) = cos^2(g t)
    field_dim = 5
    g = 1.0
    H = jc.hamiltonian(field_dim, omega_cavity=0.0, omega_atom=0.0, g=g)
    vac = fock_state(0, field_dim)
    psi0 = jc.initial_state(vac, jc.ATOM_EXCITED)
    tlist = np.linspace(0, 4 * np.pi, 60)
    psis = jc.evolve(H, psi0, tlist)
    W = jc.atomic_inversion(psis, field_dim)
    P_e = (W + 1) / 2
    expected = np.cos(g * tlist) ** 2
    assert np.allclose(P_e, expected, atol=1e-6)


def test_photon_number_conserved_total_excitation_on_resonance():
    field_dim = 8
    g = 0.7
    H = jc.hamiltonian(field_dim, omega_cavity=0.0, omega_atom=0.0, g=g)
    field0 = fock_state(2, field_dim)
    psi0 = jc.initial_state(field0, jc.ATOM_GROUND)
    tlist = np.linspace(0, 10, 40)
    psis = jc.evolve(H, psi0, tlist)
    n_field = jc.mean_photon_number(psis, field_dim)
    W = jc.atomic_inversion(psis, field_dim)
    P_e = (W + 1) / 2
    total_excitation = n_field + P_e
    assert np.allclose(total_excitation, 2.0, atol=1e-6)  # JC conserves n_field + n_atom


def test_dressed_state_splitting_matches_vacuum_rabi_frequency():
    g, omega = 0.5, 2.0
    E_minus, E_plus = jc.dressed_state_energies(n=0, omega=omega, g=g)
    assert (E_plus - E_minus) == pytest.approx(jc.vacuum_rabi_frequency(g, n=0))


# --------------------------- master equation ---------------------------

def test_cavity_decay_matches_analytic_exponential():
    field_dim = 20
    kappa = 0.3
    n0 = 5
    alpha = np.sqrt(n0)
    psi0 = coherent_state(alpha, field_dim)
    rho0 = np.outer(psi0, psi0.conj())
    H = np.zeros((field_dim, field_dim), dtype=complex)  # no free evolution needed for <n>
    L = me.cavity_decay_collapse_operator(kappa, field_dim)
    tlist = np.linspace(0, 5, 60)
    rhos = me.solve_lindblad(H, rho0, tlist, collapse_ops=[L])
    n_op = np.diag(np.arange(field_dim))
    n_t = np.array([np.real(np.trace(n_op @ r)) for r in rhos])
    expected = me.photon_number_decay(kappa, n0, tlist)
    assert np.allclose(n_t, expected, rtol=0.02, atol=0.05)


def test_lindblad_preserves_trace():
    field_dim = 10
    kappa = 0.5
    psi0 = fock_state(3, field_dim)
    rho0 = np.outer(psi0, psi0.conj())
    H = np.diag(np.arange(field_dim)).astype(complex)
    L = me.cavity_decay_collapse_operator(kappa, field_dim)
    tlist = np.linspace(0, 3, 50)
    rhos = me.solve_lindblad(H, rho0, tlist, collapse_ops=[L])
    traces = np.array([np.real(np.trace(r)) for r in rhos])
    assert np.allclose(traces, 1.0, atol=1e-6)


# --------------------------- beamsplitter / HOM ---------------------------

def test_balanced_beamsplitter_power_split_50_50():
    n_a, n_b = bs.splitting_ratio(theta=np.pi / 4, dim=12, alpha=1.5)
    total = n_a + n_b
    assert n_a / total == pytest.approx(0.5, rel=1e-3)
    assert n_b / total == pytest.approx(0.5, rel=1e-3)


def test_beamsplitter_unitary_is_unitary():
    U = bs.beamsplitter_unitary(0.37, dim=6)
    assert np.allclose(U @ U.conj().T, np.eye(U.shape[0]), atol=1e-8)


def test_hong_ou_mandel_dip_at_balanced_beamsplitter():
    p11 = bs.hong_ou_mandel_coincidence_probability(dim=6, theta=np.pi / 4)
    assert p11 == pytest.approx(0.0, abs=1e-9)


def test_hong_ou_mandel_nonzero_away_from_balance():
    p11 = bs.hong_ou_mandel_coincidence_probability(dim=6, theta=np.pi / 6)
    assert p11 > 1e-3


def test_hom_visibility_interpolates_indistinguishability():
    v_full = bs.hong_ou_mandel_visibility_vs_distinguishability(1.0)
    v_none = bs.hong_ou_mandel_visibility_vs_distinguishability(0.0)
    assert v_full == pytest.approx(0.0, abs=1e-9)
    assert v_none == pytest.approx(0.25, rel=1e-9)


# --------------------------- coherence functions (g2) ---------------------------

def test_g2_fock_state_is_zero_for_n1():
    state = fock_state(1, 10)
    assert cf.g2_zero_delay(state) == pytest.approx(0.0, abs=1e-9)
    assert cf.is_antibunched(state)


def test_g2_coherent_state_is_unity():
    state = coherent_state(2.0, 60)
    g2 = cf.g2_zero_delay(state)
    assert g2 == pytest.approx(1.0, rel=1e-2)


def test_g2_thermal_state_is_two():
    rho = thermal_state(3.0, dim=80)
    g2 = cf.g2_zero_delay(rho)
    assert g2 == pytest.approx(2.0, rel=1e-2)


def test_fano_factor_coherent_state_is_one():
    state = coherent_state(1.5, 60)
    assert cf.fano_factor(state) == pytest.approx(1.0, rel=1e-2)


# --------------------------- entanglement ---------------------------

def test_bell_state_maximally_entangled():
    psi = bell_state('phi+')
    rho = ent.density_matrix_from_state(psi)
    rho_a = ent.partial_trace(rho, [2, 2], keep=0)
    assert ent.purity(rho_a) == pytest.approx(0.5, rel=1e-6)
    assert ent.von_neumann_entropy(rho_a) == pytest.approx(1.0, rel=1e-6)  # 1 bit of entanglement


def test_product_state_zero_entanglement_entropy():
    psi_a = np.array([1, 0], dtype=complex)
    psi_b = np.array([0, 1], dtype=complex)
    psi = np.kron(psi_a, psi_b)
    rho = ent.density_matrix_from_state(psi)
    rho_a = ent.partial_trace(rho, [2, 2], keep=0)
    assert ent.von_neumann_entropy(rho_a) == pytest.approx(0.0, abs=1e-9)
    assert ent.purity(rho_a) == pytest.approx(1.0, rel=1e-9)


def test_fidelity_of_identical_pure_states_is_one():
    psi = bell_state('psi-')
    rho = ent.density_matrix_from_state(psi)
    assert ent.fidelity(rho, rho) == pytest.approx(1.0, rel=1e-3)


def test_schmidt_number_bell_state_is_two():
    psi = bell_state('phi+')
    n_schmidt = ent.schmidt_number(psi, 2, 2)
    assert n_schmidt == pytest.approx(2.0, rel=1e-6)


def test_schmidt_number_product_state_is_one():
    psi_a = np.array([1, 0], dtype=complex)
    psi_b = np.array([0, 1], dtype=complex)
    psi = np.kron(psi_a, psi_b)
    n_schmidt = ent.schmidt_number(psi, 2, 2)
    assert n_schmidt == pytest.approx(1.0, rel=1e-6)


def test_concurrence_bell_state_is_one():
    psi = bell_state('phi-')
    rho = ent.density_matrix_from_state(psi)
    assert ent.concurrence(rho) == pytest.approx(1.0, rel=1e-6)

