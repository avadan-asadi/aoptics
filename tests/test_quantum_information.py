import numpy as np
import pytest

from opticspy.qoptics import quantum_information as qi


# --------------------------- gates ---------------------------

def test_hadamard_squared_is_identity():
    assert np.allclose(qi.HADAMARD @ qi.HADAMARD, qi.I2, atol=1e-10)


def test_pauli_gates_square_to_identity():
    for P in (qi.PAULI_X, qi.PAULI_Y, qi.PAULI_Z):
        assert np.allclose(P @ P, qi.I2, atol=1e-10)


def test_hadamard_creates_superposition():
    state = qi.HADAMARD @ qi.QUBIT_0
    assert np.allclose(state, qi.QUBIT_PLUS, atol=1e-10)


def test_rotation_gate_full_turn_is_minus_identity():
    Rz = qi.rotation_gate('z', 2 * np.pi)
    assert np.allclose(Rz, -qi.I2, atol=1e-8)


def test_cnot_flips_target_when_control_is_one():
    state = np.kron(qi.QUBIT_1, qi.QUBIT_0)  # |10>
    out = qi.CNOT @ state
    expected = np.kron(qi.QUBIT_1, qi.QUBIT_1)  # |11>
    assert np.allclose(out, expected)


def test_bell_pair_is_normalized_and_entangled():
    psi = qi.bell_pair()
    assert np.linalg.norm(psi) == pytest.approx(1.0)
    expected = np.array([1, 0, 0, 1], dtype=complex) / np.sqrt(2)
    assert np.allclose(np.abs(psi), np.abs(expected), atol=1e-10)


# --------------------------- teleportation ---------------------------

def test_teleportation_recovers_arbitrary_state():
    rng = np.random.default_rng(123)
    # a handful of random single-qubit states
    for _ in range(8):
        theta = rng.uniform(0, np.pi)
        phi = rng.uniform(0, 2 * np.pi)
        psi = np.array([np.cos(theta / 2), np.exp(1j * phi) * np.sin(theta / 2)])
        bob_state, bits, fidelity = qi.teleportation_protocol(psi, rng=rng)
        assert fidelity == pytest.approx(1.0, abs=1e-8)


def test_teleportation_bits_are_valid():
    rng = np.random.default_rng(7)
    _, bits, _ = qi.teleportation_protocol(qi.QUBIT_PLUS, rng=rng)
    assert all(b in (0, 1) for b in bits)


# --------------------------- BB84 ---------------------------

def test_bb84_no_eavesdropper_zero_error_rate():
    rng = np.random.default_rng(42)
    result = qi.bb84_simulate(4000, eavesdrop=False, rng=rng)
    assert result['n_sifted'] > 1500
    assert result['qber'] < 0.01


def test_bb84_eavesdropper_introduces_expected_error_rate():
    rng = np.random.default_rng(99)
    result = qi.bb84_simulate(6000, eavesdrop=True, rng=rng)
    # intercept-resend eavesdropping on BB84 gives the well-known ~25% QBER
    assert 0.18 < result['qber'] < 0.32


def test_bb84_sifted_fraction_near_half():
    rng = np.random.default_rng(5)
    result = qi.bb84_simulate(5000, eavesdrop=False, rng=rng)
    frac = result['n_sifted'] / 5000
    assert 0.4 < frac < 0.6
