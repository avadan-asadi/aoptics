"""
opticspy.qoptics.quantum_information
=======================================
Qubit gates, the quantum teleportation protocol, and a simulation of
the BB84 quantum-key-distribution protocol (including a simple
intercept-resend eavesdropper, whose signature is the well-known 25%
induced quantum bit error rate).

Reference: Gerry & Knight, "Introductory Quantum Optics", Ch. 9
(Quantum information processing and quantum computation) and the
optical-qubit / photon-polarization framing used throughout the text.
"""

import numpy as np

I2 = np.eye(2, dtype=complex)
PAULI_X = np.array([[0, 1], [1, 0]], dtype=complex)
PAULI_Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
PAULI_Z = np.array([[1, 0], [0, -1]], dtype=complex)
HADAMARD = (1 / np.sqrt(2)) * np.array([[1, 1], [1, -1]], dtype=complex)
S_GATE = np.array([[1, 0], [0, 1j]], dtype=complex)
T_GATE = np.array([[1, 0], [0, np.exp(1j * np.pi / 4)]], dtype=complex)

QUBIT_0 = np.array([1, 0], dtype=complex)
QUBIT_1 = np.array([0, 1], dtype=complex)
QUBIT_PLUS = (QUBIT_0 + QUBIT_1) / np.sqrt(2)
QUBIT_MINUS = (QUBIT_0 - QUBIT_1) / np.sqrt(2)

CNOT = np.array([
    [1, 0, 0, 0],
    [0, 1, 0, 0],
    [0, 0, 0, 1],
    [0, 0, 1, 0],
], dtype=complex)


def rotation_gate(axis, angle):
    """Single-qubit rotation exp(-i*angle/2 * sigma_axis), axis in {'x','y','z'}."""
    from scipy.linalg import expm
    sigma = {'x': PAULI_X, 'y': PAULI_Y, 'z': PAULI_Z}[axis]
    return expm(-1j * angle / 2 * sigma)


def apply_single_qubit_gate(gate, state, qubit_index, n_qubits):
    """Apply a single-qubit `gate` to `qubit_index` (0-indexed) of an
    n_qubits combined state vector (length 2**n_qubits)."""
    ops = [I2] * n_qubits
    ops[qubit_index] = gate
    full = ops[0]
    for op in ops[1:]:
        full = np.kron(full, op)
    return full @ state


def bell_pair(gate1=HADAMARD):
    """Create the Bell pair |Phi+> = (|00>+|11>)/sqrt(2) from |00> via H (x) I then CNOT."""
    psi = np.kron(QUBIT_0, QUBIT_0)
    psi = apply_single_qubit_gate(gate1, psi, 0, 2)
    return CNOT @ psi


def teleportation_protocol(psi_to_teleport, rng=None):
    """
    Simulate the standard quantum teleportation protocol: Alice holds an
    unknown qubit `psi_to_teleport` and one half of a shared Bell pair;
    she performs a Bell-basis measurement on her two qubits and sends
    the 2 classical bits to Bob, who applies a corresponding correction
    to his half of the Bell pair, recovering `psi_to_teleport` exactly
    (up to global phase).

    Returns (bob_final_state, measurement_bits, fidelity_to_original).
    """
    rng = np.random.default_rng() if rng is None else rng

    # 3-qubit state: [Alice's unknown qubit, Alice's Bell-pair half, Bob's Bell-pair half]
    epr = bell_pair()  # qubits (1,2)
    psi = np.kron(psi_to_teleport, epr)  # qubits (0,1,2)

    # CNOT(0 -> 1), then H on qubit 0
    CNOT_01 = np.kron(CNOT, I2)
    psi = CNOT_01 @ psi
    psi = apply_single_qubit_gate(HADAMARD, psi, 0, 3)

    # measure qubits 0 and 1 in the computational basis
    probs = np.abs(psi) ** 2
    outcome = rng.choice(8, p=probs / probs.sum())
    bits = [(outcome >> 2) & 1, (outcome >> 1) & 1]  # qubit0, qubit1 measurement results

    # Bob's (qubit 2) reduced state after the projective measurement on
    # qubits 0,1 is proportional to the corresponding sub-vector of psi.
    bob_state = np.zeros(2, dtype=complex)
    for b2 in (0, 1):
        idx = (bits[0] << 2) | (bits[1] << 1) | b2
        bob_state[b2] = psi[idx]
    bob_state = bob_state / np.linalg.norm(bob_state)

    # classical correction: apply X^bits[1] Z^bits[0]
    if bits[1] == 1:
        bob_state = PAULI_X @ bob_state
    if bits[0] == 1:
        bob_state = PAULI_Z @ bob_state

    fidelity = np.abs(np.vdot(psi_to_teleport, bob_state)) ** 2
    return bob_state, bits, fidelity


def bb84_simulate(n_bits, eavesdrop=False, rng=None):
    """
    Simulate the BB84 QKD protocol over `n_bits` transmitted qubits:
    Alice picks random bits and random bases (Z or X), Bob picks random
    measurement bases; if `eavesdrop` is True, Eve intercepts every
    qubit, measures it in a random basis, and resends the (possibly
    disturbed) state before Bob receives it. Returns a dict with the
    sifted key length, the quantum bit error rate (QBER) on the sifted
    key, and the raw bit/basis arrays.
    """
    rng = np.random.default_rng() if rng is None else rng
    alice_bits = rng.integers(0, 2, n_bits)
    alice_bases = rng.integers(0, 2, n_bits)   # 0 = Z basis, 1 = X basis
    bob_bases = rng.integers(0, 2, n_bits)

    def encode(bit, basis):
        if basis == 0:
            return QUBIT_0 if bit == 0 else QUBIT_1
        return QUBIT_PLUS if bit == 0 else QUBIT_MINUS

    def measure(state, basis, rng):
        if basis == 1:
            state = HADAMARD @ state  # rotate X-basis into Z-basis for measurement
        p0 = np.abs(state[0]) ** 2
        return 0 if rng.random() < p0 else 1

    bob_bits = np.zeros(n_bits, dtype=int)
    for i in range(n_bits):
        state = encode(alice_bits[i], alice_bases[i])
        if eavesdrop:
            eve_basis = rng.integers(0, 2)
            eve_bit = measure(state, eve_basis, rng)
            state = encode(eve_bit, eve_basis)  # Eve resends her (possibly wrong) guess
        bob_bits[i] = measure(state, bob_bases[i], rng)

    sifted_mask = alice_bases == bob_bases
    alice_sifted = alice_bits[sifted_mask]
    bob_sifted = bob_bits[sifted_mask]
    n_sifted = len(alice_sifted)
    n_errors = int(np.sum(alice_sifted != bob_sifted))
    qber = n_errors / n_sifted if n_sifted > 0 else np.nan

    return dict(n_sifted=n_sifted, qber=qber, alice_bits=alice_bits, alice_bases=alice_bases,
                bob_bases=bob_bases, bob_bits=bob_bits, sifted_mask=sifted_mask)
