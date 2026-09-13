"""
Phase 4 demo -- opticspy.qoptics
===================================
(1) Jaynes-Cummings: vacuum Rabi oscillation (Fock field) vs.
    collapse-and-revival (coherent-state field).
(2) Lindblad master equation: cavity photon-number decay vs. the
    analytic exponential, and decoherence of a Schrodinger-cat state
    (loss of Wigner-function interference fringes).
(3) Hong-Ou-Mandel dip: output |1,1> coincidence probability vs.
    beamsplitter angle.
(4) g^(2)(0) photon-bunching/antibunching for Fock, coherent, and
    thermal light.

Run:  python examples/phase4_qoptics_demo.py
"""
import numpy as np
import matplotlib.pyplot as plt

from opticspy.qoptics import jaynes_cummings as jc, master_equation as me
from opticspy.qoptics import beamsplitter as bs, coherence_functions as cf
from opticspy.quantum import fock_state, coherent_state, cat_state, wigner_function

# ---- 1) Jaynes-Cummings: Fock-state Rabi vs coherent-state collapse/revival ----
field_dim = 60
g = 1.0
H0 = jc.hamiltonian(field_dim, 0.0, 0.0, g)
tlist = np.linspace(0, 40, 800)

psi0_fock = jc.initial_state(fock_state(4, field_dim), jc.ATOM_EXCITED)
W_fock = jc.atomic_inversion(jc.evolve(H0, psi0_fock, tlist), field_dim)

nbar = 16.0
alpha = np.sqrt(nbar)
psi0_coh = jc.initial_state(coherent_state(alpha, field_dim), jc.ATOM_EXCITED)
W_coh = jc.atomic_inversion(jc.evolve(H0, psi0_coh, tlist), field_dim)

t_revival = jc.collapse_and_revival_time(g, nbar)
print(f"Coherent-state JC (nbar={nbar}): expected revival time ~ {t_revival:.2f} (1/g units)")

# ---- 2) Cavity decay + cat-state decoherence ----
kappa = 0.15
n0 = 8
psi0 = coherent_state(np.sqrt(n0), 30)
rho0 = np.outer(psi0, psi0.conj())
tlist2 = np.linspace(0, 15, 60)
L = me.cavity_decay_collapse_operator(kappa, 30)
rhos = me.solve_lindblad(np.zeros((30, 30), dtype=complex), rho0, tlist2, [L])
n_op = np.diag(np.arange(30))
n_t = np.array([np.real(np.trace(n_op @ r)) for r in rhos])
n_analytic = me.photon_number_decay(kappa, n0, tlist2)

cat_dim = 30
cat = cat_state(2.0, phi=0.0, dim=cat_dim)
rho_cat0 = np.outer(cat, cat.conj())
L_cat = me.cavity_decay_collapse_operator(0.05, cat_dim)
t_cat = np.linspace(0, 6, 4)
rhos_cat = me.solve_lindblad(np.zeros((cat_dim, cat_dim), dtype=complex), rho_cat0, t_cat, [L_cat])

# ---- 3) Hong-Ou-Mandel dip ----
thetas = np.linspace(0, np.pi / 2, 60)
p11 = [bs.hong_ou_mandel_coincidence_probability(dim=6, theta=th) for th in thetas]

# ---- 4) g2(0) comparison ----
from opticspy.quantum import thermal_state
states = {
    "Fock |1>": fock_state(1, 20),
    "Fock |3>": fock_state(3, 20),
    "Coherent |alpha|^2=4": coherent_state(2.0, 40),
    "Thermal <n>=2": thermal_state(2.0, dim=60),
}
print("\ng^(2)(0) photon statistics:")
for name, s in states.items():
    print(f"  {name:24s}: g2(0) = {cf.g2_zero_delay(s):.3f}")

# ---- Figure ----
fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.5))

axes[0, 0].plot(tlist, W_fock, label='field: Fock |4>')
axes[0, 0].plot(tlist, W_coh, label=f'field: coherent, nbar={nbar}', alpha=0.8)
axes[0, 0].axvline(t_revival, color='gray', ls='--', lw=1, label='predicted revival time')
axes[0, 0].set_xlabel("g*t"); axes[0, 0].set_ylabel("atomic inversion <sigma_z>")
axes[0, 0].set_title("Jaynes-Cummings: Rabi vs. collapse & revival")
axes[0, 0].legend(fontsize=7)

axes[0, 1].plot(tlist2 * kappa, n_t, 'o', ms=3, label='Lindblad RK4')
axes[0, 1].plot(tlist2 * kappa, n_analytic, '-', label='analytic n0*exp(-kappa t)')
axes[0, 1].set_xlabel("kappa * t"); axes[0, 1].set_ylabel("<n>(t)")
axes[0, 1].set_title("Cavity photon-number decay")
axes[0, 1].legend(fontsize=8)

axes[1, 0].plot(thetas / np.pi, p11)
axes[1, 0].axvline(0.25, color='gray', ls='--', lw=1)
axes[1, 0].set_xlabel("beamsplitter angle theta / pi"); axes[1, 0].set_ylabel("P(1,1) at output")
axes[1, 0].set_title("Hong-Ou-Mandel dip (theta=pi/4: balanced BS)")

xvec = np.linspace(-5, 5, 120)
W0 = wigner_function(rhos_cat[0], xvec, xvec)
W2 = wigner_function(rhos_cat[-1], xvec, xvec)
im = axes[1, 1].contourf(xvec, xvec, W2 - 0 * W0, levels=40, cmap='RdBu_r')
axes[1, 1].set_title(f"Cat-state Wigner function after decoherence\n(t={t_cat[-1]:.1f}, kappa=0.05)")
axes[1, 1].set_xlabel("x"); axes[1, 1].set_ylabel("p")
fig.colorbar(im, ax=axes[1, 1], shrink=0.8)

fig.tight_layout()
fig.savefig("phase4_qoptics_demo.png", dpi=140)
print("\nSaved figure: phase4_qoptics_demo.png")
