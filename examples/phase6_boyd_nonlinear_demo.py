"""
Phase 6 demo -- opticspy.nonlinear_optics_boyd
=================================================
(1) Sum-frequency generation: coupled-wave amplitude growth and the
    Manley-Rowe photon-flux conservation check.
(2) Two-level-atom saturable absorption: saturated Lorentzian lineshape.
(3) Self-focusing: critical power and the Marburger collapse distance
    vs. input power.
(4) Two-photon absorption: nonlinear transmission vs. input intensity.

Run:  python examples/phase6_boyd_nonlinear_demo.py
"""
import numpy as np
import matplotlib.pyplot as plt

from opticspy.nonlinear_optics_boyd import (coupled_wave_mixing as cwm, two_level_atom as tla,
                                             self_action as sa, multiphoton_absorption as tpa)

# ---- 1) SFG coupled-wave growth + Manley-Rowe ----
omega1, omega2 = 1.5e15, 1.8e15
omega3 = omega1 + omega2
z, A1, A2, A3 = cwm.simulate_three_wave_mixing(
    A1_0=1.0, A2_0=1.0, A3_0=0.0, z_max=0.15,
    omega1=omega1, omega2=omega2, omega3=omega3,
    n1=1.5, n2=1.5, n3=1.5, deff=1e-12, delta_k=0.0, n_steps=1000)
f1, f2, f3 = cwm.manley_rowe_photon_fluxes(A1, A2, A3, omega1, omega2, omega3)
print(f"SFG: |A3| grows from {abs(A3[0]):.3e} to {abs(A3[-1]):.3e} over {z[-1]*100:.0f} cm")
print(f"Manley-Rowe invariant f1+f3 relative drift: "
      f"{(np.max(f1+f3)-np.min(f1+f3))/(f1+f3)[0]:.2e} (should be ~0)")

# ---- 2) Two-level saturable absorption ----
detunings = np.linspace(-5, 5, 400)
alpha0 = 1.0
curves = {}
for s in [0.0, 1.0, 5.0]:
    curves[s] = tla.saturated_absorption_coefficient(alpha0, intensity=s, saturation_intensity=1.0,
                                                       detuning=detunings, linewidth=1.0)

# ---- 3) Self-focusing ----
wl, n0, n2 = 800e-9, 1.0003, 5e-19 * 1e-4  # air n2, converted cm^2/W -> m^2/W
Pcr = sa.critical_power_self_focusing(wl, n0, n2)
zR = 0.6  # m, illustrative Rayleigh range
P_over_Pcr = np.linspace(1.01, 30, 300)
z_collapse = [sa.self_focusing_collapse_distance(p * Pcr, Pcr, zR) for p in P_over_Pcr]
print(f"\nSelf-focusing in air at 800 nm: Pcr = {Pcr/1e9:.2f} GW")

# ---- 4) Two-photon absorption ----
I_in = np.logspace(-2, 3, 300)  # arbitrary units
T_rel = [tpa.nonlinear_transmission_tpa(I, beta_tpa=0.05, length=1.0)[1] for I in I_in]

# ---- Figure ----
fig, axes = plt.subplots(2, 2, figsize=(11, 8.5))

axes[0, 0].plot(z * 100, np.abs(A1), label='|A1| (signal)', color='C0')
axes[0, 0].set_xlabel("z (cm)"); axes[0, 0].set_ylabel("|A1| (signal, ~undepleted)", color='C0')
ax_twin = axes[0, 0].twinx()
ax_twin.plot(z * 100, np.abs(A3), label='|A3| (sum-frequency)', color='C1')
ax_twin.set_ylabel("|A3| (sum-frequency, grows from 0)", color='C1')
axes[0, 0].set_title("Sum-frequency generation (coupled-wave)")

for s, curve in curves.items():
    axes[0, 1].plot(detunings, curve, label=f"I/Isat={s}")
axes[0, 1].set_xlabel("detuning / linewidth"); axes[0, 1].set_ylabel("absorption coefficient")
axes[0, 1].set_title("Two-level saturated absorption (power broadening)")
axes[0, 1].legend(fontsize=8)

axes[1, 0].plot(P_over_Pcr, z_collapse)
axes[1, 0].set_xlabel("P / Pcr"); axes[1, 0].set_ylabel("collapse distance (m)")
axes[1, 0].set_title("Marburger self-focusing collapse distance")
axes[1, 0].set_yscale('log')

axes[1, 1].semilogx(I_in, T_rel)
axes[1, 1].set_xlabel("input intensity (arb. units)")
axes[1, 1].set_ylabel("transmission relative to linear-only")
axes[1, 1].set_title("Two-photon-absorption nonlinear transmission")

fig.tight_layout()
fig.savefig("phase6_boyd_nonlinear_demo.png", dpi=140)
print("\nSaved figure: phase6_boyd_nonlinear_demo.png")
