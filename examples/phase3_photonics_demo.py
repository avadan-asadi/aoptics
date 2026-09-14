"""
Phase 3 demo -- aoptics.photonics
=====================================
(1) Planar waveguide: exact TE mode effective indices vs. core thickness.
(2) Optical fiber: V-number/mode-count vs wavelength, and the material
    (Sellmeier) dispersion curve of fused silica showing the ~1.27 um
    zero-dispersion wavelength.
(3) Fabry-Perot resonator: Airy transmission lineshape and finesse.
(4) Nonlinear optics: SHG phase-matching sinc^2 curve.

Run:  python examples/phase3_photonics_demo.py
"""
import numpy as np
import matplotlib.pyplot as plt

from aoptics.photonics import waveguides as wg, fibers, resonators as res, nonlinear_optics as nl

# ---- 1) Planar waveguide modes vs thickness ----
wavelength = 1.0e-6
n1, n2 = 1.50, 1.45
thicknesses = np.linspace(0.3e-6, 8e-6, 60)
mode_indices = [wg.solve_te_modes(wavelength, d, n1, n2) for d in thicknesses]
max_modes = max(len(m) for m in mode_indices)
print(f"Slab waveguide (n1={n1}, n2={n2}, lambda={wavelength*1e6} um): "
      f"up to {max_modes} TE modes for the thickest core scanned.")

# ---- 2) Fiber V-number and dispersion ----
a = 4.1e-6
n_core, n_clad = 1.4504, 1.4447
wavelengths = np.linspace(0.8e-6, 1.7e-6, 200)
V = fibers.v_number(wavelengths, a, n_core, n_clad)
single_mode_cut = wavelengths[np.argmin(np.abs(V - 2.405))]
print(f"\nFiber (a={a*1e6} um, NA={fibers.numerical_aperture(n_core,n_clad):.4f}): "
      f"single-mode for wavelength > {single_mode_cut*1e9:.0f} nm")

lam_range = np.linspace(0.5e-6, 1.7e-6, 400)
n_silica = fibers.sellmeier_index_silica(lam_range)
D_mat = np.array([fibers.material_dispersion_parameter(fibers.sellmeier_index_silica, l)
                   for l in lam_range]) * 1e6  # -> ps/(nm.km)-like units (s/m^2 * 1e6)
zdw = fibers.zero_dispersion_wavelength(fibers.sellmeier_index_silica, 1.0e-6, 1.4e-6)
print(f"Fused-silica material zero-dispersion wavelength(s): "
      f"{[f'{z*1e9:.0f} nm' for z in zdw]}")

# ---- 3) Fabry-Perot resonator ----
L_cavity = 0.05  # m
R1 = R2 = 0.95
F = res.finesse(R1, R2)
fsr = res.free_spectral_range(L_cavity)
lw = res.linewidth(L_cavity, R1, R2)
print(f"\nFabry-Perot cavity (L={L_cavity*100:.0f} cm, R={R1}): "
      f"finesse={F:.1f}, FSR={fsr/1e9:.3f} GHz, linewidth={lw/1e6:.2f} MHz")
delta = np.linspace(-3 * np.pi, 3 * np.pi, 2000)
T = res.airy_transmission(delta, R1, R2)

# ---- 4) SHG phase matching ----
Lcoh = 20e-6  # m, illustrative coherence length
delta_k_range = np.linspace(-3 * np.pi / Lcoh, 3 * np.pi / Lcoh, 400)
crystal_length = 1e-3  # 1 mm crystal
eta = nl.shg_phase_matching_efficiency_array(delta_k_range, crystal_length)

# ---- Figure ----
fig, axes = plt.subplots(2, 2, figsize=(11, 8))

for i, d in enumerate(thicknesses):
    n_eff = mode_indices[i]
    axes[0, 0].plot([d * 1e6] * len(n_eff), n_eff, 'b.', ms=2)
axes[0, 0].axhline(n1, color='gray', ls=':', lw=1)
axes[0, 0].axhline(n2, color='gray', ls=':', lw=1)
axes[0, 0].set_xlabel("core thickness (micron)"); axes[0, 0].set_ylabel("mode effective index")
axes[0, 0].set_title("Slab waveguide: TE mode indices vs. thickness")

axes[0, 1].plot(lam_range * 1e9, D_mat)
axes[0, 1].axhline(0, color='k', lw=0.6)
for z in zdw:
    axes[0, 1].axvline(z * 1e9, color='r', ls='--', lw=1)
axes[0, 1].set_xlabel("wavelength (nm)"); axes[0, 1].set_ylabel("material dispersion (arb. scaled)")
axes[0, 1].set_title("Fused-silica material dispersion (Sellmeier)")

axes[1, 0].plot(delta / np.pi, T)
axes[1, 0].set_xlabel("round-trip phase (units of pi)"); axes[1, 0].set_ylabel("transmission")
axes[1, 0].set_title(f"Fabry-Perot Airy function (F={F:.1f})")

axes[1, 1].plot(delta_k_range * Lcoh / np.pi, eta)
axes[1, 1].set_xlabel("Delta k * L_coh / pi"); axes[1, 1].set_ylabel("relative SHG efficiency")
axes[1, 1].set_title("SHG phase-matching sinc^2 curve")

fig.tight_layout()
fig.savefig("phase3_photonics_demo.png", dpi=140)
print("\nSaved figure: phase3_photonics_demo.png")

