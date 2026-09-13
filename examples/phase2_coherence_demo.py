"""
Phase 2 demo -- opticspy.coherence
====================================
(1) Temporal coherence: Lorentzian-linewidth source -> Michelson
    interferogram envelope decaying with delay, coherence length.
(2) Spatial coherence: Van Cittert-Zernike visibility of a Young's
    double slit under an extended incoherent (thermal) source.
(3) Speckle statistics: Monte Carlo random-phasor-sum simulation vs
    the analytic negative-exponential intensity distribution, and the
    contrast reduction from averaging N independent speckle patterns.

Run:  python examples/phase2_coherence_demo.py
"""
import numpy as np
import matplotlib.pyplot as plt

from opticspy.coherence import temporal, spatial, speckle

# ---- 1) Temporal coherence: a sodium-lamp-like Lorentzian line ----
wavelength0 = 589e-9
nu0 = temporal.C_LIGHT / wavelength0
fwhm = 5e11  # 500 GHz linewidth (illustrative, broad enough to see quickly)
nu = np.linspace(nu0 - 200 * fwhm, nu0 + 200 * fwhm, 4000)
S = temporal.lorentzian_spectrum(nu, nu0, fwhm)

tau_c_expected = 1.0 / (np.pi * fwhm)
L_c = temporal.coherence_length(tau_c_expected)
print(f"Lorentzian source: FWHM = {fwhm:.3e} Hz")
print(f"  coherence time  tau_c = {tau_c_expected*1e15:.2f} fs")
print(f"  coherence length L_c  = {L_c*1e6:.2f} micron")

delay = np.linspace(-8 * tau_c_expected, 8 * tau_c_expected, 4000)
interferogram = temporal.interferogram(delay, nu, S, nu0=nu0)
_, gamma_delay = temporal.complex_degree_of_coherence(nu, S, tau=delay)
tau_c_numeric = temporal.coherence_time(delay, gamma_delay)
print(f"  coherence time from numerical Wiener-Khinchin integral = {tau_c_numeric*1e15:.2f} fs "
      f"(cross-check against analytic 1/(pi*FWHM))")

# ---- 2) Spatial coherence: Young's double slit, extended thermal source ----
wavelength = 550e-9
distance = 2.0      # source-to-slits distance (m)
source_diam = 0.5e-3  # incoherent source diameter (m)
rho_c = spatial.coherence_radius_circular_source(wavelength, distance, source_diam)
print(f"\nExtended circular source (D={source_diam*1e3:.2f} mm) at {distance} m: "
      f"transverse coherence radius = {rho_c*1e3:.3f} mm")

separations = np.linspace(0, 3 * rho_c, 300)
visibility = spatial.young_double_slit_visibility(separations, wavelength, distance, source_diam)

# ---- 3) Speckle statistics ----
rng = np.random.default_rng(7)
samples = speckle.simulate_speckle_intensity(n_scatterers=300, n_samples=60_000, rng=rng)
samples = samples / np.mean(samples)
print(f"\nMonte Carlo speckle (300 scatterers, 60000 samples): "
      f"contrast = {speckle.speckle_contrast(samples):.3f} (theory: 1.0)")

N_avg = 8
avg_samples = np.mean(
    speckle.simulate_speckle_intensity(300, 60_000 * N_avg, rng=rng).reshape(-1, N_avg), axis=1)
avg_samples = avg_samples / np.mean(avg_samples)
print(f"Averaging N={N_avg} independent speckle patterns: "
      f"contrast = {speckle.speckle_contrast(avg_samples):.3f} "
      f"(theory: {speckle.summed_speckle_contrast(N_avg):.3f})")

# ---- Figure ----
fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.2))

axes[0].plot(delay * 1e15, interferogram, lw=0.8, color='C0', alpha=0.5)
period = 1.0 / nu0
delay_zoom = np.linspace(-6 * period, 6 * period, 2000)
zoom_interferogram = temporal.interferogram(delay_zoom, nu, S, nu0=nu0)
_, gamma_zoom = temporal.complex_degree_of_coherence(nu, S, tau=delay_zoom)
inset = axes[0].inset_axes([0.62, 0.55, 0.35, 0.4])
inset.plot(delay_zoom * 1e15, zoom_interferogram, lw=0.8)
inset.set_title("zoom: carrier fringes", fontsize=7)
inset.tick_params(labelsize=6)
axes[0].plot(delay * 1e15, 1 + np.abs(gamma_delay), 'r--', lw=1, label='envelope 1+|gamma|')
axes[0].plot(delay * 1e15, 1 - np.abs(gamma_delay), 'r--', lw=1)
axes[0].set_xlabel("delay (fs)"); axes[0].set_ylabel("interferogram I/I0")
axes[0].set_title("Michelson fringes, Lorentzian source")
axes[0].legend(fontsize=8)

axes[1].plot(separations * 1e3, visibility)
axes[1].axvline(rho_c * 1e3, color='gray', ls='--', lw=1, label='coherence radius')
axes[1].set_xlabel("slit separation (mm)"); axes[1].set_ylabel("fringe visibility")
axes[1].set_title("Young double slit, extended source (VCZ)")
axes[1].legend(fontsize=8)

bins = np.linspace(0, 6, 60)
axes[2].hist(samples, bins=bins, density=True, alpha=0.5, label='N=1 (Monte Carlo)')
axes[2].hist(avg_samples, bins=bins, density=True, alpha=0.5, label=f'N={N_avg} averaged')
I_plot = np.linspace(0.001, 6, 300)
axes[2].plot(I_plot, speckle.speckle_intensity_pdf(I_plot, 1.0), 'k-', lw=1.2, label='theory N=1')
axes[2].plot(I_plot, speckle.summed_speckle_pdf(I_plot, N_avg, 1.0), 'k--', lw=1.2, label=f'theory N={N_avg}')
axes[2].set_xlabel("I / <I>"); axes[2].set_ylabel("probability density")
axes[2].set_title("Speckle intensity statistics")
axes[2].legend(fontsize=7)

fig.tight_layout()
fig.savefig("phase2_coherence_demo.png", dpi=140)
print("\nSaved figure: phase2_coherence_demo.png")
