"""
Phase 5 demo -- remaining chapters: Fourier optics, lasers/detectors,
modulation, scattering, and photonic crystals.
=====================================================================
(1) Fourier optics: diffraction-limited PSF/MTF of a circular pupil.
(2) Laser rate equations: photon number turning on above threshold.
(3) Acousto-optics: Bragg regime diffraction geometry.
(4) Rayleigh scattering: 1/lambda^4 law (why the sky is blue).
(5) Photonic crystal: quarter-wave Bragg mirror reflectance spectrum.

Run:  python examples/phase5_extended_demo.py
"""
import numpy as np
import matplotlib.pyplot as plt

from aoptics.fourier_optics import transfer_functions as tf
from aoptics.sources_detectors import lasers
from aoptics.scattering import rayleigh_scattering as rs
from aoptics.photonics import photonic_crystals as pc

# ---- 1) Diffraction-limited PSF / MTF ----
N, L = 256, 1.0
aperture_radius = 0.25
pupil = tf.pupil_function_circular(N, L, aperture_radius)
psf = tf.incoherent_psf(pupil)
mtf = tf.modulation_transfer_function(pupil)

NA, wl = 0.6, 550e-9
print(f"Diffraction-limited resolution (NA={NA}): "
      f"Rayleigh={tf.rayleigh_resolution(wl, NA)*1e9:.0f} nm, "
      f"Abbe={tf.abbe_resolution(wl, NA)*1e9:.0f} nm")

# ---- 2) Laser turn-on above threshold ----
tau_2, tau_photon, gain_coeff = 1e-9, 1e-11, 1e3
R_th = lasers.threshold_pump_rate(tau_2, tau_photon, gain_coeff)
t, N_pop, phi = lasers.simulate_laser(3 * R_th, tau_2, tau_photon, gain_coeff,
                                       phi0=1.0, t_max=60 * tau_2, n_steps=3000)
phi_ss = lasers.steady_state_photon_number(3 * R_th, tau_2, tau_photon, gain_coeff)
print(f"\nLaser at 3x threshold: predicted steady-state photon number = {phi_ss:.3e}, "
      f"simulated final value = {phi[-1]:.3e}")

# ---- 3) Rayleigh scattering: blue sky ----
wavelengths = np.linspace(400e-9, 700e-9, 200)
rel_scatter = rs.relative_scattering_vs_wavelength(wavelengths, reference_wavelength=550e-9)

# ---- 4) Photonic crystal Bragg mirror ----
n_high, n_low, wl0 = 2.3, 1.45, 550e-9
layers = pc.quarter_wave_stack(n_high, n_low, wl0, n_periods=12)
wl_scan = np.linspace(400e-9, 750e-9, 300)
R_stack = pc.reflectance_spectrum(layers, wl_scan, n_substrate=1.5)
lo, hi = pc.bragg_stopband_edges(n_high, n_low, wl0)
print(f"\nQuarter-wave Bragg mirror (n_H={n_high}, n_L={n_low}, 12 periods): "
      f"predicted stopband {lo*1e9:.0f}-{hi*1e9:.0f} nm, "
      f"peak reflectance = {np.max(R_stack)*100:.2f}%")

# ---- Figure ----
fig, axes = plt.subplots(2, 2, figsize=(11, 8.5))

extent = [-L, L, -L, L]
axes[0, 0].imshow(psf[100:156, 100:156], extent=[-0.1, 0.1, -0.1, 0.1], cmap='inferno')
axes[0, 0].set_title("Diffraction-limited PSF (Airy pattern)")
axes[0, 0].set_xlabel("x (normalized)")

freqs = np.linspace(-1, 1, N)
axes[0, 1].plot(freqs, mtf[N // 2, :])
axes[0, 1].set_xlabel("spatial frequency (normalized to cutoff)")
axes[0, 1].set_ylabel("MTF")
axes[0, 1].set_title("Modulation transfer function")

axes[1, 0].plot(t / tau_2, phi)
axes[1, 0].axhline(phi_ss, color='r', ls='--', lw=1, label='predicted steady state')
axes[1, 0].set_xlabel("t / tau_2"); axes[1, 0].set_ylabel("cavity photon number")
axes[1, 0].set_title("Laser turn-on dynamics (3x threshold)")
axes[1, 0].legend(fontsize=8)

axes[1, 1].plot(wl_scan * 1e9, R_stack * 100, label='Bragg mirror R')
axes[1, 1].axvspan(lo * 1e9, hi * 1e9, alpha=0.15, color='gray', label='predicted stopband')
ax2 = axes[1, 1].twinx()
ax2.plot(wavelengths * 1e9, rel_scatter, color='tab:blue', alpha=0.5, label='Rayleigh scattering (rel.)')
axes[1, 1].set_xlabel("wavelength (nm)"); axes[1, 1].set_ylabel("reflectance (%)")
ax2.set_ylabel("relative Rayleigh scattering")
axes[1, 1].set_title("Bragg mirror spectrum & Rayleigh 1/lambda^4 law")
axes[1, 1].legend(fontsize=7, loc='upper left')

fig.tight_layout()
fig.savefig("phase5_extended_demo.png", dpi=140)
print("\nSaved figure: phase5_extended_demo.png")

