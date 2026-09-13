"""
opticspy.sources_detectors.photodetectors
============================================
Core photodetector figures of merit: responsivity, quantum efficiency,
photocurrent, and the two fundamental noise mechanisms (shot noise,
Johnson/thermal noise) that set the noise-equivalent power (NEP).

Reference: Saleh & Teich, "Fundamentals of Photonics", Ch. 18
(Photodetectors).
"""

import numpy as np

Q_ELECTRON = 1.602176634e-19  # C
H_PLANCK = 6.62607015e-34     # J*s
C_LIGHT = 299_792_458.0       # m/s
K_BOLTZMANN = 1.380649e-23    # J/K


def photon_energy(wavelength):
    return H_PLANCK * C_LIGHT / wavelength


def responsivity(quantum_efficiency, wavelength):
    """Responsivity R = eta * q * wavelength / (h*c), in A/W."""
    return quantum_efficiency * Q_ELECTRON * wavelength / (H_PLANCK * C_LIGHT)


def photocurrent(optical_power, quantum_efficiency, wavelength):
    return responsivity(quantum_efficiency, wavelength) * optical_power


def quantum_efficiency_from_absorption(absorption_coeff, thickness, surface_reflectance=0.0):
    """eta = (1 - R_surface) * (1 - exp(-alpha*d)), the standard absorbed-fraction estimate."""
    return (1 - surface_reflectance) * (1 - np.exp(-absorption_coeff * thickness))


def shot_noise_current(mean_current, bandwidth):
    """RMS shot-noise current: sqrt(2*q*I*B)."""
    return np.sqrt(2 * Q_ELECTRON * mean_current * bandwidth)


def thermal_noise_current(temperature, bandwidth, load_resistance):
    """RMS Johnson/thermal-noise current: sqrt(4*kB*T*B / R_load)."""
    return np.sqrt(4 * K_BOLTZMANN * temperature * bandwidth / load_resistance)


def total_noise_current(mean_current, bandwidth, temperature, load_resistance):
    """Combine shot and thermal noise in quadrature (independent noise sources)."""
    i_shot = shot_noise_current(mean_current, bandwidth)
    i_thermal = thermal_noise_current(temperature, bandwidth, load_resistance)
    return np.sqrt(i_shot ** 2 + i_thermal ** 2)


def noise_equivalent_power(noise_current, responsivity_A_per_W):
    """NEP = noise_current / responsivity, i.e. the optical power that gives SNR=1."""
    return noise_current / responsivity_A_per_W


def signal_to_noise_ratio(optical_power, quantum_efficiency, wavelength,
                           bandwidth, temperature, load_resistance):
    R = responsivity(quantum_efficiency, wavelength)
    I_signal = R * optical_power
    I_noise = total_noise_current(I_signal, bandwidth, temperature, load_resistance)
    return I_signal / I_noise
