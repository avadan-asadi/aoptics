"""
aoptics.nonlinear_optics_boyd.stimulated_scattering
========================================================
Stimulated Raman and stimulated Brillouin scattering: frequency shifts,
gain, and the standard exponential-gain threshold criteria
(gR*Pth*Leff/Aeff ~ 16 for SRS, gB*Pth*Leff/Aeff ~ 21 for SBS) widely
used in fiber and bulk nonlinear optics (R.G. Smith's threshold
criterion, reproduced throughout the nonlinear-fiber-optics literature).

Reference: Boyd, "Nonlinear Optics", Ch. 8 (Spontaneous Light
Scattering and Acousto-Optics) and Ch. 10 (Stimulated Raman Scattering
and Other Third-Order Parametric Processes); Agrawal, "Nonlinear Fiber
Optics", Ch. 8-9.
"""

import numpy as np


def effective_length(length, alpha):
    """Effective interaction length accounting for linear loss: Leff = (1-exp(-alpha*L))/alpha."""
    if alpha == 0:
        return length
    return (1 - np.exp(-alpha * length)) / alpha


def stokes_frequency_shift(pump_frequency, vibrational_frequency):
    """Raman Stokes line frequency: omega_S = omega_pump - omega_vibrational."""
    return pump_frequency - vibrational_frequency


def anti_stokes_frequency_shift(pump_frequency, vibrational_frequency):
    """Raman anti-Stokes line frequency: omega_AS = omega_pump + omega_vibrational."""
    return pump_frequency + vibrational_frequency


def raman_gain(pump_intensity, gain_coefficient, length):
    """Exponential Stokes gain in the undepleted-pump limit: G = exp(gR*Ip*L)."""
    return np.exp(gain_coefficient * pump_intensity * length)


def srs_threshold_power(gain_coefficient, effective_area, effective_length_value,
                         threshold_gain=16.0):
    """
    Stimulated Raman scattering threshold pump power, from the R.G.
    Smith criterion gR*Pth*Leff/Aeff ~ 16:
        Pth = threshold_gain * Aeff / (gR * Leff)
    """
    return threshold_gain * effective_area / (gain_coefficient * effective_length_value)


def sbs_threshold_power(gain_coefficient, effective_area, effective_length_value,
                         pump_linewidth=None, brillouin_linewidth=None, threshold_gain=21.0):
    """
    Stimulated Brillouin scattering threshold pump power:
        Pth = threshold_gain * Aeff / (gB * Leff)
    optionally scaled by (1 + Delta_nu_pump/Delta_nu_B) when the pump
    linewidth is not negligible compared to the (typically much
    narrower) Brillouin gain linewidth (both standard refinements of
    the basic ~21 criterion).
    """
    Pth = threshold_gain * effective_area / (gain_coefficient * effective_length_value)
    if pump_linewidth is not None and brillouin_linewidth is not None:
        Pth *= (1 + pump_linewidth / brillouin_linewidth)
    return Pth


def brillouin_frequency_shift(refractive_index, acoustic_velocity, wavelength, theta=np.pi):
    """
    Brillouin frequency shift for scattering angle theta (theta=pi is
    pure backscatter, the dominant SBS geometry in fibers/bulk media):
        nu_B = 2*n*V_a*sin(theta/2) / wavelength
    """
    return 2 * refractive_index * acoustic_velocity * np.sin(theta / 2) / wavelength


def acoustic_phonon_lifetime(brillouin_linewidth):
    """Acoustic phonon (damping) lifetime tau_B = 1/(pi*Delta_nu_B) (from the
    Brillouin gain linewidth, Lorentzian FWHM Delta_nu_B)."""
    return 1.0 / (np.pi * brillouin_linewidth)


def raman_shift_wavenumber_to_frequency(wavenumber_cm_inv, c=299_792_458.0):
    """Convert a Raman shift quoted in wavenumbers (cm^-1, as is conventional
    in the Raman-spectroscopy literature) to a frequency shift in Hz: f = c * (wavenumber in 1/m)."""
    return c * (wavenumber_cm_inv * 100.0)

