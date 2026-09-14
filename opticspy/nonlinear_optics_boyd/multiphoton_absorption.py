"""
aoptics.nonlinear_optics_boyd.multiphoton_absorption
=========================================================
Two-photon absorption (TPA): the nonlinear (intensity-squared) term
added to the Beer-Lambert propagation equation, its effective
interaction length, and the resulting intensity-dependent transmission
-- plus a simple saturable-absorber transmission model (the opposite
sign effect, used e.g. for passive mode-locking).

Reference: Boyd, "Nonlinear Optics", Ch. 12 (Optically Induced Damage
and Multiphoton Absorption).
"""

import numpy as np


def tpa_propagation_rhs(intensity, alpha_linear, beta_tpa):
    """dI/dz = -alpha*I - beta*I^2 (linear + two-photon absorption)."""
    return -alpha_linear * intensity - beta_tpa * intensity ** 2


def intensity_after_tpa(intensity_in, beta_tpa, length, alpha_linear=0.0):
    """
    Closed-form solution of dI/dz = -alpha*I - beta*I^2 for constant
    alpha (exact analytic solution of the linear + TPA Riccati equation):
        for alpha=0: I(L) = I0 / (1 + beta*I0*L)
        for alpha!=0: I(L) = I0*exp(-alpha*L) / (1 + (beta*I0/alpha)*(1-exp(-alpha*L)))
    """
    if alpha_linear == 0:
        return intensity_in / (1 + beta_tpa * intensity_in * length)
    Leff = (1 - np.exp(-alpha_linear * length)) / alpha_linear
    return intensity_in * np.exp(-alpha_linear * length) / (1 + beta_tpa * intensity_in * Leff)


def tpa_effective_length(length, alpha_linear):
    """Effective length for the TPA term under linear background loss: Leff = (1-exp(-aL))/a."""
    if alpha_linear == 0:
        return length
    return (1 - np.exp(-alpha_linear * length)) / alpha_linear


def nonlinear_transmission_tpa(intensity_in, beta_tpa, length, alpha_linear=0.0):
    I_out = intensity_after_tpa(intensity_in, beta_tpa, length, alpha_linear)
    linear_only = intensity_in * np.exp(-alpha_linear * length)
    relative_to_linear = I_out / linear_only if linear_only > 0 else np.nan
    return I_out / intensity_in, relative_to_linear


def saturable_absorber_transmission(intensity, alpha0, saturation_intensity,
                                     nonsaturable_loss=0.0):
    """
    Fast saturable absorber transmission (opposite sign nonlinearity
    from TPA -- absorption DECREASES with intensity), used e.g. to
    model passive mode-locking elements:
        T(I) = 1 - nonsaturable_loss - alpha0 / (1 + I/Isat)
    """
    return 1 - nonsaturable_loss - alpha0 / (1 + intensity / saturation_intensity)


def two_photon_absorption_cross_section_rate(beta_tpa, photon_energy, n_density):
    """
    Two-photon transition rate per molecule/atom, from the bulk TPA
    coefficient beta (m/W) and number density n_density (1/m^3):
        R_2photon = beta * I^2 / (2 * h_nu * n_density)   [returns the
    rate COEFFICIENT per unit I^2; multiply by I^2 for the actual rate]
        coefficient = beta / (2 * photon_energy * n_density)
    """
    return beta_tpa / (2 * photon_energy * n_density)

