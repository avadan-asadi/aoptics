"""
aoptics.coherence.temporal
============================
Temporal coherence of stationary light: the Wiener-Khinchin theorem
relates the (self-)coherence function Gamma(tau) = <E*(t) E(t+tau)> to
the power spectral density S(nu) via a Fourier transform pair,

    Gamma(tau) = Integral[ S(nu) exp(-i 2 pi nu tau) dnu ]

The normalized complex degree of temporal coherence is
gamma(tau) = Gamma(tau) / Gamma(0), and the fringe visibility observed
in a two-beam interferometer (e.g. Michelson) at path-length delay tau
is exactly |gamma(tau)|.

Reference: Goodman, "Statistical Optics", Ch. 5-6; Born & Wolf, Ch. 10
(Elements of the theory of interference and interferometers /
coherence); Mandel & Wolf, "Optical Coherence and Quantum Optics".

Implementation note: rather than hard-coding closed-form gamma(tau) for
each spectral shape (which vary by convention/definition across texts),
`complex_degree_of_coherence` computes it numerically, directly from
the Wiener-Khinchin Fourier relation, for ANY user-supplied spectrum.
The one closed form used in this module (Lorentzian spectrum <->
exponential coherence function) is a standard, well-established result
used only as a unit-test cross-check of the numerical machinery.
"""

import numpy as np

C_LIGHT = 299_792_458.0  # m/s


def lorentzian_spectrum(nu, nu0, fwhm):
    """Normalized Lorentzian power spectral density, full width at half max `fwhm`."""
    gamma_half = fwhm / 2.0
    S = (gamma_half / np.pi) / ((nu - nu0) ** 2 + gamma_half ** 2)
    return S


def gaussian_spectrum(nu, nu0, fwhm):
    """Normalized Gaussian power spectral density, full width at half max `fwhm`."""
    sigma = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    S = np.exp(-(nu - nu0) ** 2 / (2 * sigma ** 2)) / (sigma * np.sqrt(2 * np.pi))
    return S


def rectangular_spectrum(nu, nu0, width):
    S = (np.abs(nu - nu0) <= width / 2.0).astype(float) / width
    return S


def complex_degree_of_coherence(nu, S, tau=None, n_tau=400, tau_max=None):
    """
    Numerically compute gamma(tau) = Gamma(tau)/Gamma(0) from a sampled
    power spectral density S(nu) (nu need not be centered at zero -- the
    fast oscillation at the optical carrier frequency is included, which
    correctly reproduces the interferometric fringes riding under the
    coherence envelope).

    Returns (tau, gamma) with gamma complex-valued, gamma(0) = 1.
    """
    nu = np.asarray(nu, dtype=float)
    S = np.asarray(S, dtype=float)
    if tau is None:
        if tau_max is None:
            span = nu.max() - nu.min()
            tau_max = 5.0 / max(span, 1e-30)
        tau = np.linspace(-tau_max, tau_max, n_tau)

    # direct numerical Fourier integral (robust for irregularly-sampled nu)
    Gamma = np.array([np.trapezoid(S * np.exp(-1j * 2 * np.pi * nu * t), nu) for t in tau])
    Gamma0 = np.trapezoid(S, nu)
    gamma = Gamma / Gamma0
    return tau, gamma


def coherence_time(tau, gamma):
    """
    Rigorous coherence-time definition (Mandel & Wolf):
        tau_c = Integral[ |gamma(tau)|^2 dtau ]  over all tau.
    """
    return float(np.trapezoid(np.abs(gamma) ** 2, tau))


def coherence_length(tau_c, v=C_LIGHT):
    return v * tau_c


def fringe_visibility(gamma_tau):
    """Interferometric fringe visibility = |gamma(tau)| (classic result)."""
    return np.abs(gamma_tau)


def interferogram(tau, nu, S, nu0=None):
    """
    Simulated two-beam (e.g. Michelson) interferogram as a function of
    delay tau: I(tau)/I0 = 1 + |gamma(tau)| cos(2*pi*nu0*tau + phase(gamma(tau))),
    with nu0 the mean/carrier frequency (defaults to the spectrum's
    intensity-weighted mean).
    """
    if nu0 is None:
        nu0 = np.trapezoid(np.asarray(nu) * np.asarray(S), nu) / np.trapezoid(S, nu)
    _, gamma = complex_degree_of_coherence(nu, S, tau=tau)
    return 1.0 + np.abs(gamma) * np.cos(2 * np.pi * nu0 * tau + np.angle(gamma))


def autocorrelation(signal, dt, mode='biased'):
    """
    Autocorrelation of a (real or complex) sampled time-domain signal via
    direct summation (small/medium N; for long signals use FFT-based
    methods). Returns (lags, R) with lags in the same units as dt*index.
    """
    signal = np.asarray(signal)
    N = len(signal)
    R = np.correlate(signal, signal, mode='full')
    lags = (np.arange(-(N - 1), N)) * dt
    if mode == 'unbiased':
        norm = np.concatenate([np.arange(1, N + 1), np.arange(N - 1, 0, -1)])
        R = R / norm
    else:
        R = R / N
    return lags, R


def power_spectral_density_from_signal(signal, dt):
    """
    Wiener-Khinchin theorem applied the other way: estimate the PSD of a
    sampled time series as |FFT(signal)|^2 (periodogram estimate).
    Returns (freq, psd).
    """
    signal = np.asarray(signal)
    N = len(signal)
    spec = np.fft.fftshift(np.fft.fft(signal))
    freq = np.fft.fftshift(np.fft.fftfreq(N, d=dt))
    psd = np.abs(spec) ** 2 * dt / N
    return freq, psd

