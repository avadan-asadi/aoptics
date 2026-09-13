"""
opticspy.coherence.speckle
============================
Statistics of laser speckle: the intensity of "fully developed" speckle
(many independent unit-amplitude, uniformly-random-phase contributions,
by the central limit theorem) is a circular complex Gaussian random
variable, whose intensity I = |sum of phasors|^2 follows a negative
exponential distribution. Summing N independent, fully-developed,
same-polarization speckle patterns gives a Gamma-distributed intensity
with reduced contrast 1/sqrt(N).

Reference: Goodman, "Statistical Optics", Ch. 3 & 9 (also his dedicated
monograph "Speckle Phenomena in Optics").
"""

import numpy as np


def speckle_intensity_pdf(I, mean_I=1.0):
    """Negative-exponential PDF of fully developed (polarized) speckle intensity."""
    I = np.asarray(I, dtype=float)
    p = np.zeros_like(I)
    pos = I >= 0
    p[pos] = np.exp(-I[pos] / mean_I) / mean_I
    return p


def speckle_intensity_cdf(I, mean_I=1.0):
    I = np.asarray(I, dtype=float)
    cdf = np.zeros_like(I)
    pos = I >= 0
    cdf[pos] = 1.0 - np.exp(-I[pos] / mean_I)
    return cdf


def speckle_contrast(intensity_samples):
    """Empirical speckle contrast C = std(I)/mean(I) (1.0 for ideal fully developed speckle)."""
    I = np.asarray(intensity_samples, dtype=float)
    return float(np.std(I) / np.mean(I))


def summed_speckle_pdf(I, N, mean_I=1.0):
    """
    PDF of the intensity of the incoherent sum of N independent,
    fully-developed, equal-mean, same-polarization speckle patterns:
    a Gamma(N, mean_I/N) distribution,
        p(I) = (N/mean_I)^N * I^(N-1) * exp(-N I / mean_I) / (N-1)!
    Contrast of the sum = 1/sqrt(N).
    """
    from scipy.special import gamma as gamma_fn
    I = np.asarray(I, dtype=float)
    p = np.zeros_like(I)
    pos = I > 0
    p[pos] = ((N / mean_I) ** N) * I[pos] ** (N - 1) * np.exp(-N * I[pos] / mean_I) / gamma_fn(N)
    return p


def summed_speckle_contrast(N):
    """Theoretical contrast of the sum of N independent fully-developed speckle patterns."""
    return 1.0 / np.sqrt(N)


def simulate_speckle_intensity(n_scatterers, n_samples, rng=None, batch_size=20_000):
    """
    Monte Carlo random-phasor-sum simulation of fully developed speckle:
    for each of `n_samples` independent realizations, sum `n_scatterers`
    unit-amplitude phasors with independent uniform random phase in
    [0, 2*pi) and record the resulting intensity |sum|^2. In the limit
    of large n_scatterers this converges (by the central limit theorem)
    to the negative-exponential distribution of `speckle_intensity_pdf`.

    Processed in batches of `batch_size` samples to bound peak memory
    use for large Monte Carlo runs.
    """
    rng = np.random.default_rng() if rng is None else rng
    out = np.empty(n_samples)
    done = 0
    while done < n_samples:
        n = min(batch_size, n_samples - done)
        phases = rng.uniform(0, 2 * np.pi, size=(n, n_scatterers))
        field = np.sum(np.exp(1j * phases), axis=1)
        out[done:done + n] = np.abs(field) ** 2
        done += n
    return out


def speckle_size_free_space(wavelength, distance, aperture_diameter):
    """
    'Objective' speckle correlation width observed a distance `distance`
    from a rough surface illuminated coherently over an aperture of
    diameter `aperture_diameter` -- same diffraction-limited form as the
    Van Cittert-Zernike coherence radius: 1.22 * wavelength * distance / D.
    """
    return 1.22 * wavelength * distance / aperture_diameter


def speckle_size_imaging(wavelength, f_number, magnification=0.0):
    """
    'Subjective' (imaging-system) speckle correlation width at the image
    plane of an optical system with the given working f-number and
    magnification: 1.22 * (1 + magnification) * wavelength * f_number.
    """
    return 1.22 * (1.0 + magnification) * wavelength * f_number
