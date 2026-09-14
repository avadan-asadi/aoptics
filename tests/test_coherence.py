import numpy as np
import pytest

from aoptics.coherence import temporal, spatial, speckle


# --------------------------- temporal coherence ---------------------------

def test_lorentzian_spectrum_gives_exponential_coherence_function():
    # Classic result: Lorentzian spectrum of FWHM Dnu <-> |gamma(tau)| = exp(-pi*Dnu*|tau|)
    nu0 = 5e14
    fwhm = 1e9  # 1 GHz linewidth
    span = 200 * fwhm
    nu = np.linspace(nu0 - span, nu0 + span, 200_000)
    S = temporal.lorentzian_spectrum(nu, nu0, fwhm)
    tau_max = 5 / fwhm
    tau, gamma = temporal.complex_degree_of_coherence(nu, S, tau=np.linspace(-tau_max, tau_max, 300))
    # remove the fast carrier oscillation to compare envelopes: multiply gamma by exp(+i*2*pi*nu0*tau)
    envelope = gamma * np.exp(1j * 2 * np.pi * nu0 * tau)
    expected_envelope = np.exp(-np.pi * fwhm * np.abs(tau))
    assert np.allclose(np.abs(envelope), expected_envelope, atol=2e-2)


def test_gamma_at_zero_delay_is_unity():
    nu0 = 5e14
    nu = np.linspace(nu0 - 5e9, nu0 + 5e9, 20000)
    S = temporal.gaussian_spectrum(nu, nu0, fwhm=1e9)
    tau, gamma = temporal.complex_degree_of_coherence(nu, S, tau=np.array([0.0]))
    assert np.abs(gamma[0]) == pytest.approx(1.0, rel=1e-6)


def test_coherence_time_lorentzian_matches_known_formula():
    # For |gamma(tau)|=exp(-pi*Dnu*|tau|), tau_c = Integral|gamma|^2 dtau = 1/(pi*Dnu)
    fwhm = 1e9
    tau = np.linspace(-50 / fwhm, 50 / fwhm, 200000)
    gamma_env = np.exp(-np.pi * fwhm * np.abs(tau))
    tau_c = temporal.coherence_time(tau, gamma_env)
    expected = 1.0 / (np.pi * fwhm)
    assert tau_c == pytest.approx(expected, rel=1e-2)


def test_fringe_visibility_is_modulus_of_gamma():
    g = np.array([1 + 0j, 0.5 - 0.3j, -0.2j])
    assert np.allclose(temporal.fringe_visibility(g), np.abs(g))


def test_power_spectral_density_of_pure_tone_peaks_at_correct_frequency():
    dt = 1e-3
    N = 4096
    t = np.arange(N) * dt
    f0 = 50.0
    signal = np.cos(2 * np.pi * f0 * t)
    freq, psd = temporal.power_spectral_density_from_signal(signal, dt)
    peak_freq = freq[np.argmax(psd)]
    assert abs(peak_freq) == pytest.approx(f0, abs=1.0 / (N * dt) * 2)


# --------------------------- spatial coherence (VCZ) ---------------------------

def test_jinc_zero_value():
    assert spatial.jinc(0.0) == pytest.approx(1.0)


def test_van_cittert_zernike_circular_first_zero_at_coherence_radius():
    wavelength = 550e-9
    distance = 10.0
    D = 0.001
    rho_c = spatial.coherence_radius_circular_source(wavelength, distance, D)
    gamma_at_null = spatial.van_cittert_zernike_circular(D, wavelength, distance, rho_c)
    assert abs(gamma_at_null) < 1e-3


def test_van_cittert_zernike_slit_first_zero_at_coherence_width():
    wavelength = 550e-9
    distance = 10.0
    b = 0.002
    w_c = spatial.coherence_width_slit_source(wavelength, distance, b)
    gamma_at_null = spatial.van_cittert_zernike_slit(b, wavelength, distance, w_c)
    assert abs(gamma_at_null) < 1e-9


def test_van_cittert_zernike_zero_separation_full_coherence():
    assert spatial.van_cittert_zernike_circular(1e-3, 550e-9, 10.0, 0.0) == pytest.approx(1.0)
    assert spatial.van_cittert_zernike_slit(1e-3, 550e-9, 10.0, 0.0) == pytest.approx(1.0)


def test_young_double_slit_visibility_matches_vcz():
    wavelength, distance, D, sep = 550e-9, 5.0, 0.5e-3, 1e-3
    V = spatial.young_double_slit_visibility(sep, wavelength, distance, D, source_shape='circular')
    direct = abs(spatial.van_cittert_zernike_circular(D, wavelength, distance, sep))
    assert V == pytest.approx(direct)


# --------------------------- speckle statistics ---------------------------

def test_speckle_pdf_integrates_to_one():
    I = np.linspace(0, 30, 200000)
    p = speckle.speckle_intensity_pdf(I, mean_I=2.0)
    assert np.trapezoid(p, I) == pytest.approx(1.0, rel=1e-3)


def test_speckle_pdf_mean_matches_parameter():
    I = np.linspace(0, 50, 400000)
    mean_I = 3.0
    p = speckle.speckle_intensity_pdf(I, mean_I=mean_I)
    computed_mean = np.trapezoid(I * p, I)
    assert computed_mean == pytest.approx(mean_I, rel=1e-2)


def test_monte_carlo_speckle_matches_negative_exponential():
    rng = np.random.default_rng(42)
    samples = speckle.simulate_speckle_intensity(n_scatterers=200, n_samples=40_000, rng=rng)
    samples = samples / np.mean(samples)  # normalize to unit mean
    contrast = speckle.speckle_contrast(samples)
    assert contrast == pytest.approx(1.0, abs=0.03)
    # compare empirical histogram to the analytic negative-exponential PDF
    hist, edges = np.histogram(samples, bins=60, range=(0, 6), density=True)
    centers = 0.5 * (edges[1:] + edges[:-1])
    theoretical = speckle.speckle_intensity_pdf(centers, mean_I=1.0)
    assert np.mean(np.abs(hist - theoretical)) < 0.03


def test_summed_speckle_contrast_scaling():
    for N in [1, 4, 16]:
        assert speckle.summed_speckle_contrast(N) == pytest.approx(1.0 / np.sqrt(N))


def test_summed_speckle_pdf_reduces_to_exponential_for_N1():
    I = np.linspace(0.001, 20, 100000)
    p_sum = speckle.summed_speckle_pdf(I, N=1, mean_I=1.0)
    p_exp = speckle.speckle_intensity_pdf(I, mean_I=1.0)
    assert np.allclose(p_sum, p_exp, rtol=1e-6)


def test_summed_speckle_pdf_integrates_to_one():
    I = np.linspace(0.0001, 40, 200000)
    p = speckle.summed_speckle_pdf(I, N=5, mean_I=2.0)
    assert np.trapezoid(p, I) == pytest.approx(1.0, rel=1e-3)


def test_speckle_size_formulas_are_consistent_with_vcz_form():
    wavelength, distance, D = 632.8e-9, 1.0, 5e-3
    s1 = speckle.speckle_size_free_space(wavelength, distance, D)
    s2 = spatial.coherence_radius_circular_source(wavelength, distance, D)
    assert s1 == pytest.approx(s2)

