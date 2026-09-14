import numpy as np
import pytest

from aoptics.scattering import rayleigh_scattering as rs, metal_optics as mo


# --------------------------- Rayleigh scattering ---------------------------

def test_rayleigh_cross_section_wavelength_scaling():
    d, n = 50e-9, 1.5
    sigma_blue = rs.rayleigh_cross_section(450e-9, d, n)
    sigma_red = rs.rayleigh_cross_section(650e-9, d, n)
    ratio = sigma_blue / sigma_red
    expected_ratio = (650e-9 / 450e-9) ** 4
    assert ratio == pytest.approx(expected_ratio, rel=1e-6)


def test_rayleigh_cross_section_size_scaling():
    wl, n = 500e-9, 1.5
    sigma1 = rs.rayleigh_cross_section(wl, 20e-9, n)
    sigma2 = rs.rayleigh_cross_section(wl, 40e-9, n)
    assert sigma2 / sigma1 == pytest.approx(2 ** 6, rel=1e-6)


def test_phase_function_normalized():
    theta = np.linspace(0, np.pi, 20000)
    f = rs.rayleigh_phase_function(theta)
    integral = np.trapezoid(f * 2 * np.pi * np.sin(theta), theta)  # integrate over sphere
    assert integral == pytest.approx(1.0, rel=1e-3)


def test_polarization_degree_max_at_90deg():
    thetas = np.linspace(0.01, np.pi - 0.01, 500)
    dop = rs.rayleigh_polarization_degree(thetas)
    assert thetas[np.argmax(dop)] == pytest.approx(np.pi / 2, abs=0.01)
    assert np.max(dop) == pytest.approx(1.0, rel=1e-3)


def test_polarization_zero_forward_and_backward():
    assert rs.rayleigh_polarization_degree(0.0) == pytest.approx(0.0, abs=1e-9)
    assert rs.rayleigh_polarization_degree(np.pi) == pytest.approx(0.0, abs=1e-9)


def test_relative_scattering_vs_wavelength_blue_stronger_than_red():
    wl = np.array([450e-9, 550e-9, 650e-9])
    rel = rs.relative_scattering_vs_wavelength(wl, reference_wavelength=550e-9)
    assert rel[0] > rel[1] > rel[2]
    assert rel[1] == pytest.approx(1.0)


# --------------------------- metal optics ---------------------------

def test_normal_incidence_reflectance_lossless_limit_matches_dielectric():
    from aoptics.electromagnetic import fresnel
    n = 1.5
    R_metal_form = mo.normal_incidence_reflectance(n, 0.0, n_incident=1.0)
    R_dielectric = fresnel.normal_incidence_reflectance(1.0, n)
    assert R_metal_form == pytest.approx(R_dielectric)


def test_metal_reflectance_high_for_large_k():
    R = mo.normal_incidence_reflectance(0.2, 3.5)  # typical metal-like n,k at visible
    assert R > 0.9


def test_skin_depth_decreases_with_k():
    d1 = mo.skin_depth(500e-9, 2.0)
    d2 = mo.skin_depth(500e-9, 4.0)
    assert d2 == pytest.approx(d1 / 2)


def test_absorptance_plus_reflectance_is_one_at_normal_incidence():
    n, k = 0.2, 3.5
    R = mo.normal_incidence_reflectance(n, k)
    A = mo.absorptance(0.0, n, k, polarization='unpolarized')
    assert R + A == pytest.approx(1.0, rel=1e-6)

