import numpy as np
import pytest

from opticspy.electromagnetic import fresnel, polarization as pol, crystal_optics as co


# --------------------------- Fresnel ---------------------------

def test_normal_incidence_reflectance_glass():
    R = fresnel.normal_incidence_reflectance(1.0, 1.5)
    assert R == pytest.approx(0.04, abs=1e-3)


def test_brewster_angle_gives_zero_Rp():
    n1, n2 = 1.0, 1.5
    thetaB = fresnel.brewster_angle(n1, n2)
    rt_ = fresnel.reflectance_transmittance(thetaB, n1, n2)
    assert rt_['R_p'] == pytest.approx(0.0, abs=1e-10)


def test_brewster_angle_value():
    thetaB = fresnel.brewster_angle(1.0, 1.5)
    assert np.degrees(thetaB) == pytest.approx(56.31, abs=0.05)


def test_critical_angle_and_tir_total_reflectance():
    n1, n2 = 1.5, 1.0
    thetac = fresnel.critical_angle(n1, n2)
    assert np.degrees(thetac) == pytest.approx(41.81, abs=0.05)
    rt_ = fresnel.reflectance_transmittance(thetac + np.radians(5), n1, n2)
    assert rt_['R_s'] == pytest.approx(1.0, abs=1e-9)
    assert rt_['R_p'] == pytest.approx(1.0, abs=1e-9)
    assert rt_['T_unpolarized'] == pytest.approx(0.0, abs=1e-9)


def test_energy_conservation_below_critical_angle():
    n1, n2 = 1.0, 1.5
    for theta_deg in [0, 15, 30, 45, 60]:
        rt_ = fresnel.reflectance_transmittance(np.radians(theta_deg), n1, n2)
        assert rt_['R_s'] + rt_['T_s'] == pytest.approx(1.0, abs=1e-9)
        assert rt_['R_p'] + rt_['T_p'] == pytest.approx(1.0, abs=1e-9)


# --------------------------- Jones calculus ---------------------------

def test_polarizer_extinction_crossed():
    v = pol.linear_horizontal()
    Jv = pol.linear_polarizer(theta=np.pi / 2)  # vertical polarizer
    out = pol.apply_jones(Jv, v)
    assert pol.jones_intensity(out) == pytest.approx(0.0, abs=1e-12)


def test_malus_law_through_polarizer():
    v = pol.linear_horizontal()
    theta = np.radians(30)
    J = pol.linear_polarizer(theta)
    out = pol.apply_jones(J, v)
    assert pol.jones_intensity(out) == pytest.approx(np.cos(theta) ** 2, rel=1e-9)


def test_quarter_wave_plate_converts_linear_to_circular():
    v = pol.linear_at_angle(np.pi / 4)  # 45 deg linear
    J = pol.quarter_wave_plate(fast_axis_angle=0.0)
    out = pol.apply_jones(J, v)
    stokes = pol.jones_to_stokes(out)
    dop = pol.degree_of_polarization(stokes)
    assert dop == pytest.approx(1.0, rel=1e-9)
    assert abs(stokes[3]) == pytest.approx(stokes[0], rel=1e-6)  # |S3| = S0 -> circular
    assert stokes[1] == pytest.approx(0.0, abs=1e-9)  # S1 = 0 for circular


def test_half_wave_plate_flips_linear_polarization():
    v = pol.linear_at_angle(0.0)  # horizontal
    J = pol.half_wave_plate(fast_axis_angle=np.radians(45))
    out = pol.apply_jones(J, v)
    stokes = pol.jones_to_stokes(out)
    # HWP at 45 deg should rotate horizontal -> vertical
    assert stokes[1] == pytest.approx(-stokes[0], rel=1e-6)


def test_two_crossed_polarizers_with_45deg_between_pass_25pct():
    v = pol.linear_horizontal()
    J1 = pol.linear_polarizer(np.radians(45))
    J2 = pol.linear_polarizer(np.radians(90))
    out = J2 @ (J1 @ v)
    assert pol.jones_intensity(out) == pytest.approx(0.25, rel=1e-6)


def test_jones_to_stokes_intensity_matches():
    v = pol.elliptical(0.8, 0.6, phase=np.pi / 3)
    stokes = pol.jones_to_stokes(v)
    assert stokes[0] == pytest.approx(pol.jones_intensity(v), rel=1e-9)


# --------------------------- Mueller calculus ---------------------------

def test_mueller_polarizer_matches_jones_polarizer_on_unpolarized_light():
    theta = np.radians(20)
    unpolarized = np.array([1.0, 0.0, 0.0, 0.0])
    out = pol.apply_mueller(pol.mueller_polarizer(theta), unpolarized)
    assert out[0] == pytest.approx(0.5, rel=1e-9)  # ideal polarizer passes half of unpolarized light


def test_jones_to_mueller_matches_direct_mueller_for_polarizer():
    theta = np.radians(33)
    M_direct = pol.mueller_polarizer(theta)
    M_from_jones = pol.jones_to_mueller(pol.linear_polarizer(theta))
    assert np.allclose(M_direct, M_from_jones, atol=1e-9)


def test_jones_to_mueller_matches_direct_mueller_for_retarder():
    theta = np.radians(17)
    delta = 1.3
    M_direct = pol.mueller_retarder(delta, theta)
    M_from_jones = pol.jones_to_mueller(pol.wave_plate(delta, theta))
    assert np.allclose(M_direct, M_from_jones, atol=1e-9)


def test_mueller_matrix_preserves_intensity_for_retarder():
    stokes_in = pol.jones_to_stokes(pol.linear_at_angle(np.radians(23)))
    M = pol.mueller_retarder(1.7, np.radians(50))
    out = pol.apply_mueller(M, stokes_in)
    assert out[0] == pytest.approx(stokes_in[0], rel=1e-9)


# --------------------------- Crystal optics ---------------------------

def test_extraordinary_index_limits():
    no, ne = 1.65, 1.49
    assert co.extraordinary_index(0.0, no, ne) == pytest.approx(no, rel=1e-9)
    assert co.extraordinary_index(np.pi / 2, no, ne) == pytest.approx(ne, rel=1e-9)


def test_walk_off_zero_at_axis_and_perpendicular():
    no, ne = 1.65, 1.49
    assert co.walk_off_angle(0.0, no, ne) == pytest.approx(0.0, abs=1e-9)
    assert co.walk_off_angle(np.pi / 2, no, ne) == pytest.approx(0.0, abs=1e-9)


def test_waveplate_thickness_gives_quarter_wave_retardation():
    no, ne = 1.5443, 1.5534  # quartz
    wl = 589e-9
    d = co.waveplate_thickness(no, ne, wl, order=0, fraction=0.25)
    gamma = co.phase_retardation(no, ne, d, wl)
    assert gamma == pytest.approx(np.pi / 2, rel=1e-9)


def test_birefringence_sign_matches_crystal_type():
    d = co.COMMON_UNIAXIAL_CRYSTALS["calcite_589nm"]
    assert co.birefringence(d['no'], d['ne']) < 0
    assert not co.is_positive_uniaxial(d['no'], d['ne'])
    q = co.COMMON_UNIAXIAL_CRYSTALS["quartz_589nm"]
    assert co.is_positive_uniaxial(q['no'], q['ne'])
