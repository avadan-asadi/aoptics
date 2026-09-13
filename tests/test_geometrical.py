import numpy as np
import pytest

from opticspy.geometrical import ray_tracing as rt
from opticspy.geometrical import aberrations as ab


def test_free_space_matrix_translates_height():
    M = rt.free_space_matrix(d=10.0, n=1.0)
    y, w = M @ np.array([0.0, 1.0])  # ray starting on axis, angle 1 rad (paraxial, illustrative)
    assert y == pytest.approx(10.0)
    assert w == pytest.approx(1.0)


def test_thin_lens_matrix_focuses_parallel_ray():
    f = 100.0
    M = rt.combine_matrices(rt.thin_lens_matrix(f), rt.free_space_matrix(f))
    y_in, u_in = 5.0, 0.0  # ray parallel to axis at height 5
    y_out, w_out = M @ np.array([y_in, u_in])
    assert y_out == pytest.approx(0.0, abs=1e-9)  # must cross axis at the focal point


def test_refraction_matrix_zero_power_for_flat_surface():
    M = rt.refraction_matrix(np.inf, 1.0, 1.5)
    assert np.allclose(M, np.eye(2))


def test_thick_lens_efl_matches_lensmaker_equation():
    # Biconvex lens, thin-lens limit (t -> 0): 1/f = (n-1)(1/R1 - 1/R2)
    n_lens = 1.5
    R1, R2 = 50.0, -50.0
    t = 1e-6
    M = rt.thick_lens_matrix(R1, R2, t, n_lens, n_out=1.0, n_in=1.0)
    f_expected = 1.0 / ((n_lens - 1) * (1 / R1 - 1 / R2))
    f_from_matrix = -1.0 / M[1, 0]
    assert f_from_matrix == pytest.approx(f_expected, rel=1e-3)


def test_optical_system_cardinal_points_thin_lens_in_air():
    f = 75.0
    surf = rt.Surface(radius=np.inf, thickness=0.0, index=1.0)  # placeholder unused directly
    # Build a system directly equivalent to a thin lens using two surfaces of a
    # thin biconvex lens (n=1.5) with negligible thickness.
    R1, R2 = 2 * f * (1.5 - 1), -2 * f * (1.5 - 1)
    surfaces = [
        rt.Surface(radius=R1, thickness=1e-6, index=1.5),
        rt.Surface(radius=R2, thickness=0.0, index=1.0),
    ]
    system = rt.OpticalSystem(surfaces, object_index=1.0)
    cp = system.cardinal_points()
    assert cp['efl'] == pytest.approx(f, rel=1e-3)


def test_real_ray_trace_flat_normal_incidence_straight_through():
    surfaces = [
        rt.Surface(radius=np.inf, thickness=5.0, index=1.5, semi_diameter=10.0),
        rt.Surface(radius=np.inf, thickness=0.0, index=1.0, semi_diameter=10.0),
    ]
    system = rt.OpticalSystem(surfaces, object_index=1.0)
    path = system.trace_ray(point=(0, 0, -1), direction=(0, 0, 1))
    assert not path.vignetted and not path.tir
    # normal incidence: ray should continue straight (no bending)
    for d in path.directions:
        assert np.allclose(d, [0, 0, 1], atol=1e-9)


def test_real_ray_trace_snell_bending_direction():
    # Going from less dense (n=1) to more dense (n=1.5) medium at a flat
    # interface, an oblique ray should bend TOWARD the normal.
    surfaces = [rt.Surface(radius=np.inf, thickness=5.0, index=1.5, semi_diameter=10.0)]
    system = rt.OpticalSystem(surfaces, object_index=1.0)
    theta_i = np.radians(30.0)
    direction = (np.sin(theta_i), 0, np.cos(theta_i))
    path = system.trace_ray(point=(0, 0, -1), direction=direction)
    d_out = path.directions[-1]
    theta_t = np.arctan2(d_out[0], d_out[2])
    theta_t_expected = np.arcsin(np.sin(theta_i) / 1.5)
    assert theta_t == pytest.approx(theta_t_expected, rel=1e-6)


def test_real_ray_trace_total_internal_reflection_flag():
    # From dense (n=1.5) to rare (n=1) beyond the critical angle -> TIR
    surfaces = [rt.Surface(radius=np.inf, thickness=5.0, index=1.0, semi_diameter=10.0)]
    system = rt.OpticalSystem(surfaces, object_index=1.5)
    theta_c = np.arcsin(1.0 / 1.5)
    theta_i = theta_c + np.radians(10.0)
    direction = (np.sin(theta_i), 0, np.cos(theta_i))
    path = system.trace_ray(point=(0, 0, -1), direction=direction)
    assert path.tir


def test_real_ray_trace_vignetting_flag():
    surfaces = [rt.Surface(radius=np.inf, thickness=5.0, index=1.5, semi_diameter=1.0)]
    system = rt.OpticalSystem(surfaces, object_index=1.0)
    path = system.trace_ray(point=(0, 5, -1), direction=(0, 0, 1))  # well outside aperture
    assert path.vignetted


def test_conic_sphere_matches_spherical_intersection():
    R = 50.0
    O = np.array([0.0, 0.0, -10.0])
    d = np.array([0.0, 0.0, 1.0])
    t = rt._conic_intersect(O, d, R, k=0.0)
    hit = O + t * d
    # For an on-axis ray hitting a sphere of radius R centered at (0,0,R),
    # the near intersection with the vertex-tangent axis is at z=0.
    assert hit[2] == pytest.approx(0.0, abs=1e-9)


def test_mirror_reflection_reverses_z():
    surfaces = [rt.Surface(radius=np.inf, thickness=10.0, index=1.0, is_mirror=True, semi_diameter=10.0)]
    system = rt.OpticalSystem(surfaces, object_index=1.0)
    path = system.trace_ray(point=(0, 0, -5), direction=(0, 0, 1))
    d_out = path.directions[-1]
    assert d_out[2] == pytest.approx(-1.0, abs=1e-9)


# --------------------- Seidel aberration theory ---------------------

def test_seidel_single_surface_matches_known_marginal_invariant():
    # Single refracting surface; verify the Abbe invariant A is conserved
    # across the surface (n*A/y should be identical before/after by
    # construction of the paraxial refraction formula).
    surfaces = [rt.Surface(radius=50.0, thickness=20.0, index=1.5, semi_diameter=10.0)]
    system = rt.OpticalSystem(surfaces, object_index=1.0)
    marginal = ab._full_paraxial_trace(system, y0=5.0, u0=0.0, object_distance=0.0)
    m = marginal[0]
    A_before = m['n_before'] * (m['y'] / m['R'] + m['u_before'])
    A_after = m['n_after'] * (m['y'] / m['R'] + m['u_after'])
    assert A_before == pytest.approx(A_after, rel=1e-10)


def test_seidel_sums_zero_for_flat_plate_normal_field():
    # A flat plate (R = inf everywhere) has zero curvature -> all Seidel
    # sums that depend on curvature vanish identically for on-axis rays
    # with a flat window; spherical/coma/astigmatism/Petzval sums must be 0.
    surfaces = [
        rt.Surface(radius=np.inf, thickness=5.0, index=1.5, semi_diameter=10.0),
        rt.Surface(radius=np.inf, thickness=10.0, index=1.0, semi_diameter=10.0),
    ]
    system = rt.OpticalSystem(surfaces, object_index=1.0)
    marginal = ab._full_paraxial_trace(system, y0=2.0, u0=0.0, object_distance=0.0)
    chief = ab._full_paraxial_trace(system, y0=0.0, u0=0.05, object_distance=0.0)
    result = ab.seidel_sums(system, marginal_trace=marginal, chief_trace=chief)
    assert result.S_I == pytest.approx(0.0, abs=1e-12)
    assert result.S_II == pytest.approx(0.0, abs=1e-12)
    assert result.S_III == pytest.approx(0.0, abs=1e-12)
    assert result.S_IV == pytest.approx(0.0, abs=1e-12)


def test_petzval_curvature_single_surface_known_formula():
    R = 40.0
    n1, n2 = 1.0, 1.5
    surfaces = [rt.Surface(radius=R, thickness=10.0, index=n2)]
    system = rt.OpticalSystem(surfaces, object_index=n1)
    P = ab.petzval_curvature(system)
    expected = (1.0 / R) * (1.0 / n2 - 1.0 / n1)
    assert P == pytest.approx(expected)


def test_image_conjugate_matches_thin_lens_formula():
    f = 100.0
    surfaces = [
        rt.Surface(radius=2 * f * (1.5 - 1), thickness=1e-6, index=1.5),
        rt.Surface(radius=-2 * f * (1.5 - 1), thickness=0.0, index=1.0),
    ]
    system = rt.OpticalSystem(surfaces, object_index=1.0)
    d_o = 300.0
    v_expected = d_o * f / (d_o - f)
    m_expected = f / (f - d_o)
    v, m = system.image_conjugate(d_o)
    assert v == pytest.approx(v_expected, rel=1e-3)
    assert m == pytest.approx(m_expected, rel=1e-3)


def test_solve_marginal_ray_hits_target_height_exactly():
    surfaces = [
        rt.Surface(radius=80.0, thickness=5.0, index=1.5, semi_diameter=15.0),
        rt.Surface(radius=-80.0, thickness=50.0, index=1.0, semi_diameter=15.0),
    ]
    system = rt.OpticalSystem(surfaces, object_index=1.0)
    trace = ab.solve_marginal_ray(system, object_distance=1000.0,
                                   stop_surface_index=0, aperture_semi_diameter=10.0)
    assert trace[0]['y'] == pytest.approx(10.0, rel=1e-9)
