"""
aoptics.geometrical.aberrations
=================================
Third-order (Seidel) monochromatic aberration theory for rotationally
symmetric optical systems, computed from two paraxial rays (the marginal
ray and the chief/principal ray) traced surface-by-surface.

Reference formulas (standard "Seidel sum" formulation, consistent with
Born & Wolf Ch. 5, Welford's "Aberrations of Optical Systems", and
reproduced in numerous optical-design texts):

    Per-surface refraction invariants
        A      = n (y c + u)          Abbe (Snell) invariant, marginal ray
        A_bar  = n (ybar c + ubar)    Abbe invariant, chief ray
        H      = n * ubar * y - n * u * ybar        Lagrange invariant

    Seidel sums (summed over all surfaces)
        S_I   = -sum[ A^2 * y * Delta(u/n) ]                  spherical aberration
        S_II  = -sum[ A_bar * A * y * Delta(u/n) ]            coma
        S_III = -sum[ A_bar^2 * y * Delta(u/n) ]              astigmatism
        S_IV  = -sum[ H^2 * c * Delta(1/n) ]                  Petzval (field curvature)
        S_V   = -sum[ (A_bar^3/A) * y * Delta(u/n)
                       + (A_bar/A) * H^2 * c * Delta(1/n) ]    distortion

    with Delta(u/n) = u'/n' - u/n and Delta(1/n) = 1/n' - 1/n evaluated
    across each surface (unprimed = before, primed = after).

    Peak wavefront-aberration coefficients (standard conversion):
        W040 = S_I / 8      W131 = S_II / 2      W222 = S_III / 2
        W220 = (S_III + S_IV) / 4      W311 = S_V / 2
"""

from dataclasses import dataclass
from typing import List, Optional

import numpy as np

from .ray_tracing import OpticalSystem


def _full_paraxial_trace(system: OpticalSystem, y0: float, u0: float, object_distance: float = 0.0):
    """Paraxial trace keeping BOTH before- and after-surface angles, needed
    for Seidel sums (OpticalSystem.trace_paraxial_ray only returns the
    after-surface angle)."""
    data = []
    y = y0 + u0 * object_distance
    u = u0
    n_before = system.object_index
    for s in system.surfaces:
        c = 0.0 if np.isinf(s.radius) else 1.0 / s.radius
        n_after = -n_before if s.is_mirror else s.index
        phi = (n_after - n_before) * c
        u_after = (n_before * u - y * phi) / n_after
        data.append(dict(y=y, u_before=u, u_after=u_after,
                          n_before=n_before, n_after=n_after, c=c, R=s.radius))
        y = y + u_after * s.thickness
        u = u_after
        n_before = n_after
    return data


def solve_u0_for_height(system: OpticalSystem, y0: float, object_distance: float,
                         target_surface_index: int, target_height: float) -> float:
    """
    Solve (exactly, by linearity of the paraxial trace) for the initial
    ray angle u0 -- starting at height y0, `object_distance` before the
    first surface -- that makes the ray reach `target_height` at the
    surface indexed by `target_surface_index`.
    """
    d_base = _full_paraxial_trace(system, y0, 0.0, object_distance)
    d_unit = _full_paraxial_trace(system, 0.0, 1.0, object_distance)
    h_base = d_base[target_surface_index]['y']
    h_unit = d_unit[target_surface_index]['y']
    if abs(h_unit) < 1e-14:
        raise ValueError("Target surface height is insensitive to ray angle here "
                          "(afocal or degenerate configuration); specify u0 directly.")
    return (target_height - h_base) / h_unit


def solve_marginal_ray(system: OpticalSystem, object_distance: float,
                        stop_surface_index: int, aperture_semi_diameter: float,
                        object_height: float = 0.0):
    """Marginal ray: starts at the (usually on-axis) object point and just
    grazes the edge of the aperture stop."""
    u0 = solve_u0_for_height(system, object_height, object_distance,
                              stop_surface_index, aperture_semi_diameter)
    return _full_paraxial_trace(system, object_height, u0, object_distance)


def solve_chief_ray(system: OpticalSystem, object_distance: float,
                     stop_surface_index: int, object_height: float):
    """Chief (principal) ray: starts at the off-axis object point and
    passes through the centre of the aperture stop."""
    u0 = solve_u0_for_height(system, object_height, object_distance,
                              stop_surface_index, 0.0)
    return _full_paraxial_trace(system, object_height, u0, object_distance)


@dataclass
class SeidelResult:
    S_I: float
    S_II: float
    S_III: float
    S_IV: float
    S_V: float
    per_surface: List[dict]

    def wavefront_coefficients(self):
        """Peak wavefront-aberration coefficients W040, W131, W222, W220, W311
        (standard 1/8, 1/2 conversion, see module docstring)."""
        return dict(
            W040=self.S_I / 8.0,
            W131=self.S_II / 2.0,
            W222=self.S_III / 2.0,
            W220=(self.S_III + self.S_IV) / 4.0,
            W311=self.S_V / 2.0,
        )

    def summary(self):
        return {
            "spherical (S_I)": self.S_I,
            "coma (S_II)": self.S_II,
            "astigmatism (S_III)": self.S_III,
            "Petzval / field curvature (S_IV)": self.S_IV,
            "distortion (S_V)": self.S_V,
        }


def seidel_sums(system: OpticalSystem,
                 marginal_trace: Optional[List[dict]] = None,
                 chief_trace: Optional[List[dict]] = None,
                 object_distance: float = 0.0,
                 stop_surface_index: Optional[int] = None,
                 aperture_semi_diameter: Optional[float] = None,
                 object_height: Optional[float] = None) -> SeidelResult:
    """
    Compute the five third-order Seidel sums for `system`.

    Two ways to specify the marginal/chief rays:
      1. Pass pre-traced `marginal_trace` / `chief_trace` (as returned by
         `_full_paraxial_trace`, or `solve_marginal_ray` / `solve_chief_ray`).
      2. Pass `object_distance`, `stop_surface_index`, `aperture_semi_diameter`
         and `object_height`, and the rays are solved for automatically.
    """
    if marginal_trace is None or chief_trace is None:
        if stop_surface_index is None or aperture_semi_diameter is None or object_height is None:
            raise ValueError("Provide either explicit ray traces, or "
                              "(object_distance, stop_surface_index, "
                              "aperture_semi_diameter, object_height).")
        marginal_trace = solve_marginal_ray(system, object_distance, stop_surface_index,
                                             aperture_semi_diameter, object_height=0.0)
        chief_trace = solve_chief_ray(system, object_distance, stop_surface_index, object_height)

    S_I = S_II = S_III = S_IV = S_V = 0.0
    per_surface = []
    for j in range(len(system.surfaces)):
        m, ch = marginal_trace[j], chief_trace[j]
        y, ybar = m['y'], ch['y']
        u, uprime = m['u_before'], m['u_after']
        ubar = ch['u_before']
        n, nprime = m['n_before'], m['n_after']
        c = m['c']

        A = n * (y * c + u)
        Abar = n * (ybar * c + ubar)
        H = n * ubar * y - n * u * ybar

        d_u_n = uprime / nprime - u / n
        d_1_n = 1.0 / nprime - 1.0 / n

        sI = -(A ** 2) * y * d_u_n
        sII = -(Abar * A) * y * d_u_n
        sIII = -(Abar ** 2) * y * d_u_n
        sIV = -(H ** 2) * c * d_1_n
        sV = -((Abar ** 3 / A) * y * d_u_n + (Abar / A) * (H ** 2) * c * d_1_n) if abs(A) > 1e-13 else 0.0

        S_I += sI; S_II += sII; S_III += sIII; S_IV += sIV; S_V += sV
        per_surface.append(dict(surface=j, S_I=sI, S_II=sII, S_III=sIII,
                                 S_IV=sIV, S_V=sV, A=A, A_bar=Abar, H=H))

    return SeidelResult(S_I, S_II, S_III, S_IV, S_V, per_surface)


def petzval_curvature(system: OpticalSystem) -> float:
    """
    System Petzval curvature P = sum_j c_j * (1/n'_j - 1/n_j), independent
    of ray heights/angles -- the radius of the Petzval image surface is
    R_petzval = -1/P.
    """
    n_before = system.object_index
    P = 0.0
    for s in system.surfaces:
        c = 0.0 if np.isinf(s.radius) else 1.0 / s.radius
        n_after = -n_before if s.is_mirror else s.index
        P += c * (1.0 / n_after - 1.0 / n_before)
        n_before = n_after
    return P


def petzval_radius(system: OpticalSystem) -> float:
    P = petzval_curvature(system)
    return np.inf if abs(P) < 1e-14 else -1.0 / P

