"""
opticspy.geometrical.ray_tracing
=================================
Paraxial (Gaussian) matrix optics and exact real ray tracing through
rotationally-symmetric optical systems made of spherical / conic surfaces.

Sign convention (standard optical-design convention, e.g. Born & Wolf,
Kingslake, Smith "Modern Optical Engineering"):
    * Light nominally travels in the +z direction.
    * A surface radius R is positive if its center of curvature lies to
      the +z side of the vertex (i.e. downstream), negative otherwise.
    * Ray height y and angle u (radians, measured from the axis, positive
      counter-clockwise from +z towards +y) follow the same convention.
    * After a reflection the ray direction reverses; propagation
      distances and the refractive index carried by the paraxial
      "reduced angle" formalism are handled by flipping the sign of n
      for the reflected space (a standard book-keeping trick, see e.g.
      Born & Wolf, section on mirror systems), so a single matrix
      formalism handles lenses and mirrors together.

Paraxial part: Born & Wolf, "Principles of Optics", Ch. 4 (Geometrical
theory of optical imaging) and Saleh & Teich, "Fundamentals of Photonics",
Ch. 1 (Ray Optics), section on matrix optics.

Real ray trace part implements exact 3-D vector refraction/reflection at
spherical or general-conic surfaces (sag equation
z = (x^2+y^2)/(R (1+sqrt(1-(1+k)(x^2+y^2)/R^2)))), which reduces to a
sphere for conic constant k = 0.
"""

from dataclasses import dataclass, field
from typing import List, Optional, Sequence, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Paraxial (2x2 ABCD) matrix optics, using the "reduced angle" convention
# ray-vector = (y, w) with w = n*u  (u = true ray angle in radians)
# ---------------------------------------------------------------------------

def free_space_matrix(d, n=1.0):
    """Propagation by distance d through a medium of index n (reduced-angle form)."""
    return np.array([[1.0, d / n], [0.0, 1.0]])


def refraction_matrix(R, n1, n2):
    """Refraction at a single spherical surface of radius R, index n1 -> n2."""
    if np.isinf(R):
        power = 0.0
    else:
        power = (n2 - n1) / R
    return np.array([[1.0, 0.0], [-power, 1.0]])


def flat_interface_matrix(n1, n2):
    """Refraction at a flat (R = inf) interface. Included for readability."""
    return refraction_matrix(np.inf, n1, n2)


def thin_lens_matrix(f):
    """Thin lens of focal length f, immersed in air on both sides (n=1)."""
    return np.array([[1.0, 0.0], [-1.0 / f, 1.0]])


def mirror_matrix(R):
    """
    Thin spherical mirror of radius R (R < 0 for a concave mirror facing
    the incoming beam, following the sign convention above). Equivalent
    focal length f = R/2 (paraxial mirror formula), matrix built in
    reduced-angle form directly.
    """
    return np.array([[1.0, 0.0], [-2.0 / R, 1.0]])


def thick_lens_matrix(R1, R2, t, n_lens, n_out=1.0, n_in=1.0):
    """Thick lens: refraction at R1 (n_in->n_lens), free space t, refraction at R2 (n_lens->n_out)."""
    M1 = refraction_matrix(R1, n_in, n_lens)
    M2 = free_space_matrix(t, n_lens)
    M3 = refraction_matrix(R2, n_lens, n_out)
    return combine_matrices(M3, M2, M1)


def combine_matrices(*matrices):
    """
    Combine ABCD matrices in the optical (physical) order they are
    encountered by the ray, i.e. combine_matrices(M1, M2, M3) applies
    M1 first. Internally this multiplies M3 @ M2 @ M1, matching the
    matrix-optics convention v_out = M_total @ v_in.
    """
    M = np.eye(2)
    for m in matrices:
        M = m @ M
    return M


def apply_matrix(M, y, u, n=1.0):
    """Apply an ABCD matrix (reduced-angle form) to a ray, returning (y', u')."""
    v = M @ np.array([y, n * u])
    y2, w2 = v
    return y2, w2  # w2 is still a reduced angle; caller divides by n_out if needed


# ---------------------------------------------------------------------------
# Optical system: a stack of surfaces, both paraxial and real ray tracing
# ---------------------------------------------------------------------------

@dataclass
class Surface:
    """A single rotationally-symmetric optical surface.

    Parameters
    ----------
    radius : float
        Radius of curvature (signed). ``np.inf`` for a flat surface/plane.
    thickness : float
        Axial distance from this surface's vertex to the next surface's
        vertex (the medium *after* this surface).
    index : float
        Refractive index of the medium *after* this surface (between this
        surface and the next one).
    semi_diameter : float, optional
        Clear-aperture radius, used to flag vignetting in real ray tracing.
    conic : float
        Conic constant k (0 = sphere, -1 = paraboloid, <-1 = hyperboloid,
        -1<k<0 = prolate ellipsoid, k>0 = oblate ellipsoid).
    is_mirror : bool
        If True, the surface reflects instead of refracting (real ray
        tracing reverses the ray's z-direction; paraxial matrix uses the
        mirror_matrix formalism with n flipped in sign for the space that
        follows, which the OpticalSystem bookkeeping handles automatically).
    name : str
        Optional label.
    """
    radius: float
    thickness: float
    index: float
    semi_diameter: Optional[float] = None
    conic: float = 0.0
    is_mirror: bool = False
    name: str = ""


@dataclass
class RayPath:
    points: List[np.ndarray]
    directions: List[np.ndarray]
    vignetted: bool = False
    tir: bool = False
    stopped_at_surface: Optional[int] = None


def _normalize(v):
    v = np.asarray(v, dtype=float)
    return v / np.linalg.norm(v)


def _conic_intersect(O, d, R, k):
    """
    Intersect a ray (origin O, unit direction d, both in the surface's
    local frame with vertex at the origin and axis along +z) with the
    conicoid x^2 + y^2 + (1+k) z^2 - 2 R z = 0. Returns the ray parameter
    t for the physically relevant (near-vertex) intersection, or None.
    """
    Ox, Oy, Oz = O
    dx, dy, dz = d
    kp1 = 1.0 + k

    A = dx * dx + dy * dy + kp1 * dz * dz
    B = 2.0 * (Ox * dx + Oy * dy) + 2.0 * kp1 * Oz * dz - 2.0 * R * dz
    C = Ox * Ox + Oy * Oy + kp1 * Oz * Oz - 2.0 * R * Oz

    if np.isinf(R):
        # Flat surface: plane z = 0
        if abs(dz) < 1e-14:
            return None
        return -Oz / dz

    if abs(A) < 1e-14:
        # Degenerate to linear equation (e.g. paraboloid hit along axis)
        if abs(B) < 1e-14:
            return None
        return -C / B

    disc = B * B - 4 * A * C
    if disc < 0:
        return None
    sq = np.sqrt(disc)
    t1 = (-B + sq) / (2 * A)
    t2 = (-B - sq) / (2 * A)
    # choose the root giving the intersection point nearest the vertex
    # plane (z closest to 0) among the physically forward-going ones
    candidates = [t for t in (t1, t2) if t is not None]
    if not candidates:
        return None
    best = min(candidates, key=lambda t: abs(Oz + t * dz))
    return best


def _conic_normal(P, R, k):
    x, y, z = P
    if np.isinf(R):
        n = np.array([0.0, 0.0, -1.0])
    else:
        n = np.array([2 * x, 2 * y, 2 * (1 + k) * z - 2 * R])
    return _normalize(n)


def _refract_vector(d, N, n1, n2):
    """Vector form of Snell's law. d, N unit vectors; N oriented against d.
    Returns None on total internal reflection."""
    d = _normalize(d)
    cos_i = -np.dot(N, d)
    if cos_i < 0:
        N = -N
        cos_i = -np.dot(N, d)
    eta = n1 / n2
    sin2_t = eta ** 2 * (1.0 - cos_i ** 2)
    if sin2_t > 1.0:
        return None  # TIR
    cos_t = np.sqrt(1.0 - sin2_t)
    return _normalize(eta * d + (eta * cos_i - cos_t) * N)


def _reflect_vector(d, N):
    d = _normalize(d)
    if np.dot(N, d) > 0:
        N = -N
    return _normalize(d - 2 * np.dot(d, N) * N)


class OpticalSystem:
    """A sequential, rotationally-symmetric optical system."""

    def __init__(self, surfaces: Sequence[Surface], object_index: float = 1.0):
        self.surfaces: List[Surface] = list(surfaces)
        self.object_index = object_index

    # ---------------- paraxial ----------------

    def _surface_indices(self):
        """Return (n_before, n_after) for each surface, tracking mirror sign flips."""
        n_before = self.object_index
        pairs = []
        for s in self.surfaces:
            n_after = s.index
            pairs.append((n_before, n_after))
            n_before = n_after
        return pairs

    def paraxial_system_matrix(self) -> np.ndarray:
        """
        System matrix (reduced-angle convention) taking a ray from just
        before the first surface to just after the last surface.
        """
        M = np.eye(2)
        n_prev = self.object_index
        for i, s in enumerate(self.surfaces):
            if s.is_mirror:
                Mr = mirror_matrix(s.radius)
            else:
                Mr = refraction_matrix(s.radius, n_prev, s.index)
            M = Mr @ M
            if i < len(self.surfaces) - 1:
                Mt = free_space_matrix(s.thickness, s.index)
                M = Mt @ M
            n_prev = s.index
        return M

    def trace_paraxial_ray(self, y0: float, u0: float) -> List[Tuple[float, float]]:
        """Trace a paraxial ray (y0,u0) surface by surface.
        Returns a list of (y, u) *after* each surface (true angle u, not reduced)."""
        results = []
        y, w = y0, self.object_index * u0
        n_prev = self.object_index
        for i, s in enumerate(self.surfaces):
            if s.is_mirror:
                Mr = mirror_matrix(s.radius)
            else:
                Mr = refraction_matrix(s.radius, n_prev, s.index)
            y, w = Mr @ np.array([y, w])
            n_after = s.index
            results.append((y, w / n_after))
            if i < len(self.surfaces) - 1:
                Mt = free_space_matrix(s.thickness, s.index)
                y, w = Mt @ np.array([y, w])
            n_prev = n_after
        return results

    def effective_focal_length(self) -> float:
        M = self.paraxial_system_matrix()
        C = M[1, 0]
        if abs(C) < 1e-14:
            return np.inf
        n_first = self.object_index
        n_last = self.surfaces[-1].index
        # power phi = -C (reduced-angle form); EFL referenced to image space index
        return -n_last / C if abs(C) > 0 else np.inf

    def cardinal_points(self) -> dict:
        """
        Compute the cardinal points of the system: effective focal length,
        front/back focal distances (from first/last vertex) and the
        principal-plane offsets (from first/last vertex).
        """
        M = self.paraxial_system_matrix()
        A, B = M[0, 0], M[0, 1]
        C, D = M[1, 0], M[1, 1]
        n1 = self.object_index
        n2 = self.surfaces[-1].index
        if abs(C) < 1e-14:
            return dict(efl=np.inf, bfd=np.inf, ffd=np.inf, epl=np.inf, bpl=np.inf)
        efl = -n2 / C          # image-space effective focal length (paraxial power = -C)
        efl_obj = -n1 / C      # object-space effective focal length (equal magnitude if n1=n2)
        bfd = -A / C           # back focal distance from LAST vertex to back focal point
        ffd = -D / C           # front focal distance from FIRST vertex to front focal point (sign per convention below)
        bpl = bfd - efl        # back principal plane offset from last vertex
        fpl = ffd + efl_obj    # front principal plane offset from first vertex
        return dict(efl=efl, efl_object_space=efl_obj, bfd=bfd, ffd=-ffd,
                    back_principal_plane=bpl, front_principal_plane=fpl)

    def image_conjugate(self, object_distance: float):
        """
        Given an object at distance `object_distance` in front of the first
        surface (positive = to the left, i.e. light travels a distance
        `object_distance` before reaching surface 1), return
        (image_distance_from_last_vertex, transverse_magnification).
        """
        M_sys = self.paraxial_system_matrix()
        M_obj = free_space_matrix(object_distance, self.object_index)
        M_total_partial = M_sys @ M_obj  # ray just after last surface, from object plane
        B = M_total_partial[0, 1]
        D = M_total_partial[1, 1]
        n_img = self.surfaces[-1].index
        if abs(D) < 1e-14:
            return np.inf, np.inf
        # Imaging condition: the B element of the object-plane-to-image-plane
        # matrix must vanish (image height independent of input ray angle).
        image_distance = -B * n_img / D
        M_img = free_space_matrix(image_distance, n_img)
        M_total = M_img @ M_total_partial
        magnification = M_total[0, 0]
        return image_distance, magnification

    # ---------------- real ray tracing ----------------

    def trace_ray(self, point, direction) -> RayPath:
        """
        Exact real ray trace of a single ray through the system.
        `point` : (x, y, z) starting coordinates, z measured from the
                  vertex of the first surface.
        `direction` : (dx, dy, dz) initial propagation direction (need not
                      be normalized).
        """
        P = np.array(point, dtype=float)
        d = _normalize(direction)
        points = [P.copy()]
        directions = [d.copy()]

        z_vertex = 0.0
        n_before = self.object_index
        for i, s in enumerate(self.surfaces):
            local_O = P - np.array([0.0, 0.0, z_vertex])
            t = _conic_intersect(local_O, d, s.radius, s.conic)
            if t is None or t < -1e-9:
                return RayPath(points, directions, vignetted=True, stopped_at_surface=i)
            hit_local = local_O + t * d
            hit_global = hit_local + np.array([0.0, 0.0, z_vertex])

            if s.semi_diameter is not None:
                r_hit = np.hypot(hit_local[0], hit_local[1])
                if r_hit > s.semi_diameter:
                    points.append(hit_global)
                    directions.append(d.copy())
                    return RayPath(points, directions, vignetted=True, stopped_at_surface=i)

            N = _conic_normal(hit_local, s.radius, s.conic)

            if s.is_mirror:
                d_new = _reflect_vector(d, N)
                n_after = -n_before  # sign-flip convention for reflected space bookkeeping
            else:
                n_after = s.index
                d_new = _refract_vector(d, N, n_before, n_after)
                if d_new is None:
                    points.append(hit_global)
                    directions.append(d.copy())
                    return RayPath(points, directions, tir=True, stopped_at_surface=i)

            points.append(hit_global)
            directions.append(d_new.copy())
            P = hit_global
            d = d_new
            z_vertex += s.thickness
            n_before = s.index if not s.is_mirror else -n_before

        return RayPath(points, directions)

    def spot_diagram(self, object_point, aperture_semi_diameter, n_rays=64,
                      pattern='hexapolar', stop_distance=0.0, image_distance=None):
        """
        Trace a fan of real rays from `object_point` (x, y, z) sampled over
        a circular aperture (an approximation of the entrance pupil placed
        at `stop_distance` from the first vertex), and return their (x, y)
        coordinates at the paraxial image plane (or a user-specified
        `image_distance` measured from the last vertex).

        Returns
        -------
        xy : (N, 2) ndarray of ray intersection coordinates at the image
             plane (rows are NaN for vignetted/TIR rays).
        """
        object_point = np.asarray(object_point, dtype=float)
        pts = _pupil_sampling(aperture_semi_diameter, n_rays, pattern)

        if image_distance is None:
            obj_dist = stop_distance - object_point[2]
            image_distance, _ = self.image_conjugate(obj_dist)

        results = []
        z_stop = stop_distance
        for (px, py) in pts:
            target = np.array([px, py, z_stop])
            direction = target - object_point
            path = self.trace_ray(object_point, direction)
            if path.vignetted or path.tir:
                results.append((np.nan, np.nan))
                continue
            # propagate from last surface hit to the image plane
            z_last_vertex = sum(s.thickness for s in self.surfaces[:-1])
            P_last = path.points[-1]
            d_last = path.directions[-1]
            z_image = z_last_vertex + self.surfaces[-1].thickness + image_distance \
                if False else z_last_vertex + image_distance
            # simpler: propagate from last hit point along d_last to plane z = z_total_to_last_vertex + image_distance
            z_target = sum(s.thickness for s in self.surfaces[:-1]) + image_distance
            if abs(d_last[2]) < 1e-14:
                results.append((np.nan, np.nan))
                continue
            t = (z_target - P_last[2]) / d_last[2]
            P_img = P_last + t * d_last
            results.append((P_img[0], P_img[1]))
        return np.array(results)


def _pupil_sampling(radius, n, pattern='hexapolar'):
    pts = []
    if pattern == 'grid':
        m = int(np.ceil(np.sqrt(n)))
        xs = np.linspace(-radius, radius, m)
        for x in xs:
            for y in xs:
                if x * x + y * y <= radius * radius:
                    pts.append((x, y))
    else:  # hexapolar rings
        pts.append((0.0, 0.0))
        n_rings = max(1, int(round(np.sqrt(n / 3.0))))
        for ring in range(1, n_rings + 1):
            r = radius * ring / n_rings
            n_theta = 6 * ring
            for k in range(n_theta):
                th = 2 * np.pi * k / n_theta
                pts.append((r * np.cos(th), r * np.sin(th)))
    return pts


def ray_fan(system: OpticalSystem, object_height: float, aperture_semi_diameter: float,
            n_rays: int = 21, stop_distance: float = 0.0, image_distance=None,
            object_distance: float = 1e6):
    """
    Classic meridional ray-fan aberration plot data: transverse ray error
    at the paraxial image plane as a function of the normalized pupil
    coordinate, for an object point at height `object_height`.

    Returns (pupil_coords, transverse_error_y) both length n_rays.
    """
    object_point = np.array([0.0, object_height, -object_distance])
    pupil_frac = np.linspace(-1, 1, n_rays)
    pts_y = pupil_frac * aperture_semi_diameter

    if image_distance is None:
        image_distance, mag = system.image_conjugate(object_distance + stop_distance)
    else:
        _, mag = system.image_conjugate(object_distance + stop_distance)

    # paraxial (chief-ray) image height for reference
    paraxial_ray = system.trace_paraxial_ray(object_height, 0.0)
    ideal_y = object_height * mag

    errors = np.full(n_rays, np.nan)
    z_last_vertex = sum(s.thickness for s in system.surfaces[:-1])
    z_target = z_last_vertex + image_distance
    for i, py in enumerate(pts_y):
        target = np.array([0.0, py, stop_distance])
        direction = target - object_point
        path = system.trace_ray(object_point, direction)
        if path.vignetted or path.tir:
            continue
        P_last = path.points[-1]
        d_last = path.directions[-1]
        if abs(d_last[2]) < 1e-14:
            continue
        t = (z_target - P_last[2]) / d_last[2]
        P_img = P_last + t * d_last
        errors[i] = P_img[1] - ideal_y
    return pupil_frac, errors
