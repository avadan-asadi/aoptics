"""
Phase 1 demo -- opticspy.geometrical
=====================================
A plano-convex singlet lens: paraxial cardinal points, a real-ray spot
diagram at the paraxial image plane, a meridional ray-aberration fan,
and the third-order Seidel aberration sums (spherical, coma,
astigmatism, Petzval/field-curvature, distortion).

Run:  python examples/phase1_geometrical_demo.py
"""
import numpy as np
import matplotlib.pyplot as plt

import opticspy as op
from opticspy.geometrical import ray_tracing as rt
from opticspy.geometrical import aberrations as ab

# ---- Define a simple plano-convex singlet (BK7-like glass, n = 1.5168) ----
n_glass = 1.5168
R1 = 51.5        # convex first surface (mm)
R2 = np.inf       # flat second surface
thickness = 5.3   # center thickness (mm)
semi_dia = 12.5   # clear aperture radius (mm)

surfaces = [
    rt.Surface(radius=R1, thickness=thickness, index=n_glass, semi_diameter=semi_dia, name="S1"),
    rt.Surface(radius=R2, thickness=0.0, index=1.0, semi_diameter=semi_dia, name="S2"),
]
system = rt.OpticalSystem(surfaces, object_index=1.0)

# ---- Paraxial cardinal points ----
cp = system.cardinal_points()
print("=== Paraxial cardinal points ===")
for k, v in cp.items():
    print(f"  {k:28s} = {v: .4f} mm")

# ---- Image conjugate for an object 500 mm in front of the lens ----
obj_dist = 500.0
img_dist, mag = system.image_conjugate(obj_dist)
print(f"\nObject at {obj_dist} mm -> image at {img_dist:.4f} mm behind S2, magnification = {mag:.5f}")

# ---- Third-order Seidel aberration sums (on-axis field, full aperture stop at S1) ----
seidel = ab.seidel_sums(system, object_distance=obj_dist, stop_surface_index=0,
                         aperture_semi_diameter=semi_dia, object_height=0.5)
print("\n=== Seidel sums (on-axis marginal ray, small field for coma/astig/distortion) ===")
for name, val in seidel.summary().items():
    print(f"  {name:35s} = {val: .6e}")
print("\nPeak wavefront-aberration coefficients (waves x wavelength):")
for name, val in seidel.wavefront_coefficients().items():
    print(f"  {name} = {val: .6e}")

print(f"\nPetzval radius = {ab.petzval_radius(system):.3f} mm")

# ---- Real-ray spot diagram at the paraxial image plane ----
xy = system.spot_diagram(object_point=(0, 0, -obj_dist), aperture_semi_diameter=semi_dia,
                          n_rays=200, stop_distance=0.0, image_distance=img_dist)
rms_spot = np.sqrt(np.nanmean(xy[:, 0] ** 2 + xy[:, 1] ** 2))
print(f"\nReal-ray RMS spot radius at paraxial image plane: {rms_spot * 1000:.2f} micron "
      f"({np.sum(~np.isnan(xy[:, 0]))}/{len(xy)} rays traced successfully)")

# ---- Meridional ray-aberration fan ----
pupil, err = rt.ray_fan(system, object_height=0.0, aperture_semi_diameter=semi_dia,
                         n_rays=41, object_distance=obj_dist)

fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
axes[0].scatter(xy[:, 0] * 1000, xy[:, 1] * 1000, s=4, alpha=0.6)
axes[0].set_xlabel("x (micron)"); axes[0].set_ylabel("y (micron)")
axes[0].set_title("Real-ray spot diagram (on-axis)")
axes[0].set_aspect('equal')

axes[1].plot(pupil, err * 1000)
axes[1].axhline(0, color='k', lw=0.5)
axes[1].set_xlabel("normalized pupil coordinate")
axes[1].set_ylabel("transverse ray error (micron)")
axes[1].set_title("Meridional ray-aberration fan (spherical aberration)")

fig.tight_layout()
fig.savefig("phase1_geometrical_demo.png", dpi=140)
print("\nSaved figure: phase1_geometrical_demo.png")
