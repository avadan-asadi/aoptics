"""
Phase 1 demo -- opticspy.electromagnetic
=========================================
Fresnel reflectance vs angle for an air-glass interface, Malus's law
through a rotating analyzer, a quarter-wave plate turning linear light
circular (shown on the Poincare sphere), and quartz quarter-/half-wave
plate design from its birefringence.

Run:  python examples/phase1_polarization_demo.py
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 (registers 3d projection)

from opticspy.electromagnetic import fresnel, polarization as pol, crystal_optics as co

n1, n2 = 1.0, 1.5168  # air -> BK7 glass

# ---- 1) Fresnel reflectance vs angle of incidence ----
angles = np.linspace(0, 89.9, 400)
rt_data = fresnel.reflectance_transmittance(np.radians(angles), n1, n2)
thetaB = np.degrees(fresnel.brewster_angle(n1, n2))
print(f"Brewster angle (air->BK7): {thetaB:.3f} deg  (R_p should vanish there)")
print(f"Normal-incidence reflectance: {fresnel.normal_incidence_reflectance(n1, n2):.4f}")

# ---- 2) Malus's law: intensity through a rotating analyzer after a fixed polarizer ----
theta_analyzer = np.linspace(0, 2 * np.pi, 200)
v0 = pol.linear_horizontal()
I = [pol.jones_intensity(pol.apply_jones(pol.linear_polarizer(t), v0)) for t in theta_analyzer]

# ---- 3) Quarter-wave plate: linear -> circular, tracked on the Poincare sphere ----
v_lin = pol.linear_at_angle(np.pi / 4)
n_steps = 60
retardances = np.linspace(0, np.pi / 2, n_steps)
poincare_pts = []
for d in retardances:
    J = pol.wave_plate(d, fast_axis_angle=0.0)
    out = pol.apply_jones(J, v_lin)
    s = pol.jones_to_stokes(out)
    poincare_pts.append(pol.poincare_coordinates(s))
poincare_pts = np.array(poincare_pts)
print(f"\n45-deg linear light through a variable retarder: "
      f"S3/S0 goes from {poincare_pts[0,2]:.3f} (linear) to {poincare_pts[-1,2]:.3f} (circular) "
      "as retardance sweeps 0 -> lambda/4.")

# ---- 4) Quartz waveplate design ----
quartz = co.COMMON_UNIAXIAL_CRYSTALS["quartz_589nm"]
wl = 589e-9
d_qwp = co.waveplate_thickness(quartz['no'], quartz['ne'], wl, order=0, fraction=0.25)
d_hwp = co.waveplate_thickness(quartz['no'], quartz['ne'], wl, order=0, fraction=0.5)
print(f"\nQuartz (no={quartz['no']}, ne={quartz['ne']}) at {wl*1e9:.0f} nm:")
print(f"  zero-order quarter-wave plate thickness = {d_qwp*1e6:.3f} micron")
print(f"  zero-order half-wave plate thickness    = {d_hwp*1e6:.3f} micron")

# ---- Figure ----
fig = plt.figure(figsize=(13, 4.2))

ax1 = fig.add_subplot(1, 3, 1)
ax1.plot(angles, rt_data['R_s'], label='R_s')
ax1.plot(angles, rt_data['R_p'], label='R_p')
ax1.axvline(thetaB, color='gray', ls='--', lw=1, label=f'Brewster {thetaB:.1f} deg')
ax1.set_xlabel("angle of incidence (deg)"); ax1.set_ylabel("reflectance")
ax1.set_title("Fresnel reflectance: air -> BK7")
ax1.legend(fontsize=8)

ax2 = fig.add_subplot(1, 3, 2)
ax2.plot(np.degrees(theta_analyzer), I)
ax2.set_xlabel("analyzer angle (deg)"); ax2.set_ylabel("transmitted intensity")
ax2.set_title("Malus's law")

ax3 = fig.add_subplot(1, 3, 3, projection='3d')
u, v = np.mgrid[0:2*np.pi:30j, 0:np.pi:15j]
xs, ys, zs = np.cos(u)*np.sin(v), np.sin(u)*np.sin(v), np.cos(v)
ax3.plot_wireframe(xs, ys, zs, color='lightgray', linewidth=0.3)
ax3.plot(poincare_pts[:, 0], poincare_pts[:, 1], poincare_pts[:, 2], color='crimson', lw=2)
ax3.scatter(*poincare_pts[0], color='blue', label='linear (start)')
ax3.scatter(*poincare_pts[-1], color='green', label='circular (end)')
ax3.set_title("Poincare sphere: QWP linear->circular")
ax3.legend(fontsize=7)

fig.tight_layout()
fig.savefig("phase1_polarization_demo.png", dpi=140)
print("\nSaved figure: phase1_polarization_demo.png")
