"""
opticspy.electromagnetic.crystal_optics
========================================
Uniaxial crystal optics: the index ellipsoid, birefringence, waveplate
design, and spatial walk-off of the extraordinary ray.

Reference: Born & Wolf, "Principles of Optics", Ch. 14 (Optics of
crystals); Saleh & Teich, "Fundamentals of Photonics", Ch. 6.

Sign convention: a "positive" uniaxial crystal has n_e > n_o (e.g.
quartz); a "negative" uniaxial crystal has n_e < n_o (e.g. calcite,
BBO). theta is measured from the optic axis (z).
"""

import numpy as np

# Approximate room-temperature indices at common reference wavelengths,
# for convenient examples/tests (not a substitute for a dispersion model).
COMMON_UNIAXIAL_CRYSTALS = {
    "calcite_589nm": dict(no=1.6584, ne=1.4864),   # negative uniaxial
    "quartz_589nm": dict(no=1.5443, ne=1.5534),    # positive uniaxial
    "KDP_546nm": dict(no=1.5115, ne=1.4698),       # negative uniaxial
    "lithium_niobate_633nm": dict(no=2.286, ne=2.200),  # negative uniaxial
    "BBO_532nm": dict(no=1.6749, ne=1.5555),       # negative uniaxial
}


def birefringence(no, ne):
    """Delta n = ne - no (sign indicates positive/negative uniaxial)."""
    return ne - no


def is_positive_uniaxial(no, ne):
    return ne > no


def extraordinary_index(theta, no, ne):
    """
    Direction-dependent extraordinary index n_e(theta), theta measured
    from the optic axis: 1/n_e(theta)^2 = cos^2(theta)/no^2 + sin^2(theta)/ne^2.
    Reduces to no at theta=0 and ne at theta=90deg.
    """
    inv_n2 = np.cos(theta) ** 2 / no ** 2 + np.sin(theta) ** 2 / ne ** 2
    return 1.0 / np.sqrt(inv_n2)


def walk_off_angle(theta, no, ne):
    """
    Spatial walk-off angle rho between the Poynting vector and the wave
    vector for the extraordinary ray propagating at angle theta from the
    optic axis: tan(rho) = (ne^2 - no^2) tan(theta) / (ne^2 + no^2 tan^2(theta)).
    """
    t = np.tan(theta)
    return np.arctan2((ne ** 2 - no ** 2) * t, (ne ** 2 + no ** 2 * t ** 2))


def phase_retardation(no, ne, thickness, wavelength):
    """Phase retardation Gamma = 2*pi*thickness*|ne-no| / wavelength (radians)."""
    return 2 * np.pi * thickness * abs(ne - no) / wavelength


def waveplate_thickness(no, ne, wavelength, order=0, fraction=0.25):
    """
    Thickness for a waveplate giving retardation (order + fraction) full
    waves, e.g. fraction=0.25 -> quarter-wave plate, fraction=0.5 -> half-wave.
    A true zero-order plate uses order=0.
    """
    dn = abs(ne - no)
    if dn == 0:
        raise ValueError("no == ne: this material is not birefringent.")
    return (order + fraction) * wavelength / dn


def retarder_type(retardation_rad):
    """Classify a (reduced modulo 2*pi) retardation as quarter-/half-/full-wave, informationally."""
    frac = (retardation_rad / (2 * np.pi)) % 1.0
    labels = {0.25: "quarter-wave", 0.5: "half-wave", 0.75: "three-quarter-wave", 0.0: "full-wave (no net effect)"}
    closest = min(labels, key=lambda k: abs(k - frac))
    return labels[closest] if abs(closest - frac) < 1e-3 else f"{frac:.3f} waves"


def optical_activity_rotation(specific_rotation, thickness):
    """Rotation angle (radians) from optical activity: theta = specific_rotation * thickness
    (specific_rotation in rad per unit length, thickness in the same length unit)."""
    return specific_rotation * thickness


def ordinary_extraordinary_split_angle(theta_incident, no, ne, n_medium=1.0):
    """
    For a uniaxial crystal cut with its optic axis normal to the input
    face, compute the internal refraction angles of the o-ray (ordinary
    Snell's law) and e-ray (using the direction-dependent index) for an
    incident angle theta_incident from air/medium n_medium.
    Returns (theta_o, theta_e_approx), where the e-ray angle is obtained
    by self-consistently solving n_medium*sin(theta_incident) = n_e(theta_e)*sin(theta_e).
    """
    sin_i = n_medium * np.sin(theta_incident)
    theta_o = np.arcsin(np.clip(sin_i / no, -1, 1))

    theta_e = theta_o  # initial guess
    for _ in range(50):
        n_e_theta = extraordinary_index(theta_e, no, ne)
        sin_e = np.clip(sin_i / n_e_theta, -1, 1)
        theta_e_new = np.arcsin(sin_e)
        if abs(theta_e_new - theta_e) < 1e-12:
            theta_e = theta_e_new
            break
        theta_e = theta_e_new
    return theta_o, theta_e
