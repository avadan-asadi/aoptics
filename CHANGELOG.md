# Changelog

All notable changes to this project are documented in this file.
Versioning follows [Semantic Versioning](https://semver.org/) (MAJOR.MINOR.PATCH).

## [2.1.0]
### Added
- `opticspy.nonlinear_optics_boyd` -- nonlinear optics per Boyd's *Nonlinear Optics*:
  `susceptibility`, `coupled_wave_mixing` (SFG/DFG/OPA + Manley-Rowe), `two_level_atom`
  (saturated absorption/dispersion), `self_action` (self-focusing, Marburger collapse
  distance), `stimulated_scattering` (SRS/SBS thresholds), `multiphoton_absorption` (TPA).

## [2.0.0]
### Added
- `opticspy.fourier_optics` -- pupil function, coherent/incoherent PSF, CTF/OTF/MTF,
  resolution criteria, coherent 4f spatial-filtering system.
- `opticspy.sources_detectors` -- four-level laser rate equations, photodetector noise/NEP.
- `opticspy.modulation` -- electro-optics (Pockels effect) and acousto-optics (Bragg/Raman-Nath).
- `opticspy.scattering` -- Rayleigh scattering, metal optics (complex-index reflectance).
- `opticspy.photonics.photonic_crystals` -- multilayer transfer-matrix method, Bragg mirrors.

## [1.4.0]
### Added
- `opticspy.qoptics` -- advanced quantum optics per Gerry & Knight: field operators,
  Jaynes-Cummings model, Lindblad master equation, quantum beamsplitter / Hong-Ou-Mandel,
  g^(2) coherence functions, entanglement tools, quantum information (teleportation, BB84).

## [1.3.0]
### Added
- `opticspy.photonics` -- planar waveguides, optical fibers, resonators (Fabry-Perot +
  Gaussian-mode ABCD self-consistency), chi(2)/chi(3) nonlinear optics basics.

## [1.2.0]
### Added
- `opticspy.coherence` -- temporal coherence (Wiener-Khinchin), spatial coherence
  (Van Cittert-Zernike), speckle statistics.

## [1.1.0]
### Added
- `opticspy.electromagnetic` -- Fresnel equations, Jones/Stokes/Mueller polarization
  calculus, uniaxial crystal optics.

## [1.0.0] (baseline extended with Phase 1)
### Added
- `opticspy.geometrical` -- paraxial/real ray tracing, third-order Seidel aberration theory.
### Existing core (predates this changelog)
- `beams`, `propagation`, `interference`, `diffraction`, `wavefront`, `quantum`, `utils`,
  `visualization`, `gratings`, `talbot`, `moire`.
