# aoptics Research Library

**[English](README.md) | [فارسی](README.fa.md) | [Türkçe](README.tr.md)**

[![PyPI version](https://img.shields.io/pypi/v/aoptics.svg)](https://pypi.org/project/aoptics/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python](https://img.shields.io/pypi/pyversions/aoptics.svg)](https://pypi.org/project/aoptics/)

A comprehensive Python library for classical and quantum optics research and modeling,
built to systematically cover Born & Wolf's *Principles of Optics*, Saleh & Teich's
*Fundamentals of Photonics*, Goodman's *Statistical/Fourier Optics*, Gerry & Knight's
*Introductory Quantum Optics*, and Boyd's *Nonlinear Optics*.

**190 unit tests** validate the physics throughout (analytic cross-checks against textbook
closed-form results, energy/probability conservation, Monte Carlo agreement with theory,
and cross-validation between independent derivations of the same quantity). Every module
below has a corresponding `tests/test_*.py` file, and every numbered example script is
runnable end-to-end and produces a figure.

## Module map

### Original core (wave optics / scalar diffraction / quantum states)
- `aoptics.beams` -- Gaussian, Laguerre-Gaussian, Hermite-Gaussian, Bessel beams
- `aoptics.propagation` -- angular-spectrum, Fresnel, Fraunhofer, beam-propagation-method (BPM), ABCD Gaussian-beam propagation
- `aoptics.interference` -- two/multi-beam interference, Fabry-Perot, Michelson/Mach-Zehnder, double-slit, phase-shifting interferometry, speckle, holography (record/reconstruct)
- `aoptics.diffraction` -- gratings, Talbot carpets, Moire patterns, Fresnel zone plates
- `aoptics.wavefront` -- Zernike polynomials/decomposition, Strehl ratio, Shack-Hartmann simulation
- `aoptics.quantum` -- Fock/coherent/squeezed/thermal/cat states, Wigner function, photon statistics, Bell states, concurrence
- `aoptics.gratings`, `aoptics.talbot`, `aoptics.moire` -- focused single-topic modules
- `aoptics.utils`, `aoptics.visualization` -- plotting and numerical helpers

### Phase 1 -- Geometrical & electromagnetic optics (Born & Wolf Ch. 1,3-5,14; Saleh & Teich Ch. 1,6)
- `aoptics.geometrical.ray_tracing` -- paraxial ABCD matrix optics, cardinal points, exact real 3-D ray tracing (spherical/conic surfaces, refraction + reflection, TIR, spot diagrams, ray-fan plots)
- `aoptics.geometrical.aberrations` -- third-order Seidel aberration theory (S_I-S_V), wavefront coefficients, Petzval curvature
- `aoptics.electromagnetic.fresnel` -- Fresnel coefficients, Brewster/critical angle, TIR phase
- `aoptics.electromagnetic.polarization` -- Jones calculus, Stokes/Mueller calculus, Poincare sphere
- `aoptics.electromagnetic.crystal_optics` -- index ellipsoid, birefringence, waveplate design, spatial walk-off

### Phase 2 -- Coherence & statistical optics (Born & Wolf Ch. 10; Goodman)
- `aoptics.coherence.temporal` -- Wiener-Khinchin theorem, numerical complex degree of coherence, coherence time/length, Michelson interferograms
- `aoptics.coherence.spatial` -- Van Cittert-Zernike theorem (circular/slit/Gaussian sources), coherence radius, extended-source visibility
- `aoptics.coherence.speckle` -- negative-exponential/Gamma intensity statistics, Monte Carlo random-phasor-sum simulation, contrast reduction

### Phase 3 -- Photonics (Saleh & Teich Ch. 7-9,21)
- `aoptics.photonics.waveguides` -- planar-waveguide TE modes via exact dispersion-relation root finding, confinement factor
- `aoptics.photonics.fibers` -- V-number, NA, single-mode cutoff, Marcuse mode-field radius, Sellmeier/material dispersion, zero-dispersion wavelength
- `aoptics.photonics.resonators` -- Fabry-Perot (finesse, FSR, Airy function), stability, Kogelnik ABCD self-consistent Gaussian mode
- `aoptics.photonics.nonlinear_optics` -- SHG phase matching, SPM/Kerr effect, soliton order
- `aoptics.photonics.photonic_crystals` -- multilayer transfer-matrix method, quarter-wave Bragg mirrors, stopband width

### Phase 4 -- Advanced quantum optics (Gerry & Knight, full text)
- `aoptics.qoptics.operators` -- field operators in a truncated Fock basis, displacement/squeeze operators, quadratures
- `aoptics.qoptics.jaynes_cummings` -- full Jaynes-Cummings Hamiltonian, vacuum Rabi oscillation, coherent-state collapse & revival, dressed states
- `aoptics.qoptics.master_equation` -- general Lindblad integrator (cavity decay, cat-state decoherence)
- `aoptics.qoptics.beamsplitter` -- two-mode quantum beamsplitter, Hong-Ou-Mandel interference
- `aoptics.qoptics.coherence_functions` -- g^(2)(0), Fano factor, antibunching
- `aoptics.qoptics.entanglement` -- partial trace, von Neumann entropy, fidelity, Schmidt decomposition
- `aoptics.qoptics.quantum_information` -- qubit gates, quantum teleportation, BB84 QKD simulation (with eavesdropper QBER signature)

### Phase 5 -- Remaining chapters: Fourier optics, sources/detectors, modulation, scattering
- `aoptics.fourier_optics.transfer_functions` -- pupil function, coherent/incoherent PSF, CTF, OTF/MTF, Rayleigh/Abbe resolution
- `aoptics.fourier_optics.spatial_filtering` -- coherent 4f system, low/high-pass and phase-contrast spatial filters
- `aoptics.sources_detectors.lasers` -- four-level laser rate equations, threshold, gain clamping, slope efficiency
- `aoptics.sources_detectors.photodetectors` -- responsivity, quantum efficiency, shot/thermal noise, NEP
- `aoptics.modulation.electrooptics` -- Pockels effect, half-wave voltage, Mach-Zehnder EO modulator
- `aoptics.modulation.acoustooptics` -- Klein-Cook parameter, Bragg angle, Raman-Nath/Bragg diffraction efficiency
- `aoptics.scattering.rayleigh_scattering` -- scattering cross section, 1/lambda^4 law, sky-polarization degree
- `aoptics.scattering.metal_optics` -- complex-index reflectance, skin depth, absorptance

### Phase 6 -- Nonlinear optics per Boyd's textbook chapters (Boyd, "Nonlinear Optics")
- `aoptics.nonlinear_optics_boyd.susceptibility` -- polarization expansion, d-coefficient, Miller's rule, Kleinman-symmetry tensor check (Ch. 1)
- `aoptics.nonlinear_optics_boyd.coupled_wave_mixing` -- three-wave-mixing (SFG/DFG/OPA) coupled-amplitude ODE integration, Manley-Rowe relations (Ch. 2)
- `aoptics.nonlinear_optics_boyd.two_level_atom` -- saturated/power-broadened Lorentzian absorption & dispersion, steady-state Bloch-equation inversion (Ch. 3, 6)
- `aoptics.nonlinear_optics_boyd.self_action` -- intensity-dependent index, Marburger critical power & collapse distance, SPM phase (Ch. 4, 7)
- `aoptics.nonlinear_optics_boyd.stimulated_scattering` -- stimulated Raman/Brillouin frequency shifts, R.G. Smith threshold criteria (Ch. 8-10)
- `aoptics.nonlinear_optics_boyd.multiphoton_absorption` -- two-photon-absorption propagation/transmission, saturable-absorber model (Ch. 12)

## Known, deliberate scope limits

A few topics were intentionally left out because the available formulas carry real notational-
convention risk across sources (better to omit than to silently guess): Gaussian-Schell-model
partially-coherent beam propagation, full Mie scattering (only the small-particle Rayleigh limit
is covered), detailed mode-locked ultrafast-pulse dynamics, the full optical-Bloch-equation
coherence (u,v) components (only the population inversion and the absorption/dispersion
lineshape are provided), and exact absolute-unit parametric-gain/SHG-efficiency prefactors
(provided in normalized form instead, following the same reasoning as the Phase-3 SHG efficiency
function). Ion-trap quantum computation and semiconductor band-structure optics are out of scope
as specialized sub-fields adjacent to, but not core content of, the reference texts.

## Development

```bash
pip install -e .[dev]
pytest tests/ -q

# one demo per phase, each produces a PNG figure
python examples/phase1_geometrical_demo.py
python examples/phase1_polarization_demo.py
python examples/phase2_coherence_demo.py
python examples/phase3_photonics_demo.py
python examples/phase4_qoptics_demo.py
python examples/phase5_extended_demo.py
python examples/phase6_boyd_nonlinear_demo.py
```

## Install
```bash
pip install aoptics
```

## For maintainers: publishing to PyPI

See `docs/PUBLISHING.en.md` (also available in [Persian](docs/PUBLISHING.fa.md) and
[Turkish](docs/PUBLISHING.tr.md)) for the full guide: building the package, uploading to
PyPI/TestPyPI, and setting up automatic publishing via GitHub Actions so every new GitHub
Release is published to PyPI automatically.

## Quick Start
```python
import aoptics as op
import numpy as np

# Generate Laguerre-Gaussian beam
field = op.beams.laguerre_gaussian(p=0, l=2)
op.utils.show_field(field, title='LG(0,2) beam')

# Interference pattern
lg = op.beams.laguerre_gaussian(0, 1)
ref = op.beams.gaussian_beam()
I, phase, vis = op.interference.two_beam_interference(lg, ref)
op.utils.show_interference(I, phase)

# A geometrical-optics singlet lens with Seidel aberrations (Phase 1)
from aoptics.geometrical import ray_tracing as rt, aberrations as ab
surfaces = [rt.Surface(radius=51.5, thickness=5.3, index=1.5168, semi_diameter=12.5),
            rt.Surface(radius=np.inf, thickness=0.0, index=1.0, semi_diameter=12.5)]
system = rt.OpticalSystem(surfaces)
print(system.cardinal_points())

# A quantum Jaynes-Cummings vacuum-Rabi simulation (Phase 4)
from aoptics.qoptics import jaynes_cummings as jc
from aoptics.quantum import fock_state
H = jc.hamiltonian(field_dim=5, omega_cavity=0, omega_atom=0, g=1.0)
psi0 = jc.initial_state(fock_state(0, 5), jc.ATOM_EXCITED)
```

## License

MIT -- see [LICENSE](LICENSE).


