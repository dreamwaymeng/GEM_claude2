# GEM Hadron Framework Requirements

General framework using the Gaussian Expansion Method (GEM) for hadron calculations: mesons (qq̄), baryons (qqq), tetraquarks (4q), pentaquarks (5q), hexaquarks (6q).

## Design Principles

- [ ] Extensible architecture (2-body → 6-body)
- [ ] Arbitrary quark masses
- [ ] Modular components reusable across N-body systems
- [ ] S-wave (l=0) Gaussian basis with geometric progression
- [ ] Fully analytical matrix elements (no numerical integration)

## Physics Derivations

- [ ] N-body Jacobi coordinates (general masses)
- [ ] Jacobi transformation matrices (orthogonal rotations in mass-scaled coords)
- [ ] Color wave functions (diquark-antidiquark + molecular)
- [ ] Spin wave functions (N-particle coupling)
- [ ] Color-spin matrix elements for all pairs
- [ ] Permutation antisymmetry (identical quarks)
- [ ] N-body Cholesky decomposition for non-diagonal pairs

## Algorithm Components

- [ ] Gaussian basis generation (geometric series)
- [ ] Overlap/kinetic matrices (analytical)
- [ ] AL1 potential matrices: Coulomb, linear confinement, hyperfine (all pairs)
- [ ] Three-body force (baryons)
- [ ] Generalized eigenvalue solver
- [ ] Results analysis

## Validation

- [ ] Meson masses (π, ρ, J/ψ)
- [ ] Baryon masses (p, Δ, Λc)
- [ ] Tetraquark comparison with literature
- [ ] Pentaquark comparison with literature
