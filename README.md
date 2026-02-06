# GEM Hadron Framework

A computational framework for calculating hadron properties using the **Gaussian Expansion Method (GEM)**. This repository contains algorithm documentation for solving few-body quantum mechanical problems in QCD, including mesons, baryons, tetraquarks, and higher multiquark systems.

## Overview

The Gaussian Expansion Method is a variational approach that expands wave functions using Gaussian basis functions with geometrically distributed width parameters. Combined with Jacobi coordinates, this enables **fully analytical matrix element calculations** without numerical integration.

### Supported Systems

| System | Quarks | Jacobi Coords | Example Hadrons |
|--------|--------|---------------|-----------------|
| Meson | qq̄ | 1 | π, ρ, J/ψ, Υ |
| Baryon | qqq | 2 | p, n, Δ, Λ, Ξ, Ω |
| Tetraquark | qqq̄q̄ | 3 | X(3872), Z_c, T_cc |
| Pentaquark | qqqqq̄ | 4 | P_c states |

## Documentation Structure

### Core Algorithm Documents

| File | Description |
|------|-------------|
| [`physics.md`](physics.md) | Physical foundations: GEM method, AL1 potential model, Jacobi coordinates, Cholesky decomposition for non-diagonal pairs |
| [`algorithm_spatial_MM.md`](algorithm_spatial_MM.md) | Spatial matrix elements: Jacobi coordinate transformations, width matrix construction, Cholesky factorization algorithm |
| [`algorithm_discrete_MM.md`](algorithm_discrete_MM.md) | Discrete wave functions: spin/color/isospin construction via eigenvalue method, pairwise operators |
| [`algorithm_antisymmetrization.md`](algorithm_antisymmetrization.md) | Fermion antisymmetrization: permutation operators, identical particle handling, combined matrix elements |
| [`requirements.md`](requirements.md) | Project requirements and validation targets |

### Key Features

- **Fully analytical**: All matrix elements computed without numerical integration
- **Extensible**: Unified framework from 2-body to 6-body systems
- **Arbitrary masses**: General formulas for any quark mass configuration
- **Modular design**: Spatial and discrete sectors computed independently, combined via tensor products

## Physical Model

### AL1 Quark Potential

The quark-antiquark interaction includes three terms:

1. **Coulomb** (one-gluon exchange): `V_C = α_s/r · (λ_i·λ_j)/4`
2. **Linear confinement**: `V_conf = (-3b/4 · r + V_c) · (λ_i·λ_j)/4`
3. **Hyperfine** (spin-spin): `V_hyp ∝ exp(-τ²r²) · (λ_i·λ_j) · (s_i·s_j)`

### Quark Masses (GeV)

| u, d | s | c | b |
|------|---|---|---|
| 0.315 | 0.577 | 1.836 | 5.227 |

## Algorithm Highlights

### Cholesky Decomposition Method

For non-diagonal pair potentials (pairs not aligned with Jacobi coordinates), the Cholesky method:

1. Transforms to the Jacobi set where the target pair is diagonal
2. Constructs combined width matrix from bra and ket Gaussians
3. Applies Cholesky factorization to decouple integrals
4. Evaluates factorized integrals analytically

### Antisymmetrization via Permutations

Instead of Young tableaux, directly sum over permutations:

```
⟨Ψ_a|O|AΨ_b⟩ = (1/N!) Σ_P sign(P) ⟨Ψ_a|O|PΨ_b⟩
```

This factors into spatial and discrete parts:
- **Spatial**: Permutation induces Jacobi coordinate transformation
- **Discrete**: Permutation matrix acts on tensor product basis

### Discrete Wave Functions

Spin, color, and isospin wave functions constructed via:

1. Build group generators in N-particle tensor product space
2. Construct and diagonalize Casimir operator
3. Select eigenstates with target quantum numbers

No Clebsch-Gordan tables required—pure linear algebra.

## Getting Started

1. Review [`physics.md`](physics.md) for the theoretical foundation
2. Study [`algorithm_spatial_MM.md`](algorithm_spatial_MM.md) for spatial matrix elements
3. Study [`algorithm_discrete_MM.md`](algorithm_discrete_MM.md) for discrete wave functions
4. See [`algorithm_antisymmetrization.md`](algorithm_antisymmetrization.md) for identical fermion handling

## Validation Targets

- Meson masses: π, ρ, J/ψ, Υ
- Baryon masses: p, Δ, Λ_c
- Tetraquark and pentaquark comparison with literature

## References

- Gaussian Expansion Method: Hiyama et al., Prog. Part. Nucl. Phys. 51 (2003) 223
- AL1 Potential: Semay & Silvestre-Brac, Z. Phys. C 61 (1994) 271
- Multiquark systems: Various QCD sum rule and lattice QCD comparisons

## License

This project is for academic research purposes.
