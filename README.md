# GEM Hadron Framework

A C++ computational framework for calculating hadron properties using the **Gaussian Expansion Method (GEM)**. This repository contains architecture design and algorithm documentation for solving few-body quantum mechanical problems in QCD, including mesons, baryons, tetraquarks, and pentaquarks.

## Overview

The Gaussian Expansion Method is a variational approach that expands wave functions using Gaussian basis functions with geometrically distributed width parameters. Combined with Jacobi coordinates, this enables **fully analytical matrix element calculations** without numerical integration.

### Supported Systems

| System | Quarks | Jacobi Coords | Basis Size (N=15) | Example Hadrons |
|--------|--------|---------------|-------------------|-----------------|
| Meson | qq̄ | 1 | 15 | π, ρ, J/ψ, Υ |
| Baryon | qqq | 2 | 225 | p, n, Δ, Λ, Ξ, Ω |
| Tetraquark | qqq̄q̄ | 3 | 3,375 | X(3872), Z_c, T_cc |
| Pentaquark | qqqqq̄ | 4 | 50,625 | P_c states |

### Technology Stack

| Component | Library | Purpose |
|-----------|---------|---------|
| Language | C++17/20 | Core implementation |
| Matrix Storage | Eigen | Matrix operations, linear algebra |
| Eigenvalue Solver | LAPACK (MKL/OpenBLAS) | Primary solver with truncation |
| Iterative Solver | Spectra | Large matrices, few eigenvalues |
| Distributed Solver | ELPA | HPC cluster (pentaquark) |
| Configuration | yaml-cpp | YAML parsing |
| Parallelization | OpenMP / MPI | Shared/distributed memory |

## Documentation Structure

### Architecture

| File | Description |
|------|-------------|
| [`ARCHITECTURE.md`](ARCHITECTURE.md) | **Complete C++ software architecture**: module structure, solver design, parallelization strategy, HPC deployment |

### Algorithm Documents

| File | Description |
|------|-------------|
| [`physics.md`](physics.md) | Physical foundations: GEM method, AL1 potential model, Jacobi coordinates, Cholesky decomposition |
| [`algorithm_spatial_MM.md`](algorithm_spatial_MM.md) | Spatial matrix elements: Jacobi transformations, width matrix construction, Cholesky factorization |
| [`algorithm_discrete_MM.md`](algorithm_discrete_MM.md) | Discrete wave functions: spin/color/isospin via eigenvalue method, pairwise operators |
| [`algorithm_antisymmetrization.md`](algorithm_antisymmetrization.md) | Fermion antisymmetrization: permutation operators, identical particle handling |
| [`algorithm_eigenvalue_truncation.md`](algorithm_eigenvalue_truncation.md) | Numerical stability: eigenvalue truncation for ill-conditioned overlap matrices |
| [`requirements.md`](requirements.md) | Project requirements and validation targets |

### Key Features

- **Fully analytical**: All matrix elements computed without numerical integration
- **Extensible**: Unified framework from 2-body to 5-body systems
- **Arbitrary masses**: General formulas for any quark mass configuration
- **Modular design**: Spatial and discrete sectors computed independently
- **Numerically stable**: Eigenvalue truncation for ill-conditioned bases
- **HPC ready**: MPI+OpenMP hybrid parallelization for pentaquark scale

## Physical Model

### AL1 Quark Potential

The quark-antiquark interaction includes three terms:

1. **Coulomb** (one-gluon exchange): `V_C = α_s/r · (λ_i·λ_j)/4`
2. **Linear confinement**: `V_conf = (-3b/4 · r + V_c) · (λ_i·λ_j)/4`
3. **Hyperfine** (spin-spin): `V_hyp ∝ exp(-τ²r²) · (λ_i·λ_j) · (s_i·s_j)`

### Three-Body Force (Baryons)

```
V_3body = -C_3 / (m_1 · m_2 · m_3)
```

where C_3 = 2.02 × 10⁻³ GeV⁴.

### Quark Masses (GeV)

| u, d | s | c | b |
|------|---|---|---|
| 0.315 | 0.577 | 1.836 | 5.227 |

## Algorithm Highlights

### Cholesky Decomposition Method

For non-diagonal pair potentials (pairs not aligned with Jacobi coordinates):

1. Transform to the Jacobi set where the target pair is diagonal
2. Construct combined width matrix from bra and ket Gaussians
3. Apply Cholesky factorization to decouple integrals
4. Evaluate factorized integrals analytically

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

### Eigenvalue Truncation

For large non-orthogonal basis sets:

1. Diagonalize the overlap matrix S
2. Remove near-null directions (eigenvalues < ε · λ_max)
3. Transform H to orthonormal basis
4. Solve standard eigenvalue problem

Ensures numerical stability while preserving the variational principle.

## Deployment Strategy

| System | Hardware | Parallelization | Solver |
|--------|----------|-----------------|--------|
| Meson, Baryon, Tetraquark | Personal computer | OpenMP | LAPACK |
| Pentaquark | HPC cluster | MPI + OpenMP | ELPA |

## Getting Started

1. Review [`ARCHITECTURE.md`](ARCHITECTURE.md) for the software design
2. Review [`physics.md`](physics.md) for the theoretical foundation
3. Study algorithm documents for implementation details
4. See [`requirements.md`](requirements.md) for validation targets

## Implementation Order

1. **Phase 1**: Core framework (meson) - Eigen, LAPACK, config parsing
2. **Phase 2**: Baryon extension - spin/color coupling, antisymmetrization, three-body force
3. **Phase 3**: Tetraquark - discrete state reduction, Spectra solver
4. **Phase 4**: HPC/Pentaquark - MPI, checkpoint/restart, ELPA

## Validation Targets

- Meson masses: π (138 MeV), ρ (775 MeV), J/ψ (3097 MeV)
- Baryon masses: p (938 MeV), Δ (1232 MeV)
- Tetraquark and pentaquark comparison with literature

## References

- Gaussian Expansion Method: Hiyama et al., Prog. Part. Nucl. Phys. 51 (2003) 223
- AL1 Potential: Semay & Silvestre-Brac, Z. Phys. C 61 (1994) 271
- Multiquark systems: Various QCD sum rule and lattice QCD comparisons

## License

This project is for academic research purposes.
