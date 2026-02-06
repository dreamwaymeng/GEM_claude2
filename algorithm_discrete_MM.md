# Discrete Wave Functions and Matrix Elements via Eigenvalue Method

## 1. Overview

**Goal**: Construct spin/isospin/color wave functions and compute matrix elements using a universal eigenvalue method that works identically for SU(2) and SU(3).

**Method**:
1. Build group generators in N-particle tensor product space
2. Construct Casimir operator, diagonalize to find eigenstates
3. Compute matrix elements directly via matrix multiplication

| Sector | Group | Single-particle dim | Generators | Casimir |
|--------|-------|---------------------|------------|---------|
| Spin | SU(2) | 2 (all quarks) | σ/2 | S² |
| Isospin | SU(2) | 2 (u,d) or 1 (s,c,b) | τ/2 or 0 | I² |
| Color | SU(3) | 3 (all quarks) | λ/2 | C₂ |

---

## 2. Generator Matrices

### 2.1 SU(2): Pauli Matrices

```text
σₓ = [0, 1; 1, 0]    σᵧ = [0, -i; i, 0]    σᵤ = [1, 0; 0, -1]

Spin operators: Sₐ = σₐ/2
```

### 2.2 Isospin: Flavor-Dependent

**Quark isospin quantum numbers:**

| Quark | I | I₃ | Dimension |
|-------|---|-----|-----------|
| u | 1/2 | +1/2 | 2 |
| d | 1/2 | -1/2 | 2 |
| s, c, b | 0 | 0 | 1 (trivial) |

**Isospin operators:**
- For u, d quarks: Iₐ = τₐ/2 (same as spin Pauli matrices)
- For s, c, b quarks: Iₐ = 0 (1×1 zero matrix, trivial singlet)

**Tensor product space dimension:**
```text
dim(isospin) = Π_k dim_k

where dim_k = 2 if quark k is u or d
              1 if quark k is s, c, or b
```

**Examples:**
- Meson (uū): dim = 2 × 2 = 4
- Meson (cc̄): dim = 1 × 1 = 1 (trivial)
- Baryon (uds): dim = 2 × 2 × 1 = 4
- Baryon (udc): dim = 2 × 2 × 1 = 4
- Baryon (uuu): dim = 2 × 2 × 2 = 8

### 2.3 SU(3): Gell-Mann Matrices

```text
λ₁ = [0,1,0; 1,0,0; 0,0,0]    λ₂ = [0,-i,0; i,0,0; 0,0,0]    λ₃ = [1,0,0; 0,-1,0; 0,0,0]
λ₄ = [0,0,1; 0,0,0; 1,0,0]    λ₅ = [0,0,-i; 0,0,0; i,0,0]    λ₆ = [0,0,0; 0,0,1; 0,1,0]
λ₇ = [0,0,0; 0,0,-i; 0,i,0]   λ₈ = [1,0,0; 0,1,0; 0,0,-2]/√3
```

**Quark color operators** (representation **3**):
```text
Tₐ = λₐ/2
```

**Antiquark color operators** (representation **3̄**):
```text
T̄ₐ = -λₐ*/2 = -λₐᵀ/2
```

Properties of Gell-Mann matrices:
- Real symmetric (λ₁, λ₃, λ₄, λ₆, λ₈): λ* = λ, so T̄ₐ = -Tₐ
- Imaginary antisymmetric (λ₂, λ₅, λ₇): λ* = -λ, so T̄ₐ = +Tₐ

---

## 3. N-Particle Space Construction

### 3.1 Dimensions

| System | N | Spin (2^N) | Color (3^N) |
|--------|---|------------|-------------|
| Meson | 2 | 4 | 9 |
| Baryon | 3 | 8 | 27 |
| Tetraquark | 4 | 16 | 81 |

### 3.2 Basis Convention

Lexicographic ordering: `index = Σₖ iₖ × d^(N-1-k)` where iₖ ∈ {0, ..., d-1}

Example (2-particle spin): index 0→|↑↑⟩, 1→|↑↓⟩, 2→|↓↑⟩, 3→|↓↓⟩

### 3.3 Single-Particle Operator in N-Particle Space

```text
build_single_particle_op(Op, i, N, d):
    # Op at position i, identity elsewhere
    result = I(1)
    for k = 1 to N:
        result = kron(result, Op if k==i else eye(d))
    return result
```

### 3.4 Total Generator

```text
G_total = Σᵢ build_single_particle_op(G, i, N, d)
```

---

## 4. Casimir Operators and Eigenvalues

### 4.1 SU(2) Casimir

```text
S² = Sₓ² + Sᵧ² + Sᵤ²
Eigenvalue: S(S+1)
```

| S | Eigenvalue | Degeneracy |
|---|------------|------------|
| 0 | 0 | 1 |
| 1/2 | 3/4 | 2 |
| 1 | 2 | 3 |
| 3/2 | 15/4 | 4 |

### 4.2 SU(3) Casimir

```text
C₂ = Σₐ Tₐ²
```

| Representation | C₂ |
|----------------|-----|
| **1** (singlet) | 0 |
| **3** or **3̄** | 4/3 |
| **6** or **6̄** | 10/3 |
| **8** (octet) | 3 |

### 4.3 Tensor Product Decompositions

**Spin** (spin-1/2 particles):
- N=2: 0 ⊕ 1 (1 singlet, 1 triplet)
- N=3: 1/2 ⊕ 1/2 ⊕ 3/2 (2 doublets, 1 quartet)
- N=4: 0 ⊕ 0 ⊕ 1 ⊕ 1 ⊕ 1 ⊕ 2

**Color** (q in **3**, q̄ in **3̄**):
- Meson (qq̄): **3⊗3̄** = **1** ⊕ **8**
- Baryon (qqq): **3⊗3⊗3** = **1** ⊕ **8** ⊕ **8** ⊕ **10**
- Tetraquark (qqq̄q̄): **3⊗3⊗3̄⊗3̄** = **1** ⊕ **1** ⊕ **8**⁴ ⊕ **10** ⊕ **1̄0** ⊕ **27**

---

## 5. Wave Function Construction Algorithm

### 5.1 Spin Wave Functions (All Quarks Have Spin-1/2)

For spin, all quarks have spin-1/2, so dimension is 2 for each particle:

```text
construct_spin_wave_functions(N, target_S, target_Sz=None):
    generators = [σₓ/2, σᵧ/2, σᵤ/2]

    # Build total generators (all particles have d=2)
    G_total = [Σᵢ build_single_particle_op(G, i, N, d=2) for G in generators]

    # Build Casimir S² = Sₓ² + Sᵧ² + Sᵤ²
    C = Σₐ G_total[a] @ G_total[a]
    target_casimir = target_S * (target_S + 1)

    # Diagonalize and select eigenspace
    eigenvalues, eigenvectors = eig(C)
    eigenspace = eigenvectors[:, |eigenvalues - target_casimir| < tol]

    # Project to Sz if specified
    if target_Sz is not None:
        Sz_proj = eigenspace.H @ G_total[z] @ eigenspace
        sz_vals, sz_vecs = eig(Sz_proj)
        states = [eigenspace @ sz_vecs[:, k] for k where |sz_vals[k] - target_Sz| < tol]
    else:
        states = orthonormalize(eigenspace)

    return states
```

### 5.2 Isospin Wave Functions (Flavor-Dependent)

For isospin, only u and d quarks have I=1/2 (dim=2). Other quarks have I=0 (dim=1, trivial).

```text
construct_isospin_wave_functions(N, flavors, target_I, target_I3=None):
    # flavors: array of quark flavors, e.g., ['u', 'd', 's']
    # Determine dimension for each particle
    dims = [2 if f in ['u', 'd'] else 1 for f in flavors]
    total_dim = product(dims)

    # If all quarks have I=0 (e.g., cc̄, sss), isospin is trivially 0
    if total_dim == 1:
        return [array([1.0])]  # trivial singlet

    # Build generators with flavor-dependent dimensions
    # Use tensor product with identity for I=0 quarks
    I_total = [zeros(total_dim) for a in [x, y, z]]

    for i in 1..N:
        if flavors[i] in ['u', 'd']:
            for a in [x, y, z]:
                gen = τₐ/2  # 2×2 Pauli matrix
                I_total[a] += build_mixed_dim_op(gen, i, dims)
        # else: I=0 quarks contribute nothing (generator is 0)

    # Build Casimir I² = Iₓ² + Iᵧ² + Iᵤ²
    C = Σₐ I_total[a] @ I_total[a]
    target_casimir = target_I * (target_I + 1)

    # Diagonalize and select eigenspace
    eigenvalues, eigenvectors = eig(C)
    eigenspace = eigenvectors[:, |eigenvalues - target_casimir| < tol]

    # Project to I3 if specified
    if target_I3 is not None:
        I3_proj = eigenspace.H @ I_total[z] @ eigenspace
        i3_vals, i3_vecs = eig(I3_proj)
        states = [eigenspace @ i3_vecs[:, k] for k where |i3_vals[k] - target_I3| < tol]
    else:
        states = orthonormalize(eigenspace)

    return states

# Helper: build operator with mixed dimensions
build_mixed_dim_op(Op, i, dims):
    # dims: list of dimensions for each particle
    # Op: operator for particle i (dims[i] × dims[i])
    result = 1
    for k in 1..N:
        if k == i:
            result = kron(result, Op)
        else:
            result = kron(result, eye(dims[k]))
    return result
```

**Example isospin states:**

| System | Flavors | Isospin dim | I states |
|--------|---------|-------------|----------|
| π⁺ (ud̄) | [u, d̄] | 2×2=4 | I=0,1 |
| K⁺ (us̄) | [u, s̄] | 2×1=2 | I=1/2 only |
| J/ψ (cc̄) | [c, c̄] | 1×1=1 | I=0 only (trivial) |
| Λ (uds) | [u, d, s] | 2×2×1=4 | I=0,1 |
| Ξ (uss) | [u, s, s] | 2×1×1=2 | I=1/2 only |
| Ω (sss) | [s, s, s] | 1×1×1=1 | I=0 only (trivial) |

### 5.3 Color Wave Functions with Antiquarks

For color, all quarks have dim=3, but must use correct generators for each particle type:

```text
construct_color_wave_functions(N, particle_types, target_C2):
    # particle_types: array of 'q' or 'qbar' for each particle
    # e.g., meson: ['q', 'qbar'], baryon: ['q', 'q', 'q']

    lambda_matrices = [λ₁, λ₂, ..., λ₈]

    # Build total color generators accounting for particle types
    T_total = []
    for a in 1..8:
        T_a_total = zeros(3^N, 3^N)
        for i in 1..N:
            if particle_types[i] == 'q':
                gen = λₐ/2
            else:  # 'qbar'
                gen = -λₐᵀ/2
            T_a_total += build_single_particle_op(gen, i, N, d=3)
        T_total.append(T_a_total)

    # Build Casimir C₂ = Σₐ Tₐ²
    C = Σₐ T_total[a] @ T_total[a]

    # Diagonalize and select eigenspace with target C₂
    eigenvalues, eigenvectors = eig(C)
    eigenspace = eigenvectors[:, |eigenvalues - target_C2| < tol]
    states = orthonormalize(eigenspace)

    return states
```

### 5.4 Example: Meson Color Singlet

```text
# Meson: quark at position 1, antiquark at position 2
particle_types = ['q', 'qbar']
color_states = construct_color_wave_functions(N=2, particle_types, target_C2=0)

# Result: |singlet⟩ = (|rr̄⟩ + |gḡ⟩ + |bb̄⟩)/√3
# Note: r̄, ḡ, b̄ are anti-colors (basis states for antiquark in 3̄)
```

---

## 6. Pairwise Operators and Matrix Elements

### 6.1 Spin-Spin Operator

```text
sᵢ·sⱼ = Σₐ sᵢᵃ sⱼᵃ

build_spin_spin(i, j, N):
    result = 0
    for a in [x, y, z]:
        s_i = build_single_particle_op(σₐ/2, i, N, d=2)
        s_j = build_single_particle_op(σₐ/2, j, N, d=2)
        result += s_i @ s_j
    return result
```

**Expectation values** (spin-1/2 pairs):

| Pair spin | ⟨sᵢ·sⱼ⟩ |
|-----------|---------|
| S=0 (singlet) | -3/4 |
| S=1 (triplet) | +1/4 |

*Derivation*: ⟨sᵢ·sⱼ⟩ = (1/2)[S₁₂(S₁₂+1) - 3/4 - 3/4]

### 6.2 Color-Color Operator

The color-color operator depends on particle types (quark vs antiquark):

```text
# Quark-quark pair
build_color_color_qq(i, j, N):
    result = 0
    for a in 1..8:
        T_i = build_single_particle_op(λₐ/2, i, N, d=3)
        T_j = build_single_particle_op(λₐ/2, j, N, d=3)
        result += T_i @ T_j
    return result

# Quark-antiquark pair
build_color_color_qqbar(i, j, N):
    result = 0
    for a in 1..8:
        T_i = build_single_particle_op(λₐ/2, i, N, d=3)
        T_j = build_single_particle_op(-λₐᵀ/2, j, N, d=3)  # antiquark
        result += T_i @ T_j
    return result

# Antiquark-antiquark pair
build_color_color_qbarqbar(i, j, N):
    result = 0
    for a in 1..8:
        T_i = build_single_particle_op(-λₐᵀ/2, i, N, d=3)
        T_j = build_single_particle_op(-λₐᵀ/2, j, N, d=3)
        result += T_i @ T_j
    return result
```

**General form with particle type array**:
```text
build_color_color(i, j, N, particle_types):
    # particle_types[k] = 'q' or 'qbar'
    result = 0
    for a in 1..8:
        gen_i = λₐ/2 if particle_types[i]=='q' else -λₐᵀ/2
        gen_j = λₐ/2 if particle_types[j]=='q' else -λₐᵀ/2
        T_i = build_single_particle_op(gen_i, i, N, d=3)
        T_j = build_single_particle_op(gen_j, j, N, d=3)
        result += T_i @ T_j
    return result
```

**Expectation values** (color singlets):

| System | Pair type | ⟨Tᵢ·Tⱼ⟩ |
|--------|-----------|---------|
| Meson (qq̄) | q-q̄ | -4/3 |
| Baryon (qqq) | q-q | -2/3 |
| Tetraquark (qqq̄q̄) | q-q | depends on config |
| Tetraquark (qqq̄q̄) | q̄-q̄ | depends on config |
| Tetraquark (qqq̄q̄) | q-q̄ | depends on config |

*Derivation*: Using Tᵢ·Tⱼ = (1/2)[C₂(i+j) - C₂(i) - C₂(j)] where C₂(q) = C₂(q̄) = 4/3

### 6.3 Direct Matrix Element Computation

```text
⟨ψₐ|Ôᵢⱼ|ψᵦ⟩ = ψₐ† @ Oᵢⱼ @ ψᵦ
```

No Clebsch-Gordan tables needed — just matrix multiplication.

---

## 7. Verification Data

### 7.1 Meson Spin States

Singlet: |χ₀⟩ = (|↑↓⟩ - |↓↑⟩)/√2 → ⟨s₁·s₂⟩ = -3/4
Triplet: |χ₁⟩ = (|↑↓⟩ + |↓↑⟩)/√2 → ⟨s₁·s₂⟩ = +1/4

### 7.2 Baryon Spin States (S=1/2)

Two independent doublets exist. Using Sz=+1/2 states:

**Λ-type** (pair 1,2 in singlet):
```text
|χ_Λ⟩ = (|↑↓↑⟩ - |↓↑↑⟩)/√2
⟨s₁·s₂⟩ = -3/4,  ⟨s₁·s₃⟩ = 0,  ⟨s₂·s₃⟩ = 0
```

**Σ-type** (pair 1,2 in triplet):
```text
|χ_Σ⟩ = (2|↑↑↓⟩ - |↑↓↑⟩ - |↓↑↑⟩)/√6
⟨s₁·s₂⟩ = +1/4,  ⟨s₁·s₃⟩ = -1/2,  ⟨s₂·s₃⟩ = -1/2
```

**Sum rule**: ⟨s₁·s₂⟩ + ⟨s₁·s₃⟩ + ⟨s₂·s₃⟩ = [S(S+1) - 9/4]/2 = -3/4 for S=1/2

### 7.3 Meson Color Singlet (qq̄)

```text
particle_types = ['q', 'qbar']
|φ⟩ = (|rr̄⟩ + |gḡ⟩ + |bb̄⟩)/√3
```

Here r̄, ḡ, b̄ are the anti-color basis states for the antiquark.
Pair: ⟨T₁·T̄₂⟩ = -4/3

### 7.4 Baryon Color Singlet (qqq)

```text
particle_types = ['q', 'q', 'q']
|φ⟩ = (|rgb⟩ - |rbg⟩ + |gbr⟩ - |grb⟩ + |brg⟩ - |bgr⟩)/√6
```

Totally antisymmetric (εᵢⱼₖ). Each pair: ⟨Tᵢ·Tⱼ⟩ = -2/3

### 7.5 Tetraquark Color Singlets (qqq̄q̄)

```text
particle_types = ['q', 'q', 'qbar', 'qbar']
```

Two independent color singlets exist:
1. Diquark-antidiquark: (qq)₃̄ ⊗ (q̄q̄)₃ → singlet
2. Meson-meson type: (qq̄)₁ ⊗ (qq̄)₁ or (qq̄)₈ ⊗ (qq̄)₈ → singlet

The eigenvalue method finds both automatically. Distinguish by pair expectation values.

---

## 8. Reduced Space Method for Large N

### 8.1 Problem

Full space dimensions grow as 2^N (spin) and 3^N (color). For N≥5, this becomes expensive.

### 8.2 Solution: Highest Weight Projection

1. Find eigenstates with target (S², Sz=S) or (C₂, T₃=T₈=0)
2. Project operators to this reduced subspace
3. Compute matrix elements in reduced space

```text
# Get projection matrix
P = [eigenvectors with target quantum numbers]

# Reduce operators
O_reduced = P† @ O_full @ P

# Matrix elements in reduced space
H_reduced[a,b] = states[a]† @ O_reduced @ states[b]
```

### 8.3 Dimension Comparison

| N | Full spin | S=1/2 reduced |
|---|-----------|---------------|
| 3 | 8 | 2 |
| 4 | 16 | 0 (S=0: 2) |
| 5 | 32 | 5 |
| 6 | 64 | 5 (S=0) |

---

## 9. Implementation Summary

```text
# 1. Define particle content
flavors = ['u', 'd']  # meson (e.g., pion)
particle_types = ['q', 'qbar']
# flavors = ['u', 'd', 's'], particle_types = ['q', 'q', 'q']  # Λ baryon
# flavors = ['c', 'c'], particle_types = ['q', 'qbar']  # J/ψ

# 2. Build pairwise operators
for i,j in pairs:
    s_dot_s[i,j] = build_spin_spin(i, j, N)
    color_op[i,j] = build_color_color(i, j, N, particle_types)

# 3. Construct wave functions
spin_states = construct_spin_wave_functions(N, S, Sz)
isospin_states = construct_isospin_wave_functions(N, flavors, I, I3)
color_states = construct_color_wave_functions(N, particle_types, C2_target=0)

# 4. Compute matrix elements
for a, b in state_pairs:
    spin_ME[a,b][i,j] = spin_states[a]† @ s_dot_s[i,j] @ spin_states[b]
    color_ME[a,b][i,j] = color_states[a]† @ color_op[i,j] @ color_states[b]
    # Isospin: typically diagonal or trivial for heavy quarks

# 5. Assemble Hamiltonian
H[a,b] = Σᵢⱼ A[i,j] * spin_ME[a,b][i,j] * color_ME[a,b][i,j]
```

**Key features**:
- Spin: all quarks have dim=2
- Isospin: u,d have dim=2; s,c,b have dim=1 (trivial I=0)
- Color: all quarks have dim=3; antiquarks use T̄ = -λᵀ/2
- No Clebsch-Gordan tables needed
- Automatic handling of degeneracies

---

## 10. Permutation Matrices for Antisymmetrization

### 10.1 Motivation

For identical fermions, the total wave function must be antisymmetric under particle exchange. To compute antisymmetrized matrix elements, we need the action of permutation operators on discrete (spin, color, isospin) states.

### 10.2 Permutation Action on Tensor Product States

A permutation $P$ acts on tensor product basis states by permuting particle labels:

$$
P |i_1, i_2, \ldots, i_N\rangle = |i_{P^{-1}(1)}, i_{P^{-1}(2)}, \ldots, i_{P^{-1}(N)}\rangle
$$

**Why inverse?** The permutation moves particle $k$'s state to position $P(k)$. Equivalently, position $j$ receives the state from particle $P^{-1}(j)$.

### 10.3 Building Permutation Matrices (Uniform Dimension)

For spin (all $d=2$) or color (all $d=3$):

```text
build_permutation_matrix_uniform(perm, N, d):
    # perm: permutation array, perm[i] = P(i)
    # N: number of particles
    # d: dimension per particle

    P_inv = inverse(perm)
    total_dim = d^N
    P_matrix = zeros(total_dim, total_dim)

    FOR old_index = 0 TO total_dim - 1:
        # Convert linear index to multi-index (base d)
        multi_idx = []
        temp = old_index
        FOR k = N downto 1:
            multi_idx[k] = temp MOD d
            temp = temp DIV d

        # Apply inverse permutation
        permuted_idx = [multi_idx[P_inv[k]] for k = 1 TO N]

        # Convert back to linear index
        new_index = 0
        FOR k = 1 TO N:
            new_index = new_index * d + permuted_idx[k]

        P_matrix[new_index, old_index] = 1

    RETURN P_matrix
```

### 10.4 Building Permutation Matrices (Mixed Dimensions)

For isospin where different quarks have different dimensions:

```text
build_permutation_matrix_mixed(perm, dims):
    # perm: permutation array
    # dims: array of dimensions for each particle

    N = length(dims)
    P_inv = inverse(perm)
    total_dim = product(dims)
    P_matrix = zeros(total_dim, total_dim)

    FOR old_index = 0 TO total_dim - 1:
        # Convert to multi-index with mixed radix
        multi_idx = linear_to_multi(old_index, dims)

        # Apply inverse permutation
        # Note: dims must also be permuted for correct indexing
        permuted_idx = [multi_idx[P_inv[k]] for k = 1 TO N]
        permuted_dims = [dims[P_inv[k]] for k = 1 TO N]

        # Convert back with permuted dimensions
        new_index = multi_to_linear(permuted_idx, permuted_dims)

        P_matrix[new_index, old_index] = 1

    RETURN P_matrix

linear_to_multi(index, dims):
    N = length(dims)
    multi = zeros(N)
    FOR k = N downto 1:
        multi[k] = index MOD dims[k]
        index = index DIV dims[k]
    RETURN multi

multi_to_linear(multi, dims):
    N = length(dims)
    index = 0
    FOR k = 1 TO N:
        index = index * dims[k] + multi[k]
    RETURN index
```

### 10.5 Matrix Elements with Permutation

For discrete states $|\chi_a\rangle$, $|\chi_b\rangle$ and operator $O$:

$$
\langle \chi_a | O | P \chi_b \rangle = \chi_a^\dagger \cdot O \cdot P_{\text{matrix}} \cdot \chi_b
$$

**Overlap with permutation:**

$$
\langle \chi_a | P \chi_b \rangle = \chi_a^\dagger \cdot P_{\text{matrix}} \cdot \chi_b
$$

**Operator with permutation:**

$$
\langle \chi_a | (\mathbf{s}_i \cdot \mathbf{s}_j) | P \chi_b \rangle = \chi_a^\dagger \cdot S_{ij} \cdot P_{\text{spin}} \cdot \chi_b
$$

### 10.6 Permutation Matrix Properties

**Orthogonality:**

$$
P^T P = P P^T = I
$$

**Composition:**

$$
P_{Q \circ P} = P_Q \cdot P_P
$$

**Inverse:**

$$
P^{-1} = P^T
$$

### 10.7 Example: 2-Particle Spin Permutation

For $N=2$, $d=2$ (spin), the basis is: $|00\rangle, |01\rangle, |10\rangle, |11\rangle$

Permutation $P_{12} = [2,1]$ with $P_{12}^{-1} = [2,1]$:

| Old state | Multi-idx | Permuted | New state |
|-----------|-----------|----------|-----------|
| $\|00\rangle$ | [0,0] | [0,0] | $\|00\rangle$ |
| $\|01\rangle$ | [0,1] | [1,0] | $\|10\rangle$ |
| $\|10\rangle$ | [1,0] | [0,1] | $\|01\rangle$ |
| $\|11\rangle$ | [1,1] | [1,1] | $\|11\rangle$ |

$$
P_{12}^{\text{spin}} = \begin{pmatrix} 1 & 0 & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix}
$$

### 10.8 Example: 3-Particle Color Permutation

For $N=3$, $d=3$ (color), basis dimension is $27$.

The cyclic permutation $P_{123} = [2,3,1]$ with $P_{123}^{-1} = [3,1,2]$:

$$
P_{123} |c_1, c_2, c_3\rangle = |c_3, c_1, c_2\rangle
$$

The color singlet $|\epsilon\rangle = \frac{1}{\sqrt{6}}(|rgb\rangle - |rbg\rangle + |gbr\rangle - |grb\rangle + |brg\rangle - |bgr\rangle)$ is antisymmetric:

$$
P_{12} |\epsilon\rangle = -|\epsilon\rangle, \quad P_{123} |\epsilon\rangle = |\epsilon\rangle
$$

### 10.9 Verification

**Symmetry eigenvalues:**

| State type | $P_{12}$ eigenvalue |
|------------|---------------------|
| Spin singlet (S=0) | $-1$ |
| Spin triplet (S=1) | $+1$ |
| Color singlet (qqq) | $-1$ |
| Color octet | mixed |

**Antisymmetrization check:**

For antisymmetric state $|\psi\rangle$:

$$
\sum_P \text{sign}(P) \cdot P |\psi\rangle = N! |\psi\rangle
$$

For references to the full antisymmetrization algorithm, see `algorithm_antisymmetrization.md`.
