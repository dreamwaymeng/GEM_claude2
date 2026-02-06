# Antisymmetrization via Permutation Operators

## 1. Overview

### 1.1 The Problem

Identical fermions require antisymmetric total wave functions. Traditional approaches use Young tableaux and matching symmetry patterns between spin, color, and spatial parts. Instead, we directly antisymmetrize the total wave function using permutation operators.

### 1.2 Key Insight

For the total wave function $\Psi = \Phi_{\text{spatial}} \otimes \chi_{\text{discrete}}$, antisymmetrize via:

$$
\langle \Psi_a | O | \mathcal{A}\Psi_b \rangle = \frac{1}{\mathcal{N}} \sum_{P} \text{sign}(P) \langle \Psi_a | O | P\Psi_b \rangle
$$

where $\mathcal{A}$ is the antisymmetrizer and the sum is over relevant permutations of identical particles.

Since $\Psi$ factors as a tensor product:

$$
\langle \Psi_a | O | P\Psi_b \rangle = \langle \Phi_a | O_{\text{spatial}} | P\Phi_b \rangle \times \langle \chi_a | O_{\text{discrete}} | P\chi_b \rangle
$$

### 1.3 Advantages

1. **No Young tableaux** — Direct numerical implementation
2. **No symmetry matching** — Automatic handling via sum over permutations
3. **Basis flexibility** — Spatial and discrete bases can be constructed independently
4. **Transparent** — Easy to verify and debug

---

## 2. Permutation Representation and Generation

### 2.1 Permutation as Array

A permutation $P$ of $N$ particles is represented as an array where `perm[i] = P(i)`:

```text
# Permutation array: perm[i] = where particle i goes
# Example for N=3:
identity = [1, 2, 3]       # P(1)=1, P(2)=2, P(3)=3
P_12     = [2, 1, 3]       # P(1)=2, P(2)=1, P(3)=3 (swap 1↔2)
P_13     = [3, 2, 1]       # swap 1↔3
P_23     = [1, 3, 2]       # swap 2↔3
P_123    = [2, 3, 1]       # cyclic: 1→2→3→1
P_132    = [3, 1, 2]       # cyclic: 1→3→2→1
```

### 2.2 Permutation Sign

The sign (parity) of a permutation is $(-1)^{n}$ where $n$ is the number of inversions:

```text
PROCEDURE permutation_sign(perm):
    # Count inversions: pairs (i,j) where i < j but perm[i] > perm[j]
    inversions = 0
    N = length(perm)
    FOR i = 1 TO N-1:
        FOR j = i+1 TO N:
            IF perm[i] > perm[j]:
                inversions += 1
    RETURN (-1)^inversions
```

**Examples:**
| Permutation | Array | Inversions | Sign |
|-------------|-------|------------|------|
| Identity | [1,2,3] | 0 | +1 |
| $P_{12}$ | [2,1,3] | 1 | -1 |
| $P_{123}$ | [2,3,1] | 2 | +1 |

### 2.3 Inverse Permutation

The inverse permutation $P^{-1}$ satisfies $P^{-1}(P(i)) = i$:

```text
PROCEDURE inverse_permutation(perm):
    N = length(perm)
    inv = zeros(N)
    FOR i = 1 TO N:
        inv[perm[i]] = i
    RETURN inv
```

**Example:** If $P = [2, 3, 1]$ (i.e., $P(1)=2, P(2)=3, P(3)=1$), then $P^{-1} = [3, 1, 2]$.

### 2.4 Permutation Composition

The composition $Q \circ P$ means "first apply $P$, then apply $Q$":

```text
PROCEDURE compose_permutations(Q, P):
    # (Q ∘ P)(i) = Q(P(i))
    N = length(P)
    result = zeros(N)
    FOR i = 1 TO N:
        result[i] = Q[P[i]]
    RETURN result
```

### 2.5 Generating All Permutations

For a set of particle indices, generate all permutations using Heap's algorithm or recursion:

```text
PROCEDURE generate_permutations(particles):
    # particles: array of particle indices, e.g., [1, 2, 3]
    # Returns list of all N! permutation arrays

    IF length(particles) == 1:
        RETURN [particles]

    result = []
    FOR i = 1 TO length(particles):
        # Put particles[i] first, permute the rest
        first = particles[i]
        rest = particles without element i
        FOR sub_perm in generate_permutations(rest):
            result.append([first] + sub_perm)

    RETURN result
```

---

## 3. Relevant Permutations (Optimization)

### 3.1 Identical Particle Groups

Only permute identical particles. Different quark flavors are distinguishable and need not be permuted.

**Identification:**
- Group particles by flavor: $\{u, d, s, c, b\}$ for quarks, $\{\bar{u}, \bar{d}, \bar{s}, \bar{c}, \bar{b}\}$ for antiquarks
- Only permute within each group

### 3.2 Examples

| System | Flavors | Identical Groups | Permutations |
|--------|---------|------------------|--------------|
| $\Lambda$ (uds) | u, d, s | None | 1 (identity only) |
| $\Sigma^+$ (uud) | u, u, d | {1,2} identical | 2 (id, $P_{12}$) |
| $\Delta^{++}$ (uuu) | u, u, u | {1,2,3} identical | 6 (full $S_3$) |
| $\Xi^0$ (uss) | u, s, s | {2,3} identical | 2 (id, $P_{23}$) |
| $\Omega^-$ (sss) | s, s, s | {1,2,3} identical | 6 (full $S_3$) |
| $T_{cc}$ (cc$\bar{u}\bar{d}$) | c, c, $\bar{u}$, $\bar{d}$ | {1,2} identical | 2 (id, $P_{12}$) |
| $X$ (cc$\bar{u}\bar{u}$) | c, c, $\bar{u}$, $\bar{u}$ | {1,2}, {3,4} | 4 = 2! × 2! |

### 3.3 Algorithm for Relevant Permutations

```text
PROCEDURE relevant_permutations(flavors):
    # flavors: array of quark flavors, e.g., ['u', 'u', 'd']
    # Returns list of permutations that permute only identical particles

    N = length(flavors)

    # Group particles by flavor
    groups = {}
    FOR i = 1 TO N:
        IF flavors[i] not in groups:
            groups[flavors[i]] = []
        groups[flavors[i]].append(i)

    # Generate permutations within each group
    group_perms = []
    FOR flavor, indices in groups:
        group_perms.append(generate_permutations(indices))

    # Cartesian product of group permutations
    all_perms = []
    FOR combo in cartesian_product(group_perms):
        # Build full permutation array
        perm = [0] * N
        FOR group_perm in combo:
            FOR i, target in enumerate(group_perm):
                perm[original_index[i]] = target
        all_perms.append(perm)

    RETURN all_perms
```

### 3.4 Normalization Factor

The normalization depends on the number of identical particles in each group:

$$
\mathcal{N} = \prod_g n_g!
$$

where $n_g$ is the number of particles in identical group $g$.

**Examples:**
- $\Lambda$ (uds): $\mathcal{N} = 1! \cdot 1! \cdot 1! = 1$
- $\Sigma^+$ (uud): $\mathcal{N} = 2! \cdot 1! = 2$
- $\Delta^{++}$ (uuu): $\mathcal{N} = 3! = 6$
- $T_{cc\bar{u}\bar{u}}$: $\mathcal{N} = 2! \cdot 2! = 4$

---

## 4. Permutation Action on Discrete States

### 4.1 Tensor Product Basis

Discrete wave functions (spin, color, isospin) live in tensor product spaces:

$$
|i_1, i_2, \ldots, i_N\rangle = |i_1\rangle \otimes |i_2\rangle \otimes \cdots \otimes |i_N\rangle
$$

where $i_k \in \{1, \ldots, d_k\}$ is the basis index for particle $k$ in dimension $d_k$.

**Basis ordering:** Lexicographic with multi-index $\to$ linear index:

$$
\text{index} = \sum_{k=1}^{N} (i_k - 1) \prod_{l=k+1}^{N} d_l
$$

### 4.2 Permutation Action

A permutation $P$ acts on tensor product states by permuting particle labels:

$$
P |i_1, i_2, \ldots, i_N\rangle = |i_{P^{-1}(1)}, i_{P^{-1}(2)}, \ldots, i_{P^{-1}(N)}\rangle
$$

**Why $P^{-1}$?** The permutation moves the state of particle $k$ to position $P(k)$. Equivalently, position $j$ receives the state from particle $P^{-1}(j)$.

**Example:** For $P_{12} = [2,1,3]$ with $P^{-1} = [2,1,3]$:
$$
P_{12} |i_1, i_2, i_3\rangle = |i_2, i_1, i_3\rangle
$$

### 4.3 Building Permutation Matrices

```text
PROCEDURE build_permutation_matrix_discrete(perm, dims):
    # perm: permutation array
    # dims: dimensions for each particle, e.g., [2,2,2] for 3-quark spin
    # Returns: permutation matrix P such that P @ |ket⟩ = |permuted_ket⟩

    N = length(perm)
    P_inv = inverse_permutation(perm)

    total_dim = product(dims)
    P_matrix = zeros(total_dim, total_dim)

    FOR old_index = 0 TO total_dim - 1:
        # Convert linear index to multi-index
        multi_idx = linear_to_multi(old_index, dims)

        # Apply permutation: new position j gets state from P^{-1}(j)
        permuted_idx = [multi_idx[P_inv[k]-1] for k = 1 TO N]

        # Convert back to linear index
        new_index = multi_to_linear(permuted_idx, dims)

        P_matrix[new_index, old_index] = 1

    RETURN P_matrix

PROCEDURE linear_to_multi(index, dims):
    # Convert linear index to multi-index
    N = length(dims)
    multi = zeros(N)
    FOR k = N downto 1:
        multi[k] = index MOD dims[k]
        index = index DIV dims[k]
    RETURN multi

PROCEDURE multi_to_linear(multi, dims):
    # Convert multi-index to linear index
    N = length(dims)
    index = 0
    FOR k = 1 TO N:
        index = index * dims[k] + multi[k]
    RETURN index
```

### 4.4 Uniform Dimension Case

For spin (all $d_k = 2$) or color (all $d_k = 3$):

```text
PROCEDURE build_permutation_matrix_uniform(perm, N, d):
    # perm: permutation array
    # N: number of particles
    # d: dimension per particle (2 for spin, 3 for color)

    dims = [d] * N
    RETURN build_permutation_matrix_discrete(perm, dims)
```

### 4.5 Matrix Elements with Permutation

For discrete states $|\chi_a\rangle$, $|\chi_b\rangle$ and operator $O_{\text{discrete}}$:

$$
\langle \chi_a | O_{\text{discrete}} | P \chi_b \rangle = \chi_a^\dagger \cdot O_{\text{discrete}} \cdot P_{\text{matrix}} \cdot \chi_b
$$

For the overlap (identity operator):

$$
\langle \chi_a | P \chi_b \rangle = \chi_a^\dagger \cdot P_{\text{matrix}} \cdot \chi_b
$$

### 4.6 Examples

**2-particle spin with $P_{12}$:**

Basis: $|{\uparrow}{\uparrow}\rangle, |{\uparrow}{\downarrow}\rangle, |{\downarrow}{\uparrow}\rangle, |{\downarrow}{\downarrow}\rangle$

$P_{12}$ swaps the two spins:

$$
P_{12} = \begin{pmatrix} 1 & 0 & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix}
$$

**Spin singlet:** $|\chi_0\rangle = (|{\uparrow}{\downarrow}\rangle - |{\downarrow}{\uparrow}\rangle)/\sqrt{2} = (0, 1, -1, 0)^T/\sqrt{2}$

$$
P_{12} |\chi_0\rangle = (0, -1, 1, 0)^T/\sqrt{2} = -|\chi_0\rangle
$$

The singlet is antisymmetric under $P_{12}$ (eigenvalue $-1$).

**Spin triplet:** $|\chi_1\rangle = (|{\uparrow}{\downarrow}\rangle + |{\downarrow}{\uparrow}\rangle)/\sqrt{2} = (0, 1, 1, 0)^T/\sqrt{2}$

$$
P_{12} |\chi_1\rangle = (0, 1, 1, 0)^T/\sqrt{2} = +|\chi_1\rangle
$$

The triplet is symmetric under $P_{12}$ (eigenvalue $+1$).

---

## 5. Permutation Action on Spatial States

### 5.1 Permutation as Jacobi Coordinate Transformation

A permutation $P$ of particle labels induces a transformation of Jacobi coordinates. If the original Jacobi tree is $\mathcal{T}$, applying $P$ gives a permuted tree $P \circ \mathcal{T}$.

**Key insight:** The permutation is equivalent to a coordinate transformation between Jacobi sets.

### 5.2 Permuting a Jacobi Tree

```text
PROCEDURE apply_perm_to_tree(tree, perm):
    # tree: Jacobi clustering tree, e.g., ((1,2),3)
    # perm: permutation array
    # Returns: tree with permuted particle labels

    IF tree is a leaf (integer):
        RETURN perm[tree]  # Relabel particle

    (left, right) = tree
    new_left = apply_perm_to_tree(left, perm)
    new_right = apply_perm_to_tree(right, perm)
    RETURN (new_left, new_right)
```

**Example:** Tree $((1,2),3)$ with $P_{12} = [2,1,3]$:

$$
P_{12} \circ ((1,2),3) = ((2,1),3)
$$

This represents the same Jacobi set but with particles 1 and 2 swapped internally.

### 5.3 Transformation Matrix for Permuted Coordinates

The matrix element $\langle \Phi_a | O | P\Phi_b \rangle$ requires expressing the permuted ket in the bra's coordinate system.

```text
PROCEDURE spatial_ME_with_permutation(bra, ket, operator, perm, masses):
    # bra: (tree, widths) for bra Gaussian
    # ket: (tree, widths) for ket Gaussian
    # operator: potential type and pair
    # perm: permutation array
    # masses: particle masses

    # Step 1: Apply permutation to ket's Jacobi tree
    permuted_tree = apply_perm_to_tree(ket.tree, perm)

    # Step 2: Get transformation matrix from permuted tree to bra's tree
    #         (or to target Jacobi set for the operator)
    T_perm = build_transformation_matrix(permuted_tree, target_tree, masses)

    # Step 3: Transform ket's width matrix
    # Original: exp(-ρ^T D_ket ρ) where D_ket = diag(ket.widths)
    # After permutation: coordinates change, width matrix transforms
    D_ket = diag(ket.widths)

    # The permuted Gaussian in target coordinates:
    T_inv = inverse(T_perm)
    W_ket_permuted = T_inv^T @ D_ket @ T_inv

    # Step 4: Use existing Cholesky algorithm with permuted ket
    # (The bra transformation proceeds as usual)
    RETURN cholesky_matrix_element(bra, W_ket_permuted, operator)
```

### 5.4 Special Case: $P_{12}$ in Jacobi Set $(12)3$

For the Jacobi set $(12)3$ with coordinates $(\boldsymbol{\rho}, \boldsymbol{\lambda})$:

$$
\boldsymbol{\rho} = \mathbf{r}_1 - \mathbf{r}_2, \quad \boldsymbol{\lambda} = \frac{m_1 \mathbf{r}_1 + m_2 \mathbf{r}_2}{m_1 + m_2} - \mathbf{r}_3
$$

Under $P_{12}$ (swap particles 1 and 2):

$$
P_{12}: \boldsymbol{\rho} \to -\boldsymbol{\rho}, \quad \boldsymbol{\lambda} \to \boldsymbol{\lambda}
$$

**Transformation matrix:**

$$
\mathbf{T}_{P_{12}} = \begin{pmatrix} -1 & 0 \\ 0 & 1 \end{pmatrix}
$$

**For s-wave Gaussians:** $\exp(-\nu \rho^2)$ is unchanged by $\boldsymbol{\rho} \to -\boldsymbol{\rho}$.

Thus, for s-wave Gaussians in a Jacobi set where the permuted pair is diagonal, $P_{12}$ leaves the spatial matrix element unchanged:

$$
\langle \Phi_a | O | P_{12}\Phi_b \rangle = \langle \Phi_a | O | \Phi_b \rangle
$$

### 5.5 General Permutation Effect

For arbitrary permutations and Jacobi sets, the effect is a linear coordinate transformation. The algorithm in Section 5.3 handles all cases uniformly.

**Implementation note:** The transformation matrix between permuted coordinates is computed using the same tree-based method as for different Jacobi sets (see `algorithm_spatial_MM.md`).

---

## 6. Combined Antisymmetrized Matrix Element

### 6.1 Master Formula

For operator $O = O_{\text{spatial}} \otimes O_{\text{discrete}}$:

$$
\boxed{\langle \Psi_a | O | \mathcal{A}\Psi_b \rangle = \frac{1}{\mathcal{N}} \sum_{P} \text{sign}(P) \cdot \langle \Phi_a | O_{\text{spatial}} | P\Phi_b \rangle \cdot \langle \chi_a | O_{\text{discrete}} | P\chi_b \rangle}
$$

where:
- $\mathcal{N} = \prod_g n_g!$ (product of factorials of identical group sizes)
- Sum over all permutations of identical particles
- $\Phi$ = spatial wave function
- $\chi$ = discrete wave function (spin $\otimes$ color $\otimes$ isospin)

### 6.2 Algorithm

```text
PROCEDURE antisymmetrized_ME(
    bra_spatial, bra_discrete,
    ket_spatial, ket_discrete,
    operator_spatial, operator_discrete,
    flavors, masses
):
    # bra_spatial, ket_spatial: (tree, widths) for Gaussians
    # bra_discrete, ket_discrete: coefficient vectors
    # operator_spatial: potential type and pair
    # operator_discrete: spin/color operator matrix
    # flavors: quark flavor array
    # masses: particle masses

    # Get relevant permutations
    perms = relevant_permutations(flavors)

    # Compute normalization
    norm = 1.0
    FOR group in identical_groups(flavors):
        norm *= factorial(length(group))

    # Sum over permutations
    total = 0
    FOR perm in perms:
        sign = permutation_sign(perm)

        # Spatial matrix element with permutation
        spatial_ME = spatial_ME_with_permutation(
            bra_spatial, ket_spatial, operator_spatial, perm, masses
        )

        # Discrete matrix element with permutation
        P_discrete = build_permutation_matrix_discrete(perm, discrete_dims)
        discrete_ME = bra_discrete^† @ operator_discrete @ P_discrete @ ket_discrete

        total += sign * spatial_ME * discrete_ME

    RETURN total / norm
```

### 6.3 Overlap and Kinetic Energy

**Overlap:** $O_{\text{spatial}} = 1$, $O_{\text{discrete}} = 1$

$$
S_{ab} = \frac{1}{\mathcal{N}} \sum_P \text{sign}(P) \cdot \langle \Phi_a | P\Phi_b \rangle \cdot \langle \chi_a | P\chi_b \rangle
$$

**Kinetic energy:** $O_{\text{spatial}} = T$, $O_{\text{discrete}} = 1$

$$
T_{ab} = \frac{1}{\mathcal{N}} \sum_P \text{sign}(P) \cdot \langle \Phi_a | T | P\Phi_b \rangle \cdot \langle \chi_a | P\chi_b \rangle
$$

### 6.4 Potential Energy

For pair potential $V_{ij}$:

$$
V_{ab} = \frac{1}{\mathcal{N}} \sum_P \text{sign}(P) \cdot \langle \Phi_a | V(r_{ij}) | P\Phi_b \rangle \cdot \langle \chi_a | F_{ij} | P\chi_b \rangle
$$

where $F_{ij}$ is the discrete factor (e.g., $\mathbf{s}_i \cdot \mathbf{s}_j$ or $\boldsymbol{\lambda}_i \cdot \boldsymbol{\lambda}_j$).

---

## 7. Full Hamiltonian Construction

### 7.1 Basis Structure

The full basis is a product of spatial and discrete bases:

$$
|\Psi_{i_s, i_d}\rangle = |\Phi_{i_s}\rangle \otimes |\chi_{i_d}\rangle
$$

where:
- $i_s = 1, \ldots, N_{\text{spatial}}$ labels spatial basis functions
- $i_d = 1, \ldots, N_{\text{discrete}}$ labels discrete states

### 7.2 Algorithm

```text
PROCEDURE build_antisymmetrized_hamiltonian(
    spatial_basis,    # List of (tree, widths)
    discrete_basis,   # List of coefficient vectors
    flavors,          # Quark flavors
    masses,           # Particle masses
    potential_params  # Potential parameters
):
    N_s = length(spatial_basis)
    N_d = length(discrete_basis)
    N_total = N_s * N_d

    # Get relevant permutations (compute once)
    perms = relevant_permutations(flavors)
    norm = product(factorial(group_size) for group in identical_groups(flavors))

    # Precompute discrete permutation matrices (independent of widths)
    N_particles = length(flavors)
    P_spin = {}   # permutation -> spin permutation matrix
    P_color = {}  # permutation -> color permutation matrix

    FOR perm in perms:
        P_spin[perm] = build_permutation_matrix_uniform(perm, N_particles, d=2)
        P_color[perm] = build_permutation_matrix_uniform(perm, N_particles, d=3)

    # Precompute discrete operator matrices
    FOR pair (i,j) in all_pairs:
        spin_spin[i,j] = build_spin_spin_operator(i, j, N_particles)
        color_color[i,j] = build_color_color_operator(i, j, N_particles, particle_types)

    # Initialize Hamiltonian and overlap matrices
    H = zeros(N_total, N_total)
    S = zeros(N_total, N_total)

    # Main loop over basis pairs
    FOR i_s = 1 TO N_s:
        FOR i_d = 1 TO N_d:
            i = (i_s - 1) * N_d + i_d  # Combined index

            FOR j_s = 1 TO N_s:
                FOR j_d = 1 TO N_d:
                    j = (j_s - 1) * N_d + j_d

                    # Skip if j < i (use Hermiticity)
                    IF j < i:
                        CONTINUE

                    h_elem = 0
                    s_elem = 0

                    # Sum over permutations
                    FOR perm in perms:
                        sign = permutation_sign(perm)

                        # Discrete overlaps with permutation
                        spin_overlap = discrete_basis[i_d].spin^† @
                                       P_spin[perm] @ discrete_basis[j_d].spin
                        color_overlap = discrete_basis[i_d].color^† @
                                        P_color[perm] @ discrete_basis[j_d].color
                        discrete_overlap = spin_overlap * color_overlap

                        # Spatial overlap with permutation
                        spatial_overlap = spatial_overlap_with_perm(
                            spatial_basis[i_s], spatial_basis[j_s], perm, masses
                        )

                        # Total overlap contribution
                        s_elem += sign * spatial_overlap * discrete_overlap

                        # Kinetic energy (spatial only)
                        spatial_kinetic = kinetic_ME_with_perm(
                            spatial_basis[i_s], spatial_basis[j_s], perm, masses
                        )
                        h_elem += sign * spatial_kinetic * discrete_overlap

                        # Potential energy (sum over pairs)
                        FOR pair (k,l) in all_pairs:
                            # Spatial potential
                            spatial_V = potential_ME_with_perm(
                                spatial_basis[i_s], spatial_basis[j_s],
                                pair, potential_params, perm, masses
                            )

                            # Discrete factor for this pair
                            spin_factor = discrete_basis[i_d].spin^† @
                                          spin_spin[k,l] @ P_spin[perm] @
                                          discrete_basis[j_d].spin
                            color_factor = discrete_basis[i_d].color^† @
                                           color_color[k,l] @ P_color[perm] @
                                           discrete_basis[j_d].color

                            # Add contribution
                            h_elem += sign * spatial_V * spin_factor * color_factor * V_coeff

                    # Apply normalization and store
                    H[i,j] = h_elem / norm
                    S[i,j] = s_elem / norm

                    # Hermiticity
                    IF i != j:
                        H[j,i] = conj(H[i,j])
                        S[j,i] = conj(S[i,j])

    RETURN H, S
```

### 7.3 Optimization Notes

1. **Precompute permutation matrices** — Independent of Gaussian widths
2. **Precompute operator matrices** — $\mathbf{s}_i \cdot \mathbf{s}_j$, $\boldsymbol{\lambda}_i \cdot \boldsymbol{\lambda}_j$
3. **Cache transformation matrices** — Between Jacobi sets
4. **Exploit Hermiticity** — Only compute upper triangle
5. **Sparse structure** — Many discrete matrix elements are zero

---

## 8. Verification Test Cases

### 8.1 Meson (qq̄): No Antisymmetrization

Quark and antiquark are distinguishable. Only identity permutation.

**Check:** Results match standard (non-antisymmetrized) calculation.

### 8.2 Proton (uud): $P_{12}$ Only

Two identical u quarks → permutations: {identity, $P_{12}$}

**Spin-color structure:**
- Color singlet is antisymmetric in any qq pair
- Total wave function must be antisymmetric under $P_{12}$
- Spatial-spin combination must have definite symmetry under $P_{12}$

**Verification:**
1. Compute matrix elements with and without antisymmetrization
2. For properly constructed basis, antisymmetrized and non-antisymmetrized results should match (basis already has correct symmetry)
3. Cross-check: $\langle \chi | P_{12} | \chi \rangle$ gives symmetry eigenvalue

### 8.3 $\Delta^{++}$ (uuu): Full $S_3$

All three quarks identical → 6 permutations.

**Known properties:**
- Color singlet: fully antisymmetric (eigenvalue $-1$ under all transpositions)
- $\Delta^{++}$ has $S = 3/2$: spin fully symmetric (eigenvalue $+1$ under all transpositions)
- Therefore spatial must be fully symmetric

**Verification:**
1. All pairs equivalent after antisymmetrization
2. $\langle \mathbf{s}_1 \cdot \mathbf{s}_2 \rangle = \langle \mathbf{s}_1 \cdot \mathbf{s}_3 \rangle = \langle \mathbf{s}_2 \cdot \mathbf{s}_3 \rangle$

### 8.4 $\Lambda$ (uds): Identity Only

All three quarks distinguishable → only identity permutation.

**Check:** Antisymmetrization reduces to standard calculation (factor of 1).

### 8.5 Consistency Checks

**Hermiticity:**

$$
H_{ij} = H_{ji}^*, \quad S_{ij} = S_{ji}^*
$$

**Overlap normalization:**

For properly normalized, antisymmetrized states:

$$
\langle \Psi | \mathcal{A} | \Psi \rangle = \langle \Psi | \Psi \rangle \quad \text{(if } \Psi \text{ is already antisymmetric)}
$$

**Sum rules:**

Total spin-spin sum:

$$
\sum_{i<j} \langle \mathbf{s}_i \cdot \mathbf{s}_j \rangle = \frac{1}{2}[S(S+1) - \frac{3N}{4}]
$$

Total color-color sum (for color singlet):

$$
\sum_{i<j} \langle \boldsymbol{\lambda}_i \cdot \boldsymbol{\lambda}_j \rangle = -\frac{4N}{3} \quad \text{(for color singlet)}
$$

---

## 9. Implementation Order

### Phase 1: Permutation Utilities
1. Permutation representation (array format)
2. Sign computation (inversion count)
3. Inverse permutation
4. Permutation composition
5. Generation of all permutations
6. Relevant permutations for given flavors

### Phase 2: Discrete Permutation Matrices
1. Multi-index ↔ linear index conversion
2. Build permutation matrix for uniform dimension
3. Build permutation matrix for mixed dimensions (isospin)
4. Matrix elements with permutation

### Phase 3: Spatial Jacobi Transformation
1. Apply permutation to Jacobi tree
2. Transformation matrix for permuted coordinates
3. Special cases (P₁₂ in set (12)3, etc.)
4. Spatial matrix elements with permutation

### Phase 4: Combined Computation
1. Antisymmetrized overlap matrix element
2. Antisymmetrized kinetic energy
3. Antisymmetrized potential energy
4. Full Hamiltonian assembly

### Phase 5: Verification
1. Implement test cases (Sections 8.1-8.4)
2. Consistency checks
3. Comparison with known results

---

## 10. References

- This algorithm avoids Young tableaux by explicit summation over permutations
- For spatial matrix elements, see `algorithm_spatial_MM.md`
- For discrete wave function construction, see `algorithm_discrete_MM.md`
- For physics background, see `physics.md`
