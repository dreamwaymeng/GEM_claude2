# Spatial Matrix Elements via Cholesky Decomposition

This document describes the algorithm for computing spatial matrix elements of potentials in N-body systems using Cholesky decomposition. The method enables fully analytical calculation of matrix elements for any pair interaction, including:
- Non-diagonal pairs that depend on multiple Jacobi coordinates
- **Basis functions defined in different Jacobi coordinate sets**

---

## 1. Overview

### 1.1 The Problem

In N-body variational calculations, the wave function is expanded using Gaussian basis functions:

$$
\Psi = \sum_{\alpha} \sum_{n} c_{\alpha,n} \Phi^{(\alpha)}_n
$$

where $\alpha$ labels the **Jacobi coordinate set** and $n$ labels the basis function within that set.

**Key insight:** To achieve good variational convergence, the basis should include Gaussians defined in **multiple Jacobi sets**, not just one. This leads to matrix elements between basis functions from different Jacobi sets:

$$
\langle \Phi^{(\alpha)} | \hat{O} | \Phi^{(\beta)} \rangle
$$

where $\alpha$ and $\beta$ may be different Jacobi sets.

### 1.2 Challenges

1. **Non-diagonal pairs:** For a given Jacobi set, only one pair is "diagonal" (proportional to a single coordinate). Other pairs depend on multiple coordinates.

2. **Cross-set matrix elements:** When bra and ket are in different Jacobi sets, the combined Gaussian has a **non-diagonal width matrix** even before considering non-diagonal pairs.

**Example (3-body with Jacobi set (12)3):**
- $r_{12} = |\boldsymbol{\rho}|$ — diagonal (trivial)
- $r_{13} = |a\boldsymbol{\rho} + \boldsymbol{\lambda}|$ — non-diagonal (coupled)
- $r_{23} = |-b\boldsymbol{\rho} + \boldsymbol{\lambda}|$ — non-diagonal (coupled)

### 1.3 The Solution

The Cholesky decomposition method handles all cases through:

1. **Width Matrix Construction**: Combine Gaussian widths, transforming to a common coordinate system if needed
2. **Coordinate Transformation**: Transform to the Jacobi set where the target pair is diagonal
3. **Cholesky Factorization**: Decompose the (generally non-diagonal) width matrix to decouple integrals
4. **Analytical Integration**: Compute factorized integrals using standard formulas

---

## 2. Jacobi Coordinate Transformations

### 2.1 Jacobi Coordinate Definition

For an N-body system with particles of masses $m_1, m_2, \ldots, m_N$, a Jacobi coordinate set is defined by a binary clustering tree. Each set has $N-1$ relative coordinates $\{\boldsymbol{\rho}_1, \boldsymbol{\rho}_2, \ldots, \boldsymbol{\rho}_{N-1}\}$.

**3-body example (Jacobi set labeling (ij)k):**

The Jacobi set $(ij)k$ has coordinates:
$$
\boldsymbol{\rho} = \mathbf{r}_i - \mathbf{r}_j, \quad \boldsymbol{\lambda} = \frac{m_i \mathbf{r}_i + m_j \mathbf{r}_j}{m_i + m_j} - \mathbf{r}_k
$$

with reduced masses:
$$
\mu_\rho = \frac{m_i m_j}{m_i + m_j}, \quad \mu_\lambda = \frac{(m_i + m_j) m_k}{M}
$$

where $M = m_i + m_j + m_k$ is the total mass.

### 2.2 Transformation Between Jacobi Sets

The transformation between any two Jacobi sets is a linear transformation:

$$
\boldsymbol{\rho}^{(\beta)} = \mathbf{T}_{\alpha \to \beta} \, \boldsymbol{\rho}^{(\alpha)}
$$

where $\mathbf{T}_{\alpha \to \beta}$ is an $(N-1) \times (N-1)$ matrix that can be derived from the clustering tree structures (see Section 2.4).

**Key properties:**
- $\mathbf{T}_{\beta \to \alpha} = \mathbf{T}_{\alpha \to \beta}^{-1}$
- $\mathbf{T}_{\alpha \to \gamma} = \mathbf{T}_{\beta \to \gamma} \cdot \mathbf{T}_{\alpha \to \beta}$

### 2.3 Clustering Trees: Representation and Labeling

A Jacobi coordinate set is uniquely defined by a **binary clustering tree** that specifies how particles are grouped hierarchically.

#### 2.3.1 Tree Structure

Each clustering tree has:

- **Leaves**: Individual particles (labeled 1, 2, ..., N)
- **Internal nodes**: Clusters of particles
- **Root**: The full system (all N particles)

Each internal node represents one Jacobi coordinate: the relative position between its two children.

#### 2.3.2 Tree Notation

We use nested tuple notation: `((a, b), c)` means:
1. First cluster particles `a` and `b` → defines coordinate $\boldsymbol{\rho}_1 = \mathbf{r}_a - \mathbf{r}_b$
2. Then combine cluster `(a,b)` with particle `c` → defines coordinate $\boldsymbol{\rho}_2 = \text{CM}(a,b) - \mathbf{r}_c$

**3-body examples:**

| Tree | Jacobi Coordinates | Shorthand |
|------|-------------------|-----------|
| `((1,2),3)` | $\boldsymbol{\rho} = \mathbf{r}_1 - \mathbf{r}_2$, $\boldsymbol{\lambda} = \text{CM}_{12} - \mathbf{r}_3$ | (12)3 |
| `((1,3),2)` | $\boldsymbol{\rho} = \mathbf{r}_1 - \mathbf{r}_3$, $\boldsymbol{\lambda} = \text{CM}_{13} - \mathbf{r}_2$ | (13)2 |
| `((2,3),1)` | $\boldsymbol{\rho} = \mathbf{r}_2 - \mathbf{r}_3$, $\boldsymbol{\lambda} = \text{CM}_{23} - \mathbf{r}_1$ | (23)1 |

**4-body examples:**

| Tree | Structure | Jacobi Coordinates |
|------|-----------|-------------------|
| `(((1,2),3),4)` | Chain | $\boldsymbol{\rho}_1 = \mathbf{r}_1 - \mathbf{r}_2$, $\boldsymbol{\rho}_2 = \text{CM}_{12} - \mathbf{r}_3$, $\boldsymbol{\rho}_3 = \text{CM}_{123} - \mathbf{r}_4$ |
| `((1,2),(3,4))` | H-type | $\boldsymbol{\rho}_1 = \mathbf{r}_1 - \mathbf{r}_2$, $\boldsymbol{\rho}_2 = \mathbf{r}_3 - \mathbf{r}_4$, $\boldsymbol{\rho}_3 = \text{CM}_{12} - \text{CM}_{34}$ |
| `(((1,3),2),4)` | Chain | $\boldsymbol{\rho}_1 = \mathbf{r}_1 - \mathbf{r}_3$, $\boldsymbol{\rho}_2 = \text{CM}_{13} - \mathbf{r}_2$, $\boldsymbol{\rho}_3 = \text{CM}_{123} - \mathbf{r}_4$ |

#### 2.3.3 Data Structure for Clustering Trees

```text
# A cluster is either a leaf (particle index) or a pair of clusters
Cluster = int | (Cluster, Cluster)

# Examples:
tree_12_3 = ((1, 2), 3)      # 3-body: (12)3
tree_13_2 = ((1, 3), 2)      # 3-body: (13)2
tree_H = ((1, 2), (3, 4))    # 4-body: H-type
tree_chain = (((1, 2), 3), 4) # 4-body: chain
```

#### 2.3.4 Computing Jacobi Coordinates from a Tree

```text
PROCEDURE jacobi_from_tree(tree, masses, positions)
    # Returns list of (jacobi_coord, reduced_mass) pairs

    IF tree is a leaf (integer):
        RETURN (positions[tree], masses[tree], {tree})  # position, mass, particle set

    (left, right) = tree
    (pos_L, mass_L, particles_L) = jacobi_from_tree(left, masses, positions)
    (pos_R, mass_R, particles_R) = jacobi_from_tree(right, masses, positions)

    # Jacobi coordinate: relative position
    rho = pos_L - pos_R

    # Reduced mass for this coordinate
    mu = mass_L * mass_R / (mass_L + mass_R)

    # Center of mass of combined cluster
    cm = (mass_L * pos_L + mass_R * pos_R) / (mass_L + mass_R)
    total_mass = mass_L + mass_R
    combined_particles = particles_L ∪ particles_R

    # Store this Jacobi coordinate
    RECORD jacobi_coord(rho, mu)

    RETURN (cm, total_mass, combined_particles)
```

### 2.4 Deriving Transformation Matrices from Trees

Given two clustering trees, we can systematically derive the transformation matrix between them.

#### 2.4.1 The General Approach

The transformation between any two Jacobi sets is a linear transformation:
$$
\boldsymbol{\rho}^{(\beta)} = \mathbf{T}_{\alpha \to \beta} \boldsymbol{\rho}^{(\alpha)}
$$

To find $\mathbf{T}$:
1. Express particle positions in terms of Jacobi coordinates for set α
2. Express Jacobi coordinates for set β in terms of particle positions
3. Substitute to get β coordinates in terms of α coordinates

#### 2.4.2 Expressing Particle Positions in Jacobi Coordinates

For a clustering tree, we can recursively express each particle position as:
$$
\mathbf{r}_i = \mathbf{R}_{\text{CM}} + \sum_{k=1}^{N-1} c_{ik} \boldsymbol{\rho}_k
$$

where $c_{ik}$ are coefficients determined by the tree structure and masses.

**Algorithm to find coefficients:**

```text
PROCEDURE particle_coefficients(tree, masses, particle_index)
    # Returns coefficients c_k such that r_i = R_CM + sum_k c_k * rho_k

    coeffs = zeros(N-1)  # N-1 Jacobi coordinates
    coord_index = 0

    FUNCTION traverse(node, sign, mass_factor):
        NONLOCAL coord_index, coeffs

        IF node is a leaf:
            IF node == particle_index:
                RETURN True  # Found the particle
            ELSE:
                RETURN False

        (left, right) = node
        mass_L = total_mass(left)
        mass_R = total_mass(right)
        total = mass_L + mass_R

        found_in_left = traverse(left, sign, mass_factor * mass_R / total)
        found_in_right = traverse(right, sign, mass_factor * (-mass_L) / total)

        IF found_in_left:
            coeffs[coord_index] += sign * mass_factor * mass_R / total
        IF found_in_right:
            coeffs[coord_index] += sign * mass_factor * (-mass_L) / total

        coord_index += 1
        RETURN found_in_left OR found_in_right

    traverse(tree, 1.0, 1.0)
    RETURN coeffs
```

#### 2.4.3 Building the Transformation Matrix

```text
PROCEDURE build_transformation_matrix(tree_from, tree_to, masses)
    # Returns T such that rho_to = T @ rho_from

    N = number_of_particles(tree_from)

    # Step 1: For tree_from, find how each particle position depends on Jacobi coords
    # r_i = R_CM + sum_k C_from[i,k] * rho_from[k]
    C_from = matrix(N, N-1)
    FOR i = 1 TO N:
        C_from[i, :] = particle_coefficients(tree_from, masses, i)

    # Step 2: For tree_to, find how each Jacobi coord depends on particle positions
    # rho_to[j] = sum_i A_to[j,i] * r_i  (where sum_i A_to[j,i] = 0 for CM invariance)
    A_to = matrix(N-1, N)
    FOR j = 1 TO N-1:
        A_to[j, :] = jacobi_coord_coefficients(tree_to, masses, j)

    # Step 3: Combine: rho_to = A_to @ r = A_to @ (R_CM + C_from @ rho_from)
    #                        = A_to @ C_from @ rho_from  (since A_to @ 1 = 0)
    T = A_to @ C_from

    RETURN T
```

#### 2.4.4 Example: 3-Body Transformation Matrix Derivation

**From tree `((1,2),3)` to tree `((1,3),2)`:**

**Step 1: Particle positions in set (12)3**

$$
\mathbf{r}_1 = \mathbf{R}_{\text{CM}} + \frac{m_2}{M_{12}} \boldsymbol{\rho} + \frac{m_3}{M} \boldsymbol{\lambda}
$$
$$
\mathbf{r}_2 = \mathbf{R}_{\text{CM}} - \frac{m_1}{M_{12}} \boldsymbol{\rho} + \frac{m_3}{M} \boldsymbol{\lambda}
$$
$$
\mathbf{r}_3 = \mathbf{R}_{\text{CM}} - \frac{M_{12}}{M} \boldsymbol{\lambda}
$$

So the coefficient matrix is:
$$
\mathbf{C}_{(12)3} = \begin{pmatrix} m_2/M_{12} & m_3/M \\ -m_1/M_{12} & m_3/M \\ 0 & -M_{12}/M \end{pmatrix}
$$

**Step 2: Jacobi coordinates in set (13)2 from particle positions**

$$
\boldsymbol{\rho}' = \mathbf{r}_1 - \mathbf{r}_3
$$
$$
\boldsymbol{\lambda}' = \frac{m_1 \mathbf{r}_1 + m_3 \mathbf{r}_3}{M_{13}} - \mathbf{r}_2
$$

So the coefficient matrix is:
$$
\mathbf{A}_{(13)2} = \begin{pmatrix} 1 & 0 & -1 \\ m_1/M_{13} & -1 & m_3/M_{13} \end{pmatrix}
$$

**Step 3: Combine**

$$
\mathbf{T}_{(12)3 \to (13)2} = \mathbf{A}_{(13)2} \cdot \mathbf{C}_{(12)3}
$$

Computing:
$$
T_{11} = 1 \cdot \frac{m_2}{M_{12}} + 0 \cdot (-\frac{m_1}{M_{12}}) + (-1) \cdot 0 = \frac{m_2}{M_{12}}
$$
$$
T_{12} = 1 \cdot \frac{m_3}{M} + 0 \cdot \frac{m_3}{M} + (-1) \cdot (-\frac{M_{12}}{M}) = \frac{m_3 + M_{12}}{M} = 1
$$
$$
T_{21} = \frac{m_1}{M_{13}} \cdot \frac{m_2}{M_{12}} + (-1) \cdot (-\frac{m_1}{M_{12}}) + \frac{m_3}{M_{13}} \cdot 0 = \frac{m_1 m_2}{M_{12} M_{13}} + \frac{m_1}{M_{12}} = \frac{m_1 M}{M_{12} M_{13}}
$$
$$
T_{22} = \frac{m_1}{M_{13}} \cdot \frac{m_3}{M} + (-1) \cdot \frac{m_3}{M} + \frac{m_3}{M_{13}} \cdot (-\frac{M_{12}}{M}) = -\frac{m_3}{M_{13}}
$$

**Result:**
$$
\mathbf{T}_{(12)3 \to (13)2} = \begin{pmatrix} \frac{m_2}{M_{12}} & 1 \\ \frac{m_1 M}{M_{12} M_{13}} & -\frac{m_3}{M_{13}} \end{pmatrix}
$$

This gives the transformation matrix directly from the tree structures.

#### 2.4.5 Practical Implementation

```text
PROCEDURE get_transformation_matrix(tree_alpha, tree_beta, masses)
    # High-level interface

    IF tree_alpha == tree_beta:
        RETURN identity_matrix(len(masses) - 1)

    # Check if we have this cached
    key = (tree_alpha, tree_beta)
    IF key in T_cache:
        RETURN T_cache[key]

    # Compute from scratch
    T = build_transformation_matrix(tree_alpha, tree_beta, masses)

    # Cache both directions
    T_cache[key] = T
    T_cache[(tree_beta, tree_alpha)] = inverse(T)

    RETURN T
```

### 2.5 Identifying Which Pair is Diagonal

Given a clustering tree, we can determine which particle pair is "diagonal" (i.e., the first Jacobi coordinate).

```text
PROCEDURE diagonal_pair(tree)
    # Returns the particle pair (i,j) that is diagonal in this Jacobi set

    IF tree is a leaf:
        ERROR "Single particle has no pairs"

    (left, right) = tree

    # The first Jacobi coordinate connects the leftmost leaves of left and right subtrees
    i = leftmost_leaf(left)
    j = leftmost_leaf(right)

    RETURN (min(i,j), max(i,j))

PROCEDURE leftmost_leaf(tree)
    IF tree is a leaf:
        RETURN tree
    (left, right) = tree
    RETURN leftmost_leaf(left)
```

**Examples:**

| Tree | Diagonal Pair |
|------|---------------|
| `((1,2),3)` | (1,2) |
| `((1,3),2)` | (1,3) |
| `((2,3),1)` | (2,3) |
| `(((1,2),3),4)` | (1,2) |
| `((1,2),(3,4))` | (1,2) for $\rho_1$, (3,4) for $\rho_2$ |

### 2.6 Finding the Jacobi Set for a Given Pair

Given a particle pair $(i,j)$, find a clustering tree where this pair is diagonal:

```text
PROCEDURE jacobi_set_for_pair(pair, N)
    # Returns a tree where pair (i,j) is the first Jacobi coordinate

    (i, j) = pair
    remaining = {1, 2, ..., N} - {i, j}

    # Start with the pair
    tree = (i, j)

    # Add remaining particles one by one
    FOR k in sorted(remaining):
        tree = (tree, k)

    RETURN tree
```

**Example:** For pair (1,3) in a 4-body system:
- Start: `(1, 3)`
- Add 2: `((1, 3), 2)`
- Add 4: `(((1, 3), 2), 4)`

This tree has $r_{13}$ as the first (diagonal) Jacobi coordinate.

---

## 3. General Cholesky Algorithm

### 3.1 Basis Function Definition

A Gaussian basis function in Jacobi set $\alpha$ is:
$$
\Phi^{(\alpha)}(\boldsymbol{\rho}^{(\alpha)}) = \mathcal{N}^{(\alpha)} \exp\left( -\sum_{k=1}^{N-1} \nu_k^{(\alpha)} |\boldsymbol{\rho}_k^{(\alpha)}|^2 \right)
$$

In matrix notation:
$$
\Phi^{(\alpha)} = \mathcal{N}^{(\alpha)} \exp\left( -(\boldsymbol{\rho}^{(\alpha)})^T \mathbf{D}^{(\alpha)} \boldsymbol{\rho}^{(\alpha)} \right)
$$

where $\mathbf{D}^{(\alpha)} = \text{diag}(\nu_1^{(\alpha)}, \nu_2^{(\alpha)}, \ldots)$ is the diagonal width matrix.

### 3.2 Combined Width Matrix in Target Jacobi Set

For a matrix element $\langle \Phi^{(\alpha)} | V(r_{ij}) | \Phi^{(\beta)} \rangle$, let $\gamma$ be the Jacobi set where pair $(i,j)$ is diagonal.

**Key idea:** Transform both bra and ket width matrices directly to the target set $\gamma$, then add.

**Step 1: Transform bra width matrix to set $\gamma$**

If $\boldsymbol{\rho}^{(\gamma)} = \mathbf{T}_{\alpha \to \gamma} \boldsymbol{\rho}^{(\alpha)}$, then $\boldsymbol{\rho}^{(\alpha)} = \mathbf{T}_{\alpha \to \gamma}^{-1} \boldsymbol{\rho}^{(\gamma)}$.

The bra's Gaussian exponent transforms as:
$$
(\boldsymbol{\rho}^{(\alpha)})^T \mathbf{D}^{(\alpha)} \boldsymbol{\rho}^{(\alpha)} = (\boldsymbol{\rho}^{(\gamma)})^T \underbrace{(\mathbf{T}_{\alpha \to \gamma}^{-1})^T \mathbf{D}^{(\alpha)} \mathbf{T}_{\alpha \to \gamma}^{-1}}_{\mathbf{W}^{(\gamma)}_{\text{bra}}} \boldsymbol{\rho}^{(\gamma)}
$$

**Step 2: Transform ket width matrix to set $\gamma$**

Similarly:
$$
\mathbf{W}^{(\gamma)}_{\text{ket}} = (\mathbf{T}_{\beta \to \gamma}^{-1})^T \mathbf{D}^{(\beta)} \mathbf{T}_{\beta \to \gamma}^{-1}
$$

**Step 3: Add the transformed width matrices**

$$
\boxed{\mathbf{W}^{(\gamma)} = \mathbf{W}^{(\gamma)}_{\text{bra}} + \mathbf{W}^{(\gamma)}_{\text{ket}} = (\mathbf{T}_{\alpha \to \gamma}^{-1})^T \mathbf{D}^{(\alpha)} \mathbf{T}_{\alpha \to \gamma}^{-1} + (\mathbf{T}_{\beta \to \gamma}^{-1})^T \mathbf{D}^{(\beta)} \mathbf{T}_{\beta \to \gamma}^{-1}}
$$

**Special cases:**

- If $\alpha = \gamma$: $\mathbf{W}^{(\gamma)}_{\text{bra}} = \mathbf{D}^{(\alpha)}$ (no transformation needed for bra)
- If $\beta = \gamma$: $\mathbf{W}^{(\gamma)}_{\text{ket}} = \mathbf{D}^{(\beta)}$ (no transformation needed for ket)
- If $\alpha = \beta = \gamma$: $\mathbf{W}^{(\gamma)} = \mathbf{D}^{(\alpha)}_{\text{bra}} + \mathbf{D}^{(\alpha)}_{\text{ket}}$ (diagonal, simplest case)

**Note:** The combined width matrix $\mathbf{W}^{(\gamma)}$ is generally **non-diagonal** unless both bra and ket are in the target set $\gamma$.

### 3.3 Cholesky Factorization

Decompose the width matrix in the target coordinate system:
$$
\mathbf{W}^{(\gamma)} = \mathbf{L} \mathbf{L}^T
$$

where $\mathbf{L}$ is lower triangular.

**For 2×2 (3-body):**
$$
l_{11} = \sqrt{W_{11}}
$$
$$
l_{21} = \frac{W_{12}}{l_{11}}
$$
$$
l_{22} = \sqrt{W_{22} - l_{21}^2}
$$

**For general $(N-1) \times (N-1)$:**
$$
l_{ii} = \sqrt{W_{ii} - \sum_{k=1}^{i-1} l_{ik}^2}
$$
$$
l_{ij} = \frac{1}{l_{jj}} \left( W_{ij} - \sum_{k=1}^{j-1} l_{ik} l_{jk} \right) \quad \text{for } i > j
$$

### 3.4 The Key Insight

The Cholesky transformation $\mathbf{y} = \mathbf{L} \boldsymbol{\rho}^{(\gamma)}$ decouples the first coordinate:

$$
\rho_1^{(\gamma)} = \frac{y_1}{l_{11}}
$$

Since pair $(i,j)$ is diagonal in set $\gamma$, we have $r_{ij} = |\boldsymbol{\rho}_1^{(\gamma)}| = |y_1|/l_{11}$.

The potential $V(r_{ij})$ depends only on $y_1$, and the integral factorizes!

---

## 4. Potential Integrals

### 4.1 General Formula

After the Cholesky transformation, the matrix element factorizes:

$$
\langle \Phi^{(\alpha)} | V | \Phi^{(\beta)} \rangle = \frac{\mathcal{N}^{(\alpha)} \mathcal{N}^{(\beta)} \cdot \pi^{3(N-1)/2}}{\prod_{k=1}^{N-1} l_{kk}^3} \cdot I_V
$$

where:
- $\mathcal{N}^{(\alpha)}, \mathcal{N}^{(\beta)}$ = normalization constants for bra and ket
- $l_{kk}$ = diagonal elements of the Cholesky factor $\mathbf{L}$
- $I_V$ = potential-specific integral (depends on $l_{11}$)

### 4.2 Potential Integral Formulas

The integral $I_V$ is defined as:

$$
I_V = 4\pi \int_0^\infty r^2 \, e^{-r^2} \, V\left(\frac{r}{l_{11}}\right) \, dr
$$

**Coulomb Potential: $V(r) = 1/r$**

$$
I_V = 2\pi l_{11}
$$

**Linear Potential: $V(r) = r$**

$$
I_V = \frac{2\pi}{l_{11}}
$$

**Gaussian Potential: $V(r) = e^{-\tau^2 r^2}$**

$$
I_V = \left( \frac{\pi l_{11}^2}{l_{11}^2 + \tau^2} \right)^{3/2}
$$

**Constant Potential: $V(r) = 1$**

$$
I_V = \pi^{3/2}
$$

### 4.3 Summary Table

| Potential Type | $V(r)$ | $I_V$ |
|----------------|--------|-------|
| Coulomb | $1/r$ | $2\pi l_{11}$ |
| Linear | $r$ | $2\pi / l_{11}$ |
| Gaussian | $e^{-\tau^2 r^2}$ | $\left(\frac{\pi l_{11}^2}{l_{11}^2 + \tau^2}\right)^{3/2}$ |
| Constant | $1$ | $\pi^{3/2}$ |

---

## 5. Complete Algorithm

### 5.1 General Matrix Element Formula

For basis functions in Jacobi sets $\alpha$ (bra) and $\beta$ (ket), with potential $V(r_{ij})$ diagonal in set $\gamma$:

$$
\boxed{\langle \Phi^{(\alpha)} | V_{ij} | \Phi^{(\beta)} \rangle = \frac{\mathcal{N}^{(\alpha)} \mathcal{N}^{(\beta)} \cdot \pi^{3(N-1)/2}}{\det(\mathbf{L})^3} \cdot I_V}
$$

where:
- $\mathcal{N}^{(\alpha)} = \prod_k (2\nu_k^{(\alpha)}/\pi)^{3/4}$
- $\det(\mathbf{L}) = \prod_k l_{kk}$
- $\mathbf{L}$ comes from Cholesky decomposition of $\mathbf{W}^{(\gamma)}$

### 5.2 Special Cases

**Case A: Same Jacobi set, diagonal pair ($\alpha = \beta = \gamma$)**

Simplest case. Width matrix is diagonal:
$$
\mathbf{W} = \mathbf{D}_{\text{bra}} + \mathbf{D}_{\text{ket}} = \text{diag}(\nu_1 + \nu'_1, \nu_2 + \nu'_2, \ldots)
$$

Cholesky factor: $l_{kk} = \sqrt{\nu_k + \nu'_k}$

**Case B: Same Jacobi set, non-diagonal pair ($\alpha = \beta \neq \gamma$)**

Width matrix starts diagonal, becomes non-diagonal after transformation to set $\gamma$.

**Case C: Different Jacobi sets ($\alpha \neq \beta$)**

Width matrix is non-diagonal from the start due to cross-set combination.

### 5.3 Overlap Matrix Element

For the overlap (identity operator), no transformation to a "target" set is needed. Work in set $\alpha$:

$$
\langle \Phi^{(\alpha)} | \Phi^{(\beta)} \rangle = \frac{\mathcal{N}^{(\alpha)} \mathcal{N}^{(\beta)} \cdot \pi^{3(N-1)/2}}{\det(\mathbf{L})^3}
$$

where $\mathbf{L}$ is from Cholesky decomposition of $\mathbf{W}^{(\alpha)} = \mathbf{D}^{(\alpha)} + \mathbf{T}_{\alpha \to \beta}^T \mathbf{D}^{(\beta)} \mathbf{T}_{\alpha \to \beta}$.

### 5.4 Kinetic Energy Matrix Element

The kinetic energy operator is:
$$
T = -\frac{\hbar^2}{2} \sum_{k=1}^{N-1} \frac{1}{\mu_k} \nabla_k^2
$$

For Gaussian basis functions, use the identity:
$$
\langle \Phi | T | \Phi' \rangle = \sum_{k} \frac{\hbar^2}{2\mu_k} \cdot 6\nu_k \nu'_k \cdot \langle \phi_k | \phi'_k \rangle_{\text{modified}}
$$

where the "modified" overlap uses widths adjusted for the kinetic term.

---

## 6. Algorithm Pseudocode

### 6.1 Main Procedure (General Case)

```text
PROCEDURE ComputeMatrixElement_General(
    bra_jacobi_set,      # Jacobi set α for bra
    ket_jacobi_set,      # Jacobi set β for ket
    bra_widths,          # diagonal widths D^(α) = [ν₁^α, ν₂^α, ...]
    ket_widths,          # diagonal widths D^(β) = [ν₁^β, ν₂^β, ...]
    target_pair,         # pair (i,j) for the potential
    potential_type,      # "coulomb", "linear", "gaussian", "constant", or "overlap"
    T_matrices           # dictionary of transformation matrices between sets
)

    OUTPUT:
        matrix_element   - the computed spatial matrix element

    STEPS:

    1. DETERMINE TARGET JACOBI SET
       IF potential_type == "overlap":
           # For overlap, use bra's set as target (arbitrary choice)
           target_set = bra_jacobi_set
       ELSE:
           target_set = jacobi_set_for_pair(target_pair)

    2. TRANSFORM BRA WIDTH MATRIX TO TARGET SET
       N_coord = length(bra_widths)
       D_bra = diag(bra_widths)

       IF bra_jacobi_set == target_set:
           W_bra = D_bra
       ELSE:
           T_bra_to_target = T_matrices[bra_jacobi_set → target_set]
           T_inv = inverse(T_bra_to_target)
           W_bra = transpose(T_inv) @ D_bra @ T_inv

    3. TRANSFORM KET WIDTH MATRIX TO TARGET SET
       D_ket = diag(ket_widths)

       IF ket_jacobi_set == target_set:
           W_ket = D_ket
       ELSE:
           T_ket_to_target = T_matrices[ket_jacobi_set → target_set]
           T_inv = inverse(T_ket_to_target)
           W_ket = transpose(T_inv) @ D_ket @ T_inv

    4. COMBINE WIDTH MATRICES
       W = W_bra + W_ket

    5. CHOLESKY FACTORIZATION
       L = cholesky_lower(W)

    6. COMPUTE POTENTIAL INTEGRAL
       l11 = L[1,1]
       SWITCH potential_type:
           CASE "overlap":  I_V = π^(3/2)
           CASE "coulomb":  I_V = 2 * π * l11
           CASE "linear":   I_V = 2 * π / l11
           CASE "gaussian": I_V = (π * l11² / (l11² + τ²))^(3/2)
           CASE "constant": I_V = π^(3/2)

    7. COMPUTE NORMALIZATION
       norm_bra = product((2*bra_widths[k]/π)^(3/4) for k in 1:N_coord)
       norm_ket = product((2*ket_widths[k]/π)^(3/4) for k in 1:N_coord)

    8. COMPUTE DETERMINANT FACTOR
       det_L = product(L[k,k] for k in 1:N_coord)

    9. ASSEMBLE RESULT
       matrix_element = norm_bra * norm_ket * π^(3*N_coord/2) * I_V / det_L^3

    RETURN matrix_element
```

### 6.2 Cholesky Factorization Subroutine

```text
PROCEDURE cholesky_lower(W)

    INPUT:
        W - symmetric positive definite matrix (N × N)

    OUTPUT:
        L - lower triangular matrix such that W = L @ L^T

    STEPS:

    N = size(W, 1)
    L = zeros(N, N)

    FOR j = 1 TO N:

        # Diagonal element
        sum = 0
        FOR k = 1 TO j-1:
            sum += L[j,k]²
        L[j,j] = sqrt(W[j,j] - sum)

        # Off-diagonal elements
        FOR i = j+1 TO N:
            sum = 0
            FOR k = 1 TO j-1:
                sum += L[i,k] * L[j,k]
            L[i,j] = (W[i,j] - sum) / L[j,j]

    RETURN L
```

### 6.3 Transformation Matrix Construction (3-body)

For 3-body systems, the transformation matrix can be computed directly using the tree-based method from Section 2.4, or using these explicit formulas.

```text
PROCEDURE build_transformation_matrix_3body(m1, m2, m3, from_set, to_set)

    INPUT:
        m1, m2, m3    - particle masses
        from_set      - source Jacobi set: "12", "13", or "23"
        to_set        - target Jacobi set: "12", "13", or "23"

    OUTPUT:
        T             - transformation matrix (from_set → to_set)

    STEPS:

    IF from_set == to_set:
        RETURN identity_matrix(2)

    M = m1 + m2 + m3  # total mass

    # Use the general tree-based method (Section 2.4)
    # For 3-body, we can write explicit formulas:

    # Define mass sums
    M12 = m1 + m2
    M13 = m1 + m3
    M23 = m2 + m3

    IF (from_set, to_set) == ("12", "13"):
        # ρ' = r1 - r3 = (m2/M12)ρ + (m3/M + M12/M)λ = (m2/M12)ρ + λ
        # λ' = (m1/M13)r1 - r2 + (m3/M13)r3
        #    = (m1/M13)(m2/M12)ρ + (m1/M13)(m3/M)λ + (m1/M12)ρ - (m3/M)λ + (m3/M13)(-M12/M)λ
        T = zeros(2, 2)
        T[1,1] = m2 / M12
        T[1,2] = 1
        T[2,1] = m1 * m2 / (M12 * M13) + m1 / M12
        T[2,2] = m1 * m3 / (M13 * M) - m3 / M - m3 * M12 / (M13 * M)
        # Simplify T[2,1] = m1 * M / (M12 * M13)
        # Simplify T[2,2] = -m3 / M13

    ELSE IF (from_set, to_set) == ("12", "23"):
        T = zeros(2, 2)
        T[1,1] = -m1 / M12
        T[1,2] = 1
        T[2,1] = -m1 * m2 / (M12 * M23) - m2 / M12
        T[2,2] = -m1 / M23

    ELSE IF (from_set, to_set) == ("13", "12"):
        RETURN inverse(build_transformation_matrix_3body(m1, m2, m3, "12", "13"))

    ELSE IF (from_set, to_set) == ("13", "23"):
        T = zeros(2, 2)
        T[1,1] = m1 / M13
        T[1,2] = -1
        T[2,1] = -m3 * M / (M13 * M23)
        T[2,2] = -m2 / M23

    ELSE IF (from_set, to_set) == ("23", "12"):
        RETURN inverse(build_transformation_matrix_3body(m1, m2, m3, "12", "23"))

    ELSE IF (from_set, to_set) == ("23", "13"):
        RETURN inverse(build_transformation_matrix_3body(m1, m2, m3, "13", "23"))

    RETURN T
```

**Simplified forms for (12)→(13):**
$$
T_{11} = \frac{m_2}{M_{12}}, \quad T_{12} = 1, \quad T_{21} = \frac{m_1 M}{M_{12} M_{13}}, \quad T_{22} = -\frac{m_3}{M_{13}}
$$

### 6.4 Simplified Procedure (Same Jacobi Set)

For the common case where bra and ket are in the same Jacobi set:

```text
PROCEDURE ComputeMatrixElement_SameSet(
    bra_widths,          # [ν₁, ν₂, ...] for bra
    ket_widths,          # [ν'₁, ν'₂, ...] for ket
    target_pair,         # pair (i,j) for the potential
    potential_type,      # "coulomb", "linear", "gaussian", "constant"
    T_to_target          # transformation to target Jacobi set (identity if diagonal)
)

    1. D = diag(bra_widths + ket_widths)

    2. IF T_to_target is not identity:
           T_inv = inverse(T_to_target)
           W = transpose(T_inv) @ D @ T_inv
       ELSE:
           W = D

    3. L = cholesky_lower(W)

    4. Compute I_V using l11 = L[1,1]

    5. norm = product((2*ν/π)^(3/4) * (2*ν'/π)^(3/4))
       det_L = product(L[k,k])

    6. RETURN norm * π^(3*N_coord/2) * I_V / det_L^3
```

---

## 7. Worked Example: 3-Body Cross-Set Matrix Element

Consider $\langle \Phi^{(12)3} | V_{13} | \Phi^{(13)2} \rangle$ for equal masses $m_1 = m_2 = m_3 = m$.

**Given:**

- Bra in set (12)3 with widths $(\nu_\rho, \nu_\lambda)$
- Ket in set (13)2 with widths $(\nu'_\rho, \nu'_\lambda)$
- Potential $V_{13}$ is diagonal in set (13)2, so target set $\gamma$ = (13)2

**Step 1: Compute transformation matrix**

For equal masses, using the formulas from Section 6.3:
- $M_{12} = M_{13} = 2m$, $M = 3m$
- $T_{11} = m_2/M_{12} = 1/2$
- $T_{12} = 1$
- $T_{21} = m_1 M/(M_{12} M_{13}) = m \cdot 3m/(2m \cdot 2m) = 3/4$
- $T_{22} = -m_3/M_{13} = -m/(2m) = -1/2$

$$
\mathbf{T}_{(12)3 \to (13)2} = \begin{pmatrix} 1/2 & 1 \\ 3/4 & -1/2 \end{pmatrix}
$$

Compute the inverse:
$$
\det(\mathbf{T}) = (1/2)(-1/2) - (1)(3/4) = -1/4 - 3/4 = -1
$$
$$
\mathbf{T}^{-1} = \frac{1}{-1} \begin{pmatrix} -1/2 & -1 \\ -3/4 & 1/2 \end{pmatrix} = \begin{pmatrix} 1/2 & 1 \\ 3/4 & -1/2 \end{pmatrix}
$$

Note: For this particular case, $\mathbf{T}^{-1} = \mathbf{T}$ (the matrix is an involution).

**Step 2: Transform bra to target set (13)2**

$$
\mathbf{W}^{(13)2}_{\text{bra}} = (\mathbf{T}^{-1})^T \begin{pmatrix} \nu_\rho & 0 \\ 0 & \nu_\lambda \end{pmatrix} \mathbf{T}^{-1}
$$

Computing $(\mathbf{T}^{-1})^T = \begin{pmatrix} 1/2 & 3/4 \\ 1 & -1/2 \end{pmatrix}$:

$$
(\mathbf{T}^{-1})^T \mathbf{D} = \begin{pmatrix} \nu_\rho/2 & 3\nu_\lambda/4 \\ \nu_\rho & -\nu_\lambda/2 \end{pmatrix}
$$

$$
\mathbf{W}^{(13)2}_{\text{bra}} = \begin{pmatrix} \nu_\rho/2 & 3\nu_\lambda/4 \\ \nu_\rho & -\nu_\lambda/2 \end{pmatrix} \begin{pmatrix} 1/2 & 1 \\ 3/4 & -1/2 \end{pmatrix} = \begin{pmatrix} \frac{\nu_\rho}{4} + \frac{9\nu_\lambda}{16} & \frac{\nu_\rho}{2} - \frac{3\nu_\lambda}{8} \\ \frac{\nu_\rho}{2} - \frac{3\nu_\lambda}{8} & \nu_\rho + \frac{\nu_\lambda}{4} \end{pmatrix}
$$

**Step 3: Ket is already in target set (13)2**

No transformation needed:
$$
\mathbf{W}^{(13)2}_{\text{ket}} = \begin{pmatrix} \nu'_\rho & 0 \\ 0 & \nu'_\lambda \end{pmatrix}
$$

**Step 4: Combine width matrices**

$$
\mathbf{W}^{(13)2} = \mathbf{W}^{(13)2}_{\text{bra}} + \mathbf{W}^{(13)2}_{\text{ket}} = \begin{pmatrix} \frac{\nu_\rho}{4} + \frac{9\nu_\lambda}{16} + \nu'_\rho & \frac{\nu_\rho}{2} - \frac{3\nu_\lambda}{8} \\ \frac{\nu_\rho}{2} - \frac{3\nu_\lambda}{8} & \nu_\rho + \frac{\nu_\lambda}{4} + \nu'_\lambda \end{pmatrix}
$$

**Step 5: Cholesky factorization**

Let $W_{11} = \frac{\nu_\rho}{4} + \frac{9\nu_\lambda}{16} + \nu'_\rho$, $W_{12} = \frac{\nu_\rho}{2} - \frac{3\nu_\lambda}{8}$, $W_{22} = \nu_\rho + \frac{\nu_\lambda}{4} + \nu'_\lambda$.

$$
l_{11} = \sqrt{W_{11}}, \quad l_{21} = \frac{W_{12}}{l_{11}}, \quad l_{22} = \sqrt{W_{22} - l_{21}^2}
$$

**Step 6: Compute matrix element**

$$
\langle \Phi^{(12)3} | V_{13} | \Phi^{(13)2} \rangle = \frac{\mathcal{N}^{(12)3} \mathcal{N}^{(13)2} \cdot \pi^3}{(l_{11} l_{22})^3} \cdot I_V
$$

where $I_V = 2\pi l_{11}$ for Coulomb, $2\pi/l_{11}$ for linear, etc.

---

## 8. Algorithm Properties

### 8.1 Computational Complexity

| Operation | Complexity | Notes |
|-----------|------------|-------|
| Width matrix combination | $O((N-1)^3)$ | Matrix multiplication |
| Transformation | $O((N-1)^3)$ | Matrix multiplication |
| Cholesky | $O((N-1)^3 / 6)$ | Lower triangular factorization |
| Potential integral | $O(1)$ | Closed-form formula |
| **Per matrix element** | $O((N-1)^3)$ | Dominated by matrix operations |

For small $N$ (mesons, baryons, tetraquarks), this is essentially $O(1)$.

### 8.2 Numerical Stability

- **Cholesky decomposition** is numerically stable for positive definite matrices
- The width matrix $\mathbf{W}$ is always positive definite (sum of positive Gaussian widths)
- Cross-set combinations preserve positive definiteness
- No cancellation errors from alternating signs

**Note on solving the eigenvalue problem:** As the basis size increases, the overlap matrix $S$ may become ill-conditioned due to nearly linearly dependent Gaussians. The generalized eigenvalue problem $H\mathbf{c} = ES\mathbf{c}$ should be solved using eigenvalue truncation to ensure numerical stability. See `algorithm_eigenvalue_truncation.md` for the complete algorithm.

### 8.3 Advantages

1. **Fully analytical** — no numerical integration
2. **Unified framework** — same algorithm handles all cases (same set, cross-set, diagonal pair, non-diagonal pair)
3. **Efficient** — $O(1)$ per matrix element for few-body systems
4. **Scalable** — extends naturally to N-body systems

---

## 9. Implementation Considerations

### 9.1 Precomputing Transformation Matrices

For efficiency, precompute all transformation matrices $\mathbf{T}_{\alpha \to \beta}$ before the main calculation loop:

```text
FOR each pair (α, β) of Jacobi sets:
    T_matrices[(α, β)] = build_transformation_matrix(α, β)
    T_matrices[(β, α)] = inverse(T_matrices[(α, β)])
```

### 9.2 Caching Cholesky Factors

If the same width combination appears multiple times (e.g., for different potential types), cache the Cholesky factor $\mathbf{L}$ and determinant.

### 9.3 Basis Set Organization

Organize the basis as:
```text
basis = [
    (jacobi_set="12", widths=[ν₁, ν₂]),
    (jacobi_set="12", widths=[ν₃, ν₄]),
    (jacobi_set="13", widths=[ν₅, ν₆]),
    ...
]
```

Then loop over all pairs to build the Hamiltonian matrix.

---

## 10. Permutation of Spatial States for Antisymmetrization

### 10.1 Motivation

For identical fermions, we need to compute matrix elements of the form $\langle \Phi_a | O | P\Phi_b \rangle$ where $P$ is a permutation of particle labels. This section describes how permutations act on spatial wave functions in the Cholesky algorithm framework.

### 10.2 Permutation as Jacobi Tree Transformation

A permutation $P$ of particle labels transforms one Jacobi tree into another. For a Gaussian defined in Jacobi set with tree $\mathcal{T}$, the permuted Gaussian is naturally expressed in the permuted tree $P \circ \mathcal{T}$.

```text
PROCEDURE apply_permutation_to_tree(tree, perm):
    # tree: Jacobi clustering tree, e.g., ((1,2),3)
    # perm: permutation array where perm[i] = P(i)
    # Returns: permuted tree

    IF tree is a leaf (integer):
        RETURN perm[tree]  # Relabel particle

    (left, right) = tree
    new_left = apply_permutation_to_tree(left, perm)
    new_right = apply_permutation_to_tree(right, perm)
    RETURN (new_left, new_right)
```

**Examples for 3-body:**

| Original Tree | Permutation | Permuted Tree |
|---------------|-------------|---------------|
| ((1,2),3) | $P_{12}=[2,1,3]$ | ((2,1),3) |
| ((1,2),3) | $P_{13}=[3,2,1]$ | ((3,2),1) |
| ((1,2),3) | $P_{123}=[2,3,1]$ | ((2,3),1) |

### 10.3 Transformation Matrix for Permuted Coordinates

The permuted tree defines a different Jacobi coordinate system. The transformation between original and permuted coordinates is computed using the standard tree-based method (Section 2.4).

```text
PROCEDURE get_permutation_transformation(tree, perm, masses):
    # Returns transformation matrix T such that:
    # rho_permuted = T @ rho_original

    permuted_tree = apply_permutation_to_tree(tree, perm)
    T = build_transformation_matrix(tree, permuted_tree, masses)
    RETURN T
```

### 10.4 Matrix Element with Permutation

To compute $\langle \Phi_a | O | P\Phi_b \rangle$:

```text
PROCEDURE spatial_ME_with_permutation(
    bra,        # (tree_bra, widths_bra)
    ket,        # (tree_ket, widths_ket)
    operator,   # potential type and pair
    perm,       # permutation array
    masses      # particle masses
):
    # Step 1: Apply permutation to ket's Jacobi tree
    permuted_tree = apply_permutation_to_tree(ket.tree, perm)

    # Step 2: Create effective ket with permuted tree
    # The permuted ket has the same widths but different tree
    permuted_ket = (permuted_tree, ket.widths)

    # Step 3: Use standard matrix element algorithm
    # The algorithm handles cross-set matrix elements automatically
    RETURN ComputeMatrixElement_General(
        bra.tree, permuted_tree,
        bra.widths, ket.widths,
        operator.pair, operator.type,
        T_matrices
    )
```

### 10.5 Special Case: Transposition in Diagonal Jacobi Set

When the permutation is a transposition $P_{ij}$ of the pair that is diagonal in the ket's Jacobi set, the transformation is particularly simple.

**Example: $P_{12}$ in Jacobi set (12)3**

The coordinates are:
$$
\boldsymbol{\rho} = \mathbf{r}_1 - \mathbf{r}_2, \quad \boldsymbol{\lambda} = \frac{m_1\mathbf{r}_1 + m_2\mathbf{r}_2}{m_1+m_2} - \mathbf{r}_3
$$

Under $P_{12}$:
$$
P_{12}: \boldsymbol{\rho} \to -\boldsymbol{\rho}, \quad \boldsymbol{\lambda} \to \boldsymbol{\lambda}
$$

**Transformation matrix:**
$$
\mathbf{T}_{P_{12}} = \begin{pmatrix} -1 & 0 \\ 0 & 1 \end{pmatrix}
$$

**For s-wave Gaussians:** Since $\exp(-\nu\rho^2)$ is even in $\boldsymbol{\rho}$, the Gaussian is unchanged:
$$
P_{12} \exp(-\nu\rho^2 - \mu\lambda^2) = \exp(-\nu\rho^2 - \mu\lambda^2)
$$

Therefore, for s-wave Gaussians in this case:
$$
\langle \Phi_a | O | P_{12}\Phi_b \rangle = \langle \Phi_a | O | \Phi_b \rangle
$$

### 10.6 Width Matrix Transformation Under Permutation

For a Gaussian with diagonal width matrix $\mathbf{D}$ in tree $\mathcal{T}$, after permutation $P$:

```text
# Original Gaussian: exp(-rho^T @ D @ rho) in coordinates of tree T
# Permuted Gaussian: exp(-rho'^T @ D @ rho') in coordinates of permuted tree T'

# To express in target coordinates (e.g., for Cholesky algorithm):
T_perm = build_transformation_matrix(permuted_tree, target_tree, masses)
T_inv = inverse(T_perm)
W_permuted = T_inv^T @ D @ T_inv
```

This permuted width matrix enters the combined width matrix computation:
$$
\mathbf{W}^{(\gamma)} = \mathbf{W}^{(\gamma)}_{\text{bra}} + \mathbf{W}^{(\gamma)}_{\text{ket,permuted}}
$$

### 10.7 Algorithm Integration

The permutation handling integrates seamlessly with the existing Cholesky algorithm:

```text
PROCEDURE ComputeMatrixElement_WithPermutation(
    bra_tree, bra_widths,
    ket_tree, ket_widths,
    target_pair,
    potential_type,
    perm,
    masses
):
    # Step 1: Apply permutation to ket tree
    permuted_ket_tree = apply_permutation_to_tree(ket_tree, perm)

    # Step 2: Determine target Jacobi set
    IF potential_type == "overlap":
        target_set = bra_tree
    ELSE:
        target_set = jacobi_set_for_pair(target_pair)

    # Step 3: Transform bra to target (unchanged from standard algorithm)
    D_bra = diag(bra_widths)
    IF bra_tree == target_set:
        W_bra = D_bra
    ELSE:
        T_bra = build_transformation_matrix(bra_tree, target_set, masses)
        T_inv = inverse(T_bra)
        W_bra = T_inv^T @ D_bra @ T_inv

    # Step 4: Transform PERMUTED ket to target
    D_ket = diag(ket_widths)
    IF permuted_ket_tree == target_set:
        W_ket = D_ket
    ELSE:
        T_ket = build_transformation_matrix(permuted_ket_tree, target_set, masses)
        T_inv = inverse(T_ket)
        W_ket = T_inv^T @ D_ket @ T_inv

    # Step 5: Proceed with standard Cholesky algorithm
    W = W_bra + W_ket
    L = cholesky_lower(W)
    # ... compute potential integral and assemble result
```

### 10.8 Caching Permuted Transformations

For efficiency, precompute transformation matrices for all relevant permutations:

```text
PROCEDURE precompute_permutation_transforms(trees, perms, masses):
    T_cache = {}

    FOR tree in trees:
        FOR perm in perms:
            permuted_tree = apply_permutation_to_tree(tree, perm)
            key = (tree, perm)
            T_cache[key] = build_transformation_matrix(tree, permuted_tree, masses)

    RETURN T_cache
```

### 10.9 Verification

**Test 1: Identity permutation**

For $P = \text{identity}$, the permuted tree equals the original tree, and:
$$
\langle \Phi_a | O | P\Phi_b \rangle = \langle \Phi_a | O | \Phi_b \rangle
$$

**Test 2: Double permutation**

For any permutation $P$:
$$
\langle \Phi_a | O | P^2\Phi_b \rangle = \langle \Phi_a | O | (P \circ P)\Phi_b \rangle
$$

**Test 3: Hermiticity**

$$
\langle \Phi_a | O | P\Phi_b \rangle = \langle P^{-1}\Phi_a | O | \Phi_b \rangle^* = \langle \Phi_a | P^{-1} O^{\dagger} | \Phi_b \rangle
$$

**Test 4: S-wave transposition invariance**

For s-wave Gaussians and $P_{ij}$ diagonal in the ket's Jacobi set:
$$
\langle \Phi_a | O | P_{ij}\Phi_b \rangle = \langle \Phi_a | O | \Phi_b \rangle
$$

For the full antisymmetrization algorithm combining spatial and discrete sectors, see `algorithm_antisymmetrization.md`.

---

## 11. References

- This algorithm is based on the coordinate transformation properties of Jacobi coordinates and the standard Cholesky decomposition of linear algebra.
- The Gaussian integral identities are elementary results from quantum mechanics.
- For discrete wave function construction and permutation matrices, see `algorithm_discrete_MM.md`.
- For the complete antisymmetrization algorithm, see `algorithm_antisymmetrization.md`.
- For detailed physics context, see `physics.md`.
