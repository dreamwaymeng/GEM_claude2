# Eigenvalue Truncation for Ill-Conditioned Basis

This document describes the eigenvalue truncation algorithm for solving generalized eigenvalue problems arising from non-orthogonal Gaussian basis sets. The method handles numerical instability when the overlap matrix becomes nearly singular.

---

## 1. Overview

### 1.1 The Problem

In the Gaussian Expansion Method, the wave function is expanded in a non-orthogonal basis:

$$
\Psi = \sum_n c_n \Phi_n
$$

The variational principle leads to a **generalized eigenvalue problem**:

$$
H \mathbf{c} = E S \mathbf{c}
$$

where:
- $H$ is the Hamiltonian matrix: $H_{nm} = \langle \Phi_n | \hat{H} | \Phi_m \rangle$
- $S$ is the overlap matrix: $S_{nm} = \langle \Phi_n | \Phi_m \rangle$
- $E$ are the energy eigenvalues
- $\mathbf{c}$ are the expansion coefficients

### 1.2 Origin of Ill-Conditioning

As the basis size $N$ increases, problems arise:

1. **Near-linear dependence**: Gaussians with similar widths become nearly linearly dependent
2. **Small overlap eigenvalues**: Some eigenvalues of $S$ approach zero
3. **Large condition number**: $\kappa(S) = \lambda_{\max}/\lambda_{\min} \to \infty$

**Consequences:**
- Numerical instability in solving the eigenvalue problem
- Spurious eigenvalues from near-null space directions
- Loss of precision in computed eigenvectors
- Potential for completely wrong results

### 1.3 The Solution

Transform to an orthonormal basis by:

1. Diagonalize the overlap matrix
2. Remove directions with eigenvalues below a threshold
3. Project the Hamiltonian into the reduced space
4. Solve a standard (not generalized) eigenvalue problem

This method:
- Removes linearly dependent directions
- Preserves the variational principle for retained modes
- Provides diagnostic information about basis quality

---

## 2. Mathematical Formulation

### 2.1 Diagonalize the Overlap Matrix

The overlap matrix is symmetric and positive semi-definite, so it can be diagonalized:

$$
S = U \Lambda U^T
$$

where:
- $U$ is orthogonal: $U^T U = U U^T = I$
- $\Lambda = \text{diag}(\lambda_1, \lambda_2, \ldots, \lambda_N)$ with $\lambda_i \geq 0$

### 2.2 Identify and Truncate Small Eigenvalues

Select eigenvalues above threshold $\epsilon$:

$$
\mathcal{I} = \{i : \lambda_i > \epsilon\}
$$

Let $N_{\text{eff}} = |\mathcal{I}|$ be the effective basis dimension.

**Typical thresholds:**
- $\epsilon \sim 10^{-8}$ to $10^{-10}$ for most calculations
- Larger threshold = more stable but smaller effective basis
- Smaller threshold = larger effective basis but less stable

### 2.3 Construct Transformation Matrix

Define the truncated matrices:
- $U_{\text{trunc}}$: columns of $U$ corresponding to $\mathcal{I}$ (size $N \times N_{\text{eff}}$)
- $\Lambda_{\text{trunc}}$: diagonal matrix with $\lambda_i$ for $i \in \mathcal{I}$ (size $N_{\text{eff}} \times N_{\text{eff}}$)

The transformation to orthonormal basis:

$$
X = U_{\text{trunc}} \Lambda_{\text{trunc}}^{-1/2}
$$

**Verification:** $X^T S X = I_{N_{\text{eff}}}$

**Proof:**

First, note that $U_{\text{trunc}}^T S U_{\text{trunc}} = \Lambda_{\text{trunc}}$. This follows from:
$$
U_{\text{trunc}}^T S U_{\text{trunc}} = U_{\text{trunc}}^T (U \Lambda U^T) U_{\text{trunc}}
$$

Since the columns of $U_{\text{trunc}}$ are a subset of the orthonormal columns of $U$, we have $U^T U_{\text{trunc}} = E_{\text{select}}$ where $E_{\text{select}}$ is a selection matrix (with one 1 per column). This extracts the corresponding diagonal entries of $\Lambda$, giving $\Lambda_{\text{trunc}}$.

Therefore:
$$
X^T S X = \Lambda_{\text{trunc}}^{-1/2} U_{\text{trunc}}^T S U_{\text{trunc}} \Lambda_{\text{trunc}}^{-1/2} = \Lambda_{\text{trunc}}^{-1/2} \Lambda_{\text{trunc}} \Lambda_{\text{trunc}}^{-1/2} = I
$$

### 2.4 Transform Hamiltonian

Project the Hamiltonian into the orthonormal subspace:

$$
\tilde{H} = X^T H X
$$

This is an $N_{\text{eff}} \times N_{\text{eff}}$ matrix.

### 2.5 Solve Standard Eigenvalue Problem

The generalized eigenvalue problem reduces to a standard one:

$$
\tilde{H} \tilde{\mathbf{c}} = E \tilde{\mathbf{c}}
$$

This can be solved using standard symmetric eigenvalue solvers.

### 2.6 Back-Transform Eigenvectors

The eigenvectors in the original basis:

$$
\mathbf{c} = X \tilde{\mathbf{c}}
$$

**Important properties:**

1. The returned eigenvectors span only the $N_{\text{eff}}$-dimensional subspace. There are $N_{\text{eff}}$ eigenvalues and eigenvectors, not $N$.

2. The eigenvectors are **S-orthonormal** (not I-orthonormal):
$$
\mathbf{c}_i^T S \mathbf{c}_j = \delta_{ij}
$$

**Proof:** Since $\tilde{\mathbf{c}}_i$ are orthonormal eigenvectors of the symmetric matrix $\tilde{H}$:
$$
\mathbf{c}_i^T S \mathbf{c}_j = \tilde{\mathbf{c}}_i^T X^T S X \tilde{\mathbf{c}}_j = \tilde{\mathbf{c}}_i^T I \tilde{\mathbf{c}}_j = \delta_{ij}
$$

This is the correct normalization for generalized eigenvalue problems.

---

## 3. Algorithm

### 3.1 Main Procedure

```text
PROCEDURE solve_generalized_eigenvalue_truncated(H, S, epsilon_rel=1e-10):
    # Input:
    #   H: Hamiltonian matrix (N × N, symmetric)
    #   S: Overlap matrix (N × N, symmetric positive semi-definite)
    #   epsilon_rel: relative truncation threshold (recommended over absolute)
    #
    # Output:
    #   energies: array of eigenvalues (length N_eff)
    #   coefficients: eigenvector matrix in original basis (N × N_eff)

    N = size(H, 1)

    # Step 1: Diagonalize overlap matrix
    # S = U @ Lambda @ U^T
    eigenvalues_S, U = symmetric_eigendecomposition(S)

    # Step 1b: Compute adaptive threshold (normalization-independent)
    lambda_max = max(eigenvalues_S)
    epsilon = epsilon_rel * lambda_max

    # Step 2: Identify valid eigenvalues (above threshold)
    valid_indices = [i for i in 1..N if eigenvalues_S[i] > epsilon]
    N_eff = length(valid_indices)

    IF N_eff == 0:
        ERROR "All overlap eigenvalues below threshold - basis is degenerate"

    IF N_eff < N:
        PRINT "Truncating basis: {N} -> {N_eff} (removed {N - N_eff} near-null directions)"

    # Step 3: Construct transformation matrix
    U_trunc = U[:, valid_indices]                    # N × N_eff
    Lambda_trunc = eigenvalues_S[valid_indices]      # N_eff values

    # X = U_trunc @ diag(1/sqrt(Lambda_trunc))
    X = U_trunc @ diag(1.0 / sqrt(Lambda_trunc))     # N × N_eff

    # Step 4: Transform Hamiltonian to orthonormal basis
    H_tilde = X^T @ H @ X                            # N_eff × N_eff

    # Step 5: Solve standard eigenvalue problem
    energies, C_tilde = symmetric_eigendecomposition(H_tilde)

    # Step 6: Back-transform to original basis
    coefficients = X @ C_tilde                        # N × N_eff

    RETURN energies, coefficients
```

### 3.2 Symmetric Eigendecomposition

Use stable algorithms for symmetric eigenvalue problems:

```text
PROCEDURE symmetric_eigendecomposition(A):
    # Input: A - symmetric matrix (N × N)
    # Output: eigenvalues (sorted ascending), eigenvectors (columns)

    # Recommended implementations:
    # - LAPACK: dsyev or dsyevd (divide-and-conquer, faster for large N)
    # - NumPy: numpy.linalg.eigh
    # - Julia: LinearAlgebra.eigen(Symmetric(A))

    ENSURE A is exactly symmetric: A = (A + A^T) / 2

    eigenvalues, eigenvectors = eigh(A)  # Returns sorted eigenvalues

    RETURN eigenvalues, eigenvectors
```

---

## 4. Threshold Selection

### 4.1 Guidelines by Basis Size

| Basis Size | Recommended $\epsilon$ | Expected Truncation |
|------------|------------------------|---------------------|
| $N < 50$ | $10^{-12}$ | Minimal |
| $50 \leq N < 200$ | $10^{-10}$ | Few modes |
| $200 \leq N < 500$ | $10^{-8}$ | ~5-10% |
| $N \geq 500$ | $10^{-6}$ to $10^{-8}$ | ~10-20% |

### 4.2 Effect of Basis Normalization

The eigenvalue spectrum of $S$ depends on the normalization convention of the basis functions.

**Case 1: Normalized basis** ($\langle \Phi_n | \Phi_n \rangle = 1$)

- Diagonal elements of $S$ are all 1
- Trace constraint: $\sum_i \lambda_i = \text{tr}(S) = N$
- Eigenvalues satisfy $0 \leq \lambda_i \leq N$ (since $\lambda_i \geq 0$ and sum to $N$)
- Typical range: $\lambda_{\max} \sim O(1)$ to $O(N)$ depending on overlap structure
- For orthonormal basis: all $\lambda_i = 1$
- Absolute thresholds like $\epsilon = 10^{-10}$ are meaningful

**Case 2: Unnormalized basis**

- Diagonal elements can vary widely: $S_{nn} = \langle \Phi_n | \Phi_n \rangle$
- Eigenvalue scale depends on the normalization
- Absolute thresholds may be inappropriate

**In the GEM framework:** The spatial matrix element formulas in `algorithm_spatial_MM.md` include normalization factors:

$$
\mathcal{N}^{(\alpha)} = \prod_k \left(\frac{2\nu_k^{(\alpha)}}{\pi}\right)^{3/4}
$$

This ensures $\langle \Phi_n | \Phi_n \rangle = 1$ for each basis function. However, when combining with discrete wave functions (spin, color), the overall normalization may differ.

**Recommendation:** Always use **relative thresholds** (Section 4.3) to avoid normalization-dependent behavior.

### 4.3 Adaptive (Relative) Threshold

A more robust approach uses a relative threshold:

$$
\epsilon = \epsilon_{\text{rel}} \cdot \lambda_{\max}
$$

where $\epsilon_{\text{rel}} \sim 10^{-10}$ to $10^{-8}$.

```text
PROCEDURE compute_adaptive_threshold(S, epsilon_rel):
    eigenvalues_S = eigenvalues(S)
    lambda_max = max(eigenvalues_S)
    epsilon = epsilon_rel * lambda_max
    RETURN epsilon
```

**Advantages:**
- Automatically scales with the overall magnitude of the overlap matrix
- Independent of basis normalization convention
- More robust across different systems and basis sizes

**Note:** For normalized basis functions, $\lambda_{\max} \approx O(1)$ to $O(N)$, so relative and absolute thresholds give similar results. For unnormalized basis, relative thresholds are essential.

### 4.4 Convergence Testing

To verify threshold choice, compute eigenvalues with multiple thresholds:

```text
PROCEDURE test_threshold_stability(H, S):
    thresholds = [1e-6, 1e-8, 1e-10, 1e-12]

    FOR epsilon in thresholds:
        energies, _ = solve_generalized_eigenvalue_truncated(H, S, epsilon)
        N_eff = length(energies)
        E_0 = energies[1]  # Ground state energy
        PRINT "epsilon = {epsilon}: E_0 = {E_0:.10f}, N_eff = {N_eff}"
```

**Expected behavior:**
- Physical eigenvalues remain stable across thresholds
- Spurious eigenvalues change significantly
- $N_{\text{eff}}$ increases as threshold decreases

---

## 5. Verification and Diagnostics

### 5.1 Condition Number Analysis

```text
PROCEDURE analyze_conditioning(S, epsilon):
    eigenvalues_S = eigenvalues(S)
    eigenvalues_S = sort(eigenvalues_S, descending=True)

    lambda_max = eigenvalues_S[1]
    lambda_min_orig = eigenvalues_S[N]

    # Find minimum eigenvalue above threshold
    valid_eigenvalues = [lam for lam in eigenvalues_S if lam > epsilon]
    lambda_min_trunc = min(valid_eigenvalues)

    condition_original = lambda_max / lambda_min_orig
    condition_truncated = lambda_max / lambda_min_trunc

    PRINT "Overlap matrix diagnostics:"
    PRINT "  Maximum eigenvalue: {lambda_max}"
    PRINT "  Minimum eigenvalue: {lambda_min_orig}"
    PRINT "  Original condition number: {condition_original:.2e}"
    PRINT "  Eigenvalues below threshold: {N - length(valid_eigenvalues)}"
    PRINT "  Truncated condition number: {condition_truncated:.2e}"
    PRINT "  Improvement factor: {condition_original / condition_truncated:.2e}"

    RETURN condition_original, condition_truncated
```

### 5.2 Orthonormality Verification

After constructing $X$, verify the transformation is correct:

```text
PROCEDURE verify_orthonormality(X, S, tol=1e-10):
    # Check X^T S X = I
    I_approx = X^T @ S @ X
    I_exact = identity_matrix(size(I_approx, 1))

    error = frobenius_norm(I_approx - I_exact)

    IF error > tol:
        WARNING "Orthonormality error: {error}"
    ELSE:
        PRINT "Orthonormality verified: ||X^T S X - I|| = {error:.2e}"

    RETURN error
```

### 5.3 Eigenvalue Spectrum Plot

For visual diagnostics, plot the eigenvalue spectrum:

```text
PROCEDURE plot_eigenvalue_spectrum(S, epsilon):
    eigenvalues_S = sort(eigenvalues(S), descending=True)

    # Plot on log scale
    PLOT log10(eigenvalues_S) vs index
    DRAW horizontal line at log10(epsilon)
    MARK truncation point

    # Annotate
    LABEL "Above threshold: {N_eff} modes"
    LABEL "Below threshold: {N - N_eff} modes"
```

---

## 6. Implementation Considerations

### 6.1 Numerical Stability

1. **Symmetrize matrices before decomposition:**
   ```text
   H = (H + H^T) / 2
   S = (S + S^T) / 2
   ```

2. **Use stable eigensolvers:**
   - LAPACK `dsyev` or `dsyevd` (divide-and-conquer)
   - Avoid general eigensolvers for symmetric matrices

3. **Check for negative overlap eigenvalues:**
   ```text
   IF any(eigenvalues_S < -epsilon):
       WARNING "Overlap matrix has negative eigenvalues - check matrix construction"
   ```

### 6.2 Memory Efficiency

For large $N$:

1. **Store only truncated matrices:**
   ```text
   # Don't store full U, only U_trunc
   eigenvalues_S, U = symmetric_eigendecomposition(S)
   valid_mask = eigenvalues_S > epsilon
   U_trunc = U[:, valid_mask]
   # Delete U to free memory
   ```

2. **Compute $X$ column-by-column if needed:**
   ```text
   X = zeros(N, N_eff)
   FOR k = 1 TO N_eff:
       X[:, k] = U_trunc[:, k] / sqrt(Lambda_trunc[k])
   ```

3. **The transformed Hamiltonian is smaller:**
   - $\tilde{H}$ is $N_{\text{eff}} \times N_{\text{eff}}$
   - Significant memory savings when many modes are truncated

### 6.3 Computational Complexity

| Operation | Complexity | Notes |
|-----------|------------|-------|
| Eigendecomposition of $S$ | $O(N^3)$ | Dominant cost |
| Construct $X$ | $O(N \cdot N_{\text{eff}})$ | Column scaling |
| Transform $H$ | $O(N^2 \cdot N_{\text{eff}})$ | Two matrix multiplications |
| Eigendecomposition of $\tilde{H}$ | $O(N_{\text{eff}}^3)$ | Often much smaller |
| Back-transform eigenvectors | $O(N \cdot N_{\text{eff}}^2)$ | Matrix multiplication |

**Total:** $O(N^3)$ dominated by initial eigendecomposition.

---

## 7. Integration with GEM Framework

### 7.1 Replacing Direct Solver

The truncation algorithm wraps the standard eigenvalue solver:

```text
# Before (may be ill-conditioned):
energies, coeffs = generalized_eig(H, S)

# After (numerically stable):
energies, coeffs = solve_generalized_eigenvalue_truncated(H, S, epsilon=1e-10)
```

### 7.2 Interface with Hamiltonian Construction

```text
PROCEDURE compute_hadron_spectrum(
    spatial_basis,
    discrete_basis,
    potential_params,
    epsilon=1e-10
):
    # Build matrices using existing algorithms
    H, S = build_antisymmetrized_hamiltonian(
        spatial_basis, discrete_basis, potential_params
    )

    # Solve with truncation
    energies, coeffs = solve_generalized_eigenvalue_truncated(H, S, epsilon)

    RETURN energies, coeffs
```

### 7.3 Basis Size Optimization

The truncation method enables larger basis sets:

```text
PROCEDURE optimize_basis_size(target_eigenvalues=10, epsilon=1e-10):
    # Start with small basis
    N = 20

    WHILE True:
        basis = generate_gaussian_basis(N)
        H, S = build_hamiltonian(basis)

        # Check effective size after truncation
        eigenvalues_S = eigenvalues(S)
        N_eff = count(eigenvalues_S > epsilon)

        IF N_eff >= target_eigenvalues:
            BREAK

        N = N * 1.5  # Increase basis size

    RETURN basis, N_eff
```

---

## 8. Alternative Approaches

### 8.1 Cholesky with Diagonal Shift

Add a small shift to $S$ to make it well-conditioned:

$$
S' = S + \delta I
$$

Then solve $H \mathbf{c} = E S' \mathbf{c}$.

**Disadvantages:**
- Less principled (introduces bias)
- Doesn't provide diagnostic information
- Shifts all eigenvalues, not just small ones

### 8.2 SVD-Based Method

Use singular value decomposition instead of eigendecomposition:

$$
S = U \Sigma V^T
$$

For symmetric $S$, this is equivalent to eigendecomposition.

### 8.3 Tikhonov Regularization

Modify the eigenvalue problem:

$$
(H + \alpha S) \mathbf{c} = (E + \alpha) S \mathbf{c}
$$

**Disadvantage:** Requires choosing regularization parameter $\alpha$.

### 8.4 Gram-Schmidt Orthogonalization

Orthogonalize the basis before computing matrix elements.

**Disadvantage:** Numerically unstable for large $N$ (the problem we're trying to solve).

### 8.5 Comparison

| Method | Principled | Stable | Diagnostic | Preserves Variational |
|--------|------------|--------|------------|----------------------|
| Eigenvalue truncation | Yes | Yes | Yes | Yes (for retained modes) |
| Diagonal shift | No | Yes | No | No |
| SVD | Yes | Yes | Yes | Yes |
| Tikhonov | No | Yes | No | No |
| Gram-Schmidt | Yes | No | No | Yes |

**Recommendation:** Eigenvalue truncation is preferred for its combination of physical interpretation, numerical stability, and diagnostic capability.

---

## 9. Test Cases

### 9.1 Small Basis (N=10): No Truncation Needed

```text
TEST small_basis_no_truncation:
    # Generate well-separated Gaussians
    widths = geometric_series(0.1, 10.0, N=10)
    H, S = build_test_hamiltonian(widths)

    # Solve with and without truncation
    E_direct = generalized_eig(H, S)
    E_trunc, _ = solve_generalized_eigenvalue_truncated(H, S, epsilon=1e-12)

    ASSERT length(E_trunc) == 10  # No truncation
    ASSERT max(abs(E_direct - E_trunc)) < 1e-10  # Results match
```

### 9.2 Medium Basis (N=100): Condition Number Check

```text
TEST medium_basis_conditioning:
    widths = geometric_series(0.05, 50.0, N=100)
    H, S = build_test_hamiltonian(widths)

    kappa_orig, kappa_trunc = analyze_conditioning(S, epsilon=1e-10)

    PRINT "Condition number improvement: {kappa_orig / kappa_trunc:.2e}"
    ASSERT kappa_trunc < 1e12  # Reasonable condition number after truncation
```

### 9.3 Large Basis (N=500): Stability Verification

```text
TEST large_basis_stability:
    widths = geometric_series(0.01, 100.0, N=500)
    H, S = build_test_hamiltonian(widths)

    # Solve with multiple thresholds
    results = {}
    FOR epsilon in [1e-6, 1e-8, 1e-10]:
        E, _ = solve_generalized_eigenvalue_truncated(H, S, epsilon)
        results[epsilon] = E[1:5]  # First 5 eigenvalues

    # Check stability of lowest eigenvalues
    FOR i = 1 TO 5:
        E_values = [results[eps][i] for eps in [1e-6, 1e-8, 1e-10]]
        variation = max(E_values) - min(E_values)
        ASSERT variation < 1e-6 * abs(E_values[0])
```

### 9.4 Extreme Ill-Conditioning: Graceful Handling

```text
TEST extreme_ill_conditioning:
    # Create intentionally ill-conditioned basis
    # Many nearly identical Gaussians
    widths = [1.0 + 0.001*i for i in 1..50]
    H, S = build_test_hamiltonian(widths)

    # Should handle gracefully
    E, coeffs = solve_generalized_eigenvalue_truncated(H, S, epsilon=1e-8)

    # Effective size should be small
    ASSERT length(E) < 50
    PRINT "Reduced from 50 to {length(E)} effective basis functions"

    # Verify orthonormality
    error = verify_orthonormality(coeffs, S)
    ASSERT error < 1e-10
```

---

## 10. Physical Interpretation

### 10.1 Meaning of Truncation

Removing small eigenvalues of $S$ corresponds to:

1. **Eliminating redundant directions**: Near-null space directions don't contribute to wave function description
2. **Projecting to physical subspace**: The effective basis spans the same physical space with fewer functions
3. **Regularizing the problem**: Making it well-posed for numerical solution

### 10.2 Effect on Variational Principle

The variational principle still holds in the reduced space:

$$
E_0 \leq E_{\text{exact}}
$$

**Caveat:** With a smaller effective basis, the variational bound may be slightly higher than with the full basis (if it could be solved stably).

### 10.3 When Truncation is Significant

Significant truncation ($N_{\text{eff}} \ll N$) indicates:

1. **Basis is too large**: Consider using fewer, better-chosen Gaussians
2. **Width distribution is suboptimal**: Gaussians are too similar
3. **Physical limitations**: The system doesn't support that many independent modes

---

## 11. References

- This algorithm is a standard technique in quantum chemistry for handling non-orthogonal basis sets
- Related to canonical orthogonalization in electronic structure theory
- For spatial matrix element construction, see `algorithm_spatial_MM.md`
- For antisymmetrization, see `algorithm_antisymmetrization.md`
- For physical background, see `physics.md`
