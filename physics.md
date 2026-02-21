# Gaussian Expansion Method for Multiquark Systems

## 1. Gaussian Expansion Method (GEM)

The Gaussian Expansion Method (GEM) is a variational approach for solving few-body quantum mechanical problems. The wave function is expanded using Gaussian basis functions:

$$
\psi(\vec{r}) = \sum_{n=1}^{N} c_n \phi_n(\vec{r}), \quad \phi_n(r) = N_n e^{-\nu_n r^2}
$$

**Key features:**
- In this work, only the $l=0$ (s-wave) basis is used
- The basis with $\nu_n$ arranged in a geometric series is highly efficient
- **All matrix elements can be calculated analytically**

### 1.1 Normalization

For the s-wave Gaussian basis:

$$
N_n = \left( \frac{2\nu_n}{\pi} \right)^{3/4}
$$

### 1.2 Analytical Matrix Elements

| Matrix Element | Formula |
|----------------|---------|
| Overlap | $\langle \phi_n \| \phi_m \rangle = N_n N_m \left( \frac{\pi}{\nu_n + \nu_m} \right)^{3/2}$ |
| Kinetic | $\langle \phi_n \| T \| \phi_m \rangle = \frac{(\hbar c)^2}{2\mu} \cdot \frac{6\nu_n \nu_m}{\nu_n + \nu_m} \langle \phi_n \| \phi_m \rangle$ |
| Coulomb | $\langle \phi_n \| 1/r \| \phi_m \rangle = N_n N_m \frac{2\pi}{\nu_n + \nu_m}$ |
| Linear | $\langle \phi_n \| r \| \phi_m \rangle = N_n N_m \frac{2\pi}{(\nu_n + \nu_m)^2}$ |
| Gaussian | $\langle \phi_n \| e^{-\alpha r^2} \| \phi_m \rangle = N_n N_m \left( \frac{\pi}{\nu_n + \nu_m + \alpha} \right)^{3/2}$ |

where $(\hbar c)^2 = 0.0389$ GeV²·fm² and $\mu$ is in GeV.

---

## 2. Quark Potential Model (AL1)

The quark model includes three interaction terms between quark pairs:

### 2.1 Coulomb Potential (One-Gluon Exchange)

$$
V_{\text{Coulomb}}(r) = \frac{\alpha_s^{\text{coul}}}{r} \cdot \frac{\boldsymbol{\lambda}_i}{2} \cdot \frac{\boldsymbol{\lambda}_j}{2}
$$

### 2.2 Confinement Potential

$$
V_{\text{conf}}(r) = \left( -\frac{3b}{4} r + V_c \right) \cdot \frac{\boldsymbol{\lambda}_i}{2} \cdot \frac{\boldsymbol{\lambda}_j}{2}
$$

### 2.3 Hyperfine Potential (Spin-Spin Interaction)

$$
V_{\text{hyp}}(r) = -\frac{8\pi\alpha_s^{\text{hyp}}}{3m_i m_j} \frac{\tau^3}{\pi^{3/2}} e^{-\tau^2 r^2} \cdot \frac{\boldsymbol{\lambda}_i}{2} \cdot \frac{\boldsymbol{\lambda}_j}{2} \cdot \mathbf{s}_i \cdot \mathbf{s}_j
$$

where $\tau = (2\mu_{ij})^{b_{exp}}/a$ and $\mu_{ij} = m_i m_j/(m_i + m_j)$.

### 2.4 Model Parameters

**Quark Masses (GeV):**

| u, d | s | c | b |
|------|---|---|---|
| 0.315 | 0.577 | 1.836 | 5.227 |

**Potential Parameters:**

| Parameter | Value | Description |
|-----------|-------|-------------|
| $\alpha_s^{\text{coul}}$ | 0.380175 | Coulomb coupling |
| $\alpha_s^{\text{hyp}}$ | 1.39568 | Hyperfine coupling |
| $a$ | 1.6553 | Hyperfine smearing |
| $b_{exp}$ | 0.2204 | Hyperfine mass exponent |
| $b$ | 0.1653 | String tension (GeV²) |
| $V_c$ | 0.624075 | Confinement constant (GeV) |

### 2.5 Color Operators and Antiquarks

**Quark vs Antiquark Color Generators:**

Quarks transform under the fundamental representation **3** of SU(3)_color:
$$T_a = \frac{\lambda_a}{2}$$

Antiquarks transform under the conjugate representation **3̄**:
$$\bar{T}_a = -\frac{\lambda_a^*}{2} = -\frac{\lambda_a^T}{2}$$

For Gell-Mann matrices (Hermitian: $\lambda^\dagger = \lambda$):
- Real symmetric ($\lambda_1, \lambda_3, \lambda_4, \lambda_6, \lambda_8$): $\lambda^* = \lambda$, so $\bar{T} = -T$
- Imaginary antisymmetric ($\lambda_2, \lambda_5, \lambda_7$): $\lambda^* = -\lambda$, so $\bar{T} = +T$

**Pairwise Color Operator:**

The color-color operator depends on the particle types:

| Pair | Operator | Formula |
|------|----------|---------|
| q-q | $\sum_a T_a^{(1)} T_a^{(2)}$ | $\sum_a \frac{\lambda_a}{2} \otimes \frac{\lambda_a}{2}$ |
| q-q̄ | $\sum_a T_a^{(1)} \bar{T}_a^{(2)}$ | $-\sum_a \frac{\lambda_a}{2} \otimes \frac{\lambda_a^T}{2}$ |
| q̄-q̄ | $\sum_a \bar{T}_a^{(1)} \bar{T}_a^{(2)}$ | $\sum_a \frac{\lambda_a^T}{2} \otimes \frac{\lambda_a^T}{2}$ |

**Casimir Eigenvalues:**

| Representation | C₂ eigenvalue |
|----------------|---------------|
| **1** (singlet) | 0 |
| **3** (quark) | 4/3 |
| **3̄** (antiquark) | 4/3 |
| **6** (symmetric qq) | 10/3 |
| **3̄** (antisymmetric qq) | 4/3 |
| **8** (octet) | 3 |

### 2.6 Isospin (Flavor SU(2))

Isospin is an approximate SU(2) symmetry between u and d quarks. Only u and d quarks carry isospin; heavier quarks (s, c, b) are isospin singlets.

**Quark Isospin Quantum Numbers:**

| Quark | $I$ | $I_3$ | Hilbert space dim |
|-------|-----|-------|-------------------|
| u | 1/2 | +1/2 | 2 |
| d | 1/2 | -1/2 | 2 |
| s | 0 | 0 | 1 (trivial) |
| c | 0 | 0 | 1 (trivial) |
| b | 0 | 0 | 1 (trivial) |

**Antiquark Isospin:**

| Antiquark | $I$ | $I_3$ |
|-----------|-----|-------|
| ū | 1/2 | -1/2 |
| d̄ | 1/2 | +1/2 |
| s̄, c̄, b̄ | 0 | 0 |

**Isospin Operators:**

$$I_a = \frac{\tau_a}{2} \quad \text{(for u, d quarks)}$$
$$I_a = 0 \quad \text{(for s, c, b quarks)}$$

where $\tau_a$ are the Pauli matrices.

**Total Isospin Space Dimension:**

$$\dim(\text{isospin}) = \prod_k d_k, \quad d_k = \begin{cases} 2 & \text{if quark } k \text{ is u or d} \\ 1 & \text{if quark } k \text{ is s, c, or b} \end{cases}$$

**Examples:**

| Hadron | Quarks | Isospin dim | Possible $I$ |
|--------|--------|-------------|--------------|
| $\pi$ | $u\bar{d}$, $d\bar{u}$ | 4 | 0, 1 |
| $K$ | $u\bar{s}$, $d\bar{s}$ | 2 | 1/2 |
| $J/\psi$ | $c\bar{c}$ | 1 | 0 (trivial) |
| $\Lambda$ | $uds$ | 4 | 0, 1 |
| $\Sigma$ | $uus$, $uds$, $dds$ | 4 | 0, 1 |
| $\Xi$ | $uss$, $dss$ | 2 | 1/2 |
| $\Omega$ | $sss$ | 1 | 0 (trivial) |

**Isospin-Isospin Operator:**

For pairs where both quarks have $I=1/2$:
$$\langle \mathbf{I}_i \cdot \mathbf{I}_j \rangle = \frac{1}{2}[I_{ij}(I_{ij}+1) - \frac{3}{2}]$$

| Pair isospin | $\langle \mathbf{I}_i \cdot \mathbf{I}_j \rangle$ |
|--------------|---------------------------------------------------|
| $I_{ij}=0$ (singlet) | $-3/4$ |
| $I_{ij}=1$ (triplet) | $+1/4$ |

For pairs involving s, c, or b quarks: $\langle \mathbf{I}_i \cdot \mathbf{I}_j \rangle = 0$

### 2.7 Color and Spin Factors

**Meson ($q\bar{q}$):**

- Color singlet: $\langle \frac{\boldsymbol{\lambda}_1}{2} \cdot \frac{\bar{\boldsymbol{\lambda}}_2}{2} \rangle = -\frac{4}{3}$
- Spin: $\langle \mathbf{s}_1 \cdot \mathbf{s}_2 \rangle = -\frac{3}{4}$ (S=0) or $+\frac{1}{4}$ (S=1)

**Baryon ($qqq$):**

- Color singlet: $\langle \frac{\boldsymbol{\lambda}_i}{2} \cdot \frac{\boldsymbol{\lambda}_j}{2} \rangle = -\frac{2}{3}$ for each pair

**Tetraquark ($qq\bar{q}\bar{q}$):**

- Two color singlet configurations exist (diquark-antidiquark and meson-meson)
- Matrix elements depend on pair types (q-q, q-q̄, q̄-q̄)

**Baryon spin factors** (sum rule: $\sum_{i<j} \langle \mathbf{s}_i \cdot \mathbf{s}_j \rangle = \frac{S(S+1) - 9/4}{2}$):

| Baryon Type | S | $\langle \mathbf{s}_1 \cdot \mathbf{s}_2 \rangle$ | $\langle \mathbf{s}_1 \cdot \mathbf{s}_3 \rangle$ | $\langle \mathbf{s}_2 \cdot \mathbf{s}_3 \rangle$ |
|-------------|---|-------|-------|-------|
| Λ-type (e.g., uds, udc) | 1/2 | $-3/4$ | $0$ | $0$ |
| Σ-type (e.g., uuc, ddc) | 1/2 | $+1/4$ | $-1/2$ | $-1/2$ |
| Decuplet (e.g., Δ, Ω) | 3/2 | $+1/4$ | $+1/4$ | $+1/4$ |

- **Λ-type**: All different flavors → pair 12 in spin singlet ($S_{12}=0$)
- **Σ-type**: Identical quarks in pair 12 → Fermi statistics requires spin triplet ($S_{12}=1$)

### 2.8 Three-Body Force (Baryons Only)

A phenomenological three-body force is introduced to mimic genuine three-body effects:

$$
V_{3\text{body}} = -\frac{C_3}{m_1 m_2 m_3}
$$

where $C_3 = 2.02 \times 10^{-3}$ GeV$^4$.

**Key properties:**
- Constant potential (no spatial dependence)
- No color, spin, or isospin factor dependence
- Only applies to baryons (three-quark systems)
- Stronger effect for light quarks (inverse mass scaling)

**Matrix element:** $\langle \Phi | V_{3\text{body}} | \Phi' \rangle = V_{3\text{body}} \cdot S_{\Phi\Phi'}$

---

## 3. Jacobi Coordinates

### 3.1 Three-Body System

For particles 1, 2, 3 with masses $m_1, m_2, m_3$, choose clustering (12)3:

$$
\boldsymbol{\rho} = \mathbf{r}_1 - \mathbf{r}_2, \quad \boldsymbol{\lambda} = \frac{m_1 \mathbf{r}_1 + m_2 \mathbf{r}_2}{m_1 + m_2} - \mathbf{r}_3
$$

**Reduced masses:**

$$
\mu_\rho = \frac{m_1 m_2}{m_1 + m_2}, \quad \mu_\lambda = \frac{(m_1 + m_2) m_3}{M}
$$

**Interparticle distances:**

$$
r_{12} = |\boldsymbol{\rho}| \quad \text{(diagonal)}
$$

$$
r_{13} = \left| \frac{m_2}{m_1+m_2} \boldsymbol{\rho} + \boldsymbol{\lambda} \right| \quad \text{(non-diagonal)}
$$

$$
r_{23} = \left| -\frac{m_1}{m_1+m_2} \boldsymbol{\rho} + \boldsymbol{\lambda} \right| \quad \text{(non-diagonal)}
$$

### 3.2 N-Body Generalization

| System | Jacobi coords | Pairs | Basis size ($N_b=15$) |
|--------|---------------|-------|----------------------|
| Meson | 1 | 1 | 15 |
| Baryon | 2 | 3 | 225 |
| Tetraquark | 3 | 6 | 3,375 |
| Pentaquark | 4 | 10 | 50,625 |

---

## 4. Matrix Elements via Cholesky Decomposition

For non-diagonal pair potentials (e.g., $V_{13}$ in the (12)3 basis), we use the **Cholesky decomposition method**.

### 4.1 The Problem

The potential $V(r_{13})$ depends on both $\boldsymbol{\rho}$ and $\boldsymbol{\lambda}$:

$$
r_{13}^2 = a^2 \rho^2 + \lambda^2 + 2a \rho \lambda \cos\theta
$$

where $a = m_2/(m_1+m_2)$ and $\theta$ is the angle between the vectors.

### 4.2 The Solution: Coordinate Transformation + Cholesky

**Step 1:** Transform to Jacobi set (13)2 where $r_{13} = |\boldsymbol{\rho}_2|$ is diagonal.

The transformation from set (12)3 to set (13)2 in **mass-scaled coordinates** is an orthogonal rotation:

$$
\begin{pmatrix} \boldsymbol{\xi}_{\rho,2} \\ \boldsymbol{\xi}_{\lambda,2} \end{pmatrix} =
\begin{pmatrix} \cos\varphi & \sin\varphi \\ -\sin\varphi & \cos\varphi \end{pmatrix}
\begin{pmatrix} \boldsymbol{\xi}_{\rho,3} \\ \boldsymbol{\xi}_{\lambda,3} \end{pmatrix}
$$

where $\boldsymbol{\xi}_\rho = \sqrt{\mu_\rho} \boldsymbol{\rho}$ and $\boldsymbol{\xi}_\lambda = \sqrt{\mu_\lambda} \boldsymbol{\lambda}$.

**Transformation angle (general masses):**

$$
\cos\varphi = -\frac{m_2}{\sqrt{M_{12} M_{13}}}, \quad \sin\varphi = \frac{\sqrt{M_{12} M_{13} - m_2^2}}{\sqrt{M_{12} M_{13}}}
$$

where $M_{12} = m_1 + m_2$ and $M_{13} = m_1 + m_3$.

Note: $M_{12} M_{13} - m_2^2 = m_1^2 + m_1 m_2 + m_1 m_3 + m_2 m_3 - m_2^2 = m_1 M + m_2(m_3 - m_2)$.

**Transformation in original (non-mass-scaled) coordinates:**

$$
\begin{pmatrix} \boldsymbol{\rho}_2 \\ \boldsymbol{\lambda}_2 \end{pmatrix} = \mathbf{T}
\begin{pmatrix} \boldsymbol{\rho}_3 \\ \boldsymbol{\lambda}_3 \end{pmatrix}, \quad
\mathbf{T} = \begin{pmatrix} T_{11} & T_{12} \\ T_{21} & T_{22} \end{pmatrix}
$$

with:

$$
T_{11} = \sqrt{\frac{\mu_{\rho,3}}{\mu_{\rho,2}}} \cos\varphi, \quad
T_{12} = \sqrt{\frac{\mu_{\lambda,3}}{\mu_{\rho,2}}} \sin\varphi
$$

$$
T_{21} = -\sqrt{\frac{\mu_{\rho,3}}{\mu_{\lambda,2}}} \sin\varphi, \quad
T_{22} = \sqrt{\frac{\mu_{\lambda,3}}{\mu_{\lambda,2}}} \cos\varphi
$$

**Reduced masses:**

$$
\mu_{\rho,3} = \frac{m_1 m_2}{M_{12}}, \quad \mu_{\lambda,3} = \frac{M_{12} m_3}{M}
$$

$$
\mu_{\rho,2} = \frac{m_1 m_3}{M_{13}}, \quad \mu_{\lambda,2} = \frac{M_{13} m_2}{M}
$$

where $M = m_1 + m_2 + m_3$.

**Equal mass simplification** ($m_1 = m_2 = m_3 = m$): $\varphi = 2\pi/3$, giving:

$$
\mathbf{T}_{3\to2} = \begin{pmatrix} -1/2 & \sqrt{3}/2 \\ -\sqrt{3}/2 & -1/2 \end{pmatrix}
$$

**For V₂₃:** Transform from (12)3 to (23)1 using:

$$
\cos\varphi' = -\frac{m_1}{\sqrt{M_{12} M_{23}}}, \quad \sin\varphi' = \frac{\sqrt{M_{12} M_{23} - m_1^2}}{\sqrt{M_{12} M_{23}}}
$$

with $M_{23} = m_2 + m_3$, and analogous reduced mass formulas.

**Step 2:** Compute the width matrix in new coordinates.

$$
\mathbf{A} = (\mathbf{T}^{-1})^T \begin{pmatrix} \nu_n + \nu_{n'} & 0 \\ 0 & \nu_m + \nu_{m'} \end{pmatrix} \mathbf{T}^{-1}
$$

**Step 3:** Cholesky decomposition $\mathbf{A} = \mathbf{L} \mathbf{L}^T$ with lower triangular $\mathbf{L}$:

$$
l_{11} = \sqrt{A_{11}}, \quad l_{21} = \frac{A_{12}}{l_{11}}, \quad l_{22} = \sqrt{A_{22} - l_{21}^2}
$$

**Step 4:** The key insight — with the transformation $\mathbf{y} = \mathbf{L} \mathbf{r}$:

$$
\boxed{\boldsymbol{\rho}_2 = \frac{\mathbf{y}_\rho}{l_{11}}} \quad \text{(only rescaled, no mixing!)}
$$

This means the potential variable is **decoupled** from the other coordinate.

### 4.3 Final Formula

The matrix element factorizes into independent integrals:

$$
\boxed{\langle \Phi_{nm} | V | \Phi_{n'm'} \rangle = \frac{N_n N_m N_{n'} N_{m'} \cdot \pi^{3/2}}{(l_{11} l_{22})^3} \cdot I_V}
$$

### 4.4 Potential Integral Table

| Potential | $V(r)$ | $I_V$ |
|-----------|--------|-------|
| Coulomb | $1/r$ | $2\pi l_{11}$ |
| Linear | $r$ | $2\pi/l_{11}$ |
| Constant | $1$ | $\pi^{3/2}$ |
| Gaussian | $e^{-\tau^2 r^2}$ | $\left(\frac{\pi l_{11}^2}{l_{11}^2 + \tau^2}\right)^{3/2}$ |

**Derivation of integrals:** $I_V = 4\pi \int_0^\infty r^2 e^{-r^2} V(r/l_{11}) dr$
- Coulomb: $\int_0^\infty r e^{-r^2} dr = 1/2$, so $I_{\text{coul}} = 4\pi l_{11} \times 1/2 = 2\pi l_{11}$
- Linear: $\int_0^\infty r^3 e^{-r^2} dr = 1/2$, so $I_{\text{lin}} = (4\pi/l_{11}) \times 1/2 = 2\pi/l_{11}$

### 4.5 Algorithm

```
For each basis pair (n,m) and (n',m'):
  1. D = diag(ν_n + ν_n', ν_m + ν_m')
  2. A = (T⁻¹)ᵀ D T⁻¹
  3. l₁₁ = √A₁₁, l₂₁ = A₁₂/l₁₁, l₂₂ = √(A₂₂ - l₂₁²)
  4. I_V from table using l₁₁
  5. Result = N_bra × N_ket × π^(3/2) × I_V / (l₁₁ l₂₂)³
```

### 4.6 Advantages

1. **Fully analytical** — no numerical integration
2. **Efficient** — O(1) per matrix element
3. **Numerically stable** — Cholesky is well-conditioned
4. **Generalizable** — extends to N-body with (N-1)×(N-1) Cholesky

---

## 5. Hamiltonian Construction

### 5.1 Two-Body (Meson)

$$
H = T + V_{12}, \quad \text{solve } H\mathbf{c} = E \, S \, \mathbf{c}
$$

### 5.2 Three-Body (Baryon)

**Basis:** Product of Gaussians $\Phi_{nm}(\boldsymbol{\rho}, \boldsymbol{\lambda}) = \phi_n(\rho) \phi_m(\lambda)$

**Kinetic energy:**

$$
T = T_\rho \otimes S_\lambda + S_\rho \otimes T_\lambda
$$

**Potential energy:**

$$
V = V_{12} + V_{13} + V_{23}
$$

- $V_{12}$: diagonal in $\rho$ → $V_{12} = V_\rho \otimes S_\lambda$
- $V_{13}, V_{23}$: use Cholesky method (Section 4)

**Eigenvalue problem:**

$$
(T + V) \mathbf{c} = E \, S \, \mathbf{c}
$$

### 5.3 N-Body Systems

For tetraquarks and pentaquarks:

1. Choose a reference Jacobi tree that makes the most pairs diagonal
2. For diagonal pairs: $V_{ij} = V_k \otimes S_{\text{other}}$
3. For non-diagonal pairs: apply Cholesky method with the appropriate transformation matrix
4. The transformation matrix is an $(N-1) \times (N-1)$ orthogonal rotation in mass-scaled coordinates
