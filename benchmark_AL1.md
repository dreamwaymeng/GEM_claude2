# Benchmark Results for AL1 Potential

This document contains benchmark results for validating the GEM implementation using the AL1 quark potential model.

**Primary Reference:**
- L. Meng et al., Phys. Rev. D 108, 114016 (2023) [arXiv:2310.13354](https://arxiv.org/abs/2310.13354)

**Validation Tolerance:** All results should match within **< 1 MeV** of the reference values.

---

## 1. AL1 Potential Parameters

### 1.1 Quark Masses (GeV)

| Quark | Mass (GeV) |
|-------|------------|
| u, d  | 0.315      |
| s     | 0.577      |
| c     | 1.836      |
| b     | 5.227      |

### 1.2 Potential Parameters

| Parameter | Symbol | Value | Description |
|-----------|--------|-------|-------------|
| Coulomb coupling | $\kappa$ | 0.5069 | One-gluon-exchange strength |
| Hyperfine coupling | $\kappa'$ | 1.8609 | Spin-spin interaction strength |
| String tension | $\lambda$ | 0.1653 GeV² | Linear confinement slope |
| Smearing parameter | $a$ | 1.6553 | Hyperfine regularization |
| Hyperfine exponent | $b_{exp}$ | 0.2204 | Mass-dependent smearing |
| Confinement constant | $\Lambda$ | 0.8321 GeV | Confinement offset |

### 1.3 Hyperfine Smearing

The hyperfine regularization parameter $\tau$ is mass-dependent:

$$
\tau = \frac{(2\mu_{ij})^{b_{exp}}}{a}
$$

where $\mu_{ij} = m_i m_j / (m_i + m_j)$ is the reduced mass.

---

## 2. Meson Benchmark Masses (MeV)

All values from Table I of arXiv:2310.13354.

### 2.1 Light and Strange Mesons

| Meson | $J^P$ | Quark Content | AL1 (MeV) |
|-------|-------|---------------|-----------|
| $\pi$ | $0^-$ | $n\bar{n}$    | 138.16    |
| $\rho$| $1^-$ | $n\bar{n}$    | 767.00    |
| $K$   | $0^-$ | $n\bar{s}$    | 490.92    |
| $K^*$ | $1^-$ | $n\bar{s}$    | 903.55    |

### 2.2 Charmed Mesons

| Meson | $J^P$ | Quark Content | AL1 (MeV) |
|-------|-------|---------------|-----------|
| $D$   | $0^-$ | $c\bar{n}$    | 1862.4    |
| $D^*$ | $1^-$ | $c\bar{n}$    | 2016.1    |
| $D_s$ | $0^-$ | $c\bar{s}$    | 1962.5    |
| $D_s^*$| $1^-$| $c\bar{s}$    | 2102.0    |

### 2.3 Bottom Mesons

| Meson | $J^P$ | Quark Content | AL1 (MeV) |
|-------|-------|---------------|-----------|
| $B$   | $0^-$ | $b\bar{n}$    | 5293.5    |
| $B^*$ | $1^-$ | $b\bar{n}$    | 5350.5    |
| $B_s$ | $0^-$ | $b\bar{s}$    | 5361.0    |
| $B_s^*$| $1^-$| $b\bar{s}$    | 5417.5    |
| $B_c$ | $0^-$ | $b\bar{c}$    | 6291.6    |
| $B_c^*$| $1^-$| $b\bar{c}$    | 6343.2    |

### 2.4 Heavy Quarkonia

| Meson | $J^P$ | Quark Content | AL1 (MeV) |
|-------|-------|---------------|-----------|
| $\eta_c$ | $0^-$ | $c\bar{c}$ | 3005.3    |
| $J/\psi$ | $1^-$ | $c\bar{c}$ | 3101.3    |
| $\eta_b$ | $0^-$ | $b\bar{b}$ | 9423.7    |
| $\Upsilon$| $1^-$| $b\bar{b}$ | 9461.5    |

### 2.5 Hyperfine Splittings

| System | $\Delta M = M(1^-) - M(0^-)$ (MeV) |
|--------|-----------------------------------|
| $n\bar{n}$ | $\rho - \pi = 628.84$         |
| $n\bar{s}$ | $K^* - K = 412.63$            |
| $c\bar{n}$ | $D^* - D = 153.7$             |
| $c\bar{s}$ | $D_s^* - D_s = 139.5$         |
| $b\bar{n}$ | $B^* - B = 57.0$              |
| $b\bar{s}$ | $B_s^* - B_s = 56.5$          |
| $c\bar{c}$ | $J/\psi - \eta_c = 96.0$      |
| $b\bar{b}$ | $\Upsilon - \eta_b = 37.8$    |

---

## 3. Tetraquark Benchmark Results

All values from Tables III-IV of arXiv:2310.13354 (AL1-GEM column).

### 3.1 Dimeson Thresholds (MeV)

| Threshold | Composition | Value (MeV) |
|-----------|-------------|-------------|
| $DD^*$    | $D + D^*$   | 3878.5      |
| $\bar{B}\bar{B}^*$ | $B + B^*$ | 10644.0   |
| $D\bar{B}$ | $D + B$    | 7155.9      |
| $D\bar{B}^*$ | $D + B^*$ | 7212.9     |
| $D^*\bar{B}^*$ | $D^* + B^*$ | 7366.6   |
| $\bar{B}_s\bar{B}^*$ | $B_s + B^*$ | 10711.5 |

### 3.2 Doubly Charmed Tetraquark $T_{cc}$ ($[cc\bar{n}\bar{n}]^{I=0}$, $J^P=1^+$)

| Quantity | Value |
|----------|-------|
| Threshold ($DD^*$) | 3878.5 MeV |
| Binding Energy | **-14.0 MeV** |
| Mass | 3864.5 MeV |

### 3.3 Doubly Bottom Tetraquark $T_{bb}$ ($[bb\bar{n}\bar{n}]^{I=0}$, $J^P=1^+$)

| State | Threshold ($\bar{B}\bar{B}^*$) | Binding Energy | Mass |
|-------|-------------------------------|----------------|------|
| Ground | 10644.0 MeV | **-151.6 MeV** | 10492.4 MeV |
| Excited | 10644.0 MeV | **-0.70 MeV** | 10643.3 MeV |

### 3.4 $bc\bar{n}\bar{n}$ Tetraquarks ($I=0$)

| $J^P$ | Threshold | Binding Energy | Mass |
|-------|-----------|----------------|------|
| $0^+$ | $D\bar{B}$ = 7155.9 MeV | **-26.0 MeV** | 7129.9 MeV |
| $1^+$ | $D\bar{B}^*$ = 7212.9 MeV | **-26.5 MeV** | 7186.4 MeV |
| $2^+$ | $D^*\bar{B}^*$ = 7366.6 MeV | **-2.9 MeV** | 7363.7 MeV |

### 3.5 $bb\bar{n}\bar{s}$ Tetraquark ($J^P=1^+$)

| Quantity | Value |
|----------|-------|
| Threshold ($\bar{B}_s\bar{B}^*$) | 10711.5 MeV |
| Binding Energy | **-63.8 MeV** |
| Mass | 10647.7 MeV |

---

## 4. GEM Basis Parameters

From Section II.A of arXiv:2310.13354.

### 4.1 Gaussian Basis in Coordinate Space

The basis functions use the form:
$$
\phi_{nlm}(\bm{r}) = \sqrt{\frac{2^{l+5/2}}{\Gamma(l+\frac{3}{2})r_n^3}} \left(\frac{r}{r_n}\right)^l e^{-r^2/r_n^2} Y_{lm}(\hat{r})
$$

where $r_n$ is taken in geometric progression: $r_n = r_0 \cdot a^{n-1}$.

### 4.2 Tetraquark Basis Parameters

| Coordinate Type | $r_0$ (fm) | $r_{max}$ (fm) | $N$ |
|-----------------|------------|----------------|-----|
| $q$-$q$ or $\bar{q}$-$\bar{q}$ | 0.1 | 2 | 6 |
| $(qq)$-$(\bar{q}\bar{q})$ | 0.1 | 2 | 6 |
| $q$-$\bar{q}$ | 0.1 | 1 | 6 |
| $(q\bar{q})$-$(q\bar{q})$ | 0.1 | 5 | 6 |

Three Jacobi coordinate sets are used (diquark-type, dimeson-type 1, dimeson-type 2).

---

## 5. Validation Checklist

**Tolerance: < 1 MeV deviation from reference values**

### 5.1 Phase 1: Meson Validation

- [ ] $\pi$ mass = 138.16 MeV
- [ ] $\rho$ mass = 767.00 MeV
- [ ] $K$ mass = 490.92 MeV
- [ ] $K^*$ mass = 903.55 MeV
- [ ] $D$ mass = 1862.4 MeV
- [ ] $D^*$ mass = 2016.1 MeV
- [ ] $J/\psi$ mass = 3101.3 MeV
- [ ] $\Upsilon$ mass = 9461.5 MeV
- [ ] $\rho - \pi$ splitting = 628.84 MeV
- [ ] $J/\psi - \eta_c$ splitting = 96.0 MeV

### 5.2 Phase 2: Tetraquark Validation

- [ ] $T_{cc}$ binding energy = -14.0 MeV
- [ ] $T_{bb}$ binding energy = -151.6 MeV
- [ ] $T_{bb}$ excited binding = -0.70 MeV
- [ ] $[bc\bar{n}\bar{n}]^{I=0}_{0^+}$ binding = -26.0 MeV
- [ ] $[bc\bar{n}\bar{n}]^{I=0}_{1^+}$ binding = -26.5 MeV
- [ ] $[bc\bar{n}\bar{n}]^{I=0}_{2^+}$ binding = -2.9 MeV

---

## 6. Reference

L. Meng, Y.-K. Chen, Y. Ma, and S.-L. Zhu, "Tetraquark bound states in constituent quark models: benchmark test calculations," Phys. Rev. D 108, 114016 (2023) [arXiv:2310.13354](https://arxiv.org/abs/2310.13354)
