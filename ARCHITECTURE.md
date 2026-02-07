# GEM Program Architecture Design

## Overview

Design the software architecture for a GEM (Gaussian Expansion Method) hadron physics program in C++, balancing flexibility, efficiency, and type safety across three key areas:

1. Input settings
2. Parallelization of matrix element calculations
3. Eigenvalue solver package selection

---

## 1. Deployment Strategy

| System | Environment | Parallelization | Memory Model |
|--------|-------------|-----------------|--------------|
| Meson (N=2) | Personal Computer | OpenMP | Dense, in-memory |
| Baryon (N=3) | Personal Computer | OpenMP | Dense, in-memory |
| Tetraquark (N=4) | Personal Computer | OpenMP | Dense, in-memory |
| **Pentaquark (N=5)** | **HPC Cluster** | **MPI + OpenMP** | **Distributed** |

### Personal Computer Requirements (Meson → Tetraquark)

- **CPU**: 8-16 cores recommended (Intel/AMD with AVX2/AVX-512)
- **RAM**: 16-64 GB (tetraquark needs ~8-16 GB)
- **Compiler**: GCC 11+, Clang 14+, or MSVC 2022
- **Build**: CMake 3.20+
- **OS**: macOS, Linux, or Windows

### HPC Cluster Requirements (Pentaquark)

- **Nodes**: 10-100 nodes depending on basis size
- **Cores per node**: 32-128 (Intel Xeon, AMD EPYC)
- **RAM per node**: 128-256 GB
- **Interconnect**: InfiniBand HDR/EDR
- **Compiler**: Intel oneAPI, GCC with MPI wrappers
- **MPI**: Intel MPI, OpenMPI, or MPICH
- **Job scheduler**: SLURM (primary), PBS/Torque (alternative)

---

## 2. Language Choice Rationale

**C++ (C++17/20) with Eigen + LAPACK**

| Criterion | Rating | Notes |
|-----------|--------|-------|
| Runtime performance | Excellent | 10-100x faster than interpreted languages |
| Compile-time safety | Excellent | Template errors caught at compile time |
| Memory control | Excellent | Critical for large matrices |
| HPC integration | Excellent | Native MPI, OpenMP, LAPACK, ELPA |
| Template metaprogramming | Excellent | N-particle specialization |

### Library Roles

| Library | Role | Used For |
|---------|------|----------|
| **Eigen** | Matrix storage & operations | Matrix creation, arithmetic, transformations |
| **LAPACK** (MKL/OpenBLAS) | Eigenvalue solving | Generalized eigenvalue problem H·x = E·S·x |
| **Spectra** | Iterative eigensolving | Large matrices, few eigenvalues |
| **ELPA** | Distributed eigensolving | HPC pentaquark calculations |

**Key benefits:**
- **Eigen**: Expression templates for zero-overhead matrix operations, clean API
- **LAPACK**: Industry-standard eigensolvers, highly optimized (MKL/OpenBLAS)
- OpenMP/MPI provide mature parallelization with minimal overhead
- Template specialization for different hadron types (compile-time optimization)
- Strong type system catches physics errors at compile time
- Single codebase for both personal computer and HPC

**Note**: Eigen is used only for matrix storage and operations (addition, multiplication, etc.), not for eigenvalue solving. LAPACK provides faster and more robust eigensolvers.

---

## 3. Build System and Dependencies

### CMakeLists.txt

```cmake
cmake_minimum_required(VERSION 3.20)
project(GEM VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Build options
option(GEM_ENABLE_MPI "Enable MPI for HPC cluster (pentaquark)" OFF)
option(GEM_ENABLE_MKL "Use Intel MKL for BLAS/LAPACK (recommended)" OFF)
option(GEM_ENABLE_OPENBLAS "Use OpenBLAS for BLAS/LAPACK" OFF)
option(GEM_BUILD_TESTS "Build unit tests" ON)

# Core dependencies (always required)
find_package(Eigen3 3.4 REQUIRED)
find_package(yaml-cpp REQUIRED)
find_package(OpenMP REQUIRED)
find_package(HDF5 COMPONENTS CXX REQUIRED)

# Eigenvalue solver backends
# Priority: MKL > OpenBLAS > Eigen (fallback)

# Option 1: Intel MKL (best performance on Intel CPUs)
if(GEM_ENABLE_MKL)
    find_package(MKL CONFIG REQUIRED)
    target_compile_definitions(gem PRIVATE USE_MKL USE_LAPACK)
    target_link_libraries(gem PRIVATE MKL::MKL)
    message(STATUS "Using Intel MKL for BLAS/LAPACK")

# Option 2: OpenBLAS (good open-source alternative)
elseif(GEM_ENABLE_OPENBLAS)
    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)
    target_compile_definitions(gem PRIVATE USE_OPENBLAS USE_LAPACK)
    target_link_libraries(gem PRIVATE ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
    message(STATUS "Using OpenBLAS for BLAS/LAPACK")

# Option 3: Eigen only (fallback, slower for large matrices)
else()
    message(STATUS "Using Eigen built-in solver (consider MKL/OpenBLAS for better performance)")
endif()

# Optional: Spectra for iterative eigenvalue problems (few eigenvalues of large matrices)
find_package(Spectra)
if(Spectra_FOUND)
    target_compile_definitions(gem PRIVATE USE_SPECTRA)
    message(STATUS "Spectra found: iterative solver enabled")
endif()

# HPC/MPI support (for pentaquark on clusters)
if(GEM_ENABLE_MPI)
    find_package(MPI REQUIRED)
    target_compile_definitions(gem PRIVATE USE_MPI)
    target_link_libraries(gem PRIVATE MPI::MPI_CXX)

    # Option A: ELPA (best for large dense matrices on HPC)
    find_package(ELPA)
    if(ELPA_FOUND)
        target_compile_definitions(gem PRIVATE USE_ELPA)
        target_link_libraries(gem PRIVATE ${ELPA_LIBRARIES})
        message(STATUS "ELPA found: distributed eigensolver enabled")
    endif()

    # Option B: ScaLAPACK (standard distributed solver)
    find_package(ScaLAPACK)
    if(ScaLAPACK_FOUND)
        target_compile_definitions(gem PRIVATE USE_SCALAPACK)
        target_link_libraries(gem PRIVATE ${ScaLAPACK_LIBRARIES})
        message(STATUS "ScaLAPACK found: distributed eigensolver enabled")
    endif()
endif()

# Main library
add_library(gem
    src/config.cpp
    src/constants.cpp
    src/solver/lapack_solver.cpp
    # ... other source files
)

target_include_directories(gem PUBLIC include)
target_link_libraries(gem PUBLIC
    Eigen3::Eigen
    yaml-cpp
    OpenMP::OpenMP_CXX
    ${HDF5_LIBRARIES}
)

# Main executable
add_executable(gem_run apps/gem_run.cpp)
target_link_libraries(gem_run PRIVATE gem)

# Build type specific flags
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    target_compile_options(gem PRIVATE
        $<$<CXX_COMPILER_ID:GNU>:-O3 -march=native -ffast-math>
        $<$<CXX_COMPILER_ID:Clang>:-O3 -march=native -ffast-math>
        $<$<CXX_COMPILER_ID:Intel>:-O3 -xHost -ipo>
    )
endif()

# Tests
if(GEM_BUILD_TESTS)
    find_package(Catch2 3 REQUIRED)
    add_executable(gem_tests
        tests/test_jacobi.cpp
        tests/test_spin.cpp
        tests/test_color.cpp
        tests/test_hamiltonian.cpp
        tests/test_solver.cpp
    )
    target_link_libraries(gem_tests PRIVATE gem Catch2::Catch2WithMain)
endif()
```

### Build Commands

```bash
# Personal Computer (meson/baryon/tetraquark)
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)

# HPC Cluster (pentaquark with MPI)
module load cmake gcc openmpi mkl
cmake -B build -DCMAKE_BUILD_TYPE=Release \
      -DGEM_ENABLE_MPI=ON \
      -DGEM_ENABLE_MKL=ON
cmake --build build -j$(nproc)
```

### Dependencies

#### Core (All Systems)

| Library | Purpose | Alternative |
|---------|---------|-------------|
| **Eigen 3.4+** | Matrix storage, operations | Armadillo, Blaze |
| **yaml-cpp** | Configuration parsing | toml++, nlohmann/json |
| **OpenMP** | Shared-memory parallelization | TBB, std::execution |
| **HDF5** | Output file format | HighFive (C++ wrapper) |
| **spdlog** | Logging | - |
| **Catch2** | Unit testing | GoogleTest, doctest |

#### Eigenvalue Solver Backends (Simplified)

```
┌─────────────────────────────────────────────────────────────────┐
│                    Matrix Storage Layer                          │
│                         Eigen 3.4+                               │
│              (expression templates, zero-overhead)               │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                  Eigenvalue Solver Layer                         │
├─────────────────────────────────────────────────────────────────┤
│  Backend        │ Use Case              │ Build Flag             │
├─────────────────┼───────────────────────┼────────────────────────┤
│  LAPACK (MKL)   │ All PC systems        │ -DGEM_ENABLE_MKL=ON    │
│  LAPACK (OpenBLAS)│ All PC systems      │ -DGEM_ENABLE_OPENBLAS=ON│
│  Spectra        │ Large N, few eigs     │ (auto-detected)        │
│  ELPA           │ HPC distributed       │ -DGEM_ENABLE_MPI=ON    │
│  Eigen          │ Fallback only         │ (default if no LAPACK) │
└─────────────────────────────────────────────────────────────────┘
```

**Key Insight**: LAPACK is fast even for small matrices. Use it for everything until the matrix is too large, then switch to iterative (Spectra) or distributed (ELPA) solvers.

| Library | Purpose | When to Use |
|---------|---------|-------------|
| **Intel MKL** | Optimized LAPACK | All personal computer systems (recommended) |
| **OpenBLAS** | Open-source LAPACK | All personal computer systems (if MKL unavailable) |
| **Spectra** | Iterative (Lanczos) | N > 5000 and only need few eigenvalues |
| **ELPA** | Distributed dense | HPC pentaquark (best performance) |
| **Eigen** | Built-in solver | Fallback only (if LAPACK not installed) |

#### Solver Selection (Simplified)

| System | Matrix Size | Solver | Notes |
|--------|-------------|--------|-------|
| Meson | 15 | LAPACK | Trivial, any solver works |
| Baryon | 225 | LAPACK | Fast with LAPACK |
| Tetraquark | 3,375 | LAPACK | ~12 sec with MKL |
| Pentaquark (small) | 10,000 | LAPACK or Spectra | LAPACK if full spectrum needed |
| Pentaquark (large) | 50,000+ | **ELPA** | HPC distributed required |

#### HPC Only (Pentaquark)

| Library | Purpose | Notes |
|---------|---------|-------|
| **MPI** | Distributed-memory parallelization | OpenMPI, Intel MPI, MPICH |
| **ELPA** | Best distributed eigensolver | 2-3x faster than ScaLAPACK |
| **ScaLAPACK** | Standard distributed solver | More widely available |

---

## 4. Input Configuration Design

**Format: YAML** (via yaml-cpp)

### Configuration Structure

```yaml
# config.yaml - GEM Hadron Physics Configuration

system:
  name: "proton"
  hadron_type: "baryon"  # meson | baryon | tetraquark | pentaquark

  particles:
    - {index: 1, flavor: "u", type: "quark"}
    - {index: 2, flavor: "u", type: "quark"}
    - {index: 3, flavor: "d", type: "quark"}

  quantum_numbers:
    total_spin: 0.5
    total_isospin: 0.5
    color: "singlet"

jacobi:
  reference_tree: [[1, 2], 3]
  use_all_jacobi_sets: false

basis:
  num_gaussians: 15
  width_min: 0.05    # fm^-2
  width_max: 20.0    # fm^-2

potential:
  model: "AL1"
  quark_masses:
    u: 0.315  # GeV
    d: 0.315
    s: 0.577
    c: 1.836
    b: 5.227
  three_body_C3: 2.02e-3  # GeV^4, baryons only

solver:
  # Eigenvalue truncation for numerical stability
  truncation:
    method: "relative"
    threshold: 1.0e-10

  num_eigenvalues: 10

  # Backend selection
  # - "auto": LAPACK for PC, ELPA for HPC, Spectra for large N with few eigenvalues
  # - "lapack": Force LAPACK (MKL or OpenBLAS)
  # - "spectra": Force iterative solver (for large N, few eigenvalues)
  # - "elpa": Force ELPA (HPC only)
  backend: "auto"

  # Spectra settings (only used when backend=spectra or auto-selected)
  spectra:
    max_iterations: 1000
    tolerance: 1.0e-10

parallel:
  # Personal computer settings (OpenMP)
  num_threads: 8     # 0 for all cores
  chunk_size: 100
  schedule: "dynamic"

  # HPC cluster settings (MPI+OpenMP) - requires GEM_ENABLE_MPI=ON
  mpi:
    enabled: false
    ranks_per_node: 16
    threads_per_rank: 4

output:
  directory: "./results"
  save_eigenvectors: true
  format: "hdf5"

# HPC job configuration (only used when parallel.mpi.enabled = true)
hpc:
  scheduler: "slurm"        # slurm | pbs
  queue: "normal"
  walltime: "24:00:00"
  nodes: 4
  checkpoint_interval: 3600  # seconds
  restart_from: ""           # path to checkpoint file
```

### C++ Config Parser

```cpp
// config.hpp
#pragma once
#include <yaml-cpp/yaml.h>
#include <vector>
#include <string>
#include <optional>
#include <stdexcept>

namespace gem {

struct ParticleConfig {
    int index;
    std::string flavor;
    std::string type;
};

struct QuantumNumbers {
    double total_spin;
    double total_isospin;
    std::string color;
};

struct SystemConfig {
    std::string name;
    std::string hadron_type;
    std::vector<ParticleConfig> particles;
    QuantumNumbers quantum_numbers;

    int num_particles() const { return particles.size(); }
};

struct JacobiConfig {
    std::vector<std::vector<int>> reference_tree;
    bool use_all_jacobi_sets;
};

struct BasisConfig {
    int num_gaussians;
    double width_min;
    double width_max;
};

struct PotentialConfig {
    std::string model;
    std::map<std::string, double> quark_masses;
    double three_body_C3 = 2.02e-3;  // GeV^4, baryons only (physics.md Section 2.8)
};

struct TruncationConfig {
    std::string method = "relative";  // "relative" or "absolute"
    double threshold = 1e-10;
};

struct SpectraConfig {
    int max_iterations = 1000;
    double tolerance = 1e-10;
};

struct SolverConfig {
    TruncationConfig truncation;     // Eigenvalue truncation settings
    int num_eigenvalues = 10;
    std::string backend = "auto";    // "auto", "lapack", "spectra", "elpa"
    SpectraConfig spectra;           // Spectra-specific settings
};

struct MPIConfig {
    bool enabled;
    int ranks_per_node;
    int threads_per_rank;
};

struct ParallelConfig {
    int num_threads;
    int chunk_size;
    std::string schedule;
    MPIConfig mpi;
};

struct HPCConfig {
    std::string scheduler;
    std::string queue;
    std::string walltime;
    int nodes;
    int checkpoint_interval;
    std::string restart_from;
};

struct OutputConfig {
    std::string directory;
    bool save_eigenvectors;
    std::string format;
};

struct Config {
    SystemConfig system;
    JacobiConfig jacobi;
    BasisConfig basis;
    PotentialConfig potential;
    SolverConfig solver;
    ParallelConfig parallel;
    HPCConfig hpc;
    OutputConfig output;

    static Config from_file(const std::string& path);
    void validate() const;

    bool is_hpc_mode() const { return parallel.mpi.enabled; }
    bool is_pentaquark() const { return system.num_particles() == 5; }
};

} // namespace gem
```

---

## 5. Module Architecture

```
gem/
├── CMakeLists.txt
├── README.md
│
├── include/
│   └── gem/
│       ├── gem.hpp                 # Main header (includes all)
│       ├── config.hpp              # Configuration structures
│       ├── constants.hpp           # Physical constants, AL1 parameters
│       │
│       ├── core/
│       │   ├── types.hpp           # Particle, Hadron, BasisState
│       │   ├── quantum_numbers.hpp # Spin, Color, Isospin types
│       │   └── utils.hpp           # Utility functions
│       │
│       ├── basis/
│       │   ├── gaussian.hpp        # GaussianBasis<N>
│       │   └── combined.hpp        # CombinedBasis
│       │
│       ├── spatial/
│       │   ├── jacobi.hpp          # JacobiSystem<N>
│       │   ├── transform.hpp       # JacobiTransform
│       │   ├── cholesky.hpp        # CholeskyDecomposition
│       │   └── matrix_elements.hpp # SpatialMatrixElement
│       │
│       ├── discrete/
│       │   ├── generators.hpp      # PauliMatrices, GellMannMatrices
│       │   ├── casimir.hpp         # CasimirOperator<Group>
│       │   ├── spin.hpp            # SpinState, SpinCoupling
│       │   ├── color.hpp           # ColorState, ColorSinglet
│       │   └── isospin.hpp         # IsospinState (mixed dims: u,d=2; s,c,b=1)
│       │
│       ├── antisymmetry/
│       │   ├── permutation.hpp     # Permutation, PermutationGroup
│       │   ├── young.hpp           # YoungTableau (optional)
│       │   └── antisymmetrize.hpp  # Antisymmetrizer<N>
│       │
│       ├── hamiltonian/
│       │   ├── potential.hpp       # Potential concept, AL1Potential
│       │   ├── three_body.hpp      # ThreeBodyForce (baryons only)
│       │   ├── kinetic.hpp         # KineticEnergy
│       │   ├── builder.hpp         # HamiltonianBuilder (OpenMP)
│       │   └── matrix_element.hpp  # MatrixElementComputer<N>
│       │
│       ├── hpc/                    # HPC-only (pentaquark)
│       │   ├── mpi_builder.hpp     # MPIHamiltonianBuilder
│       │   ├── mpi_utils.hpp       # MPI helper functions
│       │   ├── distributed.hpp     # Distributed matrix operations
│       │   ├── checkpoint.hpp      # Checkpoint/restart system
│       │   └── slurm.hpp           # SLURM job script generator
│       │
│       ├── solver/
│       │   ├── truncation.hpp      # OverlapTruncation
│       │   ├── eigensolver.hpp     # Solver interface
│       │   ├── lapack_solver.hpp   # LAPACK solver (MKL/OpenBLAS) - primary
│       │   ├── spectra_solver.hpp  # Iterative solver (Spectra)
│       │   ├── elpa_solver.hpp     # Distributed solver (ELPA, HPC)
│       │   ├── eigen_solver.hpp    # Eigen fallback (if no LAPACK)
│       │   └── diagnostics.hpp     # Condition number, verification
│       │
│       └── output/
│           ├── writer.hpp          # ResultWriter interface
│           ├── hdf5_writer.hpp     # HDF5Writer
│           └── analysis.hpp        # SpectrumAnalysis
│
├── src/
│   ├── config.cpp
│   ├── constants.cpp
│   ├── core/
│   ├── basis/
│   ├── spatial/
│   ├── discrete/
│   ├── antisymmetry/
│   ├── hamiltonian/
│   ├── hpc/
│   ├── solver/
│   └── output/
│
├── apps/
│   ├── gem_run.cpp                 # Main executable
│   └── gem_submit.cpp              # HPC job submission tool
│
├── tests/
│   ├── test_jacobi.cpp
│   ├── test_spin.cpp
│   ├── test_color.cpp
│   ├── test_antisymmetry.cpp
│   ├── test_hamiltonian.cpp
│   └── test_solver.cpp
│
├── scripts/
│   ├── submit_slurm.sh             # SLURM submission template
│   └── analyze_results.py          # Post-processing (Python OK here)
│
└── examples/
    ├── meson.yaml
    ├── baryon.yaml
    ├── tetraquark.yaml
    └── pentaquark.yaml
```

---

## 6. Core Type Design

### Compile-Time Particle Count (Template Approach)

```cpp
// core/types.hpp
#pragma once
#include <Eigen/Dense>
#include <array>
#include <variant>

namespace gem {

// Compile-time hadron type for optimization
enum class HadronType { Meson = 2, Baryon = 3, Tetraquark = 4, Pentaquark = 5 };

template<int N>
constexpr int num_jacobi_coords = N - 1;

// Quark flavors
enum class Flavor { Up, Down, Strange, Charm, Bottom };

inline double get_mass(Flavor f, const std::map<std::string, double>& masses) {
    switch (f) {
        case Flavor::Up: return masses.at("u");
        case Flavor::Down: return masses.at("d");
        case Flavor::Strange: return masses.at("s");
        case Flavor::Charm: return masses.at("c");
        case Flavor::Bottom: return masses.at("b");
    }
    return 0.0;
}

// Particle representation
struct Particle {
    int index;
    Flavor flavor;
    double mass;

    bool operator==(const Particle& other) const {
        return flavor == other.flavor;  // For identical particle detection
    }
};

// Fixed-size types for compile-time optimization
template<int N>
using JacobiVector = Eigen::Matrix<double, 3 * num_jacobi_coords<N>, 1>;

template<int N>
using JacobiMatrix = Eigen::Matrix<double, 3 * num_jacobi_coords<N>, 3 * num_jacobi_coords<N>>;

// Basis state with spatial and discrete parts
template<int N>
struct BasisState {
    std::array<double, num_jacobi_coords<N>> gaussian_widths;
    int spin_state_index;
    int color_state_index;
    int isospin_state_index;
    int jacobi_set_index;
};

// Runtime polymorphism when needed
using AnyBasisState = std::variant<
    BasisState<2>,  // Meson
    BasisState<3>,  // Baryon
    BasisState<4>,  // Tetraquark
    BasisState<5>   // Pentaquark
>;

} // namespace gem
```

### Quantum Number Types

```cpp
// core/quantum_numbers.hpp
#pragma once
#include <cstdint>

namespace gem {

// Strong typing for quantum numbers (store 2*value to avoid fractions)
struct Spin {
    int twice_value;

    constexpr double value() const { return twice_value / 2.0; }
    constexpr bool is_half_integer() const { return twice_value % 2 != 0; }

    static constexpr Spin half() { return {1}; }
    static constexpr Spin one() { return {2}; }
    static constexpr Spin zero() { return {0}; }

    bool operator==(const Spin& other) const { return twice_value == other.twice_value; }
};

struct Isospin {
    int twice_value;
    constexpr double value() const { return twice_value / 2.0; }
};

enum class ColorRep {
    Singlet = 1,      // 1
    Triplet = 3,      // 3
    AntiTriplet = -3, // 3̄
    Sextet = 6,       // 6
    Octet = 8         // 8
};

// Target quantum numbers for state construction
struct TargetQuantumNumbers {
    Spin total_spin;
    Isospin total_isospin;
    ColorRep color;
};

} // namespace gem
```

---

## 7. Parallelization Strategy

### Strategy A: Personal Computer (Meson, Baryon, Tetraquark)

Use OpenMP for shared-memory parallelism on a single machine.

#### OpenMP Hamiltonian Builder

```cpp
// hamiltonian/builder.hpp
#pragma once
#include <Eigen/Dense>
#include <omp.h>
#include <vector>

namespace gem {

template<int N>
class HamiltonianBuilder {
public:
    using Matrix = Eigen::MatrixXd;
    using Basis = std::vector<BasisState<N>>;

    struct Result {
        Matrix H;
        Matrix S;
    };

    HamiltonianBuilder(const Config& config, const Basis& basis)
        : config_(config), basis_(basis) {
        int threads = config.parallel.num_threads;
        omp_set_num_threads(threads > 0 ? threads : omp_get_max_threads());
    }

    Result build() {
        const size_t dim = basis_.size();
        const size_t num_elements = dim * (dim + 1) / 2;

        // Pre-allocate flat storage for lock-free parallel writes
        std::vector<double> h_values(num_elements);
        std::vector<double> s_values(num_elements);

        #pragma omp parallel
        {
            MatrixElementComputer<N> computer(config_);

            #pragma omp for schedule(dynamic, config_.parallel.chunk_size)
            for (size_t k = 0; k < num_elements; ++k) {
                auto [i, j] = linear_to_triangular(k, dim);
                auto [h_ij, s_ij] = computer.compute(basis_[i], basis_[j]);

                // No lock needed - each k is unique
                h_values[k] = h_ij;
                s_values[k] = s_ij;
            }
        }

        // Convert to symmetric matrices
        Matrix H = Matrix::Zero(dim, dim);
        Matrix S = Matrix::Zero(dim, dim);

        for (size_t k = 0; k < num_elements; ++k) {
            auto [i, j] = linear_to_triangular(k, dim);
            H(i, j) = H(j, i) = h_values[k];
            S(i, j) = S(j, i) = s_values[k];
        }

        return {std::move(H), std::move(S)};
    }

private:
    static std::pair<size_t, size_t> linear_to_triangular(size_t k, size_t n) {
        size_t i = n - 1 - static_cast<size_t>(
            std::sqrt(2.0 * (n * (n + 1) / 2 - 1 - k) + 0.25) - 0.5
        );
        size_t j = k - i * n + i * (i + 1) / 2;
        return {i, j};
    }

    const Config& config_;
    const Basis& basis_;
};

} // namespace gem
```

### Strategy B: HPC Cluster (Pentaquark)

Use MPI for distributed-memory parallelism across multiple nodes, with OpenMP within each node (hybrid MPI+OpenMP).

```
┌─────────────────────────────────────────────────────────────────┐
│                     HPC Cluster Architecture                     │
├─────────────────────────────────────────────────────────────────┤
│                                                                  │
│  Node 0                          Node 1                          │
│  ┌─────────────────────┐         ┌─────────────────────┐        │
│  │ MPI Rank 0          │         │ MPI Rank 2          │        │
│  │  └─ OMP threads 0-7 │         │  └─ OMP threads 0-7 │        │
│  │ MPI Rank 1          │         │ MPI Rank 3          │        │
│  │  └─ OMP threads 0-7 │         │  └─ OMP threads 0-7 │        │
│  └─────────────────────┘         └─────────────────────┘        │
│           │                               │                      │
│           └───────────────┬───────────────┘                      │
│                           ▼                                      │
│                  InfiniBand Network                              │
│                           │                                      │
│           ┌───────────────┴───────────────┐                      │
│           ▼                               ▼                      │
│  Node 2                          Node N                          │
│  ┌─────────────────────┐         ┌─────────────────────┐        │
│  │ MPI Rank 4-5        │         │ MPI Rank ...        │        │
│  │  └─ OMP threads     │         │  └─ OMP threads     │        │
│  └─────────────────────┘         └─────────────────────┘        │
│                                                                  │
└─────────────────────────────────────────────────────────────────┘
```

#### MPI Hamiltonian Builder

```cpp
// hpc/mpi_builder.hpp
#pragma once
#ifdef USE_MPI

#include <mpi.h>
#include <Eigen/Dense>
#include <vector>
#include <fstream>
#include <omp.h>

namespace gem {

template<int N>
class MPIHamiltonianBuilder {
public:
    using Matrix = Eigen::MatrixXd;
    using Basis = std::vector<BasisState<N>>;

    MPIHamiltonianBuilder(const Config& config, const Basis& basis)
        : config_(config), basis_(basis) {
        MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
        MPI_Comm_size(MPI_COMM_WORLD, &size_);

        // Set OpenMP threads per MPI rank
        int threads = config.parallel.mpi.threads_per_rank;
        omp_set_num_threads(threads > 0 ? threads : 1);
    }

    std::pair<Matrix, Matrix> build() {
        const size_t dim = basis_.size();
        const size_t num_elements = dim * (dim + 1) / 2;

        // Distribute work across MPI ranks
        std::vector<size_t> my_indices;
        for (size_t k = rank_; k < num_elements; k += size_) {
            my_indices.push_back(k);
        }

        // Compute local elements with OpenMP
        std::vector<std::tuple<size_t, size_t, double, double>> local_results;
        local_results.reserve(my_indices.size());

        size_t checkpoint_count = 0;
        size_t checkpoint_interval = config_.hpc.checkpoint_interval;

        #pragma omp parallel
        {
            MatrixElementComputer<N> computer(config_);
            std::vector<std::tuple<size_t, size_t, double, double>> thread_results;

            #pragma omp for schedule(dynamic, config_.parallel.chunk_size)
            for (size_t idx = 0; idx < my_indices.size(); ++idx) {
                size_t k = my_indices[idx];
                auto [i, j] = linear_to_triangular(k, dim);
                auto [h_ij, s_ij] = computer.compute(basis_[i], basis_[j]);
                thread_results.emplace_back(i, j, h_ij, s_ij);
            }

            #pragma omp critical
            {
                local_results.insert(local_results.end(),
                                     thread_results.begin(),
                                     thread_results.end());
                checkpoint_count += thread_results.size();

                // Periodic checkpoint
                if (checkpoint_count >= checkpoint_interval) {
                    save_checkpoint(local_results);
                    checkpoint_count = 0;
                }
            }
        }

        // Final checkpoint
        save_checkpoint(local_results);

        // Gather results to rank 0
        return gather_results(local_results, dim);
    }

    // Load from checkpoint for restart
    static std::vector<std::tuple<size_t, size_t, double, double>>
    load_checkpoint(int rank) {
        std::string filename = "checkpoint_rank" + std::to_string(rank) + ".bin";
        std::ifstream ifs(filename, std::ios::binary);
        if (!ifs) return {};

        size_t n;
        ifs.read(reinterpret_cast<char*>(&n), sizeof(n));

        std::vector<std::tuple<size_t, size_t, double, double>> results(n);
        ifs.read(reinterpret_cast<char*>(results.data()),
                 n * sizeof(std::tuple<size_t, size_t, double, double>));
        return results;
    }

private:
    std::pair<Matrix, Matrix> gather_results(
        const std::vector<std::tuple<size_t, size_t, double, double>>& local,
        size_t dim) {

        // Gather sizes from all ranks
        int local_count = local.size();
        std::vector<int> counts(size_);
        MPI_Gather(&local_count, 1, MPI_INT,
                   counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

        // Prepare flat buffers for MPI communication
        std::vector<size_t> local_ij;
        std::vector<double> local_h, local_s;
        local_ij.reserve(2 * local.size());
        local_h.reserve(local.size());
        local_s.reserve(local.size());

        for (const auto& [i, j, h, s] : local) {
            local_ij.push_back(i);
            local_ij.push_back(j);
            local_h.push_back(h);
            local_s.push_back(s);
        }

        if (rank_ == 0) {
            Matrix H = Matrix::Zero(dim, dim);
            Matrix S = Matrix::Zero(dim, dim);

            // Process local results first
            for (const auto& [i, j, h, s] : local) {
                H(i, j) = H(j, i) = h;
                S(i, j) = S(j, i) = s;
            }

            // Receive from other ranks
            for (int r = 1; r < size_; ++r) {
                std::vector<size_t> recv_ij(2 * counts[r]);
                std::vector<double> recv_h(counts[r]), recv_s(counts[r]);

                MPI_Recv(recv_ij.data(), 2 * counts[r], MPI_UNSIGNED_LONG,
                         r, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(recv_h.data(), counts[r], MPI_DOUBLE,
                         r, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(recv_s.data(), counts[r], MPI_DOUBLE,
                         r, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                for (int k = 0; k < counts[r]; ++k) {
                    size_t i = recv_ij[2*k], j = recv_ij[2*k + 1];
                    H(i, j) = H(j, i) = recv_h[k];
                    S(i, j) = S(j, i) = recv_s[k];
                }
            }

            return {H, S};
        } else {
            // Send to rank 0
            MPI_Send(local_ij.data(), 2 * local_count, MPI_UNSIGNED_LONG,
                     0, 0, MPI_COMM_WORLD);
            MPI_Send(local_h.data(), local_count, MPI_DOUBLE,
                     0, 1, MPI_COMM_WORLD);
            MPI_Send(local_s.data(), local_count, MPI_DOUBLE,
                     0, 2, MPI_COMM_WORLD);

            return {Matrix(), Matrix()};
        }
    }

    void save_checkpoint(
        const std::vector<std::tuple<size_t, size_t, double, double>>& results) {
        std::string filename = "checkpoint_rank" + std::to_string(rank_) + ".bin";
        std::ofstream ofs(filename, std::ios::binary);
        size_t n = results.size();
        ofs.write(reinterpret_cast<const char*>(&n), sizeof(n));
        ofs.write(reinterpret_cast<const char*>(results.data()),
                  n * sizeof(std::tuple<size_t, size_t, double, double>));
    }

    static std::pair<size_t, size_t> linear_to_triangular(size_t k, size_t n) {
        size_t i = n - 1 - static_cast<size_t>(
            std::sqrt(2.0 * (n * (n + 1) / 2 - 1 - k) + 0.25) - 0.5
        );
        size_t j = k - i * n + i * (i + 1) / 2;
        return {i, j};
    }

    const Config& config_;
    const Basis& basis_;
    int rank_, size_;
};

} // namespace gem

#endif // USE_MPI
```

#### SLURM Job Script Generator

```cpp
// hpc/slurm.hpp
#pragma once
#include <string>
#include <sstream>
#include <fstream>

namespace gem {

inline std::string generate_slurm_script(const Config& config,
                                         const std::string& config_file) {
    int total_tasks = config.hpc.nodes * config.parallel.mpi.ranks_per_node;

    std::ostringstream ss;
    ss << "#!/bin/bash\n"
       << "#SBATCH --job-name=gem_" << config.system.name << "\n"
       << "#SBATCH --partition=" << config.hpc.queue << "\n"
       << "#SBATCH --nodes=" << config.hpc.nodes << "\n"
       << "#SBATCH --ntasks-per-node=" << config.parallel.mpi.ranks_per_node << "\n"
       << "#SBATCH --cpus-per-task=" << config.parallel.mpi.threads_per_rank << "\n"
       << "#SBATCH --time=" << config.hpc.walltime << "\n"
       << "#SBATCH --output=gem_%j.out\n"
       << "#SBATCH --error=gem_%j.err\n"
       << "\n"
       << "# Load modules (adjust for your cluster)\n"
       << "module load gcc/11\n"
       << "module load openmpi/4.1\n"
       << "module load mkl/2023\n"
       << "module load hdf5/1.12\n"
       << "\n"
       << "# Set OpenMP environment\n"
       << "export OMP_NUM_THREADS=" << config.parallel.mpi.threads_per_rank << "\n"
       << "export OMP_PROC_BIND=spread\n"
       << "export OMP_PLACES=threads\n"
       << "\n"
       << "# Prevent memory overcommit\n"
       << "export OMP_STACKSIZE=64M\n"
       << "\n"
       << "# Run calculation\n"
       << "srun --ntasks=" << total_tasks << " ./gem_run " << config_file << "\n";

    return ss.str();
}

inline void write_slurm_script(const Config& config,
                               const std::string& config_file,
                               const std::string& script_file = "submit.sh") {
    std::ofstream ofs(script_file);
    ofs << generate_slurm_script(config, config_file);
}

} // namespace gem
```

### Performance Estimates

#### Personal Computer (16 cores, OpenMP)

| System | Basis | Matrix Elements | Time | Memory |
|--------|-------|-----------------|------|--------|
| Meson | 15 | 120 | < 0.01 sec | < 1 MB |
| Baryon | 225 | 25,425 | ~0.5 sec | ~50 MB |
| Tetraquark | 3,375 | 5.7M | ~30 sec | ~1 GB |

#### HPC Cluster (Pentaquark, MPI+OpenMP)

| Configuration | Basis | Nodes | Cores | Time | Memory/Node |
|---------------|-------|-------|-------|------|-------------|
| Small | 10,000 | 4 | 256 | ~30 min | 64 GB |
| Medium | 50,000 | 16 | 1024 | ~2 hr | 128 GB |
| Large | 200,000 | 64 | 4096 | ~8 hr | 256 GB |

---

## 8. Eigenvalue Solver Design

### Simplified Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                       User Configuration                         │
│                    solver.backend: "auto"                        │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                      Solver Factory                              │
│                                                                  │
│  Personal Computer:                                              │
│    → LAPACK (MKL or OpenBLAS) for all systems                   │
│    → Spectra if N > 5000 and only few eigenvalues needed        │
│                                                                  │
│  HPC Cluster:                                                    │
│    → ELPA (best) or ScaLAPACK (widely available)                │
│                                                                  │
│  Fallback:                                                       │
│    → Eigen (only if LAPACK not installed)                       │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                    Concrete Solvers                              │
├──────────────────┬────────────────┬─────────────────────────────┤
│  LAPACKSolver    │  SpectraSolver │  ELPASolver                 │
│  (all PC cases)  │  (iterative)   │  (HPC distributed)          │
└──────────────────┴────────────────┴─────────────────────────────┘
```

### Performance (LAPACK with MKL)

| System | Matrix Size | Time | Memory |
|--------|-------------|------|--------|
| Meson | 15×15 | < 1 ms | < 1 KB |
| Baryon | 225×225 | ~5 ms | ~400 KB |
| Tetraquark | 3375×3375 | ~12 s | ~90 MB |
| Pentaquark | 50000×50000 | ~90 s* | ~20 GB |

*Pentaquark on single node; use ELPA on HPC for better scaling.

### Solver Interface

```cpp
// solver/eigensolver.hpp
#pragma once
#include <Eigen/Dense>
#include <concepts>
#include <memory>
#include <string>

namespace gem {

// Available solver backends
enum class SolverBackend {
    Auto,       // Auto-select based on matrix size
    Eigen,      // Eigen's built-in solver (fallback)
    LAPACK,     // MKL or OpenBLAS LAPACK
    Spectra,    // Iterative solver (few eigenvalues)
    ELPA,       // HPC distributed (best)
    ScaLAPACK   // HPC distributed (widely available)
};

struct EigenResult {
    Eigen::VectorXd eigenvalues;
    Eigen::MatrixXd eigenvectors;
    int num_truncated;        // Basis functions removed
    double condition_number;  // For diagnostics
    std::string solver_used;  // Which backend was used
};

template<typename T>
concept EigenSolver = requires(T solver,
                               const Eigen::MatrixXd& H,
                               const Eigen::MatrixXd& S,
                               int k) {
    { solver.solve(H, S, k) } -> std::same_as<EigenResult>;
};

// Abstract base for runtime polymorphism
class EigenSolverBase {
public:
    virtual ~EigenSolverBase() = default;
    virtual EigenResult solve(const Eigen::MatrixXd& H,
                              const Eigen::MatrixXd& S,
                              int num_eigenvalues) = 0;
    virtual std::string name() const = 0;
};

} // namespace gem
```

### Solver Factory (Simplified)

```cpp
// solver/solver_factory.hpp
#pragma once
#include "eigensolver.hpp"
#include "lapack_solver.hpp"
#ifdef USE_SPECTRA
#include "spectra_solver.hpp"
#endif
#ifdef USE_ELPA
#include "elpa_solver.hpp"
#endif

namespace gem {

class SolverFactory {
public:
    static std::unique_ptr<EigenSolverBase> create(
        const SolverConfig& config,
        size_t matrix_size) {

        SolverBackend backend = parse_backend(config.backend);

        // Auto-select: simplified logic
        if (backend == SolverBackend::Auto) {
            backend = auto_select(matrix_size, config.num_eigenvalues);
        }

        return create_backend(backend, config);
    }

private:
    static SolverBackend parse_backend(const std::string& name) {
        if (name == "auto") return SolverBackend::Auto;
        if (name == "lapack") return SolverBackend::LAPACK;
        if (name == "spectra") return SolverBackend::Spectra;
        if (name == "elpa") return SolverBackend::ELPA;
        return SolverBackend::Auto;
    }

    static SolverBackend auto_select(size_t matrix_size, int num_eigenvalues) {
        // HPC mode: use distributed solver
#ifdef USE_MPI
        int world_size;
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
        if (world_size > 1) {
#ifdef USE_ELPA
            return SolverBackend::ELPA;
#endif
        }
#endif

        // Large matrix with few eigenvalues: use iterative solver
#ifdef USE_SPECTRA
        if (matrix_size > 5000 && num_eigenvalues < 50) {
            return SolverBackend::Spectra;
        }
#endif

        // Default: LAPACK for everything
        // (fast even for small matrices, no reason to use Eigen)
        return SolverBackend::LAPACK;
    }

    static std::unique_ptr<EigenSolverBase> create_backend(
        SolverBackend backend,
        const SolverConfig& config) {

        switch (backend) {
#ifdef USE_LAPACK
            case SolverBackend::LAPACK:
                return std::make_unique<LAPACKSolver>(config.truncation.threshold);
#endif

#ifdef USE_SPECTRA
            case SolverBackend::Spectra:
                return std::make_unique<SpectraSolver>(
                    config.spectra.max_iterations,
                    config.spectra.tolerance);
#endif

#ifdef USE_ELPA
            case SolverBackend::ELPA:
                return std::make_unique<ELPASolver>(config.truncation.threshold);
#endif

            default:
                // Fallback to Eigen only if LAPACK not available
#ifdef USE_LAPACK
                return std::make_unique<LAPACKSolver>(config.truncation.threshold);
#else
                return std::make_unique<EigenDenseSolver>(config.truncation.threshold);
#endif
        }
    }
};

} // namespace gem
```

### Eigen Solver (Fallback Only)

**Note**: Only used if LAPACK (MKL/OpenBLAS) is not installed. In practice, LAPACK should always be available for scientific computing.

```cpp
// solver/eigen_solver.hpp
#pragma once
#include "eigensolver.hpp"
#include <Eigen/Eigenvalues>

namespace gem {

class EigenDenseSolver : public EigenSolverBase {
public:
    explicit EigenDenseSolver(double truncation_threshold = 1e-10)
        : threshold_(truncation_threshold) {}

    std::string name() const override { return "Eigen"; }

    EigenResult solve(const Eigen::MatrixXd& H,
                      const Eigen::MatrixXd& S,
                      int num_eigenvalues) override {
        // Step 1: Eigenvalue decomposition of overlap matrix
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> overlap_solver(S);
        const auto& eigenvalues_S = overlap_solver.eigenvalues();
        const auto& U = overlap_solver.eigenvectors();

        // Step 2: Truncation of near-zero eigenvalues
        double max_eigenvalue = eigenvalues_S.maxCoeff();
        double abs_threshold = threshold_ * max_eigenvalue;

        std::vector<int> valid_indices;
        for (int i = 0; i < eigenvalues_S.size(); ++i) {
            if (eigenvalues_S(i) > abs_threshold) {
                valid_indices.push_back(i);
            }
        }

        int n_eff = valid_indices.size();
        int n_truncated = eigenvalues_S.size() - n_eff;

        // Step 3: Build transformation matrix X = U_trunc * Lambda^{-1/2}
        Eigen::MatrixXd X(H.rows(), n_eff);
        for (int j = 0; j < n_eff; ++j) {
            int idx = valid_indices[j];
            X.col(j) = U.col(idx) / std::sqrt(eigenvalues_S(idx));
        }

        // Step 4: Transform to orthonormal basis
        Eigen::MatrixXd H_tilde = X.transpose() * H * X;

        // Step 5: Solve standard eigenvalue problem
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> hamiltonian_solver(H_tilde);

        // Step 6: Extract requested eigenvalues
        int k = std::min(num_eigenvalues, n_eff);
        Eigen::VectorXd energies = hamiltonian_solver.eigenvalues().head(k);

        // Step 7: Back-transform eigenvectors
        Eigen::MatrixXd C_tilde = hamiltonian_solver.eigenvectors().leftCols(k);
        Eigen::MatrixXd coefficients = X * C_tilde;

        double condition = max_eigenvalue / eigenvalues_S(valid_indices.front());

        return {energies, coefficients, n_truncated, condition, "Eigen"};
    }

private:
    double threshold_;
};

} // namespace gem
```

### LAPACK Solver (Primary - All Personal Computer Systems)

**This is the main solver for all personal computer calculations** (meson, baryon, tetraquark).
Uses Intel MKL or OpenBLAS backend for optimized performance.

```cpp
// solver/lapack_solver.hpp
#pragma once
#ifdef USE_LAPACK

#include "eigensolver.hpp"
#include <Eigen/Dense>

// LAPACK function declarations
extern "C" {
    // Symmetric generalized eigenvalue problem (divide-and-conquer)
    void dsygvd_(int* itype, char* jobz, char* uplo, int* n,
                 double* a, int* lda, double* b, int* ldb,
                 double* w, double* work, int* lwork,
                 int* iwork, int* liwork, int* info);

    // Symmetric eigenvalue problem (divide-and-conquer)
    void dsyevd_(char* jobz, char* uplo, int* n, double* a, int* lda,
                 double* w, double* work, int* lwork,
                 int* iwork, int* liwork, int* info);
}

namespace gem {

class LAPACKSolver : public EigenSolverBase {
public:
    explicit LAPACKSolver(double truncation_threshold = 1e-10,
                          bool use_truncation = true)
        : threshold_(truncation_threshold), use_truncation_(use_truncation) {}

    std::string name() const override { return "LAPACK"; }

    EigenResult solve(const Eigen::MatrixXd& H,
                      const Eigen::MatrixXd& S,
                      int num_eigenvalues) override {
        int n = H.rows();

        if (use_truncation_) {
            // Method 1: Eigenvalue truncation for ill-conditioned bases
            // (recommended - follows algorithm_eigenvalue_truncation.md)
            return solve_with_truncation(H, S, num_eigenvalues);
        } else {
            // Method 2: Direct generalized solver (for well-conditioned bases only)
            return solve_direct(H, S, num_eigenvalues);
        }
    }

private:
    double threshold_;
    bool use_truncation_;

    EigenResult solve_with_truncation(const Eigen::MatrixXd& H,
                                       const Eigen::MatrixXd& S,
                                       int num_eigenvalues) {
        int n = H.rows();

        // Step 1: Diagonalize overlap matrix using LAPACK dsyevd
        Eigen::MatrixXd S_copy = S;
        Eigen::VectorXd eigenvalues_S(n);

        char jobz = 'V';
        char uplo = 'U';
        int lwork = -1, liwork = -1, info;
        double work_query;
        int iwork_query;

        dsyevd_(&jobz, &uplo, &n, S_copy.data(), &n,
                eigenvalues_S.data(), &work_query, &lwork,
                &iwork_query, &liwork, &info);

        lwork = static_cast<int>(work_query);
        liwork = iwork_query;
        std::vector<double> work(lwork);
        std::vector<int> iwork(liwork);

        dsyevd_(&jobz, &uplo, &n, S_copy.data(), &n,
                eigenvalues_S.data(), work.data(), &lwork,
                iwork.data(), &liwork, &info);

        if (info != 0) {
            throw std::runtime_error("LAPACK dsyevd failed (overlap) with info = " +
                                     std::to_string(info));
        }

        // S_copy now contains eigenvectors (columns)
        Eigen::MatrixXd U = S_copy;

        // Step 2: Truncate small eigenvalues
        double max_eigenvalue = eigenvalues_S.maxCoeff();
        double abs_threshold = threshold_ * max_eigenvalue;

        std::vector<int> valid_indices;
        for (int i = 0; i < n; ++i) {
            if (eigenvalues_S(i) > abs_threshold) {
                valid_indices.push_back(i);
            }
        }

        int n_eff = valid_indices.size();
        if (n_eff == 0) {
            throw std::runtime_error("All overlap eigenvalues below threshold");
        }

        int modes_removed = n - n_eff;

        // Compute condition number
        double lambda_min = eigenvalues_S(valid_indices.front());
        double condition_number = max_eigenvalue / lambda_min;

        // Step 3: Build transformation matrix X = U_trunc * Lambda^{-1/2}
        Eigen::MatrixXd X(n, n_eff);
        for (int j = 0; j < n_eff; ++j) {
            int idx = valid_indices[j];
            X.col(j) = U.col(idx) / std::sqrt(eigenvalues_S(idx));
        }

        // Step 4: Transform Hamiltonian: H_tilde = X^T * H * X
        Eigen::MatrixXd H_tilde = X.transpose() * H * X;

        // Step 5: Solve standard eigenvalue problem using LAPACK dsyevd
        Eigen::VectorXd energies(n_eff);

        lwork = -1;
        liwork = -1;
        dsyevd_(&jobz, &uplo, &n_eff, H_tilde.data(), &n_eff,
                energies.data(), &work_query, &lwork,
                &iwork_query, &liwork, &info);

        lwork = static_cast<int>(work_query);
        liwork = iwork_query;
        work.resize(lwork);
        iwork.resize(liwork);

        dsyevd_(&jobz, &uplo, &n_eff, H_tilde.data(), &n_eff,
                energies.data(), work.data(), &lwork,
                iwork.data(), &liwork, &info);

        if (info != 0) {
            throw std::runtime_error("LAPACK dsyevd failed (H_tilde) with info = " +
                                     std::to_string(info));
        }

        // H_tilde now contains eigenvectors in reduced basis
        Eigen::MatrixXd C_tilde = H_tilde;

        // Step 6: Back-transform to original basis
        Eigen::MatrixXd coefficients = X * C_tilde;

        int k = std::min(num_eigenvalues, n_eff);
        return {
            energies.head(k),
            coefficients.leftCols(k),
            modes_removed,
            condition_number,
            "LAPACK"
        };
    }

    EigenResult solve_direct(const Eigen::MatrixXd& H,
                              const Eigen::MatrixXd& S,
                              int num_eigenvalues) {
        int n = H.rows();

        // Direct generalized solver (for well-conditioned bases only)
        Eigen::MatrixXd A = H;
        Eigen::MatrixXd B = S;
        Eigen::VectorXd W(n);

        int itype = 1;
        char jobz = 'V';
        char uplo = 'U';
        int lwork = -1, liwork = -1, info;
        double work_query;
        int iwork_query;

        dsygvd_(&itype, &jobz, &uplo, &n,
                A.data(), &n, B.data(), &n,
                W.data(), &work_query, &lwork,
                &iwork_query, &liwork, &info);

        lwork = static_cast<int>(work_query);
        liwork = iwork_query;
        std::vector<double> work(lwork);
        std::vector<int> iwork(liwork);

        dsygvd_(&itype, &jobz, &uplo, &n,
                A.data(), &n, B.data(), &n,
                W.data(), work.data(), &lwork,
                iwork.data(), &liwork, &info);

        if (info != 0) {
            throw std::runtime_error("LAPACK dsygvd failed with info = " +
                                     std::to_string(info));
        }

        int k = std::min(num_eigenvalues, n);
        return {
            W.head(k),
            A.leftCols(k),
            0,
            -1,
            "LAPACK"
        };
    }
};

} // namespace gem

#endif // USE_LAPACK
```

### Spectra Solver (Large Systems, Few Eigenvalues)

```cpp
// solver/spectra_solver.hpp
#pragma once
#ifdef USE_SPECTRA

#include "eigensolver.hpp"
#include <Spectra/SymGEigsSolver.h>
#include <Spectra/MatOp/DenseSymMatProd.h>
#include <Spectra/MatOp/DenseCholesky.h>

namespace gem {

class SpectraSolver : public EigenSolverBase {
public:
    explicit SpectraSolver(int max_iterations = 1000,
                           double tolerance = 1e-10)
        : max_iter_(max_iterations), tol_(tolerance) {}

    std::string name() const override { return "Spectra"; }

    EigenResult solve(const Eigen::MatrixXd& H,
                      const Eigen::MatrixXd& S,
                      int num_eigenvalues) override {
        using namespace Spectra;

        DenseSymMatProd<double> op_H(H);
        DenseCholesky<double> op_S(S);

        int ncv = std::min(static_cast<int>(H.rows()),
                          std::max(2 * num_eigenvalues + 1, 20));

        SymGEigsSolver<DenseSymMatProd<double>,
                       DenseCholesky<double>,
                       GEigsMode::Cholesky>
            solver(op_H, op_S, num_eigenvalues, ncv);

        solver.init();
        solver.compute(SortRule::SmallestAlge, max_iter_, tol_);

        if (solver.info() != CompInfo::Successful) {
            throw std::runtime_error("Eigenvalue computation failed");
        }

        return {
            solver.eigenvalues(),
            solver.eigenvectors(),
            0,      // No truncation in iterative solver
            -1,     // Condition number not computed
            "Spectra"
        };
    }

private:
    int max_iter_;
    double tol_;
};

} // namespace gem

#endif // USE_SPECTRA
```

### Distributed Solver (HPC Only)

```cpp
// solver/distributed_solver.hpp
#pragma once
#ifdef USE_SCALAPACK

#include "eigensolver.hpp"
#include <mpi.h>

namespace gem {

class ScaLAPACKEigenSolver : public EigenSolverBase {
public:
    ScaLAPACKEigenSolver(int block_size = 64)
        : block_size_(block_size) {
        MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
        MPI_Comm_size(MPI_COMM_WORLD, &size_);
    }

    EigenResult solve(const Eigen::MatrixXd& H,
                      const Eigen::MatrixXd& S,
                      int num_eigenvalues) override {
        // ScaLAPACK implementation for distributed eigenvalue problem
        // Uses PDSYGVX for generalized symmetric eigenvalue problem

        // This is a simplified interface - actual implementation requires:
        // 1. Setting up BLACS grid
        // 2. Distributing matrices H and S across processes
        // 3. Calling PDSYGVX
        // 4. Gathering results

        // Placeholder - actual implementation depends on ScaLAPACK bindings
        throw std::runtime_error("ScaLAPACK solver not yet implemented");
    }

private:
    int block_size_;
    int rank_, size_;
};

} // namespace gem

#endif // USE_SCALAPACK
```

---

## 9. Spatial Matrix Elements

### Jacobi Coordinate System

```cpp
// spatial/jacobi.hpp
#pragma once
#include <Eigen/Dense>
#include <array>

namespace gem {

template<int N>
class JacobiSystem {
public:
    static constexpr int num_coords = N - 1;
    static constexpr int spatial_dim = 3 * num_coords;

    using CoordVector = Eigen::Matrix<double, spatial_dim, 1>;
    using TransformMatrix = Eigen::Matrix<double, spatial_dim, spatial_dim>;
    using MassArray = std::array<double, N>;

    JacobiSystem(const MassArray& masses,
                 const std::vector<std::vector<int>>& tree)
        : masses_(masses) {
        compute_transformation(tree);
    }

    // Transform from one Jacobi set to another
    static TransformMatrix transformation(const JacobiSystem& from,
                                          const JacobiSystem& to);

    // Get transformation to center-of-mass coordinates
    const TransformMatrix& to_cm_transform() const { return T_to_cm_; }

    // Reduced masses for each Jacobi coordinate
    const std::array<double, num_coords>& reduced_masses() const {
        return reduced_masses_;
    }

private:
    void compute_transformation(const std::vector<std::vector<int>>& tree);

    MassArray masses_;
    TransformMatrix T_to_cm_;
    std::array<double, num_coords> reduced_masses_;
};

} // namespace gem
```

### Cholesky Matrix Element

```cpp
// spatial/cholesky.hpp
#pragma once
#include <Eigen/Dense>
#include <Eigen/Cholesky>
#include <cmath>

namespace gem {

class CholeskyMatrixElement {
public:
    // For diagonal width matrices (product of 1D Gaussians)
    template<int Dim>
    static double overlap(const Eigen::Matrix<double, Dim, 1>& widths_a,
                          const Eigen::Matrix<double, Dim, 1>& widths_b) {
        Eigen::Matrix<double, Dim, 1> total = widths_a + widths_b;
        double det = total.prod();

        constexpr double pi = 3.14159265358979323846;
        return std::pow(pi, Dim / 2.0) / std::sqrt(det);
    }

    // For general (non-diagonal) width matrices
    template<int Dim>
    static double overlap(const Eigen::Matrix<double, Dim, Dim>& A,
                          const Eigen::Matrix<double, Dim, Dim>& B) {
        Eigen::Matrix<double, Dim, Dim> total = A + B;

        Eigen::LLT<Eigen::Matrix<double, Dim, Dim>> llt(total);
        double log_det = 2.0 * llt.matrixL().diagonal().array().log().sum();

        constexpr double pi = 3.14159265358979323846;
        return std::pow(pi, Dim / 2.0) * std::exp(-0.5 * log_det);
    }

    // Kinetic energy matrix element
    template<int Dim>
    static double kinetic(const Eigen::Matrix<double, Dim, 1>& widths_a,
                          const Eigen::Matrix<double, Dim, 1>& widths_b,
                          const Eigen::Matrix<double, Dim, 1>& reduced_masses) {
        Eigen::Matrix<double, Dim, 1> total = widths_a + widths_b;
        double S = overlap<Dim>(widths_a, widths_b);

        double T = 0.0;
        for (int i = 0; i < Dim; ++i) {
            double a_i = widths_a(i);
            double b_i = widths_b(i);
            double mu_i = reduced_masses(i);

            // Each 3D coordinate contributes: 3 * a * b / (a + b) / mu
            T += 3.0 * a_i * b_i / total(i) / mu_i;
        }

        return T * S;
    }
};

} // namespace gem
```

---

## 10. Discrete Wave Functions

### Spin Coupling

```cpp
// discrete/spin.hpp
#pragma once
#include <Eigen/Dense>
#include <vector>
#include <map>
#include <tuple>

namespace gem {

// Clebsch-Gordan coefficient cache
class ClebschGordan {
public:
    static double coefficient(int j1_twice, int m1_twice,
                              int j2_twice, int m2_twice,
                              int J_twice, int M_twice);

private:
    static std::map<std::tuple<int,int,int,int,int,int>, double> cache_;
};

// Spin state for N spin-1/2 particles
template<int N>
class SpinStates {
public:
    using State = Eigen::VectorXd;  // Coefficients in computational basis

    // Construct all states with given total spin S and projection M
    SpinStates(Spin total_spin, int M_twice = 0);

    int num_states() const { return states_.size(); }
    const State& state(int index) const { return states_[index]; }

    // Overlap matrix between states
    Eigen::MatrixXd overlap_matrix() const;

    // Apply permutation to spin state
    State permute(const State& state, const class Permutation& perm) const;

private:
    void construct_states(Spin S, int M_twice);

    std::vector<State> states_;
    int total_spin_twice_;
};

} // namespace gem
```

### Color Wave Functions

```cpp
// discrete/color.hpp
#pragma once
#include <Eigen/Dense>
#include <vector>
#include <array>
#include <complex>

namespace gem {

// SU(3) color states
template<int N>
class ColorStates {
public:
    using State = Eigen::VectorXcd;  // Complex for general representations

    // Construct color singlet states for N quarks
    static std::vector<State> singlet_states();

    static int num_singlet_states();

    // Color matrix element <color_a | O_color | color_b>
    static std::complex<double> matrix_element(
        const State& bra,
        const State& ket,
        int quark_i,
        int quark_j);

    // SU(3) Gell-Mann matrices
    static const std::array<Eigen::Matrix3cd, 8>& gell_mann();
};

// Template specializations
template<> std::vector<ColorStates<2>::State> ColorStates<2>::singlet_states();
template<> std::vector<ColorStates<3>::State> ColorStates<3>::singlet_states();
template<> std::vector<ColorStates<4>::State> ColorStates<4>::singlet_states();
template<> std::vector<ColorStates<5>::State> ColorStates<5>::singlet_states();

} // namespace gem
```

---

## 11. Antisymmetrization

### Permutation Group

```cpp
// antisymmetry/permutation.hpp
#pragma once
#include <vector>
#include <utility>

namespace gem {

class Permutation {
public:
    explicit Permutation(std::vector<int> mapping);

    int operator()(int i) const { return mapping_[i]; }
    int sign() const { return sign_; }
    int size() const { return mapping_.size(); }

    Permutation operator*(const Permutation& other) const;
    Permutation inverse() const;
    bool is_identity() const;

    std::vector<std::pair<int, int>> transpositions() const;

private:
    std::vector<int> mapping_;
    int sign_;

    void compute_sign();
};

// Generate symmetric group S_n
std::vector<Permutation> symmetric_group(int n);

// Generate permutations of identical particles only
std::vector<Permutation> identical_particle_permutations(
    const std::vector<std::vector<int>>& identical_groups);

} // namespace gem
```

### Antisymmetrizer

```cpp
// antisymmetry/antisymmetrize.hpp
#pragma once
#include "permutation.hpp"
#include "../core/types.hpp"
#include <Eigen/Dense>

namespace gem {

template<int N>
class Antisymmetrizer {
public:
    explicit Antisymmetrizer(const std::array<Flavor, N>& flavors);

    // Compute antisymmetrized matrix element
    template<typename MatrixElementFunc>
    double compute(const BasisState<N>& bra,
                   const BasisState<N>& ket,
                   MatrixElementFunc&& me_func) const {
        double result = 0.0;

        for (const auto& perm : permutations_) {
            BasisState<N> permuted_ket = apply_permutation(ket, perm);
            double me = me_func(bra, permuted_ket);
            result += perm.sign() * me;
        }

        return result * normalization_;
    }

    int num_permutations() const { return permutations_.size(); }

private:
    BasisState<N> apply_permutation(const BasisState<N>& state,
                                    const Permutation& perm) const;

    std::vector<Permutation> permutations_;
    double normalization_;
    std::vector<std::vector<int>> identical_groups_;
};

} // namespace gem
```

---

## 12. Large Matrix Handling (Pentaquark Scale - HPC Only)

**Note**: This section applies only to pentaquark calculations running on HPC clusters.
Meson, baryon, and tetraquark calculations run on personal computers using dense in-memory storage.

### Problem Size

| System | Basis | Matrix Elements | Memory (Dense) | Environment |
|--------|-------|-----------------|----------------|-------------|
| Meson | 15 | 120 | < 1 KB | Personal Computer |
| Baryon | 225 | 25K | 400 KB | Personal Computer |
| Tetraquark | 3,375 | 5.7M | 91 MB | Personal Computer |
| **Pentaquark** | 50,625 | 1.28B | **20 GB** | **HPC Cluster** |

### Strategy 1: Matrix-Free Iterative Solver

```cpp
// solver/matrix_free.hpp
#pragma once
#ifdef USE_SPECTRA

#include <Spectra/SymEigsSolver.h>
#include <functional>
#include <omp.h>

namespace gem {

template<int N>
class MatrixFreeOperator {
public:
    MatrixFreeOperator(size_t dimension,
                       const std::vector<BasisState<N>>& basis,
                       const Config& config)
        : dim_(dimension), basis_(basis), config_(config) {}

    size_t rows() const { return dim_; }
    size_t cols() const { return dim_; }

    // y = H * x (computed on-the-fly, never store full matrix)
    void perform_op(const double* x_in, double* y_out) const {
        Eigen::Map<const Eigen::VectorXd> x(x_in, dim_);
        Eigen::Map<Eigen::VectorXd> y(y_out, dim_);
        y.setZero();

        #pragma omp parallel
        {
            Eigen::VectorXd y_local = Eigen::VectorXd::Zero(dim_);
            MatrixElementComputer<N> computer(config_);

            #pragma omp for schedule(dynamic, 100)
            for (size_t i = 0; i < dim_; ++i) {
                for (size_t j = 0; j < dim_; ++j) {
                    double h_ij = computer.compute_hamiltonian(basis_[i], basis_[j]);
                    y_local(i) += h_ij * x(j);
                }
            }

            #pragma omp critical
            y += y_local;
        }
    }

private:
    size_t dim_;
    const std::vector<BasisState<N>>& basis_;
    const Config& config_;
};

} // namespace gem

#endif
```

### Strategy 2: Block-Sparse Storage

```cpp
// solver/block_sparse.hpp
#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <map>

namespace gem {

template<int N>
class BlockSparseMatrix {
public:
    using Block = Eigen::MatrixXd;
    using BlockIndex = std::pair<int, int>;

    BlockSparseMatrix(const std::vector<int>& block_sizes)
        : block_sizes_(block_sizes) {
        offsets_.resize(block_sizes.size() + 1);
        offsets_[0] = 0;
        for (size_t i = 0; i < block_sizes.size(); ++i) {
            offsets_[i + 1] = offsets_[i] + block_sizes[i];
        }
        total_size_ = offsets_.back();
    }

    void set_block(int i, int j, Block&& block) {
        blocks_[{i, j}] = std::move(block);
        if (i != j) {
            blocks_[{j, i}] = blocks_[{i, j}].transpose();
        }
    }

    Eigen::VectorXd operator*(const Eigen::VectorXd& x) const {
        Eigen::VectorXd y = Eigen::VectorXd::Zero(total_size_);

        for (const auto& [idx, block] : blocks_) {
            auto [bi, bj] = idx;
            y.segment(offsets_[bi], block_sizes_[bi]) +=
                block * x.segment(offsets_[bj], block_sizes_[bj]);
        }

        return y;
    }

private:
    std::vector<int> block_sizes_;
    std::vector<int> offsets_;
    int total_size_;
    std::map<BlockIndex, Block> blocks_;
};

} // namespace gem
```

---

## 13. Main Program

```cpp
// apps/gem_run.cpp
#include <gem/gem.hpp>
#include <iostream>
#include <chrono>

#ifdef USE_MPI
#include <mpi.h>
#endif

int main(int argc, char* argv[]) {
#ifdef USE_MPI
    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
    int rank = 0;
#endif

    try {
        if (argc < 2) {
            if (rank == 0) {
                std::cerr << "Usage: " << argv[0] << " <config.yaml>\n";
            }
#ifdef USE_MPI
            MPI_Finalize();
#endif
            return 1;
        }

        // Load configuration
        auto config = gem::Config::from_file(argv[1]);
        config.validate();

        if (rank == 0) {
            std::cout << "GEM Calculation: " << config.system.name << "\n";
            std::cout << "Hadron type: " << config.system.hadron_type << "\n";
            std::cout << "Particles: " << config.system.num_particles() << "\n";
#ifdef USE_MPI
            std::cout << "Mode: MPI + OpenMP (HPC)\n";
#else
            std::cout << "Mode: OpenMP (Personal Computer)\n";
#endif
        }

        // Dispatch based on particle count
        gem::EigenResult result;
        auto start = std::chrono::high_resolution_clock::now();

        switch (config.system.num_particles()) {
            case 2:
                result = gem::run_calculation<2>(config);
                break;
            case 3:
                result = gem::run_calculation<3>(config);
                break;
            case 4:
                result = gem::run_calculation<4>(config);
                break;
            case 5:
                result = gem::run_calculation<5>(config);
                break;
            default:
                throw std::runtime_error("Unsupported particle count");
        }

        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        // Output results (rank 0 only)
        if (rank == 0) {
            std::cout << "\nResults:\n";
            std::cout << "  Basis dimension: " << result.eigenvectors.rows() << "\n";
            std::cout << "  Truncated states: " << result.num_truncated << "\n";
            std::cout << "  Condition number: " << result.condition_number << "\n";
            std::cout << "\nEigenvalues (MeV):\n";
            for (int i = 0; i < result.eigenvalues.size(); ++i) {
                std::cout << "  E_" << i << " = "
                          << result.eigenvalues(i) * 1000 << "\n";
            }
            std::cout << "\nComputation time: " << duration.count() << " ms\n";

            // Save results
            gem::HDF5Writer writer(config.output.directory);
            writer.write(result, config);
        }

#ifdef USE_MPI
        MPI_Finalize();
#endif
        return 0;

    } catch (const std::exception& e) {
        if (rank == 0) {
            std::cerr << "Error: " << e.what() << "\n";
        }
#ifdef USE_MPI
        MPI_Finalize();
#endif
        return 1;
    }
}
```

### Runner Template

```cpp
// runner.hpp
#pragma once
#include "config.hpp"
#include "solver/eigensolver.hpp"
#include "hamiltonian/builder.hpp"

#ifdef USE_MPI
#include "hpc/mpi_builder.hpp"
#endif

namespace gem {

template<int N>
EigenResult run_calculation(const Config& config) {
    // 1. Build basis
    auto spatial_basis = GaussianBasis<N>::generate(config.basis);
    auto spin_states = SpinStates<N>(
        Spin{static_cast<int>(config.system.quantum_numbers.total_spin * 2)});
    auto color_states = ColorStates<N>::singlet_states();

    auto combined_basis = CombinedBasis<N>::build(
        spatial_basis, spin_states, color_states
    );

    int rank = 0;
#ifdef USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    if (rank == 0) {
        std::cout << "  Spatial basis: " << spatial_basis.size() << "\n";
        std::cout << "  Spin states: " << spin_states.num_states() << "\n";
        std::cout << "  Color states: " << color_states.size() << "\n";
        std::cout << "  Combined basis: " << combined_basis.size() << "\n";
    }

    // 2. Build Hamiltonian (choose builder based on mode)
    Eigen::MatrixXd H, S;

#ifdef USE_MPI
    if (config.is_hpc_mode()) {
        MPIHamiltonianBuilder<N> builder(config, combined_basis);
        std::tie(H, S) = builder.build();
    } else
#endif
    {
        HamiltonianBuilder<N> builder(config, combined_basis);
        auto result = builder.build();
        H = std::move(result.H);
        S = std::move(result.S);
    }

    // 3. Solve eigenvalue problem (rank 0 only for MPI)
    if (rank == 0) {
        auto solver = create_solver(config.solver);
        return solver->solve(H, S, config.solver.num_eigenvalues);
    }

    return {};
}

} // namespace gem
```

---

## 14. Files to Create

| File | Purpose |
|------|---------|
| `CMakeLists.txt` | Build configuration |
| `include/gem/gem.hpp` | Main include header |
| `include/gem/config.hpp` | Configuration structures |
| `include/gem/constants.hpp` | Physical constants, AL1 parameters |
| `include/gem/core/types.hpp` | Core data types |
| `include/gem/core/quantum_numbers.hpp` | Quantum number types |
| `include/gem/basis/gaussian.hpp` | Gaussian basis generator |
| `include/gem/basis/combined.hpp` | Combined basis builder |
| `include/gem/spatial/jacobi.hpp` | Jacobi coordinates |
| `include/gem/spatial/cholesky.hpp` | Cholesky matrix elements |
| `include/gem/discrete/spin.hpp` | Spin wave functions |
| `include/gem/discrete/color.hpp` | Color wave functions |
| `include/gem/discrete/isospin.hpp` | Isospin wave functions |
| `include/gem/antisymmetry/permutation.hpp` | Permutation utilities |
| `include/gem/antisymmetry/antisymmetrize.hpp` | Antisymmetrization |
| `include/gem/hamiltonian/potential.hpp` | AL1 potential implementation |
| `include/gem/hamiltonian/three_body.hpp` | Three-body force (baryons only) |
| `include/gem/hamiltonian/kinetic.hpp` | Kinetic energy operator |
| `include/gem/hamiltonian/builder.hpp` | OpenMP Hamiltonian builder |
| `include/gem/hpc/mpi_builder.hpp` | MPI Hamiltonian builder |
| `include/gem/hpc/checkpoint.hpp` | Checkpoint/restart |
| `include/gem/hpc/slurm.hpp` | SLURM script generator |
| `include/gem/solver/eigensolver.hpp` | Solver interface |
| `include/gem/solver/lapack_solver.hpp` | LAPACK solver (primary) |
| `include/gem/solver/spectra_solver.hpp` | Spectra iterative solver |
| `include/gem/solver/elpa_solver.hpp` | ELPA distributed solver (HPC) |
| `include/gem/output/hdf5_writer.hpp` | HDF5 output |
| `apps/gem_run.cpp` | Main executable |

---

## 15. Verification Plan

1. **Unit tests**: Each module independently (Catch2)
2. **Meson validation**: π (138 MeV), ρ (775 MeV), J/ψ (3097 MeV)
3. **Baryon validation**: p (938 MeV), Δ (1232 MeV)
4. **Numerical checks**:
   - Hermiticity: `||H - H^T|| < 10^-14`
   - Orthonormality: `||X^T S X - I|| < 10^-10`
   - Sum rules: spin-spin, color-color

---

## 16. Implementation Order

### Phase 1: Core Framework (Meson)

1. CMake setup, Eigen + LAPACK integration
2. Config parsing (yaml-cpp)
3. Gaussian basis generator
4. Jacobi coordinates (N=2)
5. Spatial matrix elements
6. LAPACK eigenvalue solver (MKL or OpenBLAS)
7. Meson validation

### Phase 2: Baryon Extension

1. Template generalization (N=2,3)
2. Spin coupling and Clebsch-Gordan
3. Color singlet states
4. Three-body force (baryons only: $V_{3body} = -C_3 / (m_1 m_2 m_3)$)
5. Antisymmetrization
6. OpenMP parallelization
7. Baryon validation

### Phase 3: Tetraquark

1. N=4 template instantiation
2. Discrete state reduction
3. Iterative solver (Spectra)
4. Performance optimization
5. Tetraquark validation

### Phase 4: HPC/Pentaquark

1. MPI Hamiltonian builder
2. Hybrid MPI+OpenMP
3. Checkpoint/restart system
4. SLURM integration
5. ScaLAPACK/ELPA solver (optional)
6. Pentaquark validation

### Phase 5: Production

1. Comprehensive testing
2. Performance benchmarking
3. Documentation
4. Example configurations
