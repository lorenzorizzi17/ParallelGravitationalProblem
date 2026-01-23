# N-Body Gravitational simulation: Parallel Performance Analysis

A comprehensive study of gravitational N-body simulations comparing **sequential**, **OpenMP parallel**, and **CUDA GPU** implementations using different numerical integration schemes (Euler and Velocity Verlet).

## Project Overview

This project investigates the computational performance and numerical accuracy of solving the N-body gravitational problem across three distinct paradigms:

- **Sequential CPU**: Baseline single-threaded implementation using C++
- **Parallel CPU**: Multi-threaded OpenMP implementation with support for $k$ threads
- **GPU**: CUDA implementation with naive and tiled memory access patterns

The project focuses on:
1. **Numerical Accuracy**: Energy conservation across different integrators (implemented Euler and Verlet integration schemes)
2. **Performance**: Scalability analysis and speedup measurements

## Project Structure
The fundamental findings are reported in a Python notebook, whose only purpose is to analyze data offline and create plots and graphs (`src/report.ipynb`). The codebase can be found in the `src/` folder and is logically divided into three folders. More details on how to compile/build can be found in the local README.md

```
├── README.md                    # This file
│
├── src/                         # Source code for all implementations
│   ├── report.ipynb            # Jupyter notebook with analysis and visualizations
│   │
│   ├── CPUSequential/          # Single-threaded CPU implementation
│   │   ├── main.cpp            # Entry point to run a simulation
│   │   ├── simulation.hpp/cpp   # Core simulation data structure (header and impl. file)
│   │   ├── computeForces.cpp    # $O(N^2)$ force calculation
│   │   ├── numericalIntegrator.cpp  # Euler and Velocity Verlet integrators
│   │   ├── timer.cpp            # Performance timing utilities
│   │   └── timer.sh            # Shell script for benchmarking
│   │
│   ├── CPUParallel/            # OpenMP multi-threaded CPU implementation
│   │   ├── main.cpp            # Entry point to run a simulation
│   │   ├── simulation.hpp/cpp   # Core simulation data structure
│   │   ├── computeForces.cpp    # Parallelized force calculation
│   │   ├── numericalIntegrator.cpp  # Parallelized integrators
│   │   ├── timer.cpp            # Performance timing utilities
│   │   └── timer.sh            # Shell script for benchmarking
│   │
│   └── GPU/                     # CUDA GPU implementation
│       ├── main.cu              # Entry point and host code
│       ├── simulation.hpp/cu     # Simulation data structures
│       ├── computeForce.cu       # Naive GPU kernel for force calculation
│       ├── computeForceTiling.cu # Optimized kernel with shared memory tiling
│       ├── eulerIntegrator.cu    # Euler integration kernel
│       ├── verletIntegrator.cu   # Velocity Verlet integration kernels
│       └── timer.cu              # GPU timing utilities
```

## Key Implementation Details

**OpenMP (CPU Parallel)**
- Force calculation parallelized using OpenMP directives
- Each thread processes a subset of particle pairs (all pairs, $O(N^2)$ or reduced pairs using atomic operations)
- Synchronization barriers between force computation and position updates
- Supports variable thread counts (2, 4, 8, 16)

**CUDA (GPU)**
- One CUDA thread per particle
- Thread block configuration: 128 threads per block (4 warps)
- Global synchronization via kernel termination between phases
- Two kernels for Velocity Verlet (half-step predictor-corrector)
- Tiled memory access pattern for improved cache locality

### Physics Model

The gravitational potential energy and kinetic energy are computed as:

$$H(\vec{r}, \vec{v}) = \frac{1}{2}\sum_i m_i v_i^2 + \frac{1}{2}\sum_{i \neq j} \frac{G m_i m_j}{r_{ij}}$$

Energy conservation serves as the primary validation metric for numerical accuracy.

## Analysis & Results

The analysis is contained in [report.ipynb](src/report.ipynb), which includes:

- **Energy Conservation Plots**: Tracking $E(t)$ across implementations and integrators
- **Performance Scaling**: Execution time vs. problem size (N) with speedup curves
- **Precision Analysis**: Impact of Float vs. Double precision on energy stability
- **Throughput Metrics**: GFLOPS measurements for GPU implementations

## Compilation & Execution

Each implementation directory contains a `timer.sh` script for automated benchmarking across different particle counts.

### Sequential CPU
```bash
cd src/CPUSequential
g++ -O3 -march=native -ffast-math timer.cpp computeForces.cpp simulation.cpp numericalIntegrator.cpp -o timer.out
./timer.out
```

### Parallel CPU (OpenMP)
```bash
cd src/CPUParallel
g++ -O3 -march=native -ffast-math -fopenmp timer.cpp computeForces.cpp simulation.cpp numericalIntegrator.cpp -o timer.out
./timer.out
```

### GPU (CUDA)
```bash
cd src/GPU
nvcc -O3 main.cu simulation.cu computeForce.cu computeForceTiling.cu eulerIntegrator.cu verletIntegrator.cu timer.cu -o nbody.out
./nbody.out
```

