# N-Body Gravitational simulation: Parallel Performance Analysis

A comprehensive study of gravitational N-body simulations comparing **sequential**, **OpenMP parallel**, and **CUDA GPU** implementations using different numerical integration schemes (Euler and Velocity Verlet).

## Project Overview

This project investigates the computational performance and numerical accuracy of solving the N-body gravitational problem across three distinct paradigms:

- **Sequential CPU**: Baseline single-threaded implementation using C++
- **Parallel CPU**: Multi-threaded OpenMP implementation with support for $k$ threads
- **GPU**: CUDA implementation with naive and tiled memory access patterns

The project focuses on:
1. **Numerical Accuracy**: Energy conservation across different integrators (implemented Euler and Verlet integration schemes)
2. **Performance**: Scalability analysis and speedup measurements across the three implementations

## Project Structure
The fundamental findings are reported in a Python notebook, whose only purpose is to analyze data offline and create plots and graphs (`src/report.ipynb`). The codebase can be found in the `src/` folder and is logically divided into three folders. **More details on how to compile/build can be found in the local README.md for each directory, please have a look at those**

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

## Analysis & Results

The analysis is contained in [report.ipynb](src/report.ipynb), which includes:

- **Energy conservation plots**: Tracking $E(t)$ across implementations and integrators
- **Physical analysis**: Virial theorem, cluster's density...
- **Performance Scaling**: Execution time vs. problem size (N) with speedup curves using OpenMP or CUDA
- **Throughput Metrics**: GFLOPS measurements for GPU implementations

## Compilation & Execution
Please refer to `src/CPUSequential/README.md`


