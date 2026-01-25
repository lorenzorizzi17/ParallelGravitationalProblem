# CPU sequential

The code implements a naive $O(N^2)$ approach for force calculations with support for both Euler and Verlet integration schemes.

### Code Structure


#### `simulation.hpp`

Header file defining the `Simulation` class. The private section of the class contains
- `m_N`: Number of particles
- `m_dt`: Time step size
- `m_L`: Box length for periodic boundary conditions
- `m_x, m_y, m_z`: Particle position vectors. We use C++ template `std::vector<>`, which is essentially a heap-allocated and automatically managed array. They are slower with respect to stack-allocated native arrays, nevertheless they are needed since the size of the container $N$ is only known at runtime
- `m_vx, m_vy, m_vz`: Particle velocity vectors. Same considerations as above
- `m_fx, m_fy, m_fz`: Force vectors. Same considerations as above
- `m_mass`: Particle masses. Same considerations as above

As for the public methods:
- `Simulation(int n_particles, Real timeStep, Real length)`: Constructor with random initialization
- `Simulation(Real timeStep, Real length)`: Two-body test constructor
- `Real computeForces(bool)`: Computes gravitational forces between all particles. While doing so, it also computes the potential energy of the system (to avoid perfoming twice the same computations)
- `void integrateEuler(...)`: Euler integration scheme
- `void integrateVerlet(...)`: Verlet integration scheme

#### `simulation.cpp`
Implements the simulation constructors and initialization logic. In particular,
- `Simulation(int n_particles, Real timeStep, Real length)`: Standard constructor. It fills the vectors with $N$ elements, in particular assigning random masses between 1 and 10 and generating random positions ensuring no overlapping particles (threshold = 0.01 × box size). In addition, initializes velocities randomly in $[-1, 1]$. N.B. We used the pure C function `RAND` to generate random numbers. This would be a terrible choice in a Montecarlo simulation since `RAND` is not a professional PRNG, we here it is enough
- `Simulation(Real timeStep, Real length)`: It creates a simple two-particle system for validation. In particular,the first particle has a Large mass (200), at rest at box center while the second particle has unit mass, offset by 100 units with initial velocity parallel to an axis.

#### `computeForces.cpp`

Implements the gravitational force calculation. It's a simple naive double loop over all particle pairs leveraging Newton's third law to halve the number of computations ($0.5 N (N-1)$). The gravitational field is modified to accomodate a softening parameter. The Minimum Image Convention (MIC) is used for periodic boundary conditions. Since the function has already computed, for instance, the distance $r_{ij}$ between two particles, it's only reasonable to also compute the potential energy (and return it, optionally, according to a bool)

#### `numericalIntegrator.cpp`

Implements two integration schemes: Euler and Verlet. Both functions `integrateVerlet(...)` and `integrateEuler(...)` accepts as inputs two `std::string` defining where the user wants to save the energy (kinetic, potential, total) and the trajectories every `saveEvery` steps (another parameter of the function). The trajectory file is formatted according to the `.xyz` format, so that OVITO can easily read it

#### `main.cpp`
Entry point. Launch a simulation and runs it saving file if specified

#### `timer.cpp`
Another entry point (contains a `int main()`). Similar to `main.cpp` but contains C++ `std::chrono` utilities to accurately measure execution time. The number of particles $N$ can be specified as an `argc` 



## Compilation

To compile the main simulation:
```bash
g++ -O3 -march=native main.cpp simulation.cpp computeForces.cpp numericalIntegrator.cpp -o nbody_seq
```

To compile the timer:
```bash
g++ -O3 -march=native  timer.cpp simulation.cpp computeForces.cpp numericalIntegrator.cpp -o timer_seq
```

To run multiple times the timer, use the `timer.sh` script
