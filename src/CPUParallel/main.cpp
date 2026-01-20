#include "simulation.hpp"

int main(){
    omp_set_num_threads(8);
    // Simulation parameters
    int nParticles = 500;
    Real timeStep = 0.01;
    Real length = 1000.0f;
    int nSteps = 100000;
    // Run !
    Simulation sim(timeStep, length);
    sim.integrateEuler(nSteps, 1, "../../data/CPUParallel/energy.dat", "../../data/CPUParallel/trajectory.xyz");
}