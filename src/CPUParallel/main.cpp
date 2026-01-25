#include "simulation.hpp"

int main(){
    omp_set_num_threads(8);
    // Simulation parameters
    int nParticles = 512;
    Real timeStep = 0.01;
    Real length = 1000.0f;
    int nSteps = 300000;
    // Run !
    Simulation sim(nParticles, timeStep, length);
    sim.integrateVerlet(nSteps, 100, "../../data/CPUParallel/virial.dat", "../../data/CPUParallel/virial.xyz");
}