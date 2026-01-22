#include "simulation.hpp"

int main(){
    // Simulation parameters
    int nParticles = 512;
    Real timeStep = 0.01;
    Real length = 1000.0f;
    int nSteps = 10000;
    // Run !
    Simulation sim(nParticles, timeStep, length);
    sim.integrateVerlet(nSteps, 1, "../../data/CPUSequential/energyVerletDouble.dat", "../../data/CPUSequential/energyVerletDouble.xyz");
}