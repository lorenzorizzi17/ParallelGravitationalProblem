#include "simulation.hpp"

int main(){
    // Simulation parameters
    int nParticles = 4096;
    Real timeStep = 0.01;
    Real length = 1000.0f;
    int nSteps = 100000;
    // Run !
    Simulation sim(nParticles, timeStep, length);
    sim.integrateVerletGPU(nSteps, 100, 256, "../data/GPU/traj_N4096.dat", "../data/GPU/traj_N4096.xyz");
}