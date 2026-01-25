#include "simulation.hpp"

int main(){
    // Simulation parameters
    int nParticles = 512;
    Real timeStep = 0.01;
    Real length = 1000.0f;
    int nSteps = 10000;
    // Run !
    Simulation sim(nParticles, timeStep, length);
    sim.integrateEulerGPU(nSteps, 1, 128, "../data/GPU/energyEulerFloat.dat", "../data/GPU/energyEulerFloat.xyz");
}