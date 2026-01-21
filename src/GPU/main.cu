#include "simulation.hpp"

int main(){
    // Simulation parameters
    int nParticles = 1024;
    Real timeStep = 0.01;
    Real length = 1000.0f;
    int nSteps = 100000;
    // Run !
    Simulation sim(nParticles, timeStep, length);
    sim.integrateEulerGPU(nSteps, 100, 256, "./data.dat", "./data.xyz");
}