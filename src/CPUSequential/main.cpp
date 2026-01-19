#include "simulation.hpp"

int main(){
    // Simulation parameters
    int nParticles = 500;
    Real timeStep = 0.01;
    Real length = 1000.0f;
    int nSteps = 1000;
    // Run !
    Simulation sim(nParticles, timeStep, length);
    sim.integrateVerlet(nSteps, 1, "", "");
}