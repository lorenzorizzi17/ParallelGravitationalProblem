#include "simulation.hpp"

int main(){
    // Simulation parameters
    int nParticles = 2;
    Real timeStep = 0.05;
    Real length = 1000.0f;
    int nSteps = 100000;
    // Run !
    Simulation sim(timeStep, length);
    sim.integrateVerlet(nSteps, 1, "../../data/energy_N2_Verlet.txt", "../../data/trajectory_N2_Verlet.xyz");
}