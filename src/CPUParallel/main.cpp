#include "simulation.hpp"

int main(){
    omp_set_num_threads(2);
    // Simulation parameters
    int nParticles = 1000;
    Real timeStep = 0.01;
    Real length = 1000.0f;
    int nSteps = 1000;
    // Run !
    Simulation sim(nParticles, timeStep, length);
    sim.integrateEuler(nSteps, nSteps, "", "");
}