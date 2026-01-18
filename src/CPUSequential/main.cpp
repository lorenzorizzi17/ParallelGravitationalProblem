#include "simulation.hpp"

int main(){
    int n_particles = 128;
    Real timeStep = 0.005;
    Real length = 1000.0f;
    Simulation sim(n_particles, timeStep, length);

    int n_steps = 100000;
    sim.integrateEuler(n_steps, 1, "../../data/energy.txt", "../../data/trajectory.xyz");
}