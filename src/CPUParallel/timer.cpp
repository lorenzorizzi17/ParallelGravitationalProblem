#include <iostream>
#include <string>
#include <vector>
#include <chrono>
#include "simulation.hpp"

int main(int argc, char* argv[]) {
    omp_set_num_threads(2);
    // Valori di default
    int nParticles = 500;
    int nSteps = 1000; 
    
    if (argc > 1) nParticles = std::atoi(argv[1]);
    if (argc > 2) nSteps = std::atoi(argv[2]);

    Real dt = 0.01f;
    Real L = 1000.0f;

    Simulation sim(nParticles, dt, L);

    auto start = std::chrono::high_resolution_clock::now();
    sim.integrateEuler(nSteps, nSteps, "", ""); // Don't save to IO
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    std::cout << nParticles << " " << elapsed.count() << std::endl;

    return 0;
}