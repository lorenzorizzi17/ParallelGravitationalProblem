#ifndef SIMULATION_H
#define SIMULATION_H

#include <vector>
#include <string>

// To quickly switch precision
using Real = double; 

class Simulation {
    private:
    int m_N;        
    Real m_dt;          
    Real m_L;            
    // Fundamental data structures
    std::vector<Real> m_x, m_y, m_z;      
    std::vector<Real> m_vx, m_vy, m_vz;    
    std::vector<Real> m_fx, m_fy, m_fz;    
    std::vector<Real> m_mass;  

    public:
    // Constructor (random masses, velocities, positions); implememented in simulation.cpp
    Simulation(int n_particles, Real timeStep, Real length);

    // Crucial part, computes the forces (implemented in computeForces.cpp)
    double computeForces(bool);

    // The numerical integrators (choose which to use in main; implemented in numericalIntegrator.cpp)
    void integrateEuler(int nSteps, int saveEvery, std::string saveEnergy, std::string saveTrajectory);
    void integrateVerlet(int nSteps);

    // Methods for Analysis/Debug
    Real getTotalEnergy() const;
    void saveTrajectory(const std::string& filename, int step) const;
};

#endif