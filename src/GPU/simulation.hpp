#ifndef SIMULATION_H
#define SIMULATION_H

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

// To quickly switch precision
using Real = float; 

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

    // Pointers for device memory (living on GPU!)
    Real *d_x, *d_y, *d_z;
    Real *d_vx, *d_vy, *d_vz;
    Real *d_fx, *d_fy, *d_fz;
    Real *d_mass;
    Real *d_potEnergy;
    Real *d_kinEnergy;

    public:
    Simulation(int n_particles, Real timeStep, Real length);
    Simulation(Real timeStep, Real length);

    // The numerical integrators (choose which to use in main; implemented in numericalIntegrator.cpp)
    void integrateEulerGPU(int nSteps, int threadsPerBlock, int saveEvery, std::string saveEnergy, std::string saveTrajectory);
    void integrateVerletGPU(int nSteps, int threadsPerBlock, int saveEvery, std::string saveEnergy, std::string saveTrajectory);


};

__global__ void computeForces(const float*  x,  const float*  y,  const float*  z, const float*  mass, float* fx,  float* fy,  float* fz, int N, const float L);
__global__ void computeEnergy(const Real* x, const Real* y, const Real* z, const Real* vx, const Real* vy, const Real* vz, const Real* mass, Real* energyArray, Real* kinenergy, int N, Real L);
#endif