#include "simulation.hpp"
#include <iostream>
#include <cmath>
#include <cstdlib> // Per rand()
#include <cuda_runtime.h>

#define THRESHOLD 0.01f 


double distanceMIC(Real pos1[3], Real pos2[3], Real L) {
    double dx = *pos1 - *pos2;
    double dy = *(pos1+1) - *(pos2+1);
    double dz = *(pos1+2) - *(pos2+2);

    if (dx > L * 0.5) dx -= L;
    else if (dx < -L * 0.5) dx += L;
    if (dy > L * 0.5) dy -= L;
    else if (dy < -L * 0.5) dy += L;
    if (dz > L * 0.5) dz -= L;
    else if (dz < -L * 0.5) dz += L;

    return dx*dx + dy*dy + dz*dz;
}



Simulation::Simulation(int n_particles, Real timeStep, Real length) : m_N(n_particles), m_dt(timeStep), m_L(length) {
    #ifdef DEBUG_MODE
    std::cout << "Initializing simulation with " << m_N << " particles." << std::endl;
    #endif
    // Set masses
    for (int i = 0; i < m_N; ++i) {
        Real m = 1.0 + static_cast <Real> (rand()) /( static_cast <Real> (RAND_MAX/(9.0f)));   // RAND is a poor PRNG, but should be ok here
        m_mass.push_back(m);
    }
    // Set velocities
    for (int i = 0; i < m_N; ++i) {
        Real vx = -1.0f + static_cast <Real> (rand()) /( static_cast <Real> (RAND_MAX/(2.0f)));
        Real vy = -1.0f + static_cast <Real> (rand()) /( static_cast <Real> (RAND_MAX/(2.0f)));
        Real vz = -1.0f + static_cast <Real> (rand()) /( static_cast <Real> (RAND_MAX/(2.0f)));
        m_vx.push_back(vx);
        m_vy.push_back(vy);
        m_vz.push_back(vz);
    }
    // Set positions
    for(int i = 0; i < m_N; ++i) {
        Real temp_pos[3];
        bool overlap = true;
        while(overlap) {
            for(int k = 0; k < 3; ++k) {
                temp_pos[k] = (static_cast<double>(rand()) / RAND_MAX) * m_L;  // RANDOM ENGINE, as above, but should be ok
            }
            // Check for overlaps
            overlap = false;
            for(int j = 0; j < i; ++j) {
                Real other_pos[3] = {m_x[j], m_y[j], m_z[j]};
                double d2 = distanceMIC(temp_pos, other_pos, m_L);
                if(d2 < THRESHOLD * THRESHOLD * m_L * m_L) { // Overlap detected
                    overlap = true;
                    break;
                }
            }
        }
        m_x.push_back(temp_pos[0]);
        m_y.push_back(temp_pos[1]);
        m_z.push_back(temp_pos[2]);
    }
    #ifdef DEBUG_MODE
    std::cout << "Initialization complete." << std::endl;
    #endif
    // Initialize forces to zero
    m_fx.resize(m_N, 0.0f);
    m_fy.resize(m_N, 0.0f);
    m_fz.resize(m_N, 0.0f);
    // All done, good to go

    // GPU PART. We have to allocate device memory and transfer data
    size_t bytes = m_N * sizeof(Real);
    // Allocate memory for positions
    cudaMalloc(&d_x, bytes);
    cudaMalloc(&d_y, bytes);
    cudaMalloc(&d_z, bytes);
    // Allocate memory for velocities
    cudaMalloc(&d_vx, bytes);
    cudaMalloc(&d_vy, bytes);
    cudaMalloc(&d_vz, bytes);
    // Allocate memory for forces. This is a buffer where threads will write the force for each particle
    cudaMalloc(&d_fx, bytes);
    cudaMalloc(&d_fy, bytes);
    cudaMalloc(&d_fz, bytes);
    // Masses
    cudaMalloc(&d_mass, bytes);
    // Now transfer data from host to device
    cudaMemcpy(d_x, m_x.data(), bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(d_y, m_y.data(), bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(d_z, m_z.data(), bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(d_vx, m_vx.data(), bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(d_vy, m_vy.data(), bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(d_vz, m_vz.data(), bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(d_mass, m_mass.data(), bytes, cudaMemcpyHostToDevice);

    // Initialize forces on device to zero
    cudaMemset(d_fx, 0, bytes);
    cudaMemset(d_fy, 0, bytes);
    cudaMemset(d_fz, 0, bytes);

    #ifdef DEBUG_MODE
    std::cout << "GPU Initialization complete." << std::endl;
    #endif
}