#include "simulation.hpp"
#include <cuda_runtime.h>
#include <fstream>
#include <iostream>

#include "computeForce.cu" 


__global__ void eulerKernel(float* x, float* y, float* z, float* vx, float* vy, float* vz, const float* fx, const float* fy, const float* fz, const float* mass, int N, float dt) {
    // Index of a thread
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    // Boundary check
    if (i >= N) return;

    // From the acceleration and mass vectors, compute (local) acceleration
    float ax = fx[i] / mass[i];
    float ay = fy[i] / mass[i];
    float az = fz[i] / mass[i];

    // Load velocities
    float current_vx = vx[i];
    float current_vy = vy[i];
    float current_vz = vz[i];

    // Position update
    x[i] += current_vx * dt;
    y[i] += current_vy * dt;
    z[i] += current_vz * dt;

    //Velocity update
    vx[i] += ax * dt;
    vy[i] += ay * dt;
    vz[i] += az * dt;
}


void Simulation::integrateEulerGPU(int nSteps, int saveEvery, int threadsPerBlock, std::string saveEnergy, std::string saveTrajectory) {

    int blocksPerGrid = (m_N + threadsPerBlock - 1) / threadsPerBlock;

    // Main loop
    for (int step = 0; step < nSteps; step++) {
        // First of all, compute forces on GPU. A kernel will be launched and will fill d_fx, d_fy, d_fz
        computeForces<<<blocksPerGrid, threadsPerBlock>>>(d_x, d_y, d_z, d_mass, d_fx, d_fy, d_fz, m_N);
        // Now perform the integration step (in theory, this kernel will wait for the previous to end, so that force-computation and numerical integration are sequential) 
        eulerKernel<<<blocksPerGrid, threadsPerBlock>>>(d_x, d_y, d_z, d_vx, d_vy, d_vz, d_fx, d_fy, d_fz, d_mass, m_N, m_dt);
    }
    cudaDeviceSynchronize();
}