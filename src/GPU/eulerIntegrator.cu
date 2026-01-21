#include <cuda_runtime.h>
#include "simulation.hpp"


__global__ void eulerKernel(float* x, float* y, float* z, float* vx, float* vy, float* vz, const float* fx, const float* fy, const float* fz, const float* mass, int N, float dt, float L) {
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

    //Load positions and compute next proposal positions
    float next_x = x[i] + current_vx * dt;
    float next_y = y[i] + current_vy * dt;
    float next_z = z[i] + current_vz * dt;
    // Apply PBC
    if (next_x < 0.0f) next_x += L;
    else if (next_x >= L) next_x -= L;
    if (next_y < 0.0f) next_y += L;
    else if (next_y >= L) next_y -= L;
    if (next_z < 0.0f) next_z += L;
    else if (next_z >= L) next_z -= L;
    // Now write to memory
    x[i] = next_x;
    y[i] = next_y;
    z[i] = next_z;
    //Velocity update
    vx[i] += ax * dt;
    vy[i] += ay * dt;
    vz[i] += az * dt;
}


void Simulation::integrateEulerGPU(int nSteps, int saveEvery, int threadsPerBlock, std::string saveEnergy, std::string saveTrajectory) {
    std::ofstream energyFile; energyFile.open(saveEnergy);
    std::ofstream trajectoryFile; trajectoryFile.open(saveTrajectory);
    trajectoryFile.precision(5);
    trajectoryFile << std::scientific;

    std::vector<float> potEnergy; potEnergy.resize(m_N);
    std::vector<float> kinEnergy; kinEnergy.resize(m_N);

    int blocksPerGrid = (m_N + threadsPerBlock - 1) / threadsPerBlock;

    // Main loop
    int time = 0;
    while (time < nSteps) {
        // First of all, compute forces on GPU. A kernel will be launched and will fill d_fx, d_fy, d_fz
        computeForces<<<blocksPerGrid, threadsPerBlock>>>(d_x, d_y, d_z, d_mass, d_fx, d_fy, d_fz, m_N, m_L);
        // Now perform the integration step (in theory, this kernel will wait for the previous to end, so that force-computation and numerical integration are sequential) 
        eulerKernel<<<blocksPerGrid, threadsPerBlock>>>(d_x, d_y, d_z, d_vx, d_vy, d_vz, d_fx, d_fy, d_fz, d_mass, m_N, m_dt, m_L);
    
        //Here we decide if we want to save once every while (done on the CPU)
        bool isSavingStep = (saveEnergy != "") && ((time + 1) % saveEvery == 0);
        if (isSavingStep) {
            size_t bytes = m_N * sizeof(Real);
            cudaMemcpy(m_x.data(), d_x, bytes, cudaMemcpyDeviceToHost);
            cudaMemcpy(m_y.data(), d_y, bytes, cudaMemcpyDeviceToHost);
            cudaMemcpy(m_z.data(), d_z, bytes, cudaMemcpyDeviceToHost);
            trajectoryFile << m_N << "\n";
            trajectoryFile << "Lattice=\"" << m_L << " 0.0 0.0 " << "0.0 " << m_L << " 0.0 " << "0.0 0.0 " << m_L << "\" "<< "Properties=species:S:1:pos:R:3 "<< "Time=" << time << "\n";
            for (int i = 0; i < m_N; ++i) {
                trajectoryFile << "H " << m_x[i] << " "  << m_y[i] << " " << m_z[i] << "\n"; 
            }
            // Now the enery computation
            computeEnergy<<<blocksPerGrid, threadsPerBlock>>>(d_x, d_y, d_z, d_vx, d_vy, d_vz, d_mass, d_potEnergy, d_kinEnergy, m_N, m_L);
            cudaMemcpy(potEnergy.data(), d_potEnergy, bytes, cudaMemcpyDeviceToHost);
            cudaMemcpy(kinEnergy.data(), d_kinEnergy, bytes, cudaMemcpyDeviceToHost);
            // The sum can be computed on the CPU. It's just a O(N) computations once in a while (and the execution is asynchronous, however)
            float k = 0; float u = 0;
            for (int w = 0; w < m_N; w++){
                k += kinEnergy[w];
                u += 0.5f*potEnergy[w];
            }
            energyFile << time * m_dt << " " << k << " " << u << " " << (k + u) << "\n";
        }
        time++;
    }
    cudaDeviceSynchronize();
}
