#include <cuda_runtime.h>
#include "simulation.hpp"

// First step of the Verlet 
__global__ void verletFirstStepKernel(Real* x, Real* y, Real* z,Real* vx, Real* vy, Real* vz,const Real* fx, const Real* fy, const Real* fz, const Real* mass, int N, Real dt, Real L) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;
    // Retrieve accelerations 
    Real inv_mass = 1.f/ mass[i];
    Real ax = fx[i]* inv_mass;
    Real ay = fy[i] * inv_mass;
    Real az = fz[i] *inv_mass;
    //First half-kick
    vx[i] += 0.5f * ax * dt;
    vy[i] += 0.5f * ay * dt;
    vz[i] += 0.5f * az * dt;
    //Update position
    Real next_x = x[i] + vx[i] * dt;
    Real next_y = y[i] + vy[i] * dt;
    Real next_z = z[i] + vz[i] * dt;
    // PBC
    if (next_x < 0.0f) next_x += L; else if (next_x >= L) next_x -= L;
    if (next_y < 0.0f) next_y += L; else if (next_y >= L) next_y -= L;
    if (next_z < 0.0f) next_z += L; else if (next_z >= L) next_z -= L;
    //Write positions
    x[i] = next_x;
    y[i] = next_y;
    z[i] = next_z;
}

// Second step in Velocity Verlet
__global__ void verletSecondStepKernel(Real* vx, Real* vy, Real* vz, const Real* fx, const Real* fy, const Real* fz, const Real* mass, int N, Real dt) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;
    // update accelerations. I divided this part from the previous one to enforce synchronization!!
    Real inv_mass = 1.f / mass[i];
    Real ax = fx[i]*inv_mass;
    Real ay = fy[i]*inv_mass;
    Real az = fz[i]*inv_mass;
    // Update velocities, finally
    vx[i] += 0.5f * ax * dt;
    vy[i] += 0.5f * ay * dt;
    vz[i] += 0.5f * az * dt;
}

void Simulation::integrateVerletGPU(int nSteps, int saveEvery, int threadsPerBlock, std::string saveEnergy, std::string saveTrajectory) {
    std::ofstream energyFile; energyFile.open(saveEnergy);
    std::ofstream trajectoryFile; trajectoryFile.open(saveTrajectory);
    trajectoryFile.precision(5);
    trajectoryFile << std::scientific;
    energyFile.precision(7);
    energyFile << std::scientific;

    std::vector<Real> potEnergy; potEnergy.resize(m_N);
    std::vector<Real> kinEnergy; kinEnergy.resize(m_N);

    int blocksPerGrid = (m_N + threadsPerBlock - 1) / threadsPerBlock;
    // Main loop
    int time = 0;
    // In Velocity Verlet, we have to compute first the forces
    computeForcesTiling<<<blocksPerGrid, threadsPerBlock>>>(d_x, d_y, d_z, d_mass, d_fx, d_fy, d_fz, m_N, m_L);
    while (time < nSteps) {
        verletFirstStepKernel<<<blocksPerGrid, threadsPerBlock>>>(d_x, d_y, d_z, d_vx, d_vy, d_vz, d_fx, d_fy, d_fz, d_mass, m_N, m_dt, m_L);
        // First of all, compute forces on GPU. A kernel will be launched and will fill d_fx, d_fy, d_fz
        computeForcesTiling<<<blocksPerGrid, threadsPerBlock>>>(d_x, d_y, d_z, d_mass, d_fx, d_fy, d_fz, m_N, m_L);
        // Now perform the integration step (in theory, this kernel will wait for the previous to end, so that force-computation and numerical integration are sequential) 
        verletSecondStepKernel<<<blocksPerGrid, threadsPerBlock>>>(d_vx, d_vy, d_vz, d_fx, d_fy, d_fz, d_mass, m_N, m_dt);
    
        //Here we decide if we want to save once every while (done on the CPU). It does not need to know anything from the GPU, so it's ok if asynchronous
        bool isSavingStep = (saveEnergy != "") && ((time + 1) % saveEvery == 0);
        if (isSavingStep) {
            // If saving, copy trajectories
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