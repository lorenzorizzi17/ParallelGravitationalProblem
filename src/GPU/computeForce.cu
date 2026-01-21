#include <cuda_runtime.h>
#include "simulation.hpp"

__global__ void computeForces(const float*  x,  const float*  y,  const float*  z, const float*  mass, float* fx,  float* fy,  float* fz, int N, const float L) {
    // Identify the current thread index
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    // Bound guard
    if (i >= N) return;

    // Load positions from memory in the register
    float my_x = x[i];
    float my_y = y[i];
    float my_z = z[i];
    // Load mass from memory
    float my_mass = mass[i]; 
    // Local private accumulators for the acceleration
    float ax = 0.0f;
    float ay = 0.0f;
    float az = 0.0f;

    float half_L = 0.5f*L;
    // Loop on all other particles
    for (int j = 0; j < N; j++) {
        if (i == j) continue; // Divergence?

        // Compute distance
        float dx = x[j] - my_x;
        float dy = y[j] - my_y;
        float dz = z[j] - my_z;
        //MIC
        if (dx > half_L)       dx -= L;
        else if (dx < -half_L) dx += L;
        if (dy > half_L)       dy -= L;
        else if (dy < -half_L) dy += L;
        if (dz > half_L)       dz -= L;
        else if (dz < -half_L) dz += L;
        // Floating point operations!
        float distSqr = dx*dx + dy*dy + dz*dz + 1e-1f; 
        float invDist = rsqrtf(distSqr);
        float force = mass[j] * invDist * invDist * invDist;
        // Update acc
        ax += force * dx;
        ay += force * dy;
        az += force * dz;
    }

    // Write on global memory
    fx[i] = ax * my_mass;
    fy[i] = ay * my_mass;
    fz[i] = az * my_mass;
}

// Compute potential energy as a distinct function, since it should run only a fraction of the times as a kernel
// I'm repeting some computations, but this is to avoid a  if(saveEnergy){} in the computeForce() that may slow down performances
// In the CPU, this should not be a big problem: CPU are smart and can easily predict that the if statement is false most of the times
// GPUs cores are dumber (no branch prediction)
__global__ void computeEnergy(const Real* x, const Real* y, const Real* z, const Real* vx, const Real* vy, const Real* vz, const Real* mass, Real* potEnergyArray, Real* kinEnergyArray, int N, Real L) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;

    float my_x = x[i];
    float my_y = y[i];
    float my_z = z[i];
    
    float potential = 0.0f;
    float half_L = 0.5f * L;

    for (int j = 0; j < N; j++) {
        if (i == j) continue;

        float dx = x[j] - my_x;
        float dy = y[j] - my_y;
        float dz = z[j] - my_z;

        // MIC
        if (dx > half_L)       dx -= L;
        else if (dx < -half_L) dx += L;
        if (dy > half_L)       dy -= L;
        else if (dy < -half_L) dy += L;
        if (dz > half_L)       dz -= L;
        else if (dz < -half_L) dz += L;


        float distSqr = dx*dx + dy*dy + dz*dz + 1e-1f; 
        
        float invDist = rsqrtf(distSqr); 
        potential -= mass[j] * invDist; 
    }
    potEnergyArray[i] = potential * mass[i];
    // Kinetic part
    float kinetic = 0.5f*mass[i]*(vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i]);
    kinEnergyArray[i] = kinetic;
}