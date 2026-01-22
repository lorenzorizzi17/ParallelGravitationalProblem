#include <cuda_runtime.h>
#include "simulation.hpp"

__global__ void computeForces(const Real*  x,  const Real*  y,  const Real*  z, const Real*  mass, Real* fx,  Real* fy,  Real* fz, int N, const Real L) {
    // Identify the current thread index
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    // Bound guard
    if (i >= N) return;

    // Load positions from memory in the register
    Real my_x = x[i];
    Real my_y = y[i];
    Real my_z = z[i];
    // Load mass from memory
    Real my_mass = mass[i]; 
    // Local private accumulators for the acceleration
    Real ax = 0.0f;
    Real ay = 0.0f;
    Real az = 0.0f;

    Real half_L = 0.5f*L;
    // Loop on all other particles
    for (int j = 0; j < N; j++) {
        if (i == j) continue; // Divergence?

        // Compute distance
        Real dx = x[j] - my_x;
        Real dy = y[j] - my_y;
        Real dz = z[j] - my_z;
        //MIC
        if (dx > half_L)       dx -= L;
        else if (dx < -half_L) dx += L;
        if (dy > half_L)       dy -= L;
        else if (dy < -half_L) dy += L;
        if (dz > half_L)       dz -= L;
        else if (dz < -half_L) dz += L;
        // Floating point operations!
        Real distSqr = dx*dx + dy*dy + dz*dz + 1e-1f; 
        Real invDist = rsqrtf(distSqr);
        Real force = mass[j] * invDist * invDist * invDist;
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

    Real my_x = x[i];
    Real my_y = y[i];
    Real my_z = z[i];
    
    Real potential = 0.0f;
    Real half_L = 0.5f * L;

    for (int j = 0; j < N; j++) {
        if (i == j) continue;

        Real dx = x[j] - my_x;
        Real dy = y[j] - my_y;
        Real dz = z[j] - my_z;

        // MIC
        if (dx > half_L)       dx -= L;
        else if (dx < -half_L) dx += L;
        if (dy > half_L)       dy -= L;
        else if (dy < -half_L) dy += L;
        if (dz > half_L)       dz -= L;
        else if (dz < -half_L) dz += L;


        Real distSqr = dx*dx + dy*dy + dz*dz + 1e-1f; 
        
        Real invDist = rsqrtf(distSqr); 
        potential -= mass[j] * invDist; 
    }
    potEnergyArray[i] = potential * mass[i];
    // Kinetic part
    Real kinetic = 0.5f*mass[i]*(vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i]);
    kinEnergyArray[i] = kinetic;
}