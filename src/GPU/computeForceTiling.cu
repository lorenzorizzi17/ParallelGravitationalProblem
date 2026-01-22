#include <cuda_runtime.h>
#include "simulation.hpp"

#define BLOCK_SIZE 128 

__global__ void computeForcesTiling(const Real* x, const Real*  y, const Real*  z, const Real*  mass, Real*  fx, Real*  fy, Real*  fz, int N, const Real L) {
    // idx of the thread
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    // Load the data relevant to the i-th particle in the thred register
    Real my_x, my_y, my_z, my_mass;
    Real ax = 0.0f;
    Real ay = 0.0f;
    Real az = 0.0f;

    // Controllo se il thread 'i' è valido per caricare i propri dati
    if (i < N) {
        my_x = x[i];
        my_y = y[i];
        my_z = z[i];
        my_mass = mass[i];
    } else {
        my_x = 0.0f; 
        my_y = 0.0f; 
        my_z = 0.0f; 
        my_mass = 0.0f;
    }
    // Now the thread will upload those values in the shared region (shared between threads in a block)
    // Create a pointer to that region ...
    __shared__ Real sh_x[BLOCK_SIZE];
    __shared__ Real sh_y[BLOCK_SIZE];
    __shared__ Real sh_z[BLOCK_SIZE];
    __shared__ Real sh_mass[BLOCK_SIZE];

    Real half_L = 0.5f * L;
    
    int numTiles = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;

    // Now fill the shared vectors. Each block will update a tile at a time
    for (int tile = 0; tile < numTiles; tile++) {

        // Upload the particle i in the shared region
        int t_idx = tile * BLOCK_SIZE + threadIdx.x;

        if (t_idx < N) {
            sh_x[threadIdx.x] = x[t_idx];
            sh_y[threadIdx.x] = y[t_idx];
            sh_z[threadIdx.x] = z[t_idx];
            sh_mass[threadIdx.x] = mass[t_idx];
        } else {
            // padding
            sh_x[threadIdx.x] = 0.0f;
            sh_y[threadIdx.x] = 0.0f;
            sh_z[threadIdx.x] = 0.0f;
            sh_mass[threadIdx.x] = 0.0f; 
        }
        // Now, fundamental, wait for all threads in a block to complete the shared upload
        __syncthreads();

        // And finally we can compute forces
        if (i < N) {
            // For all threads in the block (hence all the particles in a block)
            for (int k = 0; k < BLOCK_SIZE; k++) {
                
                //if ((tile * BLOCK_SIZE + k) == i) continue; 
                // Retrieve the position from shared memory and compure fistance
                Real dx = sh_x[k] - my_x;
                Real dy = sh_y[k] - my_y;
                Real dz = sh_z[k] - my_z;
                // MIC (Minimum Image Convention)
                if (dx > half_L)       dx -= L;
                else if (dx < -half_L) dx += L;
                if (dy > half_L)       dy -= L;
                else if (dy < -half_L) dy += L;
                if (dz > half_L)       dz -= L;
                else if (dz < -half_L) dz += L;
                // 
                Real distSqr = dx*dx + dy*dy + dz*dz + 1e-1f; 
                Real invDist = rsqrtf(distSqr);
                Real invDist3 = invDist * invDist * invDist;
                
                // Mass is read from shared memory
                Real force = sh_mass[k] * invDist3; 

                ax += force * dx;
                ay += force * dy;
                az += force * dz;
            }
        }
        __syncthreads();
    }

    // 5. Scrittura finale in memoria globale
    if (i < N) {
        fx[i] = ax * my_mass;
        fy[i] = ay * my_mass;
        fz[i] = az * my_mass;
    }
}

// Compute potential energy as a distinct function, since it should run only a fraction of the times as a kernel
// I'm repeting some computations, but this is to avoid a  if(saveEnergy){} in the computeForce() that may slow down performances
// In the CPU, this should not be a big problem: CPU are smart and can easily predict that the if statement is false most of the times
// GPUs cores are dumber (no branch prediction)


// GPU NAIVE, N = 1024: 0.38
// GPU NAIVE, N = 4096: 6 s