#include <cuda_runtime.h>
#include "simulation.hpp"

#define BLOCK_SIZE 256

// Compute force vector using tiling. The idea is the following:
// We know each thread handles one particle. When this kernel is called, each thread will load into the shared memory
// the positions of a portion of all particles (a tile, indeed) for all the other threads in the block to see. When the computation is finished, another
// shared loading procedure is performed and so on
__global__ void computeForcesTiling(const Real* x, const Real*  y, const Real*  z, const Real*  mass, Real*  fx, Real*  fy, Real*  fz, int N, const Real L) {
    // Index of the thread
    int i = blockIdx.x * blockDim.x + threadIdx.x;
   
    Real my_x, my_y, my_z, my_mass;
    Real ax = 0.0f;
    Real ay = 0.0f;
    Real az = 0.0f;

    // Bound check (if the n. of particles is not a power of 2)
    if (i < N) {
        my_x = x[i];
        my_y = y[i];
        my_z = z[i];
        my_mass = mass[i];
    } else { //padding
        my_x = 0.0f; 
        my_y = 0.0f; 
        my_z = 0.0f; 
        my_mass = 0.0f;
    }

    // Now the thread will upload those values in the shared region (shared between threads in a block)
    // First create a pointer to that region ...
    __shared__ Real sh_x[BLOCK_SIZE];
    __shared__ Real sh_y[BLOCK_SIZE];
    __shared__ Real sh_z[BLOCK_SIZE];
    __shared__ Real sh_mass[BLOCK_SIZE];

    Real half_L = 0.5f * L;
    int numTiles = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;

    // Now fill the shared vectors!
    for (int tile = 0; tile < numTiles; tile++) { // Iterate on all tiles
        int t_idx = tile * BLOCK_SIZE + threadIdx.x;
        // Given a tile whose index is tile, upload the particle of index t_idx from the VRAM to the shared region
        // Note that this is a good case of memory coalescing. All threads in this warp will request data in continguous memory regions!
        // For example, here we'd have
        // sh_x[0] = x[tile*BLOCK_SIZE]
        // sh_x[1] = x[tile*BLOCK_SIZE + 1]
        // sh_x[2] = x[tile*BLOCK_SIZE + 2]
        // ...
        if (t_idx < N) {
            sh_x[threadIdx.x] = x[t_idx];
            sh_y[threadIdx.x] = y[t_idx];
            sh_z[threadIdx.x] = z[t_idx];
            sh_mass[threadIdx.x] = mass[t_idx];
        } else {// padding
            sh_x[threadIdx.x] = 0.0f;
            sh_y[threadIdx.x] = 0.0f;
            sh_z[threadIdx.x] = 0.0f;
            sh_mass[threadIdx.x] = 0.0f; 
        }
        // Now, fundamental, wait for all threads in a block to complete the shared upload
        __syncthreads();

        // And finally we can compute (partial) forces
        if (i < N) {
            // For all threads in the block (hence all the particles in a block), compute the force considering only the particles whose position is available in the shared memory
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
                Real distSqr = dx*dx + dy*dy + dz*dz + 1e-2f; 
                Real invDist = rsqrtf(distSqr);
                Real invDist3 = invDist * invDist * invDist;
                // Mass is read from shared memory
                Real force = sh_mass[k] * invDist3; 
                // update acceleration (living in the thread registers)
                ax += force * dx;
                ay += force * dy;
                az += force * dz;
            }
        }
        // Now we have completed a step inside a tile. Wait for all threads in the block to have finished and start again
        __syncthreads();
    }

    // All tiles are done. We can now write from the registers (ax, ay, az) to the memory
    if (i < N) {
        fx[i] = ax * my_mass;
        fy[i] = ay * my_mass;
        fz[i] = az * my_mass;
    }
}
