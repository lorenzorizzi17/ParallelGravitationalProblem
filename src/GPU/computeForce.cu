#include <cuda_runtime.h>

__global__ void computeForces(const float*  x,  const float*  y,  const float*  z, const float*  mass, float* fx,  float* fy,  float* fz, int N) {
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

    // Loop on all other particles
    for (int j = 0; j < N; j++) {
        if (i == j) continue;
        // Compute distance
        float dx = x[j] - my_x;
        float dy = y[j] - my_y;
        float dz = z[j] - my_z;
        // Floating point operations!
        float distSqr = dx*dx + dy*dy + dz*dz + 1e-1f; 
        float invDist = rsqrtf(distSqr);
        float force = mass[j] * invDist * invDist * invDist;

        ax += force * dx;
        ay += force * dy;
        az += force * dz;
    }

    // Write on global memory
    fx[i] = ax * my_mass;
    fy[i] = ay * my_mass;
    fz[i] = az * my_mass;
}