# GPU - CUDA

The code implemented here is the CUDA version of the same code meant to paralleliza the $N$-body gravitational problem leveraging a GPU (we used a Jetson Nano to conduct all benchmarks).

Again, the structure of the code is similar to the one reported in `CPUSequential/README.md`, so we are only going to highlight the differences:

### simulation.hpp
Together with the host allocated private vectors (`m_x, m_y, ...`), we now need to instantiate a set of specular pointers (`d_x, d_y, ...`) that will point to the device memory. Instead of using `std::vector`, we used native arrays (that decay to pointers to the first element in C) according to CUDA requirement. Those pointers will indicate the region of the device memory where the relevant arrays (positions, velocities, forces and masses) are located.

Another difference is that the public methods `computeForce()` and `integrateVerlet/Euler()` have now to be defined as free functions, since a class method cannot be a kernel. However, to mantain at least the same structure as before, we choose to:
- The kernel `computeForce()` was declared in `simulation.hpp` as a free function. Its definition is in `computeForce.cu`
- The functions `integrateVerlet/integrateEuler()` are normal host functions, public methods of `class Simulation`. They are defined in `eulerIntegrator.cu` and `verletIntegrator.cu`. This was done so that the pre-processing phase (creation of `std::ofstream`, ...) can still be performed on the CPU side. Inside their definition, the force computation is triggered with the kernel `computeForce()` and a brand new kernel is created (`eulerKernel`, `verletKernel`) to perform the numerical integration. 

Let us make an example. When the function `integrateEuler()` is called, the GPU still doesn't receive any job. The function will be called from the host side and will initialize the necessary IO stuff. Then, the main loop starts:
```cpp
while (time < nSteps) {
    computeForces<<<blocksPerGrid, threadsPerBlock>>>(d_x, d_y, d_z, d_mass, d_fx, d_fy, d_fz, m_N, m_L);
    eulerKernel<<<blocksPerGrid, threadsPerBlock>>>(d_x, d_y, d_z, d_vx, d_vy, d_vz, d_fx, d_fy, d_fz, d_mass, m_N, m_dt, m_L);
    if (saveStep){...}
}
```
We have to remember that the execution between CPU and GPU is asynchronous. This implies that the CPU will run through the while loop and instantiate a collection of kernels that will be sent to the GPU. The GPU will receive those kernels and will start performing the jobs in order (computeForce($t=0$), eulerKernel($t=0$), computeForce($t=1$), eulerKernel($t=1$), ...). The CPU won't wait for the GPU to complete all its in-pending kernels since the `saveStep` evaluation only depends on the host variable `time`. For this reason, it's essential to include a `cudaDeviceSynchronize()` at the end of the while loop if we want to measure the real GPU time. 

If we want to save the energy or the particles trajectories every `saveEvery` step, then the CPU will eventually enter the `if` structure. Here, we transfer the data from device to host memory (synchronization is implicit in `cudaMemcpy`) and we then write those data to I/O. This has to be done on the CPU side, since the GPU has no direct access to the mass storage. An alternative route to avoid the host-device synchronization and data transfer would have been to allocate device memory to hold buffers whose purpose is to store the trajectories and energies of particles every `saveEvery` timesteps. Only at the end of the simulation, those buffers can be sent to the CPU and written to disk. However, trajectory files can get quite big and the VRAM only has $2$ GB, so this is risky. When evaluating the performances of our GPU-code, we simply switched off the saving process

### computeForce()
This is the crucial part of the GPU implementation. The function `computeForce()` is promoted to a CUDA kernel (__global__) and operates under a one-thread-per-particle policy. We launch a 1D grid of threads such that the total thread count covers at least $N$ particles.

Each thread will iterate through all other $N-1$ particles in the system, precisely as we did with the OpenMP implementation. In fact, unlike the CPU sequential version, this kernel give up on Newton's third law and each thread computes the total net force on its particle independently. This redundancy doubles the total floating-point operations but is necessary to avoid race conditions and the heavy overhead of atomic operations required to write to shared force accumulators


### computeForceTiling()

In the `computeForce()` implementation, each thread has to access the global memory a lot of times to access particles position. In a GPU application, this can be problematic since data transfer is usually a critical bottleneck that can significantly slow down performances (threads will wait in idle). However, in this specific case, memory is accessed in a fast and coalesced way, since a warp (= 32 threads) will usually request the position of the same particle and those positions are contiguous in memory, enabling cache hitting. For example, imagine a warp is currently running. Each thread of the warp will want to compute the force acting on its particle (from index 0 to 31, as an example). When the loop over all other $N$ particles start, each thread will request the position of the same particle (say particle 54). A single memory transaction is issued, instead of $32$ different access. But there's more: at the next iteration of the loop, threads will request the position of particles 55. Fortunately, particles position are continguous in the VRAM (since they are contained in an array), hence the GPU can efficiently cache those data in the L2 memory and greatly reduce the global memory access

However, we can try to reduce further the memory bandwidth bottleneck implementing the tiled algorithm, which leverages shared Memory to reduce global memory traffic. In particular:
1. The interaction loop is split into tiles. In each iteration, all threads in a block collaborate to load a subset of particles (`BLOCK_SIZE`) from global memory into shared memory. Since threads in a warp access contiguous memory addresses, this triggers coalesced memory Access, maximizing bandwidth efficiency.
2.  We employ `__syncthreads()` barriers to ensure the entire tile is populated and visible to the block before computation begins.
3.  Once loaded, the data in shared memory is accessed repeatedly with extremely low latency
4.  The kernel handles cases where $N$ is not a multiple of `BLOCK_SIZE` by loading "phantom particles" (mass = 0) to preserve execution uniformity without affecting physical results. Once the force computation for a tile are finished, we repeat until all tiles are considered