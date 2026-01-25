# CPU parallel

The code implemented here is an attempt to parallelize the $N$-body gravitational problem leveraging the multi-threaded execution offered by most CPUs.  We will use OpenMP. The structure of the code is similar to the one reported in `CPUSequential/README.md`, so we are only going to highlight the differences:

### Force computation, version 1
When the program calls one of the numerical integrator functions (either `integrateVerlet()` or `integrateEuler()`), a `#pragma omp parallel` section is encountered. When this happens, the process spawns a collection of $k$ threads that will run concurrently. Note that a set of three variables (`time`, `kineticEnergy`, `potentialEnergy`) is shared between all threads. 

Now all $k$ threads will run the method `Real computeForce(...)`: inside this function, a `#pragma omp for` directive tells each thread to only consider a portion of the whole particle datasets. Since the load is essentially uniform, we used a static scheduling. For example, if $k = 8$ and $N = 80$, then the first thread will compute the forces associated to the particles whose index is between $0$ and $9$. Threads will access (in read-only mode) the whole vector of particles position, but this is not an issue. At the end of the computation, each thread writes on the correct index of the shared force vector: again, provided the index is different, this is not a race condition. Each thread will also compute its local potential energy (return value of the function). At the end of the main loop, we will need to sum up all of those partial energy to obtain the global one. 

There is a difference in how we handle the Verlet and the Euler integration scheme. In general, we first have to compute the forces updating the shared vector `m_fx, m_fy, m_fz`; only then, we can perform the numerical integration using the updated acceleration values. Hence, we need to insert a barrier at the end of the `computeForce()` method (which is implicit, actually). In the Verlet scheme, we must be more careful, since the split of the velocity step requires two barriers. At the end of each time-step, one thread (the master) will be in charge of updating the shared variables `time, kineticEnergy, potentialEnergy`.

### Force computation, version 2
In the previous implementation, each thread has to compute all of the $N-1$ interactions with all the other particles, hence giving up on Newton's third law used in the sequential implementation (for a total of $N(N-1)$ computations). This is necessary since we want to avoid race conditions at any cost.

A possible strategy to recover Newton's third law is to use atomic operations (and dynamic scheduling). This is done in `computeForce2()`.