#include "simulation.hpp"
#include <omp.h>
#include <cmath>

// This function will be called by a single thread, hence potentialEnergy is a private variable, no need for a reduction
// The #pragma omp for inside will split the work of computing forces among threads (#pragma omp parallel was called in the caller function)
Real Simulation::computeForces(bool calcEnergy) {
    Real potentialEnergy = 0.0;
    Real inv_L = 1.0 / m_L;

    // Pragma clause; schedule static (uniform computations load). This thread will compute forces (and pot energy) for a portion of particles
    #pragma omp for schedule(static)
    for (int i = 0; i < m_N; ++i) {
        // Extract the position of particle i
        Real xi = m_x[i];
        Real yi = m_y[i];
        Real zi = m_z[i];
        // prepare forces
        Real fxi = 0.0;
        Real fyi = 0.0;
        Real fzi = 0.0;
        Real ui_local = 0.0;

        // Nested loop; here j is private by default. Note that here we do not leverage Newton's third law, to avoid race conditions!!
        for (int j = 0; j < m_N; ++j) {
            if (i == j) continue; 
            Real dx = m_x[j] - xi;
            Real dy = m_y[j] - yi;
            Real dz = m_z[j] - zi;
            // Minimum Image Convention (PBC)
            dx -= m_L * std::round(dx * inv_L);
            dy -= m_L * std::round(dy * inv_L);
            dz -= m_L * std::round(dz * inv_L);
            Real r2 = dx*dx + dy*dy + dz*dz + 1e-2f; 
            Real dist = std::sqrt(r2);
            Real f_mag = (m_mass[i] * m_mass[j]) / (r2 * dist);
            fxi += f_mag * dx;
            fyi += f_mag * dy;
            fzi += f_mag * dz;
            // energy computation
            if (calcEnergy) {
                ui_local -= 0.5 * (m_mass[i] * m_mass[j]) / dist;
            }
        }
        // Write forces on the global arrays (no race condition! i is different for each thread here)
        m_fx[i] = fxi;
        m_fy[i] = fyi;
        m_fz[i] = fzi;
        // Partial potential energy
        potentialEnergy += ui_local;
    }    
    return potentialEnergy;
}


// This function will be called by a single thread, hence potentialEnergy is a private variable, no need for a reduction
// The #pragma omp for inside will split the work of computing forces among threads (#pragma omp parallel was called in the caller function)
Real Simulation::computeForces2(bool calcEnergy) {
    Real potentialEnergy = 0.0;
    Real inv_L = 1.0 / m_L;

    // Pragma clause; schedule dynamic now! The first particles are more computationally heavy, so dynamic scheduling should help load balancing (hopefully)
    #pragma omp for schedule(dynamic)
    for (int i = 0; i < m_N; ++i) {
        // Extract the position of particle i
        Real xi = m_x[i];
        Real yi = m_y[i];
        Real zi = m_z[i];
        // prepare forces
        Real fxi = 0.0;
        Real fyi = 0.0;
        Real fzi = 0.0;
        Real ui_local = 0.0;

        // Nested loop; here j is private by default. Note that here we do not leverage Newton's third law, to avoid race conditions!!
        for (int j = i+1; j < m_N; ++j) {
            Real dx = m_x[j] - xi;
            Real dy = m_y[j] - yi;
            Real dz = m_z[j] - zi;
            // Minimum Image Convention (PBC)
            dx -= m_L * std::round(dx * inv_L);
            dy -= m_L * std::round(dy * inv_L);
            dz -= m_L * std::round(dz * inv_L);
            Real r2 = dx*dx + dy*dy + dz*dz + 1e-2f; 
            Real dist = std::sqrt(r2);
            Real f_mag = (m_mass[i] * m_mass[j]) / (r2 * dist);
            fxi += f_mag * dx;
            fyi += f_mag * dy;
            fzi += f_mag * dz;
            // In this version, we update also forces on j (with atomic operations to avoid race conditions)
            #pragma omp atomic
            m_fx[j] += -f_mag * dx;
            #pragma omp atomic
            m_fy[j] += -f_mag * dy;
            #pragma omp atomic
            m_fz[j] += -f_mag * dz;
            // Energy computation
            if (calcEnergy) {
                ui_local -= 0.5 * (m_mass[i] * m_mass[j]) / dist;
            }
        }
        // Write forces on the global arrays (with atomic operations)
        #pragma omp atomic
        m_fx[i] += fxi;
        #pragma omp atomic
        m_fy[i] += fyi;
        #pragma omp atomic
        m_fz[i] += fzi;

        // Partial potential energy
        potentialEnergy += ui_local;
    }    
    return potentialEnergy;
}