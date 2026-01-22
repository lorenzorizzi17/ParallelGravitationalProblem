#include "simulation.hpp"
#include <omp.h>
#include <cmath>

// This function will be called by a single thread, hence potentialEnergy is a private variable, no need for a reduction
double Simulation::computeForces(bool calcEnergy) {
    Real potentialEnergy = 0.0;
    Real inv_L = 1.0 / m_L;
    // pragma clause; schedule static (uniform computations)
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
        // Nested loop; here j is private by default
        for (int j = 0; j < m_N; ++j) {
            if (i == j) continue; 
            Real dx = m_x[j] - xi;
            Real dy = m_y[j] - yi;
            Real dz = m_z[j] - zi;
            // Minimum Image Convention (PBC)
            dx -= m_L * std::round(dx * inv_L);
            dy -= m_L * std::round(dy * inv_L);
            dz -= m_L * std::round(dz * inv_L);
            Real r2 = dx*dx + dy*dy + dz*dz + 0.01; 
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

        // Write forces on the global arrays (no race condition!)
        m_fx[i] = fxi;
        m_fy[i] = fyi;
        m_fz[i] = fzi;
        
        // Accumulo per la reduction finale
        potentialEnergy += ui_local;
    }    
    return potentialEnergy;
}