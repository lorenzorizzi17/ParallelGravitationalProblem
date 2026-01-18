// Here, we will define the computeForces function for the sequential CPU implementation (a naive O(N^2) approach)
#include <cmath>
#include "simulation.hpp"


double Simulation::computeForces(bool saveEnergy){
    // Reset forces to zero 
    std::fill(m_fx.begin(), m_fx.end(), 0.0f);
    std::fill(m_fy.begin(), m_fy.end(), 0.0f);
    std::fill(m_fz.begin(), m_fz.end(), 0.0f);

    double inv_m_L = 1.0 / m_L; // Cache inverse length
    Real U = 0.0;               // Potential energy
    // Compute pairwise forces in a naive O(N^2) loop
    for(int i = 0; i < m_N; ++i) {
        for(int j = i + 1; j < m_N; ++j) {
            Real dx = m_x[j] - m_x[i];
            Real dy = m_y[j] - m_y[i];
            Real dz = m_z[j] - m_z[i];
            // Apply MIC
            dx -= m_L * std::round(dx * inv_m_L);
            dy -= m_L * std::round(dy * inv_m_L);
            dz -= m_L * std::round(dz * inv_m_L);

            Real r2 = dx*dx + dy*dy + dz*dz + 1e-1f; // Add small term to avoid division by zero
            Real r = std::sqrt(r2);                  // COmpute square root once
            Real f_mag = (m_mass[i] * m_mass[j]) / (r2 * r); // Gravitational force magnitude; G=1 in our units
            // Update forces on both i and j
            m_fx[i] += f_mag * dx;
            m_fy[i] += f_mag * dy;
            m_fz[i] += f_mag * dz;
            m_fx[j] -= f_mag * dx;
            m_fy[j] -= f_mag * dy;
            m_fz[j] -= f_mag * dz;

            // Compute here the potential energy
            if (saveEnergy){
                U -= (m_mass[i] * m_mass[j]) / r; // G=1 in our units
            }
        }
    }    
    return U;
}

