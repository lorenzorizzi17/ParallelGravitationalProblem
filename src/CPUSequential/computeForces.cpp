// Here, we will define the computeForces function for the sequential CPU implementation (a naive O(N^2) approach)
#include <cmath>
#include "simulation.hpp"


void Simulation::computeForces(){
    // Reset forces to zero 
    std::fill(m_fx.begin(), m_fx.end(), 0.0f);
    std::fill(m_fy.begin(), m_fy.end(), 0.0f);
    std::fill(m_fz.begin(), m_fz.end(), 0.0f);

    // Compute pairwise forces in a naive O(N^2) loop
    for(int i = 0; i < m_N; ++i) {
        for(int j = i + 1; j < m_N; ++j) {
            Real dx = m_x[j] - m_x[i];
            Real dy = m_y[j] - m_y[i];
            Real dz = m_z[j] - m_z[i];

            // Apply MIC
            if (dx >  m_L * 0.5f) dx -= m_L;
            else if (dx < -m_L * 0.5f) dx += m_L;
            if (dy >  m_L * 0.5f) dy -= m_L;
            else if (dy < -m_L * 0.5f) dy += m_L;
            if (dz >  m_L * 0.5f) dz -= m_L;
            else if (dz < -m_L * 0.5f) dz += m_L;

            Real r2 = dx*dx + dy*dy + dz*dz + 1e-10f; // Add small term to avoid division by zero
            // Gravitational force magnitude
            Real f_mag = (m_mass[i] * m_mass[j]) / (r2 * std::sqrt(r2)); // G=1 in our units
            // Update forces on both i and j
            m_fx[i] += f_mag * dx;
            m_fy[i] += f_mag * dy;
            m_fz[i] += f_mag * dz;
            m_fx[j] -= f_mag * dx;
            m_fy[j] -= f_mag * dy;
            m_fz[j] -= f_mag * dz;
        }
    }    
}

