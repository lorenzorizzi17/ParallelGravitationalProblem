#include <random>
#include "simulation.hpp"
#define THRESHOLD 0.01f

// Distance using minimum image criterion
double distanceMIC(Real pos1[3], Real pos2[3], Real L) {
    double dx = *pos1 - *pos2;
    double dy = *(pos1+1) - *(pos2+1);
    double dz = *(pos1+2) - *(pos2+2);

    if (dx > L * 0.5) dx -= L;
    else if (dx < -L * 0.5) dx += L;
    if (dy > L * 0.5) dy -= L;
    else if (dy < -L * 0.5) dy += L;
    if (dz > L * 0.5) dz -= L;
    else if (dz < -L * 0.5) dz += L;

    return dx*dx + dy*dy + dz*dz;
}

// Basic constructor. As requested in the guidelines,
// --> Each mass is a random number between 1 and 10
// --> Particles position are random; avoid particles dangerously close to each other (defined as THRESHOLD*LinearSize)
// --> Initial velocities are randomly distributed in [-1,1]
Simulation::Simulation(int n_particles, Real timeStep, Real length) : m_N(n_particles), m_dt(timeStep), m_L(length) {
    // Set masses
    for (int i = 0; i < m_N; ++i) {
        Real m = 1.0 + static_cast <Real> (rand()) /( static_cast <Real> (RAND_MAX/(9.0f)));   // RAND is a poor PRNG, but should be ok here
        m_mass.push_back(m);
    }
    // Set velocities
    for (int i = 0; i < m_N; ++i) {
        Real vx = -1.0f + static_cast <Real> (rand()) /( static_cast <Real> (RAND_MAX/(2.0f)));
        Real vy = -1.0f + static_cast <Real> (rand()) /( static_cast <Real> (RAND_MAX/(2.0f)));
        Real vz = -1.0f + static_cast <Real> (rand()) /( static_cast <Real> (RAND_MAX/(2.0f)));
        m_vx.push_back(vx);
        m_vy.push_back(vy);
        m_vz.push_back(vz);
    }
    // Set positions
    for(int i = 0; i < m_N; ++i) {
        Real temp_pos[3];
        bool overlap = true;
        while(overlap) {
            for(int k = 0; k < 3; ++k) {
                temp_pos[k] = (static_cast<double>(rand()) / RAND_MAX) * m_L;  // RANDOM ENGINE, as above, but should be ok
            }
            // Check for overlaps
            overlap = false;
            for(int j = 0; j < i; ++j) {
                Real other_pos[3] = {m_x[j], m_y[j], m_z[j]};
                double d2 = distanceMIC(temp_pos, other_pos, m_L);
                if(d2 < THRESHOLD * THRESHOLD * m_L * m_L) { // Overlap detected
                    overlap = true;
                    break;
                }
            }
        }
        m_x.push_back(temp_pos[0]);
        m_y.push_back(temp_pos[1]);
        m_z.push_back(temp_pos[2]);
    }
    // Initialize forces to zero
    m_fx.resize(m_N, 0.0f);
    m_fy.resize(m_N, 0.0f);
    m_fz.resize(m_N, 0.0f);
    // All done, good to go
}


Real Simulation::getTotalEnergy() const {
    Real kineticEnergy = 0.0f;
    Real potentialEnergy = 0.0f;

    // Kinetic energy
    for(int i = 0; i < m_N; ++i) {
        Real v2 = m_vx[i]*m_vx[i] + m_vy[i]*m_vy[i] + m_vz[i]*m_vz[i];
        kineticEnergy += 0.5f * m_mass[i] * v2;
    }

    // Potential energy (naive O(N^2) approach, again)
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

            Real r = std::sqrt(dx*dx + dy*dy + dz*dz + 1e-10f); // Avoid division by zero
            potentialEnergy -= (m_mass[i] * m_mass[j]) / r; // G=1 in our units
        }
    }

    return kineticEnergy + potentialEnergy;
}