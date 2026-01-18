#include <fstream>
#include "simulation.hpp"

void Simulation::integrateEuler(int n_steps, int saveEvery, std::string saveEnergy, std::string saveTrajectory) {
    // Full Euler integrator

    std::ofstream energyFile; energyFile.open(saveEnergy);
    std::ofstream trajectoryFile; trajectoryFile.open(saveTrajectory);
    trajectoryFile.precision(5);
    trajectoryFile << std::scientific;

    int time = 0;
    while (time < n_steps) {
        // First, compute forces
        this->computeForces();
        // Update positions and velocities using Euler's method (O(N) complexity)
        for (int i = 0; i < m_N; ++i) {
            // Update velocities
            m_vx[i] += (m_fx[i] / m_mass[i]) * m_dt;
            m_vy[i] += (m_fy[i] / m_mass[i]) * m_dt;
            m_vz[i] += (m_fz[i] / m_mass[i]) * m_dt;

            // Update positions
            m_x[i] += m_vx[i] * m_dt;
            m_y[i] += m_vy[i] * m_dt;
            m_z[i] += m_vz[i] * m_dt;

            // Apply MIC
            if (m_x[i] < 0) m_x[i] += m_L;
            else if (m_x[i] >= m_L) m_x[i] -= m_L;

            if (m_y[i] < 0) m_y[i] += m_L;
            else if (m_y[i] >= m_L) m_y[i] -= m_L;

            if (m_z[i] < 0) m_z[i] += m_L;
            else if (m_z[i] >= m_L) m_z[i] -= m_L;
        }
        time++;

        // Save energy and or trajectory if needed
        if ((saveEnergy!="") && (time % saveEvery == 0)) {
            Real totalEnergy = this->getTotalEnergy();
            energyFile << time << " " << totalEnergy << std::endl;
        }
        if ((saveTrajectory!="") && (time % saveEvery == 0)) {
            trajectoryFile << m_N << "\n";
            trajectoryFile << "Lattice=\"" << m_L << " 0.0 0.0 " << "0.0 " << m_L << " 0.0 " << "0.0 0.0 " << m_L << "\" "<< "Properties=species:S:1:pos:R:3 "<< "Time=" << time << "\n";
            for (int i = 0; i < m_N; ++i) {
                trajectoryFile << "H " << m_x[i] << " "  << m_y[i] << " " << m_z[i] << "\n"; 
            }
        }
    }
    energyFile.close();
    trajectoryFile.close();
}