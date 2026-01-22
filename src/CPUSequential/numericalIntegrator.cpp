#include <fstream>
#include "simulation.hpp"

// EULER INTEGRATOR
void Simulation::integrateEuler(int n_steps, int saveEvery, std::string saveEnergy, std::string saveTrajectory) {
    // Output files (energy, trajectories)
    std::ofstream energyFile; energyFile.open(saveEnergy);
    std::ofstream trajectoryFile; trajectoryFile.open(saveTrajectory);
    trajectoryFile.precision(5);
    trajectoryFile << std::scientific;
    energyFile.precision(7);  
    energyFile << std::scientific;

    int time = 0;
    // Main loop
    while (time < n_steps) {
        // Decides whether to save energy and trajectory this step
        bool isSavingStep = (saveEnergy != "") && ((time + 1) % saveEvery == 0);
        // First of all, compute forces (O(N^2)) and potential energy all in one (when needed, otherwise U = 0)
        double potentialEnergy = this->computeForces(isSavingStep); 

        // Now, update positions and velocities using Euler's method (O(N) complexity)
        for (int i = 0; i < m_N; ++i) {
            // Update positions
            m_x[i] += m_vx[i] * m_dt;
            m_y[i] += m_vy[i] * m_dt;
            m_z[i] += m_vz[i] * m_dt;

            // Update velocities
            m_vx[i] += (m_fx[i] / m_mass[i]) * m_dt;
            m_vy[i] += (m_fy[i] / m_mass[i]) * m_dt;
            m_vz[i] += (m_fz[i] / m_mass[i]) * m_dt;

            // Apply PBC
            if (m_x[i] < 0) m_x[i] += m_L;
            else if (m_x[i] >= m_L) m_x[i] -= m_L;
            if (m_y[i] < 0) m_y[i] += m_L;
            else if (m_y[i] >= m_L) m_y[i] -= m_L;
            if (m_z[i] < 0) m_z[i] += m_L;
            else if (m_z[i] >= m_L) m_z[i] -= m_L;
        }

        // Save kinetic energy and potential energy when needed
        if (isSavingStep) {
            // Potential energy is already computed (so we don't loop again O(N^2))
            Real kineticEnergy = 0.0f;
            for (int i = 0; i < m_N; ++i) {
                Real v2 = m_vx[i]*m_vx[i] + m_vy[i]*m_vy[i] + m_vz[i]*m_vz[i];
                kineticEnergy += 0.5f * m_mass[i] * v2;
            }
            energyFile << time * m_dt << " " << kineticEnergy << " " << potentialEnergy << " " << (kineticEnergy + potentialEnergy) << "\n";
        }
        if (isSavingStep) {
            trajectoryFile << m_N << "\n";
            trajectoryFile << "Lattice=\"" << m_L << " 0.0 0.0 " << "0.0 " << m_L << " 0.0 " << "0.0 0.0 " << m_L << "\" "<< "Properties=species:S:1:pos:R:3 "<< "Time=" << time << "\n";
            for (int i = 0; i < m_N; ++i) {
                trajectoryFile << "H " << m_x[i] << " "  << m_y[i] << " " << m_z[i] << "\n"; 
            }
        }
        // Update time
        time++;
    }
    energyFile.close();
    trajectoryFile.close();
}


// VERLET INTEGRATOR, WIP
void Simulation::integrateVerlet(int n_steps, int saveEvery, std::string saveEnergy, std::string saveTrajectory) {
    // Output files (energy, trajectories)
    std::ofstream energyFile; energyFile.open(saveEnergy);
    std::ofstream trajectoryFile; trajectoryFile.open(saveTrajectory);
    trajectoryFile.precision(5);
    trajectoryFile << std::scientific;
    energyFile.precision(7);  
    energyFile << std::scientific;

    int time = 0;
    this->computeForces(false); 
    // Main loop
    while (time < n_steps) {
        // Decides whether to save energy and trajectory AT THE END of this step
        bool isSavingStep = (saveEnergy != "") && ((time + 1) % saveEvery == 0);

        Real half_dt = 0.5 * m_dt;

        // first half-kick
        for (int i = 0; i < m_N; ++i) {
            Real inv_m = 1.0 / m_mass[i];
            // Update velocities by half dt
            m_vx[i] += (m_fx[i] * inv_m) * half_dt;
            m_vy[i] += (m_fy[i] * inv_m) * half_dt;
            m_vz[i] += (m_fz[i] * inv_m) * half_dt;
            // Update positions by full dt with the newly computed velocities
            m_x[i] += m_vx[i] * m_dt;
            m_y[i] += m_vy[i] * m_dt;
            m_z[i] += m_vz[i] * m_dt;
            // PBC
            if (m_x[i] < 0) m_x[i] += m_L; else if (m_x[i] >= m_L) m_x[i] -= m_L;
            if (m_y[i] < 0) m_y[i] += m_L; else if (m_y[i] >= m_L) m_y[i] -= m_L;
            if (m_z[i] < 0) m_z[i] += m_L; else if (m_z[i] >= m_L) m_z[i] -= m_L;
        }
        //Recompute forces (O(N^2))
        double potentialEnergy = this->computeForces(isSavingStep); 
        // Second half-kick
        for (int i = 0; i < m_N; ++i) {
            Real inv_m = 1.0 / m_mass[i];
            m_vx[i] += (m_fx[i] * inv_m) * half_dt;
            m_vy[i] += (m_fy[i] * inv_m) * half_dt;
            m_vz[i] += (m_fz[i] * inv_m) * half_dt;
        }

        // Update time
        time++;

        if (isSavingStep) {
            Real kineticEnergy = 0.0f;
            for (int i = 0; i < m_N; ++i) {
                Real v2 = m_vx[i]*m_vx[i] + m_vy[i]*m_vy[i] + m_vz[i]*m_vz[i];
                kineticEnergy += 0.5f * m_mass[i] * v2;
            }
            // Use scientific notation for precision
            energyFile << time * m_dt << " " << kineticEnergy << " " << potentialEnergy << " " << (kineticEnergy + potentialEnergy) << "\n";
        }

        if (isSavingStep && (saveTrajectory != "")) {
            trajectoryFile << m_N << "\n";
            trajectoryFile << "Lattice=\"" << m_L << " 0.0 0.0 " << "0.0 " << m_L << " 0.0 " << "0.0 0.0 " << m_L << "\" "
                           << "Properties=species:S:1:pos:R:3 " << "Time=" << time << "\n";
            for (int i = 0; i < m_N; ++i) {
                trajectoryFile << "H " << m_x[i] << " "  << m_y[i] << " " << m_z[i] << "\n"; 
            }
        }
    }
    energyFile.close();
    trajectoryFile.close();
}