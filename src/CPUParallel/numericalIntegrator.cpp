#include <iostream>
#include <fstream>

#include "simulation.hpp"

void Simulation::integrateEuler(int n_steps, int saveEvery, std::string saveEnergy, std::string saveTrajectory) {
    // Output files (energy, trajectories)
    std::ofstream energyFile; energyFile.open(saveEnergy);
    std::ofstream trajectoryFile; trajectoryFile.open(saveTrajectory);
    trajectoryFile.precision(5);
    trajectoryFile << std::scientific;

    int time = 0;
    double kineticEnergy = 0.0;
    double potentialEnergy = 0.0;

    // Parallel region. All threads will execute in parallel every single time-step (a final barrier ensures synchronization in time)
    #pragma omp parallel 
    {
        while (time < n_steps) {
            // Determine if we need to save this step
            bool isSavingStep = (saveEnergy != "") && ((time + 1) % saveEvery == 0);
            // Compute forces (O(N^2)) and potential energy in a parallel fashion. A single thread will compute the forces for a portion of particles
            // (updating a portion of the vectors m_fx, m_fy, m_fz and compute a partial potential energy; then a reduction will sum all partial energies.
            double partialPotentialEnergy = this->computeForces(isSavingStep); 
            // An implicit barrier is here, to ensure all threads have computed the forces before proceeding
            // Now we can update positions and velocities in parallel (again, static scheduling, uniform load)
            #pragma omp for schedule(static)
            for (int i = 0; i < m_N; ++i) {
                m_x[i] += m_vx[i] * m_dt;
                m_y[i] += m_vy[i] * m_dt;
                m_z[i] += m_vz[i] * m_dt;

                Real inv_m = 1.0 / m_mass[i];
                m_vx[i] += m_fx[i] * inv_m * m_dt;
                m_vy[i] += m_fy[i] * inv_m * m_dt;
                m_vz[i] += m_fz[i] * inv_m * m_dt;

                if (m_x[i] < 0) m_x[i] += m_L; else if (m_x[i] >= m_L) m_x[i] -= m_L;
                if (m_y[i] < 0) m_y[i] += m_L; else if (m_y[i] >= m_L) m_y[i] -= m_L;
                if (m_z[i] < 0) m_z[i] += m_L; else if (m_z[i] >= m_L) m_z[i] -= m_L;
            }
            // Here again an implicit barrier

            // If it is a saving step, compute kinetic energy and sum the partial potential energies
            if (isSavingStep) {
                #pragma omp for schedule(static) reduction(+:kineticEnergy)
                for (int i = 0; i < m_N; ++i) {
                    Real v2 = m_vx[i]*m_vx[i] + m_vy[i]*m_vy[i] + m_vz[i]*m_vz[i];
                    kineticEnergy += 0.5 * m_mass[i] * v2;
                }
                // An atomic update should be fine here
                #pragma omp atomic
                potentialEnergy += partialPotentialEnergy;
            }
            // Ensure all threads have finished computing energies before the master writes to file. If not savingStep, just a barrier to sync threads
            #pragma omp barrier

            // And now the master thread will reset the energies and save to file if needed
            #pragma omp master
            {
                time++;
                if (isSavingStep) {
                    energyFile << (time + 1) * m_dt << " " << kineticEnergy << " " << potentialEnergy << " " << (kineticEnergy + potentialEnergy) << "\n";
                    trajectoryFile << m_N << "\n";
                    trajectoryFile << "Lattice=\"" << m_L << " 0.0 0.0 0.0 " << m_L << " 0.0 0.0 0.0 " << m_L << "\" " << "Properties=species:S:1:pos:R:3 Time=" << (time + 1) << "\n";
                    for (int i = 0; i < m_N; ++i) {
                        trajectoryFile << "H " << m_x[i] << " " << m_y[i] << " " << m_z[i] << "\n"; 
                    }
                }
                kineticEnergy = 0.0;
                potentialEnergy = 0.0;
            }
            #pragma omp barrier
        }
    }
}


void Simulation::integrateVerlet(int n_steps, int saveEvery, std::string saveEnergy, std::string saveTrajectory) {
    // Output files (energy, trajectories)
    std::ofstream energyFile; energyFile.open(saveEnergy);
    std::ofstream trajectoryFile; trajectoryFile.open(saveTrajectory);
    trajectoryFile.precision(5);
    trajectoryFile << std::scientific;

    int time = 0;
    double kineticEnergy = 0.0;
    double potentialEnergy = 0.0;

    // Parallel region. All threads will execute in parallel every single time-step (a final barrier ensures synchronization in time)
    #pragma omp parallel 
    {
        // Initial force computation before starting the time-stepping
        double partialPotentialEnergy = this->computeForces(false); // An implicit barrier is here, to ensure all threads have computed the forces before proceeding
        while (time < n_steps) {
            // Determine if we need to save this step
            bool isSavingStep = (saveEnergy != "") && ((time + 1) % saveEvery == 0);
            
            // Verlet integration first half step velocity + full step position. 
            // Each thread will take care of a portion of particles, updating their positions and velocities (verlet)
            #pragma omp for schedule(static)
            for (int i = 0; i < m_N; ++i) {
                Real inv_m = 1.0 / m_mass[i];
                m_vx[i] += 0.5 * m_fx[i] * inv_m * m_dt;
                m_vy[i] += 0.5 * m_fy[i] * inv_m * m_dt;
                m_vz[i] += 0.5 * m_fz[i] * inv_m * m_dt;

                m_x[i] += m_vx[i] * m_dt;
                m_y[i] += m_vy[i] * m_dt;
                m_z[i] += m_vz[i] * m_dt;

                if (m_x[i] < 0) m_x[i] += m_L; else if (m_x[i] >= m_L) m_x[i] -= m_L;
                if (m_y[i] < 0) m_y[i] += m_L; else if (m_y[i] >= m_L) m_y[i] -= m_L;
                if (m_z[i] < 0) m_z[i] += m_L; else if (m_z[i] >= m_L) m_z[i] -= m_L;
            } // Here a barrier, meaning that all threads have finished updating positions before computing forces

            // Now we need to compute forces (O(N^2)) and potential energy in a parallel fashion. A single thread will compute the forces for a portion of particles
            // (updating a portion of the vectors m_fx, m_fy, m_fz and compute a partial potential energy; then a reduction will sum all partial energies.
            double partialPotentialEnergy = this->computeForces(isSavingStep); // An implicit barrier is here, to ensure all threads have computed the forces before proceeding

            // Verlet integration second half step velocity update
            #pragma omp for schedule(static)
            for (int i = 0; i < m_N; ++i) {
                Real inv_m = 1.0 / m_mass[i];
                m_vx[i] += 0.5 * m_fx[i] * inv_m * m_dt;
                m_vy[i] += 0.5 * m_fy[i] * inv_m * m_dt;
                m_vz[i] += 0.5 * m_fz[i] * inv_m * m_dt;
            } // Here again an implicit barrier

            // If it is a saving step, compute kinetic energy and sum the partial potential energies
            if (isSavingStep) {
                #pragma omp for schedule(static) reduction(+:kineticEnergy)
                for (int i = 0; i < m_N; ++i) {
                    Real v2 = m_vx[i]*m_vx[i] + m_vy[i]*m_vy[i] + m_vz[i]*m_vz[i];
                    kineticEnergy += 0.5 * m_mass[i] * v2;
                }
                // An atomic update should be fine here
                #pragma omp atomic
                potentialEnergy += partialPotentialEnergy;
            }
            // Ensure all threads have finished computing energies before the master writes to file. If not savingStep, just a barrier to sync threads
            #pragma omp barrier

            // And now the master thread will reset the energies and save to file if needed
            #pragma omp master
            {
                time++;
                if (isSavingStep) {
                    energyFile << (time + 1) * m_dt << " " << kineticEnergy << " " << potentialEnergy << " " << (kineticEnergy + potentialEnergy) << "\n";
                    trajectoryFile << m_N << "\n";
                    trajectoryFile << "Lattice=\"" << m_L << " 0.0 0.0 0.0 " << m_L << " 0.0 0.0 0.0 " << m_L << "\" " << "Properties=species:S:1:pos:R:3 Time=" << (time + 1) << "\n";
                    for (int i = 0; i < m_N; ++i) {
                        trajectoryFile << "H " << m_x[i] << " " << m_y[i] << " " << m_z[i] << "\n"; 
                    }
                }
                kineticEnergy = 0.0;
                potentialEnergy = 0.0;
            }
            #pragma omp barrier
        }
    }
}