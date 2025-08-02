#include "wcsph.h"

#include <iostream>

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>

#include "../pressureSolver/io.h"
#include "../pressureSolver/kernel.h"
#include "../pressureSolver/sampling.h"

void showProgressBar(int progress, int total, int barWidth = 50) {
	float percentage = (float)progress / total;

	std::cout << "\r[";

	int pos = barWidth * percentage;
	for (int i = 0; i < barWidth; ++i) {
		if (i < pos) std::cout << "=";
		else if (i == pos) std::cout << ">";
		else std::cout << " ";
	}

	// Show percentage completed
	std::cout << "] " << int(percentage * 100.0) << " %";

	// Flush the output to update the display
	std::cout.flush();
}

namespace pressureSolver {
	namespace wcsph {

		Simulation::Simulation(std::vector<Eigen::Vector3d> &positions, const double volume, const double freq, const double fps, const double rho0, const double particle_radius) :
			m_positions(positions),
			m_freq(freq),
			m_fps(fps),
			m_t(0),
			m_t_frame(0),
			m_t_frame_step(1.0 / fps),
			m_delta_t(1.0 / freq),
			m_fluid_particle_mass(rho0 * volume/positions.size()),
			m_rho0(rho0),
			m_particle_radius(particle_radius),
			m_particle_diameter(2.0 * particle_radius),
			m_fluid_sampling_distance(m_particle_diameter),
			m_boundary_sampling_distance(0.8 * m_particle_diameter),
			m_smoothing_length(1.2 * m_particle_diameter),
			m_compact_support(2.0 * m_smoothing_length),
			m_denom_eps(0.01 * m_smoothing_length * m_smoothing_length),
			m_nsearch(m_compact_support) {
			//initialize Kernel Parameters
			kernel::setAlpha(m_smoothing_length);
			m_kernel_0 = kernel::W(Eigen::Vector3d{ 0,0,0 });
			//resize particle data
			boundary = true;
			//intialize densities of particles		
			sampleFluidAndBoundaries();
			// computeDensities();

			std::string filename = "../res/boundary_volume" + std::to_string(0) + ".vtk";
			pressureSolver::write_particles_to_vtk(filename, m_boundary_positions, m_boundary_volumes);

			std::cout << "Simulation created" << std::endl;
		}

		Simulation::~Simulation() {
			std::cout << "Simulation destroyed" << std::endl;
		}

		void Simulation::run() {
			std::cout << "Simulation running" << std::endl;
			int steps = 30000;
			for (int i = 0; i < steps; i++) {
				m_nsearch.find_neighbors(); // find neighbors
				computeDensities(); 
				computeAccelerations();
				integrateParticles();
				m_t += m_delta_t;
				if (m_t >= m_t_frame) {
					showProgressBar(i, steps);
					removeOutsideParticles();
					exportParticlesToVTK();
					m_t_frame += m_t_frame_step;
					m_frame += 1;
				}				
			}
		}


		void Simulation::sampleFluidAndBoundaries() {
			std::cout << "Loading mesh" << std::endl;
			const std::vector<pressureSolver::TriMesh> meshes = pressureSolver::read_tri_meshes_from_obj("../res/boxEx2_6b.obj"); //needs to be parametrized into constructor
			const pressureSolver::TriMesh& box = meshes[0];

			std::cout << "Calculating boundary box" << std::endl;
			for (const auto& v : box.vertices) {
				m_viewbox_min = m_viewbox_min.cwiseMin(v);
				m_viewbox_max = m_viewbox_max.cwiseMax(v);
			}

			std::cout << "Sampling fluid and boundary particles" << std::endl;
			sampling::triangle_mesh(m_boundary_positions, box.vertices, box.triangles, m_particle_radius);
			m_nsearch.add_point_set(m_positions.front().data(), m_positions.size(), true, true);
			m_nsearch.add_point_set(m_boundary_positions.front().data(), m_boundary_positions.size(), false, true);

			m_velocities.resize(m_positions.size(), Eigen::Vector3d{ 0.0,0.0,0.0 });

			m_densities.resize(m_positions.size(), 0.0);
			m_accelerations.resize(m_positions.size(), Eigen::Vector3d{ 0.0,0.0,0.0 });
			m_boundary_masses.resize(m_boundary_positions.size(), 0.0);
			m_boundary_volumes.resize(m_boundary_positions.size(), 0.0);
			m_pressures.resize(m_positions.size(), 0.0);

			m_nsearch.find_neighbors();

			auto const& d = m_nsearch.point_set(BOUNDARY_NEIGHBORHOOD);
			
			#pragma omp parallel for
			for (int i = 0; i < d.n_points(); i++) {// iterate over all particles stores in the nsearch point set
				auto t = d.n_neighbors(BOUNDARY_NEIGHBORHOOD, i); // get number of neighbors of particle i
				double tmp_kern_sum = m_kernel_0;
				for (int j = 0; j < t; j++) { //iterate over all neighbors of particle i
					auto n = d.neighbor(BOUNDARY_NEIGHBORHOOD, i, j);
					tmp_kern_sum += kernel::W(m_boundary_positions[i] - m_boundary_positions[n]);
				}
				m_boundary_volumes[i] = 1.0 / tmp_kern_sum;
				m_boundary_masses[i] = m_boundary_volumes[i] * m_rho0;
			}
			std::cout << "Finished boundary correction" << std::endl;
		}

		void Simulation::computeDensities() {
			auto const& d = m_nsearch.point_set(FLUID_NEIGHBORHOOD);
			#pragma omp parallel for
			for (int i = 0; i < d.n_points(); i++) {// iterate over all particles stores in the nsearch point set
				auto f = d.n_neighbors(FLUID_NEIGHBORHOOD, i); // get number of neighbors of particle i
				auto b = d.n_neighbors(BOUNDARY_NEIGHBORHOOD, i);
				double tmp = m_kernel_0; //particle is its own neighbor
				for (int j = 0; j < f; j++) { //iterate over all fluid neighbors of particle i
					auto n = d.neighbor(FLUID_NEIGHBORHOOD, i, j);
					tmp += kernel::W(m_positions[i] - m_positions[n]);
				}
				tmp *= m_fluid_particle_mass; //weight sum of kernel values by mass
				if (boundary) {
					for (int k = 0; k < b; k++) { //iterate over all boundary neighbors of particle i
						auto n = d.neighbor(BOUNDARY_NEIGHBORHOOD, i, k);
						tmp += kernel::W(m_positions[i] - m_boundary_positions[n]) * m_boundary_masses[n];
					}
				}
				// #pragma omp critical
				// as only one thread can access the densities vector at a time, remove the critical section
				m_densities[i] = tmp;
				m_pressures[i] = std::max(0.0, m_B * (m_densities[i] - m_rho0*1.02));
			}
			//std::cout << "Computing densities" << std::endl;
		}

		void Simulation::computeAccelerations() {
			auto const& d = m_nsearch.point_set(FLUID_NEIGHBORHOOD);
			#pragma omp parallel for
			for (int i = 0; i < d.n_points(); i++) {// iterate over all particles stores in the nsearch point set
				auto f = d.n_neighbors(FLUID_NEIGHBORHOOD, i);
				// get number of neighbors of particle i
				auto b = d.n_neighbors(BOUNDARY_NEIGHBORHOOD, i);
				Eigen::Vector3d tmp_acc_fluid = { 0,0,0 };
				Eigen::Vector3d tmp_acc_viscosity = { 0,0,0 };
				double p_d = m_pressures[i] / (m_densities[i] * m_densities[i]);
				for (int j = 0; j < f; j++) { //iterate over all fluid neighbors of particle i
					//fluid part of fluid acceleration
					auto n = d.neighbor(FLUID_NEIGHBORHOOD, i, j);
					auto p_n = m_pressures[n] / (m_densities[n] * m_densities[n]);
					auto p_vec = m_positions[i] - m_positions[n];
					auto kernel_grad = kernel::W_grad(p_vec);
					tmp_acc_fluid -= m_fluid_particle_mass * (p_d + p_n) * kernel_grad;
					//fluid part of viscosity acceleration
					// auto norm = (m_positions[i] - m_positions[n]).norm(); // avoid use of norm * norm, use squaredNorm instead
					tmp_acc_viscosity += (m_fluid_particle_mass / m_densities[n]) * (m_velocities[i] - m_velocities[n]) * (p_vec.transpose() * kernel_grad) / (p_vec.squaredNorm() + m_denom_eps);
				}
				tmp_acc_viscosity = tmp_acc_viscosity * 2 * m_fluid_vis;
				if(boundary){
					for (int k = 0; k < b; k++) { //iterate over all boundary neighbors of particle i
						//boundary part of fluid acceleration
						auto n = d.neighbor(BOUNDARY_NEIGHBORHOOD, i, k);
						auto p_vec = m_positions[i] - m_boundary_positions[n];
						auto kernel_grad = kernel::W_grad(p_vec);
						tmp_acc_fluid -= m_boundary_masses[k] * p_d * kernel_grad;
						//boundary part of viscosity acceleration
						tmp_acc_viscosity += 2 * m_boundary_vis * m_boundary_volumes[k] * m_velocities[i] * (p_vec.transpose() * kernel_grad) / (p_vec.squaredNorm() + m_denom_eps);
					}
				}
				//std::cout << "tmp_acc_fluid: " << tmp_acc_fluid << std::endl;
				// #pragma omp critical
				// as only one thread can access the densities vector at a time, remove the critical section
				m_accelerations[i] = tmp_acc_fluid + tmp_acc_viscosity + m_gravity;

			}
			//std::cout << "Computing accelerations" << std::endl;
		}

		void Simulation::integrateParticles() {
			bool cfl_condition = true;
			//#pragma omp parallel for
			//collect list of items to throw away to enable paralellism
			double cfl_velocity = 0.5 * m_particle_diameter / m_delta_t;
			double squared_cfl_velocity = cfl_velocity * cfl_velocity;
			std::vector<bool> to_erase(m_positions.size(), false);

			#pragma omp parallel for
			for (int i = m_positions.size() - 1; i >=0; i--)
			{
				//#pragma omp critical
				m_velocities[i] = m_velocities[i] + m_delta_t * m_accelerations[i];
				double squared_velocity = m_velocities[i].squaredNorm();

				if (squared_velocity > squared_cfl_velocity)
				{
					//std::cout << "CFL condition not met" << std::endl;
					m_velocities[i] = m_velocities[i] * (cfl_velocity / (sqrt(squared_velocity)*1.1));// evtl define cfl condition and clamp max velocity
				}
				m_positions[i] = m_positions[i] + m_delta_t * m_velocities[i];
			}
		}

		void Simulation::exportParticlesToVTK() 
		{
			pressureSolver::write_particles_to_vtk("../res/fluid_new_velo" + std::to_string(m_frame) + ".vtk", m_positions, m_densities);
			//std::cout << "Exporting particles to VTK" << std::endl;
		}

		void Simulation::removeOutsideParticles() {
			std::vector<bool> to_erase(m_positions.size(), false); // local variable to store which particles to erase

			#pragma omp parallel for
			for (int i = m_positions.size() - 1; i >= 0; i--)
			{
				if (m_positions[i].x() < m_viewbox_min.x() || m_positions[i].x() > m_viewbox_max.x() ||
					m_positions[i].y() < m_viewbox_min.y() || m_positions[i].y() > m_viewbox_max.y() ||
					m_positions[i].z() < m_viewbox_min.z() || m_positions[i].z() > m_viewbox_max.z())
				{
					to_erase[i] = true;
				}
			}

			for (int i = m_positions.size() - 1; i >= 0; i--)
			{
				if (to_erase[i])
				{
					m_positions.erase(m_positions.begin() + i);
					m_velocities.erase(m_velocities.begin() + i);
					m_densities.erase(m_densities.begin() + i);
					m_accelerations.erase(m_accelerations.begin() + i);
					m_pressures.erase(m_pressures.begin() + i);
				}
			}

			m_nsearch.resize_point_set(FLUID_NEIGHBORHOOD, m_positions.front().data(), m_positions.size());
		}

	}
}

