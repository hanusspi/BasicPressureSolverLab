#include "wcsphSolver.h"
#include "util.h"

#include <iostream>

#define PI 3.14159265358979323846

namespace wcsph
{

	void PositionBasedFluids::initialize(
		std::vector<Eigen::Vector3d>& fluid_positions,
		std::vector<Eigen::Vector3d>& boundary_positions,
		std::vector<Eigen::Vector3d>& velocities,
		std::vector<Eigen::Vector3d>& accelerations,
		std::vector<double>& densities,
		std::vector<double>& boundary_masses,
		std::vector<double>& boundary_volumes,
		std::vector<Eigen::Vector3d>& external_forces,
		std::vector<bool>& is_scripted,
		double fluid_mass,
		double particle_diameter,
		double dt,
		double stiffness
	)
	{
		m_fluid_positions = &fluid_positions;
		m_boundary_positions = &boundary_positions;
		m_velocities = &velocities;
		m_accelerations = &accelerations;
		m_densities = &densities;
		m_boundary_masses = &boundary_masses;
		m_boundary_volumes = &boundary_volumes;
		m_external_forces = &external_forces;
		m_is_scripted = &is_scripted;
		m_fluid_mass = fluid_mass;
		m_stiffness = stiffness;
		m_B = stiffness;
		m_particle_diameter = particle_diameter;
		m_kernel_0 = pressureSolver::kernel::W(Eigen::Vector3d{ 0,0,0 });
		m_epsilon = 0.0000001;
		m_dt = dt;
		m_rho_0 = 1000.0;
		std::cout << "wcsph initialized" << std::endl;
	}
	void PositionBasedFluids::step(pressureSolver::NeighborhoodSearch& nsearch)
	{
		m_pressures.resize(m_fluid_positions->size(),0.0);
		//nsearch.run();
		compute_densities(nsearch);
		
		compute_accelerations(nsearch);
		//std::cout << m_accelerations->at(0).transpose() << std::endl;
		advect_particles();

		//std::cout << "end wcsph step" << std::endl;
	}

	void PositionBasedFluids::compute_densities(pressureSolver::NeighborhoodSearch& nserach)
	{
#pragma omp parallel for
		for (int i = 0; i < m_fluid_positions->size(); i++)
		{
			
			double tmp_kern_sum = m_kernel_0;
			nserach.for_each_neighbor(0, 0, i,
				[&](const int n)
				{
					tmp_kern_sum += pressureSolver::kernel::W(m_fluid_positions->at(i) - m_fluid_positions->at(n));
				}
			);
			tmp_kern_sum *= m_fluid_mass;
			nserach.for_each_neighbor(0, 1, i,
				[&](const int n)
				{
					tmp_kern_sum += m_boundary_masses->at(n) * pressureSolver::kernel::W(m_fluid_positions->at(i) - m_boundary_positions->at(n));
				}
			);
			m_densities->at(i) = tmp_kern_sum;
			m_pressures[i] = std::max(0.0, m_B * (m_densities->at(i) - m_rho_0 * 1.02));
		}
	}

	void PositionBasedFluids::compute_accelerations(pressureSolver::NeighborhoodSearch& nserach)
	{
		double drag_factor = 1.225 * 0.47 * PI * m_particle_diameter / 2.0 * m_particle_diameter / 2.0;
#pragma omp parallel for
		for (int i = 0; i < m_fluid_positions->size(); i++)
		{
			if (m_is_scripted->at(i))
				continue;
			Eigen::Vector3d tmp_acc_fluid = { 0,0,0 };
			Eigen::Vector3d tmp_acc_viscosity = { 0,0,0 };
			double p_d = m_pressures[i] / (m_densities->at(i) * m_densities->at(i));
			nserach.for_each_neighbor(0, 0, i,
				[&](const int n)
				{
					auto p_n = m_pressures[n] / (m_densities->at(n) * m_densities->at(n));
					auto p_vec = m_fluid_positions->at(i) - m_fluid_positions->at(n);
					auto kernel_grad = pressureSolver::kernel::W_grad(p_vec);
					tmp_acc_fluid -= m_fluid_mass * (p_d + p_n) * kernel_grad;
					//fluid part of viscosity acceleration
					// auto norm = (m_positions[i] - m_positions[n]).norm(); // avoid use of norm * norm, use squaredNorm instead
					tmp_acc_viscosity += (m_fluid_mass / m_densities->at(n)) * (m_velocities->at(i) - m_velocities->at(n)) * (p_vec.transpose() * kernel_grad) / (p_vec.squaredNorm() + m_epsilon);
				}
			);
			tmp_acc_viscosity = tmp_acc_viscosity * 2 * m_fluid_vis;
			nserach.for_each_neighbor(0, 1, i,
				[&](const int n)
				{
					auto p_vec = m_fluid_positions->at(i) - m_boundary_positions->at(n);
					auto kernel_grad = pressureSolver::kernel::W_grad(p_vec);
					tmp_acc_fluid -= m_boundary_masses->at(n) * p_d * kernel_grad;
					//boundary part of viscosity acceleration
					tmp_acc_viscosity += 2 * m_boundary_vis * m_boundary_volumes->at(n) * m_velocities->at(i) * (p_vec.transpose() * kernel_grad) / (p_vec.squaredNorm() + m_epsilon);
				}
			);
			double speed = m_velocities->at(i).norm();
			double drag = drag_factor * speed * speed;
			Eigen::Vector3d drag_acceleration = drag / m_fluid_mass * m_velocities->at(i);
		m_accelerations->at(i) = tmp_acc_fluid + tmp_acc_viscosity + m_accelerations->at(i);
		
		}
		//std::cout << "Computing accelerations" << std::endl;
	}

	//computes positions in m_fluid_positions_tmp
	void PositionBasedFluids::advect_particles()
	{
		double cfl_velocity = 0.4 * m_particle_diameter / m_dt;
		double squared_cfl_velocity = cfl_velocity * cfl_velocity;
#pragma omp parallel for
		for (int i = m_fluid_positions->size() - 1; i >= 0; i--)
		{
			//#pragma omp critical
			if (m_is_scripted->at(i))
			{
				m_fluid_positions->at(i) = m_fluid_positions->at(i) + m_dt * m_velocities->at(i);
				continue;
			}
			m_velocities->at(i) = m_velocities->at(i) + m_dt * m_accelerations->at(i);
			double squared_velocity = m_velocities->at(i).squaredNorm();

			if (squared_velocity > squared_cfl_velocity)
			{
				//std::cout << "CFL condition not met" << std::endl;
				m_velocities->at(i) = m_velocities->at(i) * (cfl_velocity / (sqrt(squared_velocity)));// evtl define cfl condition and clamp max velocity
			}
			m_fluid_positions->at(i) = m_fluid_positions->at(i) + m_dt * m_velocities->at(i);
		}
		//std::cout << "Advecting particles" << std::endl;
	}

}
