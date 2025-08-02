#include "pbf.h"

#include <iostream>

#define PI 3.14159265358979323846

namespace pbf
{


	void PositionBasedFluids::initialize(
		std::vector<Eigen::Vector3d> &fluid_positions, 
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
		int num_iterations
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
		m_num_iterations = num_iterations;
		m_dt = dt;
		m_particle_diameter = particle_diameter;
		m_kernel_0 = pressureSolver::kernel::W(Eigen::Vector3d{ 0,0,0 });
		m_epsilon = 0.00001;
		m_fluid_positions_tmp = *m_fluid_positions;
		std::cout << "PBF initialized" << std::endl;
	}
	void PositionBasedFluids::step(pressureSolver::NeighborhoodSearch& nsearch)
	{
		m_fluid_positions_tmp.resize(m_fluid_positions->size());
		m_S_i.resize(m_fluid_positions->size());
		m_lambda_i.resize(m_fluid_positions->size());
		m_delta_x_i.resize(m_fluid_positions->size());
		//std::cout << "finished resize" << std::endl;
		//std::cout << "fluid particles" << m_fluid_positions->size() << std::endl;
		//std::cout << m_fluid_positions_tmp[100].transpose() << std::endl;
		//nsearch.run();
		
		compute_densities(nsearch);
		
		compute_accelerations(nsearch);
		
		advect_particles();

		nsearch.run();
		
		for (int iteration = 0; iteration < m_num_iterations; iteration++) {

			compute_densities_tmp(nsearch);
			//std::cout << "updated dens" << std::endl;
			//std::cout << m_densities->at(100) << std::endl;
			
			compute_S_i(nsearch);
			//std::cout << "updated S" << std::endl;
			//std::cout << m_S_i[100] << std::endl;

			compute_lambda_i(nsearch);
			//std::cout << "updated lambda" << std::endl;
			//std::cout << m_lambda_i[100] << std::endl;

			compute_delta_x_i(nsearch);
			//std::cout << "updated delta x" << std::endl;
			//std::cout << m_delta_x_i[100].transpose() << std::endl;
			//std::cout << m_fluid_positions_tmp[100].transpose() << std::endl;
			//
			//std::cout << "updated dens" << std::endl;
			//std::cout << m_densities->at(100) << std::endl;
			//std::cout << "updated S" << std::endl;
			//std::cout << m_S_i[100] << std::endl;
			//std::cout << "updated lambda" << std::endl;
			//std::cout << m_lambda_i[100] << std::endl;
			//std::cout << "updated delta x" << std::endl;
			//std::cout << m_delta_x_i[100].transpose() << std::endl;
		}

		update_velocities();

		
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
		}
		//std::cout << "Computing densities" << std::endl;
	}

	void PositionBasedFluids::compute_densities_tmp(pressureSolver::NeighborhoodSearch& nsearch)
	{

#pragma omp parallel for
		for (int i = 0; i < m_fluid_positions->size(); i++)
		{
			//bool debug = (i == 100);
			double tmp_kern_sum = m_kernel_0;
			nsearch.for_each_neighbor(0, 0, i,
				[&](const int n)
				{
					//if (debug) {
					//	std::cout << "n: " << n << ", max n: " << m_fluid_positions_tmp.size() << std::endl;
					//	std::cout << m_fluid_positions_tmp[i].transpose() << ", n: " << m_fluid_positions_tmp[n].transpose() << std::endl;
					//	std::cout << "distance: " << (m_fluid_positions_tmp[i] - m_fluid_positions_tmp[n]).norm() << " compact support: " << m_particle_diameter*1.2 << std::endl;
					//}
					tmp_kern_sum += pressureSolver::kernel::W(m_fluid_positions_tmp[i] - m_fluid_positions_tmp[n]);
				}
			);
			tmp_kern_sum *= m_fluid_mass;
			nsearch.for_each_neighbor(0, 1, i,
				[&](const int n)
				{
					tmp_kern_sum += m_boundary_masses->at(n) * pressureSolver::kernel::W(m_fluid_positions_tmp[i] - m_boundary_positions->at(n));
				}
			);
			m_densities->at(i) = tmp_kern_sum;
		}
		//std::cout << "Computing densities tmp" << std::endl;
	}

	void PositionBasedFluids::compute_accelerations(pressureSolver::NeighborhoodSearch& nsearch)
	{
		//compute accelerations: compute viscosity (boundary and fluid part like sph), add gravity
		double drag_factor = 1.225 * 0.47 * PI * m_particle_diameter / 2.0 * m_particle_diameter / 2.0;
#pragma omp parallel for
		for (int i = 0; i < m_fluid_positions->size(); i++)
		{
			if (m_is_scripted->at(i))
				continue;
			Eigen::Vector3d tmp_acc_viscosity = { 0,0,0 };
			nsearch.for_each_neighbor(0, 0, i,
				[&](const int n)
				{
					Eigen::Vector3d p_vec = m_fluid_positions->at(i) - m_fluid_positions->at(n);
					auto kernel_grad = pressureSolver::kernel::W_grad(p_vec);
					tmp_acc_viscosity += (m_fluid_mass / m_densities->at(n))
						* (m_velocities->at(i) - m_velocities->at(n))
						* (p_vec.transpose() * kernel_grad) / (p_vec.squaredNorm() + m_denom_eps);
				}
			);
			tmp_acc_viscosity = tmp_acc_viscosity * 2 * m_fluid_vis;
			// loop over boundary neighbourhood
			nsearch.for_each_neighbor(0, 1, i,
				[&](const int n)
				{
					auto p_vec = m_fluid_positions->at(i) - m_boundary_positions->at(n);
					auto kernel_grad = pressureSolver::kernel::W_grad(p_vec);
					tmp_acc_viscosity += 2 * m_boundary_vis * m_boundary_volumes->at(n) * m_velocities->at(i)
						* (p_vec.transpose() * kernel_grad) / (p_vec.squaredNorm() + m_denom_eps);
				}
			);
			double speed = m_velocities->at(i).norm();
			double drag = drag_factor * speed * speed;
			Eigen::Vector3d drag_acceleration = drag / m_fluid_mass * m_velocities->at(i);
			m_accelerations->at(i) += tmp_acc_viscosity - drag_acceleration;
		}
		//std::cout << "Computing accelerations" << std::endl;
	}

	//computes positions in m_fluid_positions_tmp
	void PositionBasedFluids::advect_particles()
	{
		
		//ensure cfl and clamp velocity
		double cfl_velocity = 0.4 * m_particle_diameter / m_dt;
		double squared_cfl_velocity = cfl_velocity * cfl_velocity;
#pragma omp parallel for
		for (int i = 0; i < m_fluid_positions->size(); i++)
		{
			if (m_is_scripted->at(i))
				continue;
			if(m_accelerations->at(i).hasNaN())
				std::cout << "nan acceleration" << std::endl;
			m_velocities->at(i) += m_accelerations->at(i) * m_dt;
			double squared_velocity = m_velocities->at(i).squaredNorm();
			if (squared_velocity > squared_cfl_velocity)
			{
				m_velocities->at(i) = m_velocities->at(i).normalized() * (cfl_velocity / (sqrt(squared_velocity)));
			}
			
			m_fluid_positions_tmp[i] = m_velocities->at(i) * m_dt + m_fluid_positions->at(i);
		}
		//std::cout << "Advecting particles" << std::endl;
	}

	void PositionBasedFluids::compute_S_i(pressureSolver::NeighborhoodSearch& nsearch)
	{
		double m_fluid_mass_inverse = 1 / m_fluid_mass;
		double m_rho_0_inverse = 1 / m_rho_0;
		double factor = m_fluid_mass / m_rho_0;

#pragma omp parallel for
		for (int i = 0; i < m_fluid_positions->size(); i++)
		{
			if (m_is_scripted->at(i))
				continue;
			//bool debug = (i == 100);
			if (m_densities->at(i) < m_rho_0)
				m_S_i[i] = 10;
			else {
				Eigen::Vector3d tmp_fluid_kernel_grad = { 0,0,0 };
				double part_2 = 0;
				Eigen::Vector3d tmp_boundary_kernel_grad = { 0,0,0 };
				Eigen::Vector3d tmp_boundary_kernel_grad_withV = { 0,0,0 };
				// loop over fluid neighborhood
				nsearch.for_each_neighbor(0, 0, i,
					[&](const int n)
					{
						//if(debug)
						//	std::cout << i << ", " << n << m_fluid_positions_tmp[i].transpose() << ", " << m_fluid_positions_tmp[n].transpose() << std::endl;
						Eigen::Vector3d grad = pressureSolver::kernel::W_grad(m_fluid_positions_tmp[i] - m_fluid_positions_tmp[n]);
						tmp_fluid_kernel_grad += factor * grad;
						part_2 += m_fluid_mass_inverse * (-1 * factor * grad).squaredNorm();
					}
				);
				// loop over boundary neighbourhood
				nsearch.for_each_neighbor(0, 1, i,
					[&](const int n)
					{
						Eigen::Vector3d gradient = pressureSolver::kernel::W_grad(m_fluid_positions_tmp[i] - m_boundary_positions->at(n));
						tmp_boundary_kernel_grad_withV += m_boundary_volumes->at(n) * gradient;
					}
				);
				m_S_i[i] = m_fluid_mass_inverse * (tmp_fluid_kernel_grad + tmp_boundary_kernel_grad_withV).squaredNorm() + part_2;
				//if (debug)
				//	std::cout << "S_i: " << m_S_i[i] << std::endl;
			}
		}
		//std::cout << "Computing S_i" << std::endl;
	}

	void PositionBasedFluids::compute_lambda_i(pressureSolver::NeighborhoodSearch& nsearch)
	{
		double m_rho_0_inverse = 1 / m_rho_0;
#pragma omp parallel for
		for (int i = 0; i < m_fluid_positions->size(); i++)
		{
			if (m_is_scripted->at(i))
				continue;
			double tmp_C_i = std::max(0.0, m_densities->at(i) * m_rho_0_inverse - 1.0);
			m_lambda_i[i] = -1 * tmp_C_i / (m_S_i[i] + m_denom_eps);
		}
		//std::cout << "Computing lambda_i" << std::endl;
	}
	//updates positions in m_fluid_positions_tmp
	void PositionBasedFluids::compute_delta_x_i(pressureSolver::NeighborhoodSearch& nsearch)
	{
		double rho_0_inverse = 1 / m_rho_0;
		double fluid_mass_inverse = 1 / m_fluid_mass;
#pragma omp parallel for
		for (int i = 0; i < m_fluid_positions->size(); i++)
		{
			if (m_is_scripted->at(i))
				continue;
			Eigen::Vector3d tmp_fluid_contribution = { 0,0,0 };
			Eigen::Vector3d tmp_boundary_contribution = { 0,0,0 };
			auto lambda_i = m_lambda_i[i];
			// loop over fluid neighborhood
			nsearch.for_each_neighbor(0, 0, i,
				[&](const int j)
				{
					tmp_fluid_contribution += (lambda_i + m_lambda_i[j])
						* pressureSolver::kernel::W_grad(m_fluid_positions_tmp[i] - m_fluid_positions_tmp[j]);
				}
			);
			//tmp_fluid_contribution *= rho_0_inverse;
			// loop over neighbourhood
			nsearch.for_each_neighbor(0, 1, i,
				[&](const int k)
				{
					tmp_boundary_contribution += m_boundary_masses->at(k) * lambda_i * fluid_mass_inverse
						* pressureSolver::kernel::W_grad(m_fluid_positions_tmp[i] - m_boundary_positions->at(k));
				}
			);
			//if(m_fluid_positions->size()<1500)
			//	std::cout << i << ", " << tmp_fluid_contribution.transpose() << ", " << m_fluid_positions->size() << std::endl;
			//tmp_boundary_contribution *= (rho_0_inverse * fluid_mass_inverse);
			m_delta_x_i[i] = rho_0_inverse * (tmp_fluid_contribution + tmp_boundary_contribution);
			m_fluid_positions_tmp[i] += m_delta_x_i[i];
		}

		//std::cout << "Computing delta x_i" << std::endl;
	}


	void PositionBasedFluids::update_velocities()
	{

		//std::cout << "old Fluid pos" << std::endl;
		//std::cout << m_fluid_positions->at(100).transpose() << std::endl;
#pragma omp parallel for
		for (int i = 0; i < m_fluid_positions->size(); i++)
		{
			if (m_is_scripted->at(i)) {
				//m_velocities->at(i) += m_accelerations->at(i) * m_dt;
				//m_velocities->at(i) = Eigen::Vector3d(0, 0, -9.81) * m_dt;
				m_fluid_positions->at(i) += m_velocities->at(i) * m_dt;
			}
			else {
				m_velocities->at(i) = (m_fluid_positions_tmp[i] - m_fluid_positions->at(i)) / m_dt;
				//double squared_velocity = m_velocities->at(i).squaredNorm();
				//max_vel = std::max(m_velocities->at(i).norm(), max_vel);
				// std::cout<<m_velocities->at(i)<<std::endl;
				/*if (squared_velocity > squared_cfl_velocity)
					m_velocities->at(i) = m_velocities->at(i) * (cfl_velocity / (sqrt(squared_velocity) * 1.1));*/
				m_fluid_positions->at(i) = m_fluid_positions_tmp[i];
			}


		}



		//std::cout << "Updating velocities" << std::endl;
	}

	//copy m_fluid_positions_tmp to m_fluid_positions
	void PositionBasedFluids::update_position()
	{
		//std::cout << "Updating position" << std::endl;
	}

}
