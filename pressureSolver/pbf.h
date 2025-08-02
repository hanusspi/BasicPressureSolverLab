#pragma once
#include "NeighborhoodSearch.h"
#include "kernel.h"
#include <vector>
#include <Eigen/Dense>

namespace pbf
{
	class PositionBasedFluids
	{
	public:
		PositionBasedFluids() = default;
		~PositionBasedFluids() = default;
		void initialize(std::vector<Eigen::Vector3d>& fluid_positions,
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
			int num_iterations);
		void step(pressureSolver::NeighborhoodSearch& nsearch);
	private:
		void compute_densities(pressureSolver::NeighborhoodSearch& nserach);
		void compute_densities_tmp(pressureSolver::NeighborhoodSearch& nserach);
		void compute_accelerations(pressureSolver::NeighborhoodSearch& nsearch);
		void advect_particles();
		void compute_S_i(pressureSolver::NeighborhoodSearch& nsearch);
		void compute_lambda_i(pressureSolver::NeighborhoodSearch& nsearch);
		void compute_delta_x_i(pressureSolver::NeighborhoodSearch& nsearch);
		void update_position();
		void update_velocities();
		
		int m_num_iterations;
		double m_dt;
		double m_rho_0 = 1000.0;
		double m_epsilon;
		double m_k;
		double m_fluid_mass;
		double m_kernel_0;
		double m_particle_diameter;
		double m_denom_eps;
		double m_smoothing_length;
		const double m_fluid_vis = 0.0025;
		const double m_boundary_vis = 0.0001;
		Eigen::Vector3d m_ext_forces = Eigen::Vector3d(0.0, 0.0, -9.81);
		std::vector<Eigen::Vector3d>* m_fluid_positions;
		std::vector<Eigen::Vector3d> m_fluid_positions_tmp;
		std::vector<Eigen::Vector3d>* m_boundary_positions;
		std::vector<Eigen::Vector3d>* m_velocities;
		std::vector<Eigen::Vector3d>* m_accelerations;
		std::vector<double>* m_boundary_volumes;
		std::vector<double>* m_boundary_masses;
		std::vector<double>* m_densities;
		std::vector<double> m_S_i;
		std::vector<double> m_lambda_i;
		std::vector<Eigen::Vector3d> m_delta_x_i;
		std::vector<Eigen::Vector3d>* m_external_forces;
		std::vector<bool>* m_is_scripted;
	};
}
