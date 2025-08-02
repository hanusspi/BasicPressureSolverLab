#pragma once
#include "NeighborhoodSearch.h"
#include "kernel.h"
#include <vector>
#include <Eigen/Dense>

namespace wcsph
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
			std::vector<bool>& is_scritped,			
			double fluid_mass,
			double particle_diameter,
			double dt,
			double stiffness);
		void step(pressureSolver::NeighborhoodSearch& nsearch);
	private:
		void compute_densities(pressureSolver::NeighborhoodSearch& nserach);
		void compute_accelerations(pressureSolver::NeighborhoodSearch& nserach);
		void advect_particles();

		double m_particle_diameter;
		double m_stiffness;
		double m_dt;
		double m_rho_0;
		double m_epsilon;
		double m_k;
		double m_fluid_mass;
		double m_kernel_0;
		double m_fluid_vis = 0.0010;
		double m_boundary_vis = 0.0001;
		double m_B = 1000.0;
		Eigen::Vector3d m_gravity = { 0.0,0.0,-9.81 }; // still need to make parametrized
		std::vector<Eigen::Vector3d> *m_fluid_positions;
		std::vector<Eigen::Vector3d> *m_boundary_positions;
		std::vector<Eigen::Vector3d> *m_velocities;
		std::vector<Eigen::Vector3d> *m_accelerations;
		std::vector<double> *m_boundary_volumes;
		std::vector<double> *m_boundary_masses;
		std::vector<Eigen::Vector3d>* m_external_forces;
		std::vector<double> *m_densities;
		std::vector<bool>* m_is_scripted;

		std::vector<double> m_pressures;
	};
}
