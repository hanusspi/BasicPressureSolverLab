#pragma once
#include "pbf.h"
#include "wcsphSolver.h"

namespace pressureSolver
{
	// DECLARATIONS ================================================================================
	class Solver
	{
		//only select one
		static constexpr bool USE_WCSPH = true;
		static constexpr bool USE_PBF = false;

	public:
		/* Fields */
		pbf::PositionBasedFluids pSolver;
		wcsph::PositionBasedFluids wSolver;


		double timing_cns = 0.0;
		double timing_tns = 0.0;
		int n_sets = 0;

		/* Methods */
		Solver() = default;
		~Solver() = default;

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
			double num_iterations);
		void step(pressureSolver::NeighborhoodSearch& nsearch);
	};

	// DEFINITIONS ================================================================================

	inline void Solver::initialize(std::vector<Eigen::Vector3d>& fluid_positions,
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
		double params)
	{
		if constexpr (USE_PBF) {
			pSolver.initialize(fluid_positions,
				boundary_positions,
				velocities,
				accelerations,
				densities,
				boundary_masses,
				boundary_volumes,
				external_forces,
				is_scripted,
				fluid_mass,
				particle_diameter,
				dt,
				params);
		}
		if constexpr (USE_WCSPH) {
			wSolver.initialize(fluid_positions,
				boundary_positions,
				velocities,
				accelerations,
				densities,
				boundary_masses,
				boundary_volumes,
				external_forces,
				is_scripted,
				fluid_mass,
				particle_diameter,
				dt,
				params);
		}
	}

	inline void Solver::step(pressureSolver::NeighborhoodSearch& nsearch)
	{
		if constexpr (USE_PBF) {
			pSolver.step(nsearch);
		}
		if constexpr (USE_WCSPH) {
			wSolver.step(nsearch);
		}
	}

}