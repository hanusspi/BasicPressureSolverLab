#pragma once
#include "pbf.h"
#include "wcsphSolver.h"

namespace pressureSolver
{
	// Enum for solver types
	enum class SolverType {
		WCSPH,
		PBF
	};

	// DECLARATIONS ================================================================================
	class Solver
	{
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

		void setSolverType(SolverType type);

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
			double params);

		void step(pressureSolver::NeighborhoodSearch& nsearch);

	private:
		SolverType m_solver_type = SolverType::WCSPH; // Default to WCSPH
	};

	// DEFINITIONS ================================================================================

	inline void Solver::setSolverType(SolverType type) {
		m_solver_type = type;
	}

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
		switch (m_solver_type) {
		case SolverType::PBF:
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
				static_cast<int>(params)); // PBF expects num_iterations as int
			break;
		case SolverType::WCSPH:
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
				params); // WCSPH expects stiffness as double
			break;
		}
	}

	inline void Solver::step(pressureSolver::NeighborhoodSearch& nsearch)
	{
		switch (m_solver_type) {
		case SolverType::PBF:
			pSolver.step(nsearch);
			break;
		case SolverType::WCSPH:
			wSolver.step(nsearch);
			break;
		}
	}
}