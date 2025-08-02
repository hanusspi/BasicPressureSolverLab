#pragma once

#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include <algorithm>

#include <CompactNSearch/CompactNSearch>

#include "../pressureSolver/kernel.h"

namespace pressureSolver {
	namespace wcsph {
		class Simulation {
		public:
			Simulation(std::vector<Eigen::Vector3d> &positions, 
				const double volume,
				const double freq,
				const double fps,
				const double rho0,
				const double particle_radius);
			~Simulation();

			void run();
		private:
			void sampleFluidAndBoundaries();
			void computeDensities();
			void computeAccelerations();
			void integrateParticles();
			void exportParticlesToVTK();
			void removeOutsideParticles();

			double m_freq;
			double m_fps;
			double m_t;
			double m_t_frame;
			double m_t_frame_step; 
			double m_delta_t; 
			double m_kernel_0;
			int num_fluid_particles;
			int num_boundary_particles;
			int m_frame = 0;
			bool boundary;
			Eigen::Vector3d m_gravity = { 0.0,0.0,-9.81 };
			const double m_fluid_particle_mass;
			const double m_rho0;
			const double m_particle_radius;
			const double m_particle_diameter;
			const double m_fluid_sampling_distance;
			const double m_boundary_sampling_distance;
			const double m_smoothing_length;
			const double m_compact_support;
			const double m_B = 1000.0;
			const double m_fluid_vis = 0.0025;
			const double m_boundary_vis = 0.0;
			double offset = 0.5;
			double m_denom_eps;
			Eigen::Vector3d m_viewbox_min = Eigen::Vector3d::Constant(std::numeric_limits<double>::max()); 
			Eigen::Vector3d m_viewbox_max = Eigen::Vector3d::Constant(std::numeric_limits<double>::min());

			std::vector<Eigen::Vector3d> m_positions;
			std::vector<Eigen::Vector3d> m_velocities;
			std::vector<Eigen::Vector3d> m_accelerations;
			std::vector<Eigen::Vector3d> m_boundary_positions;
			std::vector<double> m_boundary_volumes;
			std::vector<double> m_boundary_masses;
			std::vector<double> m_densities;
			std::vector<double> m_pressures;

			CompactNSearch::NeighborhoodSearch m_nsearch;

			enum NeighborhoodId { 
				FLUID_NEIGHBORHOOD = 0, 
				BOUNDARY_NEIGHBORHOOD = 1 
			};

		};
	}
}