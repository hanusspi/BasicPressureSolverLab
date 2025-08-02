#pragma once
#include <string>
#include "io.h"
#include "neighborhoodSearch.h"
#include "sampling.h"

namespace pressureSolver {

	namespace scene {

		struct EmitterSetting {
			double velocity;
			Eigen::Vector3d rotation;
			Eigen::Vector3d positon;
			std::vector<Eigen::Vector3d> pattern;
			double time_step;
			double emitter_t;
			double radius;
			double scripted_time;
			double alpha;
			double boundaryThickness;
			double cylinderLength;
			double activation_time;
			double deactivation_time;
		};

		class ParticleEmitter {
		public:
			ParticleEmitter() = default;
			~ParticleEmitter() = default;
			void init(std::vector<Eigen::Vector3d>* fluid_positions,
				std::vector<Eigen::Vector3d>* scripted_positions,
				std::vector<Eigen::Vector3d>* boundary_positions,
				std::vector<Eigen::Vector3d>* velocities,
				std::vector<bool>* is_scripted,
				std::vector<double>* descript,
				double* t,
				double dt,
				double velocity,
				Eigen::Vector3d rotation,
				Eigen::Vector3d positon,
				double radius,
				double particle_diameter,
				double alpha,
				double scripted_time,
				double boundaryThickness,
				double cylinderLength,
				double activation_time,
				double deactivation_time);
			void emit(NeighborhoodSearch& nsearch);
		private:
			Eigen::Vector3d m_positon;
			Eigen::Vector3d m_rotation;
			double m_velocity;
			double m_radius;
			double m_particle_diameter;
			double m_activation_time;
			double m_deactivation_time;
			double m_alpha; // parameter to slow down emission for increased stability
			double m_scripted_time;
			int m_scirpted_steps;
			std::vector<Eigen::Vector3d> *m_fluid_positions;
			std::vector<Eigen::Vector3d>* m_scripted_positions;
			std::vector<Eigen::Vector3d> *m_boundary_positions;
			std::vector<Eigen::Vector3d> *m_velocities;
			double* m_t;
			double m_dt;
			std::vector<Eigen::Vector3d> m_newParticles;
			std::vector<Eigen::Vector3d>* m_boundaryParticles;
			double m_boundaryThickness;
			double m_cylinderLength;
			Eigen::Vector3d m_particle_velocity;
			double m_t_prev;
			double m_requiredTime;
			std::vector<bool> *m_is_scripted;
			std::vector<double>* m_descript;
		};

		/*class Scene {
		public:
			Scene() = default;
			~Scene() = default;
			void init(
				double* t,
				double dt,
				double radius,
				double rho_0,

				std::vector<Eigen::AlignedBox3d> fluidboxes,
				std::string boundary_path,
				std::string export_path,
				std::vector<EmitterSetting> emitters,

				Eigen::AlignedBox3d* viewbox,

				std::vector<Eigen::Vector3d>* fluid_positions,
				std::vector<Eigen::Vector3d>* velocities,
				std::vector<Eigen::Vector3d>* accelerations,
				std::vector<Eigen::Vector3d>* forces,
				std::vector<double>* densities,
				std::vector<bool>* is_scripted,
				std::vector<double>* descript,

				std::vector<Eigen::Vector3d>* boundary_positions,
				std::vector<double>* boundary_masses,
				std::vector<double>* boundary_volumes,
				learnSPH::NeighborhoodSearch& nsearch
			);
			void step();
		private:

			//Scene Data
			std::vector<ParticleEmitter> m_emitters;
			std::string m_boundary_path;
			std::vector<Eigen::AlignedBox3d> m_fluidboxes;
			Eigen::AlignedBox3d* m_viewbox;

			//Boundary Particle Data
			std::vector<Eigen::Vector3d> *m_boundary_positions;
			std::vector<double> *m_boundary_masses;
			std::vector<double>* m_boundary_volumes;
		
			//Fluid Particle Data
			std::vector<Eigen::Vector3d> *m_fluid_positions;
			std::vector<Eigen::Vector3d> *m_velocities;
			std::vector<Eigen::Vector3d> *m_accelerations;
			std::vector<Eigen::Vector3d> *m_forces;
			std::vector<double> *m_densities;
			std::vector<bool> *m_is_scripted;
			std::vector<double> *m_descript;
			//Simulation Data
			int m_frame = 0;
			double* m_t;
			double m_dt;
			double m_velocity;
			double m_radius;
			double m_diameter;
			double m_fluid_sampling_distance;
			double m_boundary_sampling_distance;
			double m_smoothing_length;
			double m_compact_support;
			double m_kernel_0;
			double m_fluid_mass;
			double m_rho_0;
			std::string m_export_path;
			NeighborhoodSearch *m_nsearch;

			void resize_data();
		};*/
	
	}
}	