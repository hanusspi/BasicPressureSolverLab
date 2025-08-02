#pragma once
#include <Eigen/Dense>
#include <vector>
#include <string>
#include "Solver.h"
#include "NeighborhoodSearch.h"
#include "scene.h"

namespace pressureSolver
{
	enum NeighborhoodId {
		FLUID_NEIGHBORHOOD = 0,
		BOUNDARY_NEIGHBORHOOD = 1
	};

	class Simulator
	{
	public:
		Simulator();
		~Simulator();

		void setParams(double particle_radius, std::string export_path, double rho_0, int freq, int fps, std::vector<double> params, double beta, double gamma);
		void setScene(std::vector<Eigen::AlignedBox3d> boxes, std::vector<scene::EmitterSetting> emitters, std::string boundary_path);
		void run(int iterations);
	private:

		void compute_external_forces();

		Solver m_solver;
		//scene::Scene m_scene;
		double m_particle_radius;
		NeighborhoodSearch m_nsearch;
		double m_particle_diameter;
		double m_fluid_sampling_distance;
		double m_boundary_sampling_distance;
		double m_smoothing_length;
		double m_compact_support;
		double m_kernel_0;
		double m_fluid_mass;
		std::string m_boundary_path;
		std::string m_export_path;
		Eigen::AlignedBox3d m_fluidbox;
		Eigen::AlignedBox3d m_viewbox;
		double m_rho_0;
		int m_freq;
		int m_fps;
		int m_iterations;
		double m_dt;
		int m_frame = 0;
		double m_t = 0.0;
		double m_t_frame = 0;
		double m_t_frame_step;
		double m_gamma;
		double m_beta;
		std::vector<Eigen::Vector3d> m_fluid_positions;
		std::vector<Eigen::Vector3d> m_scirpted_positions;
		std::vector<Eigen::Vector3d> m_boundary_positions;
		std::vector<Eigen::Vector3d> m_velocities;
		std::vector<Eigen::Vector3d> m_accelerations;
		std::vector<Eigen::Vector3d> m_forces;
		std::vector<Eigen::Vector3d> m_n_adh;
		std::vector<double> m_descript;
		std::vector<double> m_densities;
		std::vector<double> m_boundary_volumes;
		std::vector<double> m_boundary_masses;
		std::vector<double> m_params;
		Eigen::Vector3d m_gravity = Eigen::Vector3d{ 0,0,-9.81 };
		std::vector<scene::ParticleEmitter> m_emitters;
		std::vector<scene::EmitterSetting> m_emitter_settings;
		std::vector<bool> m_is_scripted;
		std::vector<Eigen::Vector3d> m_pattern;
		const bool log = false;
	};
}