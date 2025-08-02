#pragma once
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <string>
#include "NeighborhoodSearch.h"
#include "../pressureSolver/kernel.h"
#include "../pressureSolver/sampling.h"
#include "../pressureSolver/io.h"
#include "../pressureSolver/kernel.h"
#include "../pressureSolver/marching_cubes_lut.h"
#include "../pressureSolver/NeighborhoodSearch.h"

namespace util {

	void sampleFluidAndBoundaries(Eigen::AlignedBox3d fluidBox,
		Eigen::AlignedBox3d& viewBox,
		std::string obj_path,
		double particle_radius,
		std::vector<Eigen::Vector3d>& boundary_positions,
		std::vector<Eigen::Vector3d>& fluid_positions,
		std::vector<Eigen::Vector3d>& velocities,
		std::vector<double>& densities,
		std::vector<Eigen::Vector3d>& accelerations,
		std::vector<double>& boundary_masses,
		std::vector<double>& boundary_volumes,
		pressureSolver::NeighborhoodSearch& nsearch,
		double kernel_0,
		double rho0,
		double fluid_sampling_dist,
		double boundary_sampling_dist);
	void showProgressBar(int progress, int total, int barWidth = 50);
	void removeOutsideParticles(Eigen::AlignedBox3d* viewbox,
		pressureSolver::NeighborhoodSearch* nsearch,
		std::vector<Eigen::Vector3d>* fluid_positions,
		std::vector<Eigen::Vector3d>* velocities,
		std::vector<double>* densities,
		std::vector<Eigen::Vector3d>* accelerations,
		std::vector<Eigen::Vector3d>* forces,
		std::vector<Eigen::Vector3d>* n_adh,
		std::vector<double>* descript,
		std::vector<bool>* is_scripted
	);
	void exportParticlesToVTK(std::string filename, std::vector<Eigen::Vector3d>& positions, std::vector<double>& densities);
	void exportParticlesToVTK(std::string filename, std::vector<Eigen::Vector3d>& positions, std::vector<double>& densities, std::vector<Eigen::Vector3d>& velocities);
	void exportParticlesToVtk(std::string path,
		const std::vector<Eigen::Vector3d>& particles,
		const std::vector<double>& particle_scalar_data,
		const std::vector<std::vector<Eigen::Vector3d>>& particle_vector_data_arrays,
		const std::vector<std::string>& vector_data_names);
	void exportMesthToVTK(std::string filename,
		const std::vector<Eigen::Vector3d>& positions,
		//std::vector<double>& densities,
		pressureSolver::NeighborhoodSearch& nsearch,
		//double c_,
		double compact_support_,
		double dist_,
		double fluid_particle_mass_,
		Eigen::Vector3d origin_,
		Eigen::AlignedBox3d viewbox
	);

}