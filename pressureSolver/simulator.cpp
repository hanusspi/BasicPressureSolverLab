#include "simulator.h"
#include "util.h"
#include "kernel.h"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <algorithm>

#define PI 3.14159265358979323846

namespace pressureSolver
{
	Simulator::Simulator() : m_solver(), m_nsearch()
	{
	}
	Simulator::~Simulator()
	{
	}
	void Simulator::setParams(double particle_radius, std::string export_path, double rho_0, int freq, int fps, std::vector<double> params, double beta, double gamma)
	{
		m_particle_radius = particle_radius;
		//m_boundary_path = boundary_path;
		m_export_path = export_path;
		//m_fluidbox = fluidbox;
		m_rho_0 = rho_0;
		m_freq = freq;
		m_fps = fps;
		m_params = params;
		m_dt = 1.0 / freq;
		m_t_frame_step = 1.0 / fps;
		m_particle_diameter = 2.0 * m_particle_radius;
		m_fluid_sampling_distance = m_particle_diameter;
		m_boundary_sampling_distance = 0.6 * m_particle_diameter;
		m_smoothing_length = 1.2 * m_particle_diameter;
		m_compact_support = 2.0 * m_smoothing_length;
		kernel::setAlpha(m_smoothing_length);
		m_kernel_0 = kernel::W(Eigen::Vector3d{ 0,0,0 });
		m_nsearch.set_search_radius(m_compact_support);
		m_gamma = gamma;
		m_beta = beta;
		kernel::c = m_compact_support;
		kernel::a_coh = 32/(PI*pow(m_compact_support,9));
		kernel::a_adh = 0.007 / pow(m_compact_support, 3.25);
		kernel::sub_coh = pow(kernel::c, 6) / 64;
	}

	void Simulator::setScene(std::vector<Eigen::AlignedBox3d> boxes, std::vector<scene::EmitterSetting> emitters, std::string boundary_path)
	{
		std::cout << "Loading Mesh" << std::endl;
		const std::vector<pressureSolver::TriMesh> meshes = pressureSolver::read_tri_meshes_from_obj(boundary_path); //needs to be parametrized into constructor
		for (auto& mesh : meshes)
		{
			const pressureSolver::TriMesh& box = mesh;
			for (const auto& v : box.vertices) {
				m_viewbox = m_viewbox.extend(v);
			}
			pressureSolver::sampling::triangle_mesh(m_boundary_positions, box.vertices, box.triangles, m_boundary_sampling_distance);
		}
		double volume = 0;
		for (auto& m_fluidboxes : boxes)
		{
			pressureSolver::sampling::fluid_box(m_fluid_positions, m_fluidboxes.min(), m_fluidboxes.max(), m_fluid_sampling_distance);
			volume += m_fluidboxes.diagonal().prod();
		}
		volume = (volume == 0) ? 4.0/3.0 * PI * m_particle_radius * m_particle_radius * m_particle_radius : volume / m_fluid_positions.size();
		m_fluid_mass = m_rho_0 * volume;
		for (auto& emitter : emitters)
		{
			m_emitters.push_back(scene::ParticleEmitter());
			m_emitters.back().init(&m_fluid_positions, &m_scirpted_positions, &m_boundary_positions, &m_velocities, &m_is_scripted, &m_descript, &m_t, m_dt, emitter.velocity, emitter.rotation, emitter.positon, emitter.radius, m_particle_diameter, emitter.alpha, emitter.scripted_time, emitter.boundaryThickness, emitter.cylinderLength, emitter.activation_time, emitter.deactivation_time);
			m_emitter_settings.push_back(emitter);
		}
		m_boundary_masses.resize(m_boundary_positions.size());
		m_boundary_volumes.resize(m_boundary_positions.size());
		m_fluid_positions.push_back(m_viewbox.min() - Eigen::Vector3d::Constant(0.1));
		m_nsearch.add_point_set(m_fluid_positions.front().data(), m_fluid_positions.size());
		m_nsearch.add_point_set(m_boundary_positions.front().data(), m_boundary_positions.size());
		m_nsearch.set_active_search(0, 0);
		m_nsearch.set_active_search(0, 1);
		m_nsearch.set_active_search(1, 0);
		m_nsearch.set_active_search(1, 1);
		m_velocities.resize(m_fluid_positions.size(), Eigen::Vector3d{ 0.0,0.0,0.0 });
		m_accelerations.resize(m_fluid_positions.size(), Eigen::Vector3d{ 0.0,0.0,0.0 });
		m_forces.resize(m_fluid_positions.size(), Eigen::Vector3d{ 0.0,0.0,0.0 });
		m_densities.resize(m_fluid_positions.size(), 0.0);
		m_boundary_masses.resize(m_boundary_positions.size(), 0.0);
		m_boundary_volumes.resize(m_boundary_positions.size(), 0.0);
		m_is_scripted.resize(m_fluid_positions.size(), false);
		m_descript.resize(m_fluid_positions.size(), 0);
		m_n_adh.resize(m_fluid_positions.size(), Eigen::Vector3d{ 0.0,0.0,0.0 });
		m_fluid_positions.reserve(m_fluid_positions.size() + 100000);


		m_nsearch.run();
		std::cout << "init fluid particles: " << m_fluid_positions.size() << std::endl;


#pragma omp parallel for
		for (int i = 0; i < m_boundary_positions.size(); i++) {

			double tmp_kern_sum = m_kernel_0;
			m_nsearch.for_each_neighbor(1, 1, i,
				[&](const int n)
				{
					tmp_kern_sum += pressureSolver::kernel::W(m_boundary_positions[i] - m_boundary_positions[n]);
				}
			);
			m_boundary_volumes[i] = 1.0 / tmp_kern_sum;
			m_boundary_masses[i] = m_boundary_volumes[i] * m_rho_0;
		}
		m_solver.initialize(m_fluid_positions, 
			m_boundary_positions, 
			m_velocities, 
			m_accelerations, 
			m_densities, 
			m_boundary_masses, 
			m_boundary_volumes, 
			m_forces, 
			m_is_scripted, 
			m_fluid_mass, 
			m_particle_diameter, 
			m_dt, m_params[0]);
		
		for (auto& setting : m_emitter_settings) {
			Eigen::Vector3d normal = setting.rotation.normalized();

			// Find arbitrary vector not parallel to normal to create coordinate system
			Eigen::Vector3d tangent;
			if (std::abs(normal.dot(Eigen::Vector3d::UnitX())) < 0.9) {
				// Use cross product with X axis if normal is not too close to X
				tangent = normal.cross(Eigen::Vector3d::UnitX()).normalized();
			}
			else {
				// Use Y axis if normal is close to X
				tangent = normal.cross(Eigen::Vector3d::UnitY()).normalized();
			}

			// Complete orthonormal basis
			Eigen::Vector3d bitangent = normal.cross(tangent);
			
			// Create rotation matrix
			Eigen::Matrix4d transform = Eigen::Matrix4d::Identity();
			transform.block<3, 1>(0, 0) = tangent;
			transform.block<3, 1>(0, 1) = bitangent;
			transform.block<3, 1>(0, 2) = normal;

			// Add translation
			transform.block<3, 1>(0, 3) = setting.positon;

			// Hexagonal close packing implementation
			Eigen::Matrix4d worldTransform = transform;

			setting.time_step = (m_particle_diameter * setting.alpha / std::abs(setting.velocity)) ;
			double m_radius = setting.radius;
			double rowHeight = m_particle_diameter * 0.866;
			setting.emitter_t = setting.activation_time;
			double cfl_velocity = 0.4 * m_particle_diameter / m_dt;
			double squared_cfl_velocity = cfl_velocity * cfl_velocity;
			double squared_velocity = setting.velocity * setting.velocity;
			if (squared_velocity > squared_cfl_velocity)
			{
				setting.velocity = setting.velocity * (cfl_velocity / (sqrt(squared_velocity)));
			}
			int numRows = static_cast<int>((2.0 * m_radius) / rowHeight);
			double tmp = m_particle_diameter * 1.2;
			std::cout << "cfl velocity" << cfl_velocity << std::endl;
			for (int row = -numRows / 2; row <= numRows / 2; ++row) {
				double y = row * rowHeight;
				double rowWidth = 2.0 * std::sqrt(m_radius * m_radius - y * y);
				double xOffset = (row % 2 == 0) ? 0.0 : tmp / 2.0;

				int particlesInRow = static_cast<int>(rowWidth / tmp);

				for (int col = -particlesInRow / 2; col <= particlesInRow / 2; ++col) {
					double x = col * tmp + xOffset;

					if (x * x + y * y <= m_radius * m_radius) {
						// Create position in local space (x,y in disk, z=0)
						Eigen::Vector4d localPos(x, y, 0.0 + tmp, 1.0);

						// Transform to world space
						Eigen::Vector4d worldPos = worldTransform * localPos;

						setting.pattern.push_back(worldPos.head<3>());
					}
				}
			}
		}
		std::cout << "finished init, nr fluid par: " << m_fluid_positions.size() << std::endl;
		//m_forces = std::vector<Eigen::Vector3d>(m_fluid_positions.size(), Eigen::Vector3d{ 0.0,0.0,-9.81 });
	}

	void Simulator::run(int iterations)
	{
		std::cout << m_fluid_mass << std::endl;
		std::cout << kernel::W(Eigen::Vector3d{ 0,0,m_particle_radius }) << std::endl;
		static bool fileInitialized = false;
		std::string filename = m_export_path + "log.csv";
		std::ofstream file;
		if (log) {

			if (!fileInitialized) {
				// Open in write mode and add header for the first time
				file.open(filename, std::ios::out);
				file << "Step,PotentialEnergy,AccumulatedOverdensity,MaxDens,1pDens,5pDens,AvgDensity,#Particles\n";
				fileInitialized = true;
			}
			else {
				// Open in append mode for subsequent writes
				file.open(filename, std::ios::app);
			}
		}
		double frame_t = 0.0;
		double height_0 = -1.0;
		auto start = std::chrono::high_resolution_clock::now();
		std::vector<Eigen::Vector3d> boundaryVector(m_boundary_positions.size(),
			Eigen::Vector3d::Zero()
			);

		util::exportParticlesToVTK(m_export_path + "boundary_particles" + std::to_string(m_frame) + ".vtk", m_boundary_positions, m_boundary_masses, boundaryVector);
		std::vector<double> doubleArray(m_is_scripted.size());
		m_solver.step(m_nsearch);
		for (int i = 0; i < iterations; i++)
		{
			double energy = 0.0;
			double accumulatedOverdensity = 0.0;
			m_t += m_dt;
			
			for(auto& setting : m_emitter_settings) {
				if (m_t>= setting.emitter_t && m_t>=setting.activation_time && m_t<=setting.deactivation_time) {
					setting.emitter_t += setting.time_step;
					std::cout << "Adding particles: " << m_t << ", " << setting.emitter_t << ", " << m_fluid_positions.size() << std::endl;
					
					for (int i = 0; i < setting.pattern.size(); i++) {
						m_fluid_positions.push_back(setting.pattern[i]);
						m_velocities.push_back(setting.rotation.normalized() * abs(setting.velocity));
						m_accelerations.push_back(Eigen::Vector3d{ 0,0,0 });
						m_forces.push_back(Eigen::Vector3d{ 0,0,0 });
						m_densities.push_back(1000.0);
						m_is_scripted.push_back(true);
						m_descript.push_back(m_t + setting.scripted_time * m_dt);
					}
				
				}
			}
			m_nsearch.resize_point_set(0, m_fluid_positions.front().data(), m_fluid_positions.size());
			m_nsearch.run();
			m_accelerations.resize(m_fluid_positions.size());
			m_forces.resize(m_fluid_positions.size());
			m_densities.resize(m_fluid_positions.size());
			m_n_adh.resize(m_fluid_positions.size());
			
			for (int i = 0; i < m_fluid_positions.size(); i++)
			{
				if (m_descript[i] <= m_t) {
					m_is_scripted[i] = false;
				}
			}

			compute_external_forces();
			m_solver.step(m_nsearch);

			if (m_t >= m_t_frame) {
				util::showProgressBar(i, iterations);
				util::removeOutsideParticles(&m_viewbox, &m_nsearch, &m_fluid_positions, &m_velocities, &m_densities, &m_accelerations, &m_forces, &m_n_adh, &m_descript, &m_is_scripted);

				if (log)
				{
					double localMaxDensity = -std::numeric_limits<double>::infinity();

					double localSumDensity = 0.0;

						double myMaxDensity = -std::numeric_limits<double>::infinity();

#pragma omp for
						for (int i = 0; i < m_fluid_positions.size(); i++)
						{
							energy += m_fluid_mass * (
								(m_fluid_positions[i].z() - height_0) * 9.81
								+ 0.5 * m_velocities[i].squaredNorm()
								);
							accumulatedOverdensity += std::max(0.0, (m_densities[i] - m_rho_0));

							if (m_densities[i] > myMaxDensity)
								myMaxDensity = m_densities[i];

							localSumDensity += m_densities[i];

#pragma omp critical
							{
								if (myMaxDensity > localMaxDensity)
								{
									localMaxDensity = myMaxDensity;
								}
							}
						}
					

					// Compute average density:
					double averageDensity = localSumDensity / static_cast<double>(m_densities.size());

					// Create a copy of m_densities so we can sort it without altering the original:
					std::vector<double> densCopy = m_densities;
					std::sort(densCopy.begin(), densCopy.end());

					// For 1% "max" we want the 99th percentile; for 5% "max" we want the 95th percentile
					// Index calculation: idx = percentile * (N - 1).
					int N = static_cast<int>(densCopy.size());
					int idx99 = static_cast<int>(std::floor(0.99 * (N - 1)));
					int idx95 = static_cast<int>(std::floor(0.95 * (N - 1)));
					double density1p = densCopy[idx99];
					double density5p = densCopy[idx95];

					// Now log everything including the new statistics:
					if (file.is_open())
					{
						file << m_t << ","
							<< energy << ","
							<< accumulatedOverdensity << ","
							<< localMaxDensity << ","
							<< density1p << ","
							<< density5p << ","
							<< averageDensity << ","
							<< m_fluid_positions.size()
							<< "\n";
					}
					else
					{
						std::cerr << "Error: Unable to open file " << filename << "\n";
					}
				}
				
				doubleArray.resize(m_is_scripted.size());
				for (size_t i = 0; i < m_is_scripted.size(); ++i) {
					doubleArray[i] = static_cast<double>(m_is_scripted[i]);
				}
				std::stringstream ss;
				ss << std::setw(5) << std::setfill('0') << m_frame;
				std::string frame_str = ss.str();
				if (m_fluid_positions.size() > 400) {
					std::vector<std::vector<Eigen::Vector3d>> outVector = { m_velocities, m_accelerations };
					std::vector<std::string> vector_names = { "velocity", "acceleration" };
					util::exportParticlesToVTK(m_export_path + "fluid_particles" + frame_str + ".vtk", m_fluid_positions, m_densities, m_velocities);
					//util::exportMesthToVTK(m_export_path + "mesh_new_velo" + std::to_string(m_frame) + ".vtk", as_const(m_fluid_positions), /*m_densities,*/ /*-0.55,*/m_nsearch, m_compact_support, m_particle_radius, m_fluid_mass, m_fluidbox.min(), m_viewbox);
					
				}
				m_t_frame += m_t_frame_step;
				m_frame += 1;
			}
		}
		if (log) {
			if (file.is_open())
				file.close();
		}
		
		auto end = std::chrono::high_resolution_clock::now();
		std::cout << "Simulation took " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << " seconds" << std::endl;
	}

	void Simulator::compute_external_forces()
	{
		//gravity
		for (int i = 0; i < m_fluid_positions.size(); i++) {
			m_nsearch.for_each_neighbor(0, 0, i,
				[&](const int k)
				{
					m_n_adh[i] += (m_fluid_mass / m_densities[k]) * kernel::W_grad(m_fluid_positions[i] - m_fluid_positions[k]);
				}
			);
			m_n_adh[i] = m_n_adh[i] * m_compact_support;
		}
		for (int i = 0; i < m_fluid_positions.size(); i++)
		{
			m_accelerations[i] = Eigen::Vector3d{ 0.0,0.0,0.0 };
			m_accelerations[i] += m_gravity;
			//compute cohesion forces
			double fluid_mass_sqrd = m_fluid_mass * m_fluid_mass;
			m_nsearch.for_each_neighbor(0, 0, i,
				[&](const int k)
				{
					Eigen::Vector3d f_cohesion = -m_gamma * fluid_mass_sqrd * kernel::W_cohesion(m_fluid_positions[i] - m_fluid_positions[k]) * ((m_fluid_positions[i]-m_fluid_positions[k])/(m_fluid_positions[i] - m_fluid_positions[k]).norm());
					Eigen::Vector3d f_curvature = -m_gamma * m_fluid_mass * (m_n_adh[i] - m_n_adh[k]);
					m_accelerations[i] += ((2 * m_rho_0 / (m_densities[i] + m_densities[k])) * (f_cohesion + f_curvature)) / m_fluid_mass;
				}
			);
			////compute adhesion forces
			m_nsearch.for_each_neighbor(0, 1, i,
				[&](const int k)
				{
					m_accelerations[i] += -m_beta * /*m_fluid_mass **/ m_boundary_masses[k] * pressureSolver::kernel::W_adhesion(m_fluid_positions[i] - m_boundary_positions[k]) * (m_fluid_positions[i] - m_boundary_positions[k]) / (m_fluid_positions[i] - m_boundary_positions[k]).norm(); //divide by fluid mass to get acceleration
				}
			);

		}
	}
}