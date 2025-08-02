#include "scene.h"
#include "kernel.h"

namespace pressureSolver {

	namespace scene {

#define PI 3.14159265358979323846
		void ParticleEmitter::init(
			std::vector<Eigen::Vector3d>* fluid_positions,
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
			double deactivation_time
		)
		{
			m_fluid_positions = fluid_positions;
			m_scripted_positions = scripted_positions;
			m_boundary_positions = boundary_positions;
			m_velocities = velocities;
			m_is_scripted = is_scripted;
			m_descript = descript;
			m_t = t;
			m_dt = dt;
			m_radius = radius;
			m_particle_diameter = particle_diameter;
			m_alpha = alpha;
			m_scripted_time = scripted_time;
			m_boundaryThickness = boundaryThickness;
			m_cylinderLength = cylinderLength;
			m_activation_time = activation_time;
			m_deactivation_time = deactivation_time;
			m_boundaryParticles = boundary_positions;
			m_positon = Eigen::Vector3d(0, 0, 0);
			m_rotation = Eigen::Vector3d(0, 0, 0);
			m_t_prev = 0;
			m_requiredTime = 0;
			m_newParticles.clear();
			m_particle_velocity = Eigen::Vector3d(0, 0, 0);
			m_scirpted_steps = 0;
			m_velocity = velocity;
			m_rotation = rotation;
			m_positon = positon;
			//create raw positions
			

			Eigen::Vector3d normal = m_rotation.normalized();

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
			transform.block<3, 1>(0, 3) = m_positon;

			// Hexagonal close packing implementation
			Eigen::Matrix4d worldTransform = transform;

			// Hexagonal close packing implementation
			double rowHeight = m_particle_diameter * 0.866;
			int numRows = static_cast<int>((2.0 * m_radius) / rowHeight);

			for (int row = -numRows / 2; row <= numRows / 2; ++row) {
				double y = row * rowHeight;
				double rowWidth = 2.0 * std::sqrt(m_radius * m_radius - y * y);
				double xOffset = (row % 2 == 0) ? 0.0 : m_particle_diameter / 2.0;

				int particlesInRow = static_cast<int>(rowWidth / m_particle_diameter);

				for (int col = -particlesInRow / 2; col <= particlesInRow / 2; ++col) {
					double x = col * m_particle_diameter + xOffset;

					if (x * x + y * y <= m_radius * m_radius) {
						// Create position in local space (x,y in disk, z=0)
						Eigen::Vector4d localPos(x, y, 0.0, 1.0);

						// Transform to world space
						Eigen::Vector4d worldPos = worldTransform * localPos;

						m_newParticles.push_back(worldPos.head<3>());
					}
				}
			}


			// Parameters for boundary cylinder
			double boundaryRadius = m_radius + 2* m_boundaryThickness;
			double angleStep = m_particle_diameter / boundaryRadius;  // Angular step based on particle size
			double heightStep = m_particle_diameter * 0.866;         // Vertical spacing (similar to emission grid)
			int numLayers = static_cast<int>(m_cylinderLength / heightStep);
			auto oldBoundarySize = m_boundaryParticles->size();
			// Generate rings of particles
			for (int layer = 0; layer <= numLayers; ++layer) {
				double z = layer * heightStep;

				// Calculate number of particles needed for this ring
				int numParticles = static_cast<int>(2.0 * PI * boundaryRadius / m_particle_diameter);

				for (int i = 0; i < numParticles; ++i) {
					double angle = i * (2.0 * PI / numParticles);

					// Create particle position in local space
					double x = boundaryRadius * std::cos(angle);
					double y = boundaryRadius * std::sin(angle);

					Eigen::Vector4d localPos(x, y, z, 1.0);
					Eigen::Vector4d worldPos = worldTransform * localPos;

					m_boundaryParticles->push_back(worldPos.head<3>());
				}

				// Add particles for the back wall (only on first layer)
				if (layer == 0) {
					// Sample a disk for the back wall
					double radialStep = m_particle_diameter * 0.866;
					int numRings = static_cast<int>(boundaryRadius / radialStep);

					// Center particle
					Eigen::Vector4d centerPos(0, 0, 0, 1.0);
					m_boundaryParticles->push_back((worldTransform * centerPos).head<3>());

					// Generate concentric rings
					for (int ring = 1; ring <= numRings; ++ring) {
						double ringRadius = ring * radialStep;
						int particlesInRing = static_cast<int>(2.0 * PI * ringRadius / m_particle_diameter);

						for (int i = 0; i < particlesInRing; ++i) {
							double angle = i * (2.0 * PI / particlesInRing);
							double x = ringRadius * std::cos(angle);
							double y = ringRadius * std::sin(angle);

							Eigen::Vector4d localPos(x, y, -m_particle_diameter, 1.0);
							Eigen::Vector4d worldPos = worldTransform * localPos;
							m_boundaryParticles->push_back(worldPos.head<3>());
						}
					}
				}
			}
			auto newBoundarySize = m_boundaryParticles->size();
			std::cout << "Boundary particles: " << newBoundarySize - oldBoundarySize << std::endl;
			m_particle_velocity = m_rotation.normalized() * m_velocity;
			m_scirpted_steps = std::min(0, static_cast<int>(m_scripted_time / m_dt));
			double requiredDistance = m_particle_diameter * m_alpha;
			double m_requiredTime = requiredDistance / m_velocity;


			//transform spawn positions 
			// 
			//create boundary positions
			// 
			//std::cout << "ParticleEmitter initialized" << std::endl;
		}

		void ParticleEmitter::emit(NeighborhoodSearch& nsearch)
		{
			//determine if particles should be emitted (based on timestepsize, velocity and alpha)
			if (*m_t >= m_activation_time && *m_t <= m_deactivation_time)
			{
				
				if (*m_t - m_t_prev >= m_requiredTime)
				{
					m_t_prev = *m_t;
					m_fluid_positions->insert(m_fluid_positions->end(), m_newParticles.begin(), m_newParticles.end());
					m_velocities->insert(m_velocities->end(), m_newParticles.size(), m_particle_velocity);
					m_is_scripted->insert(m_is_scripted->end(), m_newParticles.size(), true);
					m_descript->insert(m_descript->end(), m_newParticles.size(), (*m_t + m_scripted_time * m_dt));
					nsearch.resize_point_set(0, m_fluid_positions->front().data(), m_fluid_positions->size());
				}
				//create particles and set scripted_property and velocity
				
				//add particles

				//track recently added particles and activate after a certain time
				
			}
			//std::cout << "ParticleEmitter emitted" << std::endl;
		}

		/*void Scene::setSceneParams(
			std::vector<Eigen::AlignedBox3d> fluidboxes,
			std::string boundary_path,
			std::string export_path,
			std::vector<ParticleEmitter> emitters)
		{
			m_fluidboxes = fluidboxes;
			m_boundary_path = boundary_path;
			m_export_path = export_path;
			m_emitters = emitters;
		}
		*/
		/*void Scene::setSolverParams(
			double t,
			double dt,
			double radius,
			double rho_0
		)
		{
			m_t = t;
			m_dt = dt;
			m_radius = radius;
			m_rho_0 = rho_0;
			m_diameter = 2 * radius;
			m_fluid_sampling_distance = m_diameter;
			m_boundary_sampling_distance = 0.8 * m_diameter;
			m_smoothing_length = 1.2 * m_diameter;
			m_compact_support = 2.0 * m_smoothing_length;
			kernel::setAlpha(m_smoothing_length);
			m_kernel_0 = kernel::W(Eigen::Vector3d{ 0,0,0 });
			m_nsearch.set_search_radius(m_compact_support);
			m_nsearch.set_search_radius(m_compact_support);
			//std::cout << "Solver initialized" << std::endl;


		}*/

		/*void Scene::init(double* t,
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
			std::vector<int>* descript,

			std::vector<Eigen::Vector3d>* boundary_positions,
			std::vector<double>* boundary_masses,
			std::vector<double>* boundary_volumes,
			learnSPH::NeighborhoodSearch& nsearch
			) {

			m_t = t;
			m_dt = dt;
			m_radius = radius;
			m_rho_0 = rho_0;
			m_diameter = 2 * radius;
			m_fluid_sampling_distance = m_diameter;
			m_boundary_sampling_distance = 0.8 * m_diameter;
			m_smoothing_length = 1.2 * m_diameter;
			m_compact_support = 2.0 * m_smoothing_length;
			kernel::setAlpha(m_smoothing_length);
			m_kernel_0 = kernel::W(Eigen::Vector3d{ 0,0,0 });
			m_viewbox = viewbox;
			//m_nsearch->set_search_radius(m_compact_support);
#
			m_fluidboxes = fluidboxes;
			m_boundary_path = boundary_path;
			m_export_path = export_path;
			
			m_fluid_positions = fluid_positions;
			m_velocities = velocities;
			m_accelerations = accelerations;
			m_forces = forces;
			m_densities = densities;
			m_is_scripted = is_scripted;
			m_descript = descript;

			m_boundary_positions = boundary_positions;
			m_boundary_masses = boundary_masses;
			m_boundary_volumes = boundary_volumes;

			//create boundary particles
			std::cout << "Loading Mesh" << std::endl;
			const std::vector<learnSPH::TriMesh> meshes = learnSPH::read_tri_meshes_from_obj(m_boundary_path); //needs to be parametrized into constructor
			for (auto& mesh : meshes)
			{
				const learnSPH::TriMesh& box = mesh;
				for (const auto& v : box.vertices) {
					*m_viewbox = m_viewbox->extend(v);
				}
				learnSPH::sampling::triangle_mesh(*m_boundary_positions, box.vertices, box.triangles, m_boundary_sampling_distance);
			}
			for (auto& emitter : emitters)
			{
				m_emitters.push_back(ParticleEmitter());
				//m_emitters.back().init(m_fluid_positions, m_boundary_positions, m_velocities, m_is_scripted, m_descript, m_t, m_dt, m_velocity, m_radius, m_diameter, emitter.alpha, emitter.scripted_time, emitter.boundaryThickness, emitter.cylinderLength, emitter.activation_time, emitter.deactivation_time);
				//emitter.init(m_fluid_positions, m_boundary_positions, m_velocities, m_is_scripted, m_descript, m_t, m_dt, m_velocity, m_radius, m_diameter, m_alpha, emitter.getScriptedTime(), m_boundaryThickness, m_cylinderLength, emitter.getActivationTime(), emitter.getDeactivationTime());
			}
			//create fluid particles
			for (auto& m_fluidboxes : m_fluidboxes)
			{
				learnSPH::sampling::fluid_box(*m_fluid_positions, m_fluidboxes.min(), m_fluidboxes.max(), m_fluid_sampling_distance);
			}
			//add residual fluid particle for nsearch init
			m_fluid_positions->push_back(m_viewbox->min() - Eigen::Vector3d::Constant(-0.1));
			nsearch.add_point_set(m_fluid_positions->front().data(), m_fluid_positions->size());
			nsearch.add_point_set(m_boundary_positions->front().data(), m_boundary_positions->size());
			nsearch.set_active_search(0, 0);
			nsearch.set_active_search(0, 1);
			nsearch.set_active_search(1, 0);
			nsearch.set_active_search(1, 1);
			resize_data();
			nsearch.run();
#pragma omp parallel for
			for (int i = 0; i < m_boundary_positions->size(); i++) {// iterate over all particles stores in the nsearch point set

				double tmp_kern_sum = m_kernel_0;
				nsearch.for_each_neighbor(1, 1, i,
					[&](const int n)
					{
						tmp_kern_sum += learnSPH::kernel::W(m_boundary_positions->at(i) - m_boundary_positions->at(n));
					}
				);
				m_boundary_volumes->at(i) = 1.0 / tmp_kern_sum;
				m_boundary_masses->at(i) = m_boundary_volumes->at(i) * m_rho_0;
			}
		}

		void Scene::resize_data() {
			m_velocities->resize(m_fluid_positions->size(), Eigen::Vector3d{ 0.0,0.0,0.0 });
			m_accelerations->resize(m_fluid_positions->size(), Eigen::Vector3d{ 0.0,0.0,0.0 });
			m_forces->resize(m_fluid_positions->size(), Eigen::Vector3d{ 0.0,0.0,0.0 });
			m_densities->resize(m_fluid_positions->size(), 0.0);
			m_boundary_masses->resize(m_boundary_positions->size(), 0.0);
			m_boundary_volumes->resize(m_boundary_positions->size(), 0.0);
			m_is_scripted->resize(m_fluid_positions->size(), false);
			m_descript->resize(m_fluid_positions->size(), 0);
		}

		void Scene::step()
		{
			for (auto& emitter : m_emitters)
			{
				//emitter.emit();
			}
			//std::cout << "Scene acted" << std::endl;
		}*/
	}
}