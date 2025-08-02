
#include "util.h"
#include "definitions.h"
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <cmath>

namespace util
{
#define PI 3.14159265358979323846

	struct Vector2iHash {
		std::size_t operator()(const Eigen::Vector2i& v) const noexcept {
			// Combine the hashes of the x and y components
			std::size_t h1 = std::hash<int>{}(v.x());
			std::size_t h2 = std::hash<int>{}(v.y());
			return h1 ^ (h2 << 1);  // XOR and bit-shift
		}
	};

	// Custom equality comparison for Eigen::Vector2i
	struct Vector2iEqual {
		bool operator()(const Eigen::Vector2i& v1, const Eigen::Vector2i& v2) const noexcept {
			return v1.x() == v2.x() && v1.y() == v2.y();
		}
	};

	Eigen::Vector3d grid2coordinate(uint i, uint j, uint k, double dist, Eigen::Vector3d grid_origin)
	{
		//returns the world coordinate of the grid index i,j,k
		return Eigen::Vector3d(i, j, k) * dist + grid_origin;
	}

	void sampleFluidAndBoundaries(
		Eigen::AlignedBox3d fluidBox,
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
		double boundary_sampling_dist

	) {
		std::cout << "Loading mesh" << std::endl;
		const std::vector<pressureSolver::TriMesh> meshes = pressureSolver::read_tri_meshes_from_obj(obj_path); //needs to be parametrized into constructor
		const pressureSolver::TriMesh& box = meshes[0];

		std::cout << "Calculating boundary box" << std::endl;
		for (const auto& v : box.vertices) {
			viewBox = viewBox.extend(v);

		}

		pressureSolver::sampling::fluid_box(fluid_positions, fluidBox.min(), fluidBox.max(), fluid_sampling_dist);

		std::cout << "Sampling fluid and boundary particles" << std::endl;
		pressureSolver::sampling::triangle_mesh(boundary_positions, box.vertices, box.triangles, boundary_sampling_dist);
		nsearch.add_point_set(fluid_positions.front().data(), fluid_positions.size());
		nsearch.add_point_set(boundary_positions.front().data(), boundary_positions.size());
		nsearch.set_active_search(0, 0);
		nsearch.set_active_search(0, 1);
		nsearch.set_active_search(1, 0);
		nsearch.set_active_search(1, 1);

		velocities.resize(fluid_positions.size(), Eigen::Vector3d{ 0.0,0.0,0.0 });

		densities.resize(fluid_positions.size(), 0.0);
		accelerations.resize(fluid_positions.size(), Eigen::Vector3d{ 0.0,0.0,0.0 });
		boundary_masses.resize(boundary_positions.size(), 0.0);
		boundary_volumes.resize(boundary_positions.size(), 0.0);
		std::cout << "boundary particles: " << boundary_positions.size() << std::endl;
		std::cout << "fluid particles: " << fluid_positions.size() << std::endl;
		nsearch.run();
		std::cout << "Running neighborhood search" << std::endl;

#pragma omp parallel for
		for (int i = 0; i < boundary_positions.size(); i++) {// iterate over all particles stores in the nsearch point set
			
			double tmp_kern_sum = kernel_0;
			nsearch.for_each_neighbor(1, 1, i,
				[&](const int n)
				{
					tmp_kern_sum += pressureSolver::kernel::W(boundary_positions[i] - boundary_positions[n]);
				}
			);
			boundary_volumes[i] = 1.0 / tmp_kern_sum;
			boundary_masses[i] = boundary_volumes[i] * rho0;
		}
		std::cout << "Finished boundary correction" << std::endl;
	}

	void showProgressBar(int progress, int total, int barWidth)
	{
		std::cout << "[";
		int pos = barWidth * progress / total;
		for (int i = 0; i < barWidth; ++i) {
			if (i < pos) std::cout << "=";
			else if (i == pos) std::cout << ">";
			else std::cout << " ";
		}
		std::cout << "] " << int((double(progress) * 100.0) / double(total) ) << " %\r";
		std::cout.flush();
	}

	void exportParticlesToVTK(std::string filename, std::vector<Eigen::Vector3d>& positions, std::vector<double>& densities)
	{
		pressureSolver::write_particles_to_vtk(filename, positions, densities);
	}

	void exportParticlesToVTK(std::string filename, std::vector<Eigen::Vector3d>& positions, std::vector<double>& densities, std::vector<Eigen::Vector3d>& velocities)
	{
		pressureSolver::write_particles_to_vtk(filename, positions, densities, velocities);
	}

	void exportParticlesToVtk(std::string path,
		const std::vector<Eigen::Vector3d>& particles,
		const std::vector<double>& particle_scalar_data,
		const std::vector<std::vector<Eigen::Vector3d>>& particle_vector_data_arrays,
		const std::vector<std::string>& vector_data_names)
	{
		pressureSolver::write_particles_to_vtk(path, particles, particle_scalar_data, particle_vector_data_arrays, vector_data_names);
	}

	void removeOutsideParticles(
		Eigen::AlignedBox3d* viewbox,
		pressureSolver::NeighborhoodSearch* nsearch,
		std::vector<Eigen::Vector3d>* fluid_positions,
		std::vector<Eigen::Vector3d>* velocities,
		std::vector<double>* densities,
		std::vector<Eigen::Vector3d>* accelerations,
		std::vector<Eigen::Vector3d>* forces,
		std::vector<Eigen::Vector3d>* n_adh,
		std::vector<double>* descript,
		std::vector<bool>* is_scripted
		)
	{
		//std::cout << "Removing outside particles" << std::endl;
		std::vector<bool> to_erase(fluid_positions->size(), false); // local variable to store which particles to erase
#pragma omp parallel for
		for (int i = fluid_positions->size() - 1; i >= 0; i--)
		{
			if (!viewbox->contains((*fluid_positions)[i]))
			{
				to_erase[i] = true;
			}
		}
		for (int i = fluid_positions->size() - 1; i >= 0; i--)
		{
			if (to_erase[i])
			{
				//std::cout << "Erasing particle: " << fluid_positions->at(i).transpose() << ", " << densities->at(i) <<", " << accelerations->at(i).transpose() << std::endl;
				fluid_positions->erase(fluid_positions->begin() + i);
				velocities->erase(velocities->begin() + i);
				densities->erase(densities->begin() + i);
				accelerations->erase(accelerations->begin() + i);
				forces->erase(forces->begin() + i);
				n_adh->erase(n_adh->begin() + i);
				descript->erase(descript->begin() + i);
				is_scripted->erase(is_scripted->begin() + i);
			}
		}
		nsearch->resize_point_set(0, fluid_positions->front().data(), fluid_positions->size());
		nsearch->run();
	}

	class LaplacianSmoother {
	private:
		struct Edge {
			int v1, v2;
			double length;
	
			Edge(int v1, int v2, double length) : v1(v1), v2(v2), length(length) {}
	
			bool operator<(const Edge& other) const {
				return length > other.length;  
			}
		};
	
	public:
		static void processMesh(
			std::vector<Eigen::Vector3d>& vertices,
			std::vector<std::array<int, 3>>& triangles,
			std::vector<Eigen::Vector3d>& normals,
			int smoothing_iterations = 5,
			double smoothing_lambda = 0.5,
			double edge_length_threshold = 0.01 
		) {
			smoothMesh(vertices, triangles, normals, smoothing_iterations, smoothing_lambda);
			decimateMesh(vertices, triangles, normals, edge_length_threshold);
		}
	
	private:
		static void smoothMesh(
			std::vector<Eigen::Vector3d>& vertices,
			const std::vector<std::array<int, 3>>& triangles,
			std::vector<Eigen::Vector3d>& normals,
			int iterations,
			double lambda
		) {
			std::vector<std::vector<int>> adjacency(vertices.size());
			buildAdjacencyList(triangles, adjacency);
	
			std::vector<Eigen::Vector3d> new_positions(vertices.size());
			std::vector<Eigen::Vector3d> new_normals(vertices.size());
	
			for (int iter = 0; iter < iterations; ++iter) {
				for (size_t i = 0; i < vertices.size(); ++i) {
					if (adjacency[i].empty()) {
						new_positions[i] = vertices[i];
						new_normals[i] = normals[i];
						continue;
					}
	
					Eigen::Vector3d avg_pos = Eigen::Vector3d::Zero();
					Eigen::Vector3d avg_normal = Eigen::Vector3d::Zero();
					double total_weight = 0.0;
					//calculate centroid
					for (int neighbor : adjacency[i]) {
						double weight = 1.0 / (vertices[i] - vertices[neighbor]).norm();
						avg_pos += weight * vertices[neighbor];
						avg_normal += weight * normals[neighbor];
						total_weight += weight;
					}
	
					avg_pos /= total_weight;
					avg_normal /= total_weight;
					avg_normal.normalize();
					//update position and normal
					new_positions[i] = vertices[i] + lambda * (avg_pos - vertices[i]);
					new_normals[i] = normals[i] + lambda * (avg_normal - normals[i]);
					new_normals[i].normalize();
				}
	
				vertices = new_positions;
				normals = new_normals;
			}
		}
	
		static void decimateMesh(
			std::vector<Eigen::Vector3d>& vertices,
			std::vector<std::array<int, 3>>& triangles,
			std::vector<Eigen::Vector3d>& normals,
			double length_threshold
		) {
			std::priority_queue<Edge> edge_heap;
	
			for (const auto& triangle : triangles) {
				for (int i = 0; i < 3; ++i) {
					int v1 = triangle[i];
					int v2 = triangle[(i + 1) % 3];
					if (v1 < v2) {  
						double length = (vertices[v1] - vertices[v2]).norm();
						if (length < length_threshold) {
							edge_heap.push(Edge(v1, v2, length));
						}
					}
				}
			}
	
			std::vector<bool> valid_vertex(vertices.size(), true);
			std::vector<bool> valid_triangle(triangles.size(), true);
	
			while (!edge_heap.empty()) {
				Edge edge = edge_heap.top();
				edge_heap.pop();
	
				if (!valid_vertex[edge.v1] || !valid_vertex[edge.v2]) continue;
	
				vertices[edge.v1] = (vertices[edge.v1] + vertices[edge.v2]) * 0.5;
				normals[edge.v1] = (normals[edge.v1] + normals[edge.v2]).normalized();
				valid_vertex[edge.v2] = false;
	
				for (size_t i = 0; i < triangles.size(); ++i) {
					if (!valid_triangle[i]) continue;
	
					auto& triangle = triangles[i];
					for (int j = 0; j < 3; ++j) {
						if (triangle[j] == edge.v2) {
							triangle[j] = edge.v1;
						}
					}

					if (triangle[0] == triangle[1] ||
						triangle[1] == triangle[2] ||
						triangle[2] == triangle[0]) {
						valid_triangle[i] = false;
					}
				}
			}
			compactMesh(vertices, triangles, normals, valid_vertex, valid_triangle);
		}
	
		static void compactMesh(
			std::vector<Eigen::Vector3d>& vertices,
			std::vector<std::array<int, 3>>& triangles,
			std::vector<Eigen::Vector3d>& normals,
			const std::vector<bool>& valid_vertex,
			const std::vector<bool>& valid_triangle
		) {
			std::vector<int> vertex_map(vertices.size(), -1);
			int new_vertex_count = 0;
	
			for (size_t i = 0; i < vertices.size(); ++i) {
				if (valid_vertex[i]) {
					vertex_map[i] = new_vertex_count++;
				}
			}
	
			std::vector<Eigen::Vector3d> new_vertices;
			std::vector<Eigen::Vector3d> new_normals;
			new_vertices.reserve(new_vertex_count);
			new_normals.reserve(new_vertex_count);
	
			for (size_t i = 0; i < vertices.size(); ++i) {
				if (valid_vertex[i]) {
					new_vertices.push_back(vertices[i]);
					new_normals.push_back(normals[i]);
				}
			}
	
			std::vector<std::array<int, 3>> new_triangles;
			for (size_t i = 0; i < triangles.size(); ++i) {
				if (valid_triangle[i]) {
					std::array<int, 3> new_triangle;
					for (int j = 0; j < 3; ++j) {
						new_triangle[j] = vertex_map[triangles[i][j]];
					}
					new_triangles.push_back(new_triangle);
				}
			}
	
			vertices = std::move(new_vertices);
			triangles = std::move(new_triangles);
			normals = std::move(new_normals);
		}
	
		static void buildAdjacencyList(
			const std::vector<std::array<int, 3>>& triangles,
			std::vector<std::vector<int>>& adjacency
		) {
			for (const auto& triangle : triangles) {
				for (int i = 0; i < 3; ++i) {
					int v1 = triangle[i];
					int v2 = triangle[(i + 1) % 3];
	
					adjacency[v1].push_back(v2);
					adjacency[v2].push_back(v1);
				}
			}
	
			for (auto& neighbors : adjacency) {
				std::unordered_set<int> unique(neighbors.begin(), neighbors.end());
				neighbors.assign(unique.begin(), unique.end());
			}
		}
	};



class MeshSmoother {
public:
	struct SmoothingParams {
		int iterations = 5;                  // Number of smoothing iterations
		double max_curvature = 2.0;          // Maximum allowed mean curvature
		double min_feature_size = 0.1;       // Minimum feature size to preserve
		double tension_weight = 0.5;         // Weight for surface tension effect (0-1)
		double volume_preserve = 0.8;        // Volume preservation strength (0-1)
		double edge_length_threshold = 0.01; // Minimum edge length for decimation
	};

	static void processMesh(
		std::vector<Eigen::Vector3d>& vertices,
		std::vector<std::array<int, 3>>& triangles,
		std::vector<Eigen::Vector3d>& normals,
		const SmoothingParams& params = SmoothingParams()
	) {

		std::vector<std::vector<int>> adjacency(vertices.size());
		std::vector<std::vector<int>> vertex_triangles(vertices.size());
		buildConnectivity(triangles, adjacency, vertex_triangles);

		double initial_volume = calculateMeshVolume(vertices, triangles);
		std::vector<double> vertex_areas = calculateVertexAreas(vertices, triangles, vertex_triangles);

		std::vector<Eigen::Vector3d> new_positions(vertices.size());
		std::vector<Eigen::Vector3d> new_normals(vertices.size());

#pragma omp parallel for
		for (int iter = 0; iter < params.iterations; ++iter) {
			std::vector<Eigen::Vector3d> curvature_normals(vertices.size(), Eigen::Vector3d::Zero());
			std::vector<double> mean_curvatures(vertices.size(), 0.0);

			calculateMeanCurvature(vertices, adjacency, vertex_areas, curvature_normals, mean_curvatures);

			for (size_t i = 0; i < vertices.size(); ++i) {
				if (adjacency[i].empty()) continue;

				Eigen::Vector3d tension_force = calculateSurfaceTension(
					vertices[i], adjacency[i], vertices, normals[i], params.tension_weight
				);

				Eigen::Vector3d move_dir = curvature_normals[i];
				double curvature_magnitude = mean_curvatures[i];

				if (curvature_magnitude > params.max_curvature) {
					move_dir *= params.max_curvature / curvature_magnitude;
				}

				move_dir += tension_force * params.min_feature_size;

				Eigen::Vector3d new_pos = vertices[i] + move_dir * params.min_feature_size;

				new_positions[i] = new_pos;
				//maybe take out of loop
				Eigen::Vector3d avg_normal = Eigen::Vector3d::Zero();
				for (int neighbor : adjacency[i]) {
					double weight = 1.0 / (vertices[i] - vertices[neighbor]).norm();
					avg_normal += weight * normals[neighbor];
				}
				new_normals[i] = avg_normal.normalized();
			}
			//maybe take out of loop
			if (params.volume_preserve > 0) {
				double current_volume = calculateMeshVolume(new_positions, triangles);
				double scale = std::pow(initial_volume / current_volume, 1.0 / 3.0);

				Eigen::Vector3d centroid = Eigen::Vector3d::Zero();
				for (const auto& pos : new_positions) centroid += pos;
				centroid /= new_positions.size();

				for (auto& pos : new_positions) {
					pos = centroid + (pos - centroid) *
						(1.0 + (scale - 1.0) * params.volume_preserve);
				}
			}

			vertices = new_positions;
			normals = new_normals;
		}

		decimateMesh(vertices, triangles, normals, params.edge_length_threshold);
	}

private:
	static void buildConnectivity(
		const std::vector<std::array<int, 3>>& triangles,
		std::vector<std::vector<int>>& adjacency,
		std::vector<std::vector<int>>& vertex_triangles
	) {
		for (size_t i = 0; i < triangles.size(); ++i) {
			const auto& tri = triangles[i];
			for (int j = 0; j < 3; ++j) {
				vertex_triangles[tri[j]].push_back(i);
				adjacency[tri[j]].push_back(tri[(j + 1) % 3]);
				adjacency[tri[j]].push_back(tri[(j + 2) % 3]);
			}
		}

		for (auto& neighbors : adjacency) {
			std::unordered_set<int> unique(neighbors.begin(), neighbors.end());
			neighbors.assign(unique.begin(), unique.end());
		}
	}

	static void calculateMeanCurvature(
		const std::vector<Eigen::Vector3d>& vertices,
		const std::vector<std::vector<int>>& adjacency,
		const std::vector<double>& vertex_areas,
		std::vector<Eigen::Vector3d>& curvature_normals,
		std::vector<double>& mean_curvatures
	) {
		for (size_t i = 0; i < vertices.size(); ++i) {
			if (adjacency[i].empty()) continue;

			Eigen::Vector3d laplacian = Eigen::Vector3d::Zero();
			for (int j : adjacency[i]) {
				laplacian += vertices[j] - vertices[i];
			}
			laplacian /= adjacency[i].size();

			curvature_normals[i] = laplacian;
			mean_curvatures[i] = laplacian.norm() / vertex_areas[i];
		}
	}


	static Eigen::Vector3d calculateSurfaceTension(
		const Eigen::Vector3d& vertex,
		const std::vector<int>& neighbors,
		const std::vector<Eigen::Vector3d>& vertices,
		const Eigen::Vector3d& normal,
		double tension_weight
	) {
		Eigen::Vector3d centroid = Eigen::Vector3d::Zero();
		for (int neighbor : neighbors) {
			centroid += vertices[neighbor];
		}
		centroid /= neighbors.size();

		Eigen::Vector3d to_center = centroid - vertex;
		double normal_component = to_center.dot(normal);
		Eigen::Vector3d tangential = to_center - normal_component * normal;
		return (tension_weight * tangential + (1.0 - tension_weight) * normal_component * normal);

	}

	static std::vector<double> calculateVertexAreas(
		const std::vector<Eigen::Vector3d>& vertices,
		const std::vector<std::array<int, 3>>& triangles,
		const std::vector<std::vector<int>>& vertex_triangles
	) {
		std::vector<double> areas(vertices.size(), 0.0);

		for (size_t i = 0; i < vertices.size(); ++i) {
			double area = 0.0;
			for (int tri_idx : vertex_triangles[i]) {
				const auto& tri = triangles[tri_idx];
				Eigen::Vector3d v1 = vertices[tri[1]] - vertices[tri[0]];
				Eigen::Vector3d v2 = vertices[tri[2]] - vertices[tri[0]];
				area += v1.cross(v2).norm() / 6.0; // (1/3 of triangle area)
			}
			areas[i] = area;
		}

		return areas;
	}

	static double calculateMeshVolume(
		const std::vector<Eigen::Vector3d>& vertices,
		const std::vector<std::array<int, 3>>& triangles
	) {
		double volume = 0.0;
		for (const auto& triangle : triangles) {
			const auto& v1 = vertices[triangle[0]];
			const auto& v2 = vertices[triangle[1]];
			const auto& v3 = vertices[triangle[2]];
			volume += v1.dot(v2.cross(v3)) / 6.0;
		}
		return std::abs(volume);
	}

	static void decimateMesh(
		std::vector<Eigen::Vector3d>& vertices,
		std::vector<std::array<int, 3>>& triangles,
		std::vector<Eigen::Vector3d>& normals,
		double length_threshold
	) {
		std::priority_queue<std::pair<double, std::pair<int, int>>> edge_heap;

		for (const auto& triangle : triangles) {
			for (int i = 0; i < 3; ++i) {
				int v1 = triangle[i];
				int v2 = triangle[(i + 1) % 3];
				if (v1 < v2) {
					double length = (vertices[v1] - vertices[v2]).norm();
					if (length < length_threshold) {
						edge_heap.push({ -length, {v1, v2} });
					}
				}
			}
		}

		std::vector<bool> valid_vertex(vertices.size(), true);
		std::vector<bool> valid_triangle(triangles.size(), true);

		while (!edge_heap.empty()) {
			auto edge = edge_heap.top().second;
			edge_heap.pop();

			if (!valid_vertex[edge.first] || !valid_vertex[edge.second]) continue;


			vertices[edge.first] = (vertices[edge.first] + vertices[edge.second]) * 0.5;
			normals[edge.first] = (normals[edge.first] + normals[edge.second]).normalized();
			valid_vertex[edge.second] = false;


			for (size_t i = 0; i < triangles.size(); ++i) {
				if (!valid_triangle[i]) continue;

				auto& triangle = triangles[i];
				for (int j = 0; j < 3; ++j) {
					if (triangle[j] == edge.second) {
						triangle[j] = edge.first;
					}
				}

				if (triangle[0] == triangle[1] ||
					triangle[1] == triangle[2] ||
					triangle[2] == triangle[0]) {
					valid_triangle[i] = false;
				}
			}
		}


		compactMesh(vertices, triangles, normals, valid_vertex, valid_triangle);
	}

	static void compactMesh(
		std::vector<Eigen::Vector3d>& vertices,
		std::vector<std::array<int, 3>>& triangles,
		std::vector<Eigen::Vector3d>& normals,
		const std::vector<bool>& valid_vertex,
		const std::vector<bool>& valid_triangle
	) {
		std::vector<int> vertex_map(vertices.size(), -1);
		int new_vertex_count = 0;

		for (size_t i = 0; i < vertices.size(); ++i) {
			if (valid_vertex[i]) {
				vertex_map[i] = new_vertex_count++;
			}
		}

		std::vector<Eigen::Vector3d> new_vertices;
		std::vector<Eigen::Vector3d> new_normals;
		new_vertices.reserve(new_vertex_count);
		new_normals.reserve(new_vertex_count);

		for (size_t i = 0; i < vertices.size(); ++i) {
			if (valid_vertex[i]) {
				new_vertices.push_back(vertices[i]);
				new_normals.push_back(normals[i]);
			}
		}

		std::vector<std::array<int, 3>> new_triangles;
		for (size_t i = 0; i < triangles.size(); ++i) {
			if (valid_triangle[i]) {
				std::array<int, 3> new_triangle;
				for (int j = 0; j < 3; ++j) {
					new_triangle[j] = vertex_map[triangles[i][j]];
				}
				new_triangles.push_back(new_triangle);
			}
		}

		vertices = std::move(new_vertices);
		triangles = std::move(new_triangles);
		normals = std::move(new_normals);
	}
};


	void exportMesthToVTK(
		std::string filename, 
		const std::vector<Eigen::Vector3d>& positions, 
		//std::vector<double>& densities,
		//double c_,
		pressureSolver::NeighborhoodSearch& nsearch,
		double compact_support_, 
		double dist_,
		double fluid_particle_mass_,
		Eigen::Vector3d origin_,
		Eigen::AlignedBox3d viewbox
		) {
		//double c = c_;
		double compact_support = compact_support_;
		double dist = dist_;
		double fluid_particle_mass = fluid_particle_mass_;
		Eigen::Vector3d pos_origin = origin_;
		double padding = 2 * compact_support_;
		//Eigen::Vector3d grid_origin = origin_ - Eigen::Vector3d::Constant(padding);

		Eigen::AlignedBox3d aabb;
		//aabb in world coordinate
		for (auto& pos : positions)
		{
			aabb.extend(pos);
		}

		aabb.intersection(viewbox);

		aabb.max() += Eigen::Vector3d::Constant(padding);
		//lowest frontest leftest point with padding in world space
		//This point is supposed to have the index 0,0,0 in the grid
		aabb.min() -= Eigen::Vector3d::Constant(padding);

		//Offset determines the position of the grid index 0,0,0 in world space
		//Eigen::Vector3d grid_offset = aabb.min() - pos_origin;
		Eigen::Vector3d grid_origin = aabb.min();



		Eigen::Vector3d diagonal = aabb.diagonal();

		Eigen::Vector3i ijk_max = (diagonal / dist).array().ceil().cast<int>();

		uint m_max_i = ijk_max[0];
		uint m_max_j = ijk_max[1];
		uint m_max_k = ijk_max[2];

		//std::cout << "Compute phi" << std::endl;
		std::unordered_map<uint, Eigen::Vector4d> phis;

		Eigen::Vector3d dif = Eigen::Vector3d(compact_support, compact_support, compact_support);
		// sample phi
		uint max_phi = m_max_i * m_max_j * m_max_k;
		double reverse_mass = 1.0 / fluid_particle_mass;

		//#pragma omp parallel for
		for (int p = 0; p < positions.size(); p++)
		{
			Eigen::Vector3d pos = positions[p];
			double norm_density = pressureSolver::kernel::W(Eigen::Vector3d{0,0,0});
			nsearch.for_each_neighbor(0, 0, p,
				[&](const int n)
				{

					norm_density += pressureSolver::kernel::W(positions[p] - positions[n]);
				}
			);

			//double norm_density =  densities[p]*reverse_mass;
			double reverse_density = 1.0 / norm_density;

			double tmp = pressureSolver::kernel::W(Eigen::Vector3d(0.0, 0.0, 0.0));
			double tmp2 = pressureSolver::kernel::h;

			Eigen::Vector3d bottom = pos - dif; //world coordinate
			Eigen::Vector3d bottom_grid = bottom - grid_origin;
			Eigen::Vector3d top = pos + dif;
			//determine the index of the point in the grid, needing to turn world coordinate into grid coordinate
			Eigen::Vector3i ijk_low = ((bottom - aabb.min()) / dist).array().ceil().cast<int>().cwiseMax(0);
			Eigen::Vector3i ijk_high = ((top - aabb.min()) / dist).array().ceil().cast<int>().cwiseMax(0);
			for (uint i = ijk_low[0]; i <= ijk_high[0]; i++)
			{
				for (uint j = ijk_low[1]; j <= ijk_high[1]; j++)
				{
					for (uint k = ijk_low[2]; k <= ijk_high[2]; k++)
					{
						Eigen::Vector3d grid_pos = grid2coordinate(i, j, k, dist, grid_origin);
						Eigen::Vector3d dist = grid_pos - pos;
						if (dist.squaredNorm() < (compact_support * compact_support)) {
							double w = pressureSolver::kernel::W(dist);
							Eigen::Vector3d w_grad = pressureSolver::kernel::W_grad(dist);
							//Eigen::Vector4d packet = Eigen::Vector4d{ w_grad[0], w_grad[1], w_grad[2], reverse_density * w -0.55};
							Eigen::Vector4d c{ 0,0,0,-0.55 };
							//#pragma omp critical
							if (V(i, j, k) >= 0 && V(i, j, k) < max_phi)
							{
								if (phis.find(V(i, j, k)) == phis.end()) {
									phis[V(i, j, k)] = Eigen::Vector4d{ w_grad[0], w_grad[1], w_grad[2], reverse_density * w - 0.55};
								}
								else {
									phis[V(i, j, k)] += Eigen::Vector4d{ w_grad[0], w_grad[1], w_grad[2], reverse_density * w };
								}
							}
						}
					}
				}
			}
		}
		//std::cout << "identify voxels" << std::endl;
		uint ijk_off = m_max_i * m_max_j * m_max_k + 1;
		std::vector<uint> voxels_vec(8 * phis.size(), ijk_off);
		std::vector<std::pair<uint, Eigen::Vector4d>> phi_vec(phis.begin(), phis.end());
		Eigen::Vector3i x = { 1, 0, 0 };
		Eigen::Vector3i y = { 0, 1, 0 };
		Eigen::Vector3i z = { 0, 0, 1 };
		Eigen::Vector3i o0 = { -1,-1,-1 };
		Eigen::Vector3i o1 = { 0,-1,-1 };
		Eigen::Vector3i o2 = { 0,0,-1 };
		Eigen::Vector3i o3 = { -1,0,-1 };
		Eigen::Vector3i o4 = { -1,-1,0 };
		Eigen::Vector3i o5 = { 0,-1,0 };
		Eigen::Vector3i o6 = { 0,0,0 };
		Eigen::Vector3i o7 = { -1,0,0 };
//#pragma omp parallel for
		for (int i = 0; i < phi_vec.size(); i++) {
			auto point = phi_vec[i];
			uint vertex_id = point.first;
			double phi = point.second[3];
			//std::cout << "phi: " << phi << std::endl;
			if (phi < 0.0) {
				Eigen::Vector3i ijk = { (int)(vertex_id / (m_max_j * m_max_k)), (int)((vertex_id % (m_max_j * m_max_k)) / m_max_k), (int)(vertex_id % m_max_k) };
				//check suroundings
				uint i0 = V((ijk + o0)[0], (ijk + o0)[1], (ijk + o0)[2]);
				uint i1 = V((ijk + o1)[0], (ijk + o1)[1], (ijk + o1)[2]);
				uint i2 = V((ijk + o2)[0], (ijk + o2)[1], (ijk + o2)[2]);
				uint i3 = V((ijk + o3)[0], (ijk + o3)[1], (ijk + o3)[2]);
				uint i4 = V((ijk + o4)[0], (ijk + o4)[1], (ijk + o4)[2]);
				uint i5 = V((ijk + o5)[0], (ijk + o5)[1], (ijk + o5)[2]);
				uint i6 = V((ijk + o6)[0], (ijk + o6)[1], (ijk + o6)[2]);
				uint i7 = V((ijk + o7)[0], (ijk + o7)[1], (ijk + o7)[2]);
				//check six surounding vertices
				//since we only check for phi < 0, we test if neighboring phi is larger or equal to zero
				//int voxelgrid = (m_phi_grad.find(V((ijk + x)[0], (ijk + x)[1], (ijk + x)[2])) != m_phi_grad.end() && m_phi_grad.find(V((ijk + x)[0], (ijk + x)[1], (ijk + x)[2]))->second[3] >= 0) ? 0b01100110 : 0b00000000;
				//voxelgrid |= (m_phi_grad.find(V((ijk - x)[0], (ijk - x)[1], (ijk - x)[2])) != m_phi_grad.end() && m_phi_grad.find(V((ijk - x)[0], (ijk - x)[1], (ijk - x)[2]))->second[3] >= 0) ? 0b10011001 : 0b00000000;
				//voxelgrid |= (m_phi_grad.find(V((ijk + y)[0], (ijk + y)[1], (ijk + y)[2])) != m_phi_grad.end() && m_phi_grad.find(V((ijk + y)[0], (ijk + y)[1], (ijk + y)[2]))->second[3] >= 0) ? 0b00110011 : 0b00000000;
				//voxelgrid |= (m_phi_grad.find(V((ijk - y)[0], (ijk - y)[1], (ijk - y)[2])) != m_phi_grad.end() && m_phi_grad.find(V((ijk - y)[0], (ijk - y)[1], (ijk - y)[2]))->second[3] >= 0) ? 0b11001100 : 0b00000000;
				//voxelgrid |= (m_phi_grad.find(V((ijk + z)[0], (ijk + z)[1], (ijk + z)[2])) != m_phi_grad.end() && m_phi_grad.find(V((ijk + z)[0], (ijk + z)[1], (ijk + z)[2]))->second[3] >= 0) ? 0b00001111 : 0b00000000;
				//voxelgrid |= (m_phi_grad.find(V((ijk - z)[0], (ijk - z)[1], (ijk - z)[2])) != m_phi_grad.end() && m_phi_grad.find(V((ijk - z)[0], (ijk - z)[1], (ijk - z)[2]))->second[3] >= 0) ? 0b11110000 : 0b00000000;
				voxels_vec[8 * i] = i7;
				voxels_vec[8 * i + 1] = i6;
				voxels_vec[8 * i + 2] = i5;
				voxels_vec[8 * i + 3] = i4;
				voxels_vec[8 * i + 4] = i3;
				voxels_vec[8 * i + 5] = i2;
				voxels_vec[8 * i + 6] = i1;
				voxels_vec[8 * i + 7] = i0;
			}
		}
		std::unordered_set<uint> voxels(voxels_vec.begin(), voxels_vec.end());
		voxels.erase(voxels.find(ijk_off));
		//std::cout << voxels.size() << ";" << std::endl;
		//identify relevant edges, vertices and prebuild triangles
		std::vector<std::array<int, 3>> triangles;
		//idea to parallize, didnt work yet
		//std::vector<uint> voxel_vec(voxels.begin(), voxels.end());
		std::unordered_map<int, int> edge_map;
		std::vector<Eigen::Vector3d> vertices;
		std::vector<Eigen::Vector3d> normals;
		//also tried a parallel appraoch, but seemed to be slower, so its inlined
		//#pragma omp parallel for
		for (auto& voxel : voxels) {
			Eigen::Vector3i ijk = { (int)(voxel / (m_max_j * m_max_k)), (int)((voxel % (m_max_j * m_max_k)) / m_max_k), (int)(voxel % m_max_k) };
			int MARCHING_CUBES_TABLE_idx = 0;
			//determine marching cubes table index
//#pragma unroll 8
			for (int j = 0; j < 8; j++) {
				Eigen::Vector3i v = Eigen::Vector3i(CELL_VERTICES[j][0], CELL_VERTICES[j][1], CELL_VERTICES[j][2]);
				Eigen::Vector3i adjusted = ijk + v;
				auto tmp_v = V(adjusted[0], adjusted[1], adjusted[2]);
				auto tmp = phis.find(V(adjusted[0], adjusted[1], adjusted[2]));
				MARCHING_CUBES_TABLE_idx += (tmp != phis.end() && tmp->second[3] >= 0.0) ? 1 << j : 0;
			}
			auto tmp_cell = MARCHING_CUBES_TABLE[MARCHING_CUBES_TABLE_idx];
			//iterate over all triangles
			for (auto& t : tmp_cell) {
				//run only if triangle is relevant
				if (t[0] != -1) {
					std::array<int, 3> triangle;
//#pragma unroll 3
					//one triangle ha three vertices. Each vertex is defined by an edge
					//Determine two nodes that define edge and if the edge hasnt been found yet, add it to the edge map
					//Elsewise just add index to edge map
					for (int i = 0; i < 3; i++) {
						std::array<int, 3> c1 = CELL_VERTICES[CELL_EDGES[t[i]][0]];
						Eigen::Vector3i c1_v = Eigen::Vector3i(c1[0], c1[1], c1[2]);
						int v1 = V((ijk + c1_v)[0], (ijk + c1_v)[1], (ijk + c1_v)[2]);
						int edge = v1 * 3 + CELL_EDGES_DIRECTION[t[i]];
						if (edge_map.find(edge) == edge_map.end()) {
							{
								std::array<int, 3> c2 = CELL_VERTICES[CELL_EDGES[t[i]][1]];
								Eigen::Vector3i c2_v = Eigen::Vector3i(c2[0], c2[1], c2[2]);
								int v2 = V((ijk + c2_v)[0], (ijk + c2_v)[1], (ijk + c2_v)[2]);
								Eigen::Vector3d p1 = grid2coordinate((ijk + c1_v)[0], (ijk + c1_v)[1], (ijk + c1_v)[2], dist, grid_origin);
								Eigen::Vector3d p2 = grid2coordinate((ijk + c2_v)[0], (ijk + c2_v)[1], (ijk + c2_v)[2], dist, grid_origin);
								Eigen::Vector4d phi1 = phis.find(v1) == phis.end() ? Eigen::Vector4d{ 0.0, 0.0, 0.0, -0.55 } : phis[v1];
								Eigen::Vector4d phi2 = phis.find(v2) == phis.end() ? Eigen::Vector4d{ 0.0, 0.0, 0.0, -0.55 } : phis[v2];
								Eigen::Vector3d n1 = Eigen::Vector3d(phi1[0], phi1[1], phi1[2]);
								Eigen::Vector3d n2 = Eigen::Vector3d(phi2[0], phi2[1], phi2[2]);
								double alpha = (0.0 - phi1[3]) / (phi2[3] - phi1[3]);
								vertices.push_back((1 - alpha) * p1 + alpha * p2);
								normals.push_back((((1 - alpha) * n1.normalized() + alpha * n2.normalized())));
								edge_map[edge] = vertices.size() - 1;
							}
						}
						triangle[i] = edge_map[edge];
					}
					triangles.push_back(triangle);
				}
			}
		}
		MeshSmoother::SmoothingParams params;
		params.iterations = 15;
		params.max_curvature = 1.5;         // Lower = smoother
		params.min_feature_size = dist*4;   // Minimum size to preserve
		params.tension_weight = 1.0;        // Higher = rounder droplets
		params.volume_preserve = 0.7;       // Higher = better volume preservation
		params.edge_length_threshold = dist/3.0;
		
		//uncomment, if mesh smoothing should be activated (takes alot longer depending on settings)
		//MeshSmoother::processMesh(vertices, triangles, normals, params);
		//LaplacianSmoother::processMesh(vertices, triangles, normals);
		pressureSolver::write_tri_mesh_to_vtk(filename, vertices, triangles, normals);
	}


}