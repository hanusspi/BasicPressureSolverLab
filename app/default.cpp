#include "default.h"

#include "../pressureSolver/io.h"
#include "../pressureSolver/sampling.h"

void DefaultMain::run() {
	std::cout << "Welcome to the learnSPH framework!!" << std::endl;
	std::cout << "Generating a sample scene...";

	// Load a obj surface mesh
	const std::vector<pressureSolver::TriMesh> meshes = pressureSolver::read_tri_meshes_from_obj("../res/box.obj");
	const pressureSolver::TriMesh& box = meshes[0];

	// Sample the mesh with particles
	const double sampling_distance = 0.05;
	std::vector<Eigen::Vector3d> particles;
	pressureSolver::sampling::triangle_mesh(particles, box.vertices, box.triangles, sampling_distance);

	// Initialize data vectors
	std::vector<double> particles_scalar_data(particles.size());
	std::vector<Eigen::Vector3d> particles_vector_data(particles.size());

	// Scalar data will be the particle_id
	for (int i = 0; i < (int)particles.size(); i++) {
		particles_scalar_data[i] = i;
	}

	// Simulation loop
	for (int time_step = 0; time_step < 100; time_step++) {

		for (int particle_i = 0; particle_i < (int)particles.size(); particle_i++) {
			// Move particles a bit down in the Z direction
			particles[particle_i][2] -= 0.025;

			// Clamp the Z coord to the floor
			particles[particle_i][2] = std::max(particles[particle_i][2], -1.0);

			// Vector data is going to be the position
			particles_vector_data[particle_i] = particles[particle_i] - Eigen::Vector3d(-1, -1, -1);
		}

		// Save output
		const std::string filename = "../res/example_" + std::to_string(time_step) + ".vtk";
		pressureSolver::write_particles_to_vtk(filename, particles, particles_scalar_data, particles_vector_data);
	}

	std::cout << "completed!" << std::endl;
	std::cout << "The scene files have been saved in the folder `<build_folder>/res`. \nYou can visualize them with Blender by using the Blender Sequence Loaded addon." << std::endl;

}