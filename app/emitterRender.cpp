#include "emitterRender.h"
#include "../pressureSolver/simulator.h"
#include <string> // std::string
#include <vector> // std::vector
#include <Eigen/Dense> // Eigen::Vector3d
#include <filesystem> 

void EmitterRender::run() {
	double particle_radius = 0.005;
	std::string boundary_path = "../res/boxEx2_6b.obj";
	Eigen::AlignedBox3d box1(Eigen::Vector3d(-1 + 2 * particle_radius, -1 + 2 * particle_radius, -0.5 + 2 * particle_radius), Eigen::Vector3d(-0.5 - 2.0 * particle_radius, -0.5 - 2.0 * particle_radius, -0.3));
	Eigen::AlignedBox3d box(Eigen::Vector3d(-0.075 + 2 * particle_radius, -0.4 + 2 * particle_radius, -0.5 + 2 * particle_radius), Eigen::Vector3d(0.075 - 2 * particle_radius, -0.15, 0.0));

	double rho_0 = 1000.0;
	int freq = 1000;
	int fps = 60;
	int iterations = 3 * freq; 
	std::vector<double> params = { 5 };
	std::string export_path = "../res/fountainMeshSmoothedAndDecimated/";
	try {
		if (std::filesystem::create_directories(export_path)) {
			std::cout << "Directory created successfully." << std::endl;
		}
		else {
			std::cout << "Directory already exists or failed to create." << std::endl;
		}
	}
	catch (const std::filesystem::filesystem_error& e) {
		std::cerr << "Filesystem error: " << e.what() << std::endl;
	}


	pressureSolver::scene::EmitterSetting emitter2;
	emitter2.radius = 0.06;
	emitter2.scripted_time = 15;
	emitter2.velocity = 4.0;
	emitter2.rotation = Eigen::Vector3d(0.0, 0.0, 1.0);
	emitter2.positon = Eigen::Vector3d(0.0, 0.4, -0.8);
	emitter2.alpha = 1.5;
	emitter2.boundaryThickness = particle_radius;
	emitter2.cylinderLength = 7 * particle_radius;
	emitter2.activation_time = 0.0;
	emitter2.deactivation_time = 2.0;


	double beta = 0.1;
	double gamma = 0.1;
	pressureSolver::Simulator sim;
	sim.setParams(particle_radius, export_path, rho_0, freq, fps, params, beta, gamma);
	sim.setScene({box}, {}, boundary_path);

	sim.run(iterations);

}
