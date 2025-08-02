#include "singleDamBreak.h"
#include "../pressureSolver/simulator.h"
#include <string> // std::string
#include <vector> // std::vector
#include <Eigen/Dense> // Eigen::Vector3d
#include <filesystem> 

void SingleDamBreak::run() {
	double particle_radius = 0.005;
	std::string boundary_path = "../res/boxEx2_6b.obj";
	Eigen::AlignedBox3d box1(Eigen::Vector3d(-0.075 + 2 * particle_radius, -0.4 + 2 * particle_radius, -0.5 + 2 * particle_radius), Eigen::Vector3d(0.075 - 2 * particle_radius, -0.15, 0.0));
	std::vector<Eigen::AlignedBox3d> boxes = { box1 };
	double rho_0 = 1000.0;
	int freq = 3000;
	int fps = 60;
	int iterations = 2 * freq; 
	//choose params depending on solver
	std::vector<double> params = { 1000 }; //wcsph
	//std::vector<double> params = { 5 }; //pbf
	std::string export_path = "../res/SingleDamBreak/";
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

	double beta = 0.001;
	double gamma = 0.001;
	pressureSolver::Simulator sim;
	sim.setParams(particle_radius, export_path, rho_0, freq, fps, params, beta, gamma);
	sim.setScene(boxes, {}, boundary_path);

	sim.run(iterations);

}
