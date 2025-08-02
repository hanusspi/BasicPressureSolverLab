#include "doubleDamBreak.h"
#include "../pressureSolver/simulator.h"
#include <string> // std::string
#include <vector> // std::vector
#include <Eigen/Dense> // Eigen::Vector3d
#include <filesystem> 

void DoubleDamBreak::run() {
	double particle_radius = 0.005;
	Eigen::AlignedBox3d box1(Eigen::Vector3d(-0.075 + 2 * particle_radius, -0.1 + 2 * particle_radius, -0.5 + 2 * particle_radius), Eigen::Vector3d(0.075 - 2 * particle_radius, 0.0, -0.4));
	Eigen::AlignedBox3d box2(Eigen::Vector3d(-0.055 +2*particle_radius, -0.055 + 2 * particle_radius, -0.5 + 2 * particle_radius), Eigen::Vector3d(0.055, 0.055, -0.39));
	std::vector<Eigen::AlignedBox3d> boxes = { box1, box2 };
	std::string boundary_path = "../res/box.obj";

	double rho_0 = 1000.0;
	int freq = 5000;
	int fps = 60;
	int iterations = 4 * freq; 
	//choose params depending on solver
	std::vector<double> params = { 1000 }; //wcsph
	//std::vector<double> params = { 5 }; //pbf
	std::string export_path = "../res/DoubleDamBreak/";
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


	double beta = 0.5;
	double gamma = 0.5;
	pressureSolver::Simulator sim;
	sim.setParams(particle_radius, export_path, rho_0, freq, fps, params, beta, gamma);
	sim.setScene(boxes, {}, boundary_path);

	sim.run(iterations);

}
