#include "runSimulator.h"
#include "../pressureSolver/simulator.h"
#include <string> // std::string
#include <vector> // std::vector
#include <Eigen/Dense> // Eigen::Vector3d

void RunSimulator::run() {


	
	double particle_radius = 0.005;
	std::string boundary_path = "../res/boxEx2_6b.obj";
	Eigen::AlignedBox3d box(Eigen::Vector3d(-0.075 + 2 * particle_radius, -0.4 + 2 * particle_radius, -0.5 + 2 * particle_radius), Eigen::Vector3d(0.075- 2 * particle_radius, -0.15, 0.0));
	double rho_0 = 1000.0;
	int freq = 2500;
	int fps = 120;
	std::vector<double> params = { 10 };
	int iterations = 1*freq;
	std::string export_path = "../res/c/c/";
	pressureSolver::Simulator sim;
}