#include "runSim.h"

#include "../pressureSolver/io.h"
#include "../pressureSolver/kernel.h"
#include "../pressureSolver/wcsph.h"
#include "../pressureSolver/sampling.h"
#include <CompactNSearch/CompactNSearch> 
#include <string> 

using namespace pressureSolver;

void RunSim::run() {
	
	double particle_radius = 0.005;
	std::vector<Eigen::Vector3d> positions;	
	double dist = 2 * particle_radius;
	Eigen::Vector3d bottom{ -0.075 + 2 * particle_radius, -0.4 + 2 * particle_radius,-.5 + 2 * particle_radius };
	Eigen::Vector3d top{ 0.075 - 2 * particle_radius, -0.15 - 2 * particle_radius,0.0- 2 * particle_radius };
	sampling::fluid_box(positions, bottom, top, dist);
	double rho_0 = 1000.0;
	double volume = positions.size() * dist * dist * dist;
	
	wcsph::Simulation sim(positions, volume, 5000, 60, rho_0, particle_radius);
	sim.run();

}