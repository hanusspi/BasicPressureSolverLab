#pragma once
#include <string>
#include <vector>
#include <Eigen/Dense>
#include "../pressureSolver/Solver.h"
#include "../pressureSolver/scene.h"

struct EmitterConfig {
    double velocity;
    Eigen::Vector3d rotation;
    Eigen::Vector3d position;
    double radius;
    double scripted_time;
    double alpha;
    double boundary_thickness;
    double cylinder_length;
    double activation_time;
    double deactivation_time;
};

struct SimConfig {
    double particle_radius;
    std::string boundary_path;
    std::vector<Eigen::AlignedBox3d> fluid_boxes;
    double rho_0;
    int freq;
    int fps;
    int iterations;
    std::string export_path;
    pressureSolver::SolverType solver_type;
    double solver_param;
    double beta;
    double gamma;
    std::vector<EmitterConfig> emitters;
};

class ConfigReader {
public:
    static SimConfig readConfig(const std::string& filepath);
    static void runSimulation(const SimConfig& config);
};