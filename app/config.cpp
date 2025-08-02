#include "config.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <filesystem>
#include "../pressureSolver/simulator.h"

// Simple JSON parser for our specific use case
// For production, consider using nlohmann/json
class SimpleJsonParser {
public:
    static SimConfig parseConfig(const std::string& json_content) {
        SimConfig config;

        // Parse particle_radius
        config.particle_radius = parseDouble(json_content, "particle_radius");

        // Parse boundary_path
        config.boundary_path = parseString(json_content, "boundary_path");

        // Parse rho_0
        config.rho_0 = parseDouble(json_content, "rho_0");

        // Parse freq, fps, iterations
        config.freq = parseInt(json_content, "freq");
        config.fps = parseInt(json_content, "fps");
        config.iterations = parseInt(json_content, "iterations");

        // Parse export_path
        config.export_path = parseString(json_content, "export_path");

        // Parse solver type
        std::string solver_str = parseString(json_content, "solver_type");
        config.solver_type = (solver_str == "PBF") ? pressureSolver::SolverType::PBF : pressureSolver::SolverType::WCSPH;

        // Parse solver param, beta, gamma
        config.solver_param = parseDouble(json_content, "solver_param");
        config.beta = parseDouble(json_content, "beta");
        config.gamma = parseDouble(json_content, "gamma");

        // Parse fluid boxes
        config.fluid_boxes = parseFluidBoxes(json_content);

        // Parse emitters
        config.emitters = parseEmitters(json_content);

        return config;
    }

private:
    static double parseDouble(const std::string& json, const std::string& key) {
        auto pos = json.find("\"" + key + "\"");
        if (pos == std::string::npos) return 0.0;

        pos = json.find(":", pos);
        pos = json.find_first_of("0123456789.-", pos);

        std::string value;
        while (pos < json.length() && (std::isdigit(json[pos]) || json[pos] == '.' || json[pos] == '-')) {
            value += json[pos++];
        }
        return std::stod(value);
    }

    static int parseInt(const std::string& json, const std::string& key) {
        return static_cast<int>(parseDouble(json, key));
    }

    static std::string parseString(const std::string& json, const std::string& key) {
        auto pos = json.find("\"" + key + "\"");
        if (pos == std::string::npos) return "";

        pos = json.find(":", pos);
        pos = json.find("\"", pos);
        if (pos == std::string::npos) return "";

        auto start = ++pos;
        auto end = json.find("\"", start);
        return json.substr(start, end - start);
    }

    static Eigen::Vector3d parseVector3d(const std::string& json, size_t& pos) {
        // Find opening bracket
        pos = json.find("[", pos);
        double x = 0, y = 0, z = 0;

        // Parse x
        pos = json.find_first_of("0123456789.-", pos);
        std::string num;
        while (pos < json.length() && (std::isdigit(json[pos]) || json[pos] == '.' || json[pos] == '-')) {
            num += json[pos++];
        }
        x = std::stod(num);

        // Parse y
        pos = json.find(",", pos);
        pos = json.find_first_of("0123456789.-", pos);
        num.clear();
        while (pos < json.length() && (std::isdigit(json[pos]) || json[pos] == '.' || json[pos] == '-')) {
            num += json[pos++];
        }
        y = std::stod(num);

        // Parse z
        pos = json.find(",", pos);
        pos = json.find_first_of("0123456789.-", pos);
        num.clear();
        while (pos < json.length() && (std::isdigit(json[pos]) || json[pos] == '.' || json[pos] == '-')) {
            num += json[pos++];
        }
        z = std::stod(num);

        pos = json.find("]", pos);
        return Eigen::Vector3d(x, y, z);
    }

    static std::vector<Eigen::AlignedBox3d> parseFluidBoxes(const std::string& json) {
        std::vector<Eigen::AlignedBox3d> boxes;
        auto pos = json.find("\"fluid_boxes\"");
        if (pos == std::string::npos) return boxes;

        pos = json.find("[", pos);
        while (true) {
            pos = json.find("{", pos);
            if (pos == std::string::npos) break;

            auto min_pos = json.find("\"min\"", pos);
            auto max_pos = json.find("\"max\"", pos);

            if (min_pos != std::string::npos && max_pos != std::string::npos) {
                size_t temp_pos = min_pos;
                Eigen::Vector3d min_point = parseVector3d(json, temp_pos);
                temp_pos = max_pos;
                Eigen::Vector3d max_point = parseVector3d(json, temp_pos);
                boxes.emplace_back(min_point, max_point);
            }

            pos = json.find("}", pos) + 1;
            if (json.find(",", pos) < json.find("]", pos)) {
                pos = json.find(",", pos);
            }
            else {
                break;
            }
        }
        return boxes;
    }

    static std::vector<EmitterConfig> parseEmitters(const std::string& json) {
        std::vector<EmitterConfig> emitters;
        auto pos = json.find("\"emitters\"");
        if (pos == std::string::npos) return emitters;

        pos = json.find("[", pos);
        while (true) {
            pos = json.find("{", pos);
            if (pos == std::string::npos) break;

            EmitterConfig emitter;
            auto end_pos = json.find("}", pos);
            std::string emitter_json = json.substr(pos, end_pos - pos);

            emitter.velocity = parseDouble(emitter_json, "velocity");
            emitter.radius = parseDouble(emitter_json, "radius");
            emitter.scripted_time = parseDouble(emitter_json, "scripted_time");
            emitter.alpha = parseDouble(emitter_json, "alpha");
            emitter.boundary_thickness = parseDouble(emitter_json, "boundary_thickness");
            emitter.cylinder_length = parseDouble(emitter_json, "cylinder_length");
            emitter.activation_time = parseDouble(emitter_json, "activation_time");
            emitter.deactivation_time = parseDouble(emitter_json, "deactivation_time");

            // Parse rotation and position vectors
            size_t rot_pos = 0, pos_pos = 0;
            emitter.rotation = parseVector3d(emitter_json, rot_pos);
            emitter.position = parseVector3d(emitter_json, pos_pos);

            emitters.push_back(emitter);

            pos = end_pos + 1;
            if (json.find(",", pos) < json.find("]", pos)) {
                pos = json.find(",", pos);
            }
            else {
                break;
            }
        }
        return emitters;
    }
};

SimConfig ConfigReader::readConfig(const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open config file: " + filepath);
    }

    std::ostringstream buffer;
    buffer << file.rdbuf();
    std::string json_content = buffer.str();

    return SimpleJsonParser::parseConfig(json_content);
}

void ConfigReader::runSimulation(const SimConfig& config) {
    // Create export directory
    try {
        if (std::filesystem::create_directories(config.export_path)) {
            std::cout << "Created export directory: " << config.export_path << std::endl;
        }
    }
    catch (const std::filesystem::filesystem_error& e) {
        std::cerr << "Filesystem error: " << e.what() << std::endl;
    }

    // Apply particle radius offset to fluid boxes
    std::vector<Eigen::AlignedBox3d> boxes;
    double offset = 2 * config.particle_radius;
    for (const auto& box : config.fluid_boxes) {
        Eigen::Vector3d min_offset = box.min() + Eigen::Vector3d::Constant(offset);
        Eigen::Vector3d max_offset = box.max() - Eigen::Vector3d::Constant(offset);
        boxes.emplace_back(min_offset, max_offset);
    }

    // Convert EmitterConfig to EmitterSetting
    std::vector<pressureSolver::scene::EmitterSetting> emitters;
    for (const auto& emitter : config.emitters) {
        pressureSolver::scene::EmitterSetting setting;
        setting.velocity = emitter.velocity;
        setting.rotation = emitter.rotation;
        setting.positon = emitter.position; // Note: keeping the typo from original code
        setting.radius = emitter.radius;
        setting.scripted_time = emitter.scripted_time;
        setting.alpha = emitter.alpha;
        setting.boundaryThickness = emitter.boundary_thickness;
        setting.cylinderLength = emitter.cylinder_length;
        setting.activation_time = emitter.activation_time;
        setting.deactivation_time = emitter.deactivation_time;
        emitters.push_back(setting);
    }

    // Create and run simulation
    pressureSolver::Simulator sim;
    std::vector<double> params = { config.solver_param };

    sim.setParams(config.solver_type, config.particle_radius, config.export_path,
        config.rho_0, config.freq, config.fps, params, config.beta, config.gamma);
    sim.setScene(boxes, emitters, config.boundary_path);
    sim.run(config.iterations);
}