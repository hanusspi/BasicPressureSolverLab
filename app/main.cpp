#include "config.h"
#include <nfd.h>
#include <iostream>
#include <filesystem>

namespace fs = std::filesystem;

int main()
{
    if (NFD_Init() != NFD_OKAY) {
        std::cerr << "Error: Could not initialize file dialog." << std::endl;
        return 1;
    }

    fs::path relativePath = "../configs";
    fs::path absolutePath = fs::absolute(relativePath);

    // Create file filter
    nfdfilteritem_t filterItem[1] = { { "JSON", "json" } };
    // Open file dialog
    nfdchar_t* outPath;
    nfdresult_t result = NFD_OpenDialog(&outPath, filterItem, 1, absolutePath.string().c_str());
    if (result == NFD_OKAY) {
        std::cout << "Selected config: " << outPath << std::endl;
        try {
            // Read and run simulation
            SimConfig config = ConfigReader::readConfig(outPath);
            std::cout << "Solver: " << (config.solver_type == pressureSolver::SolverType::WCSPH ? "WCSPH" : "PBF") << std::endl;
            std::cout << "Starting simulation..." << std::endl;
            ConfigReader::runSimulation(config);
            std::cout << "Simulation completed!" << std::endl;
        }
        catch (const std::exception& e) {
            std::cerr << "Error: " << e.what() << std::endl;
        }
        // Free the path
        NFD_FreePath(outPath);
    }
    else if (result == NFD_CANCEL) {
        std::cout << "User cancelled." << std::endl;
    }
    else {
        std::cerr << "Error: " << NFD_GetError() << std::endl;
    }
    NFD_Quit();
    return 0;
}


//"../configs/damBreak.json"