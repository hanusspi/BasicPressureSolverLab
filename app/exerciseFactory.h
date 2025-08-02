#pragma once
#include <memory>  // For std::unique_ptr
#include <stdexcept>  // For std::invalid_argument
#include "default.h"
#include "exercise.h"
#include "runSim.h"
#include "runSimulator.h"
#include "EmitterRender.h"
#include "doubleDamBreak.h"
#include "singleDamBreak.h"


// Factory to create Task objects
class ExerciseFactory {
public:
    static std::unique_ptr<Exercise> createTask(const std::string& type) {
        if (type == "default_main") {
            return std::make_unique<DefaultMain>();
        }
        else if (type == "run_sim") {
			return std::make_unique<RunSim>();
        }
		else if (type == "simulator") {
			return std::make_unique<RunSimulator>();
		}
        else if (type == "emitter_render") {
            return std::make_unique<EmitterRender>();
        }
		else if (type == "double_dam_break") {
			return std::make_unique<DoubleDamBreak>();
		}
		else if (type == "single_dam_break") {
			return std::make_unique<SingleDamBreak>();
		}
        else {
            throw std::invalid_argument("Unknown task type: " + type);
        }

    }
};
