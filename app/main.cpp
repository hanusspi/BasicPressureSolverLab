#include "exerciseFactory.h"

/*
* Selectable Scenes:
double_dam_break
single_dam_break

*/

//Select solver in Solver.h
//need to set parameters depending on solver in the chosen scene file
//to get logging data, logging needs to be activated in simulator.h

int main()
{
	std::string exercise = "single_dam_break";
	try {
		auto task = ExerciseFactory::createTask(exercise);
		task->run();  // Execute the task
	}
	catch (const std::invalid_argument& e) {
		std::cerr << "Error: " << e.what() << std::endl;
	}
	
	return 0;
}