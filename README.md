# General

Project is a result of a Lab of the RWTH. 

Goal of the project is the implementation of a lagrangian fluid simuatlor using wcsph and pbf.

# Setup

Build using CMake. No further external depencies required. 

Boundary Geometry is stores in the res folder. 

# Usage

In main select, which scene is supposed to be run.

In Solver.h WCSPH or PBF can be selected. Dependent on the chosen solver the input parameters for the scene need to be adapted.

I am aware that this is not very userfriendly and may update in the future

# Visualizing the results

The simulator does not offer any direct visualization as of now. 

The results will be exported as .vtk files, that can be visualized using paraview





