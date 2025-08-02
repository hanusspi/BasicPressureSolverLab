# General

Project is a result of a Lab of the RWTH. 

Goal of the project is the implementation of a lagrangian fluid simuatlor using wcsph and pbf.

# Setup

Build using CMake. No further external depencies required. 

Boundary Geometry is stores in the res folder. 

# Usage

Simply run the pressure Solver App. This leads you to the dialoge where the config file can be selected.

The options for the config files are currently a bit limited. All possibilities can be seen in the existing example files.

It is possible to add multiple fluid boxes or emitters. Currently only one boundary object is supported though. 

# Visualizing the results

The simulator does not offer any direct visualization as of now. 

The results will be exported as .vtk files, that can be visualized using paraview







