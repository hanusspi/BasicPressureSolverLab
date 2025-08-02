# External dependencies
## OpenMP
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

## Eigen
set( EIGEN3_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/extern/Eigen )
if( NOT EIGEN3_INCLUDE_DIR )
    message( FATAL_ERROR "Please point the environment variable EIGEN3_INCLUDE_DIR to the include directory of your Eigen3 installation.")
endif()

set(SOURCE_FILES
	NeighborhoodSearch.h
	Solver.h
	pbf.h
	pbf.cpp
	util.h
	util.cpp
	simulator.h
	simulator.cpp
	io.h
	io.cpp
	sampling.h
	sampling.cpp
	kernel.cpp
	kernel.h
	wcsph.h
	wcsph.cpp
	wcsphSolver.h
	wcsphSolver.cpp
	marching_cubes_lut.h
	definitions.h
	scene.cpp
	scene.h
)
    
# Create a library version so it can be linked with other modules
add_library(pressureSolver ${SOURCE_FILES})
target_link_libraries(pressureSolver CompactNSearch)
target_include_directories(pressureSolver PUBLIC ${EIGEN3_INCLUDE_DIR} )
target_include_directories(pressureSolver PUBLIC ${PROJECT_SOURCE_DIR}/extern ) # So you can include header files under pressureSolver/extern
