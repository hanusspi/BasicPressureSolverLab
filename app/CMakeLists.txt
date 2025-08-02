cmake_minimum_required(VERSION 3.5)
set(CMAKE_CXX_STANDARD 17)

# Local
set(SOURCE_FILES
    exerciseFactory.h
    exercise.h
    default.h
    default.cpp
    runSim.cpp
    runSim.h
    main.cpp
    runSimulator.h
    runSimulator.cpp
    emitterRender.h
    emitterRender.cpp
    singleDamBreak.h
    singleDamBreak.cpp
    DoubleDamBreak.h
    DoubleDamBreak.cpp
    )
    
add_executable(pressureSolver_app ${SOURCE_FILES})
add_dependencies(pressureSolver_app pressureSolver)
target_link_libraries(pressureSolver_app pressureSolver)
