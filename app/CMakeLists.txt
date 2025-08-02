cmake_minimum_required(VERSION 3.14)  # FetchContent requires 3.14+
set(CMAKE_CXX_STANDARD 17)

# Fetch nativefiledialog-extended automatically
include(FetchContent)

FetchContent_Declare(
    nativefiledialog-extended
    GIT_REPOSITORY https://github.com/btzy/nativefiledialog-extended.git
    GIT_TAG v1.1.1  # Use latest stable version
)

FetchContent_MakeAvailable(nativefiledialog-extended)

# Local source files
set(SOURCE_FILES
    config.h
    config.cpp
    main.cpp
)
    
add_executable(pressureSolver_app ${SOURCE_FILES})
add_dependencies(pressureSolver_app pressureSolver)
target_link_libraries(pressureSolver_app pressureSolver nfd)