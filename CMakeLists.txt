cmake_minimum_required(VERSION 3.16)
project(finite-element-solver)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find required packages
find_package(Qt5 COMPONENTS Core Widgets PrintSupport REQUIRED)
find_package(CGAL REQUIRED)
find_package(Armadillo REQUIRED)

# Add SuperLU (needed for Armadillo's sparse solver)
find_library(SUPERLU_LIBRARY 
             NAMES superlu
             PATHS /opt/homebrew/lib
                   /opt/homebrew/Cellar/superlu/7.0.0/lib
                   /usr/local/lib)
find_path(SUPERLU_INCLUDE_DIR 
          NAMES slu_ddefs.h supermatrix.h
          PATHS /opt/homebrew/Cellar/superlu/7.0.0/include
                /opt/homebrew/include
                /usr/local/include)

# Add muparser (you may need to install this separately)
find_library(MUPARSER_LIBRARY muparser)
find_path(MUPARSER_INCLUDE_DIR muParser.h)

# Define that we want to use SuperLU with Armadillo
add_compile_definitions(ARMA_USE_SUPERLU)
add_compile_definitions(ARMA_SUPERLU_INCLUDE_DIR=${SUPERLU_INCLUDE_DIR})

# Include directories
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/src/armadillo_bits
    ${MUPARSER_INCLUDE_DIR}
    ${ARMADILLO_INCLUDE_DIRS}
    ${CGAL_INCLUDE_DIRS}
    ${SUPERLU_INCLUDE_DIR}
)

# Set Qt automoc
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)

# Source files
set(SOURCES
    src/main.cpp
    src/simulation.cpp
    src/qcustomplot.cpp
    src/superlu_wrapper.c
)

# Header files
set(HEADERS
    src/Simulation.h
    src/qcustomplot.h
    src/armadillo_bits/superlu_wrapper.hpp
    src/superlu_wrapper.h
)

# Create executable
add_executable(finite-element-solver ${SOURCES} ${HEADERS})

# Link libraries
target_link_libraries(finite-element-solver
    Qt5::Core
    Qt5::Widgets
    Qt5::PrintSupport
    ${MUPARSER_LIBRARY}
    ${ARMADILLO_LIBRARIES}
    ${CGAL_LIBRARIES}
    ${SUPERLU_LIBRARY}
)

# Copy data files to build directory
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/data.csv DESTINATION ${CMAKE_CURRENT_BINARY_DIR}) 