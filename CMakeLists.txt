cmake_minimum_required(VERSION 3.24)
project(Assignment_3)

set(CMAKE_CXX_STANDARD 17)

#add_executable(Assignment_3 main.cpp complex_vector_utils.cpp _defs.hpp forces_and_integrators.cpp)
add_executable(Assignment_3 part_2/matrix.hpp part_2/matrix.cpp part_2/LP_solvers.cpp part_2/vector_utils.cpp part_2/main.cpp)