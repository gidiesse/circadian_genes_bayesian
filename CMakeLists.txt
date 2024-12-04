cmake_minimum_required(VERSION 3.26)
project(Bayesian_Statistics_project)

set(CMAKE_CXX_STANDARD 14)

# Set the path to Armadillo headers + Link to Armadillo library
include_directories(C:/msys64/mingw64/include)
link_directories(C:/msys64/mingw64/lib)

add_executable(Bayesian_Statistics_project MCMC.cpp)

target_link_libraries(Bayesian_Statistics_project armadillo)
