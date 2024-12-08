cmake_minimum_required(VERSION 3.26)
project(Bayesian_Statistics_project)

set(CMAKE_CXX_STANDARD 14)

include_directories(C:/msys64/mingw64/include)
link_directories(C:/msys64/mingw64/lib)

find_package(Armadillo REQUIRED)

add_executable(Bayesian_Statistics_project MCMC.cpp)

target_link_libraries(Bayesian_Statistics_project armadillo C:/msys64/mingw64/lib/libopenblas.a C:/msys64/mingw64/lib/liblapack.a)
