cmake_minimum_required(VERSION 3.26)
project(circadian_genes)

set(CMAKE_CXX_STANDARD 17)

find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})
add_executable(circadian_genes test.cpp)


#add_executable(circadian_genes
#        ../../circadian_genes/circadian_genes_bayesian/MCMC.cpp)

target_link_libraries(circadian_genes ${ARMADILLO_LIBRARIES})
