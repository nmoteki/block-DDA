cmake_minimum_required (VERSION 3.7.0)
set(CMAKE_CXX_COMPILER "g++-6")
set(CMAKE_CXX_FLAGS "-fopenmp -O3 -DNDEBUG")

add_library(gmie gmie.cpp)
add_library(gmie_coeff  gmie_coeff.cpp)
