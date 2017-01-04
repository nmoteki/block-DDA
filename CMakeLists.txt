cmake_minimum_required (VERSION 3.7.0)
set(CMAKE_CXX_COMPILER "g++-6")
set(CMAKE_CXX_FLAGS "-O3 -fopenmp -DNDEBUG -DEIGEN_NO_DEBUG -DEIGEN_FFTW_DEFAULT")
include_directories(/Users/moteki/eigen_3_2_10)
include_directories(/usr/local/include)
link_directories(/usr/local/lib)
link_libraries(fftw3)
link_libraries(fftw3_omp)
link_libraries(fftw3_threads)


project(block_DDA_main)

add_subdirectory(gmie)
add_subdirectory(mvp_fft)
add_subdirectory(bl_krylov_mvp_fft)

add_executable(block_DDA_main block_DDA_main.cpp)
add_library(Material Material.cpp)
add_library(Incident_Wave Incident_Wave.cpp)
add_library(Target Target.cpp)
add_library(Target_manager Target_manager.cpp)
add_library(Solver Solver.cpp)
add_library(Result_manager Result_manager.cpp)

target_link_libraries(block_DDA_main
Material
Incident_Wave
Target
Target_manager
Solver
Result_manager
gmie
gmie_coeff
bl_cocg_rq_mvp_fft
bl_cocg_rq_jacobi_mvp_fft
bl_bicgstab_mvp_fft
bl_bicgstab_jacobi_mvp_fft
bl_bicggr_mvp_fft
bl_bicggr_jacobi_mvp_fft
qr_reduced
application_function
BT_fft
BT_pad
BT_reconstruct
MBT_fft_init
MBT_fft_mvp)
