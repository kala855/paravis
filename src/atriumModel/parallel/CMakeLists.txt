# Archivo para realizar la generación de las Makefiles de manera automática

cmake_minimum_required (VERSION 2.8)
project(cardiac_cell)
find_package(Armadillo REQUIRED)
find_package(ArrayFire REQUIRED)
include_directories(${ArrayFire_INCLUDE_DIRS})
include_directories(${ARMADILLO_INCLUDE_DIRS})
find_package(CUDA REQUIRED)
cuda_add_executable(cardiac_cell simulation2D-CN.cc )
LIST(APPEND CUDA_NVCC_FLAGS --compiler-options -fno-strict-aliasing -lineinfo -use_fast_math -Xptxas -dlcm=cg)
LIST(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_20,code=sm_20)
LIST(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_30,code=sm_30)
LIST(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_35,code=sm_35)

target_link_libraries(cardiac_cell cell armadillo ${ArrayFire_CUDA_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT} ${CUDA_LIBRARIES})
add_subdirectory(utilities)

