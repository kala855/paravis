# Archivo para realizar la generación de las Makefiles de manera automática

cmake_minimum_required (VERSION 2.8)
project(cardiac_cell)
find_package(Armadillo REQUIRED)
find_package(ArrayFire REQUIRED)
include_directories(${ArrayFire_INCLUDE_DIRS})
include_directories(${ARMADILLO_INCLUDE_DIRS})
find_package(CUDA REQUIRED)
SET(CUDA_SEPARABLE_COMPILATION ON)
#LIST(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_20,code=sm_20)
#LIST(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_30,code=sm_30)
#LIST(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_50,code=sm_50)
cuda_add_executable(cardiac_cell cell.cuh cell.cu simulation2D-CN.cu)
target_link_libraries(cardiac_cell armadillo ${ArrayFire_CUDA_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT} ${CUDA_LIBRARIES})
