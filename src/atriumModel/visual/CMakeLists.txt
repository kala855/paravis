cmake_minimum_required(VERSION 2.8.8)
project(visual)


find_package(Armadillo REQUIRED)
find_package(ArrayFire REQUIRED)

include_directories(${ArrayFire_INCLUDE_DIRS})
include_directories(${ARMADILLO_INCLUDE_DIRS})

find_package(CUDA REQUIRED)
set(CUDA_SEPARABLE_COMPILATION ON)

set(USE_CATALYST ON CACHE BOOL "Link the simulator with Catalyst")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "/home/john/Downloads/paraviewBuild")
if(USE_CATALYST)
    find_package(ParaView 5.0 REQUIRED COMPONENTS vtkPVPythonCatalyst)
    include("${PARAVIEW_USE_FILE}")
    set(paravisAdaptor_SRCS paravisAdaptor.cxx)
    add_library(paravisAdaptor ${paravisAdaptor_SRCS})
    target_link_libraries(paravisAdaptor vtkPVPythonCatalyst vtkParallelMPI)
    add_definitions("-DUSE_CATALYST")
else()
    find_package(MPI REQUIRED)
    include_directories(${MPI_CXX_INCLUDE_PATH})
endif()
add_executable(visual cell.cuh cell.cu simulation2D-CN.cu)

if(USE_CATALYST)
    target_link_libraries(visual LINK_PRIVATE paravisAdaptor armadillo ${ArrayFire_CUDA_LIBRARIES}
        ${CMAKE_THREAD_LIBS_INIT} ${CUDA_LIBRARIES})
    include(vtkModuleMacros)
    include(vtkMPI)
    vtk_mpi_link(visual)
else()
    target_link_libraries(visual LINK_PRIVATE ${MPI_LIBRARIES})
endif()

