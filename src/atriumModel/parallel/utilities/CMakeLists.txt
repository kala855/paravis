include_directories(${cardiac_model_SOURCE_DIR}/utilities)
link_directories(${cardiac_model_BINARY_DIR}/utilities)
find_package(CUDA REQUIRED)
cuda_add_library(cell cell.cc)
