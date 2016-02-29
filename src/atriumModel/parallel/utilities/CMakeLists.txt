include_directories(${cardiac_model_SOURCE_DIR}/utilities)
link_directories(${cardiac_model_BINARY_DIR}/utilities)
add_library(cell cell.cc)
