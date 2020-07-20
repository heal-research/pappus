project(pappus_lib)

find_package(Eigen3 REQUIRED)

add_library(pappus INTERFACE)

target_include_directories(pappus INTERFACE ${CMAKE_CURRENT_SOURCE_DIR} ${EIGEN3_INCLUDE_DIR})
target_compile_features(pappus INTERFACE cxx_std_17)
