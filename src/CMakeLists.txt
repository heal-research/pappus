project(pappus_lib)

find_package(Eigen3 REQUIRED)

add_library(pappus
    SHARED
    affine.cpp
)

target_include_directories(pappus PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${EIGEN3_INCLUDE_DIR})
target_compile_features(pappus PRIVATE cxx_std_17)
