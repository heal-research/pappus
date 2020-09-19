project(pappus_lib)

find_package(Eigen3 REQUIRED)

add_library(pappus
    STATIC
    affine/approx.cpp
    affine/common.cpp
    affine/linear.cpp
    interval.cpp
)

target_include_directories(pappus PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${EIGEN3_INCLUDE_DIR})
target_compile_features(pappus PRIVATE cxx_std_17)
target_link_libraries(pappus PRIVATE crlibm)
target_compile_definitions(pappus PRIVATE "$<$<BOOL:${USE_CRLIBM}>:USE_CRLIBM>")
