project(pappus_lib)

find_package(Eigen3 REQUIRED)

add_library(pappus
    STATIC
    affine/approx.cpp
    affine/common.cpp
    affine/linear.cpp
    interval/interval.cpp
    interval/iterator.cpp
)

target_include_directories(pappus PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${EIGEN3_INCLUDE_DIR})
target_compile_features(pappus PRIVATE cxx_std_17)
target_link_libraries(pappus PRIVATE crlibm)
target_compile_definitions(pappus PRIVATE "$<$<BOOL:${USE_CRLIBM}>:USE_CRLIBM>")
target_compile_options(pappus PRIVATE "$<$<CONFIG:Debug>:-g;>$<$<CONFIG:Release>:-O3;-g;-march=znver2>;-fPIC")

if(CMAKE_EXPORT_COMPILE_COMMANDS)
    set(CMAKE_CXX_STANDARD_INCLUDE_DIRECTORIES ${CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES})
endif()
