project(pappus_test)

set(PAPPUS_TARGET_ARCH "native" CACHE STRING "Target architecture for the test binary (-march)")
find_package(Catch2 3 REQUIRED)
find_package(nanobench REQUIRED)
find_package(aaflib REQUIRED)
include(Catch)

add_executable(paptest
    affine.cpp
    interval.cpp
    benchmarks.cpp
    util.cpp
    )
target_link_libraries(paptest pappus aaflib::aaf nanobench::nanobench Catch2::Catch2WithMain)
target_compile_features(paptest PRIVATE cxx_std_23)
target_compile_options(paptest PRIVATE "$<$<CONFIG:Debug>:-g;>$<$<CONFIG:Release>:-O3;-g;-march=${PAPPUS_TARGET_ARCH}>;-fPIC")

catch_discover_tests(paptest)

if(CMAKE_EXPORT_COMPILE_COMMANDS)
    set(CMAKE_CXX_STANDARD_INCLUDE_DIRECTORIES ${CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES})
endif()
