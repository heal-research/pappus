cmake_minimum_required(VERSION 3.14)

project(pappusTests LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

include(../cmake/project-is-top-level.cmake)
include(../cmake/folders.cmake)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE INTERNAL "")
if(CMAKE_EXPORT_COMPILE_COMMANDS)
    set(CMAKE_CXX_STANDARD_INCLUDE_DIRECTORIES ${CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES})
endif()

# ---- Dependencies ----

if(PROJECT_IS_TOP_LEVEL)
  find_package(pappus REQUIRED)
  enable_testing()
endif()

# "native" is silently dropped by nix's sandbox (NIX_ENFORCE_NO_NATIVE),
# falling back to a near-scalar baseline and making SIMD benchmarks
# meaningless without anyone noticing (no warning surfaces past `nix build`
# logs). Default to a concrete AVX2 target instead; override for other ISAs.
set(PAPPUS_TARGET_ARCH "x86-64-v3" CACHE STRING "Target architecture for the test binary (-march)")
find_package(Catch2 3 REQUIRED)
find_package(nanobench REQUIRED)
find_package(aaflib REQUIRED)
include(Catch)

# ---- Tests ----

add_executable(paptest
    source/affine.cpp
    source/interval.cpp
    source/benchmarks.cpp
    source/util.cpp
    source/fuzz.cpp
)
target_link_libraries(paptest PRIVATE
    pappus::pappus
    aaflib::aaf
    nanobench::nanobench
    Catch2::Catch2WithMain
)
target_compile_features(paptest PRIVATE cxx_std_23)
target_compile_options(paptest PRIVATE "$<$<CONFIG:Debug>:-g;>$<$<CONFIG:Release>:-O3;-g;>-march=${PAPPUS_TARGET_ARCH};-fPIC")

catch_discover_tests(paptest)

# ---- End-of-file commands ----

add_folders(Test)
