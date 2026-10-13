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
# logs). Use a concrete architecture for each native runner; CI overrides this
# cache variable for x86-64 and AArch64 explicitly.
set(PAPPUS_TARGET_ARCH "" CACHE STRING "Target architecture for the test binary (-march); empty uses the compiler default")
find_package(Catch2 3 REQUIRED)
find_package(nanobench REQUIRED)
find_package(aaflib REQUIRED)
option(PAPPUS_ENABLE_MPFR_SOUNDNESS "Build MPFR-based interval soundness tests." ON)
if(PAPPUS_ENABLE_MPFR_SOUNDNESS)
  find_path(MPFR_INCLUDE_DIR mpfr.h REQUIRED)
  find_library(MPFR_LIBRARY mpfr REQUIRED)
  find_library(GMP_LIBRARY gmp REQUIRED)
endif()
include(Catch)

# ---- Tests ----

add_executable(paptest
    source/affine.cpp
    source/interval.cpp
    source/benchmarks.cpp
    source/util.cpp
    source/fuzz.cpp
)
if(PAPPUS_ENABLE_MPFR_SOUNDNESS)
  target_sources(paptest PRIVATE source/soundness.cpp)
  target_include_directories(paptest SYSTEM PRIVATE ${MPFR_INCLUDE_DIR})
  target_link_libraries(paptest PRIVATE ${MPFR_LIBRARY} ${GMP_LIBRARY})
endif()
target_link_libraries(paptest PRIVATE
    pappus::pappus
    aaflib::aaf
    nanobench::nanobench
    Catch2::Catch2WithMain
)
target_compile_features(paptest PRIVATE cxx_std_23)
target_compile_options(paptest PRIVATE "$<$<CONFIG:Debug>:-g;>$<$<CONFIG:Release>:-O3;-g;>-fPIC")
if(PAPPUS_TARGET_ARCH)
  target_compile_options(paptest PRIVATE "-march=${PAPPUS_TARGET_ARCH}")
endif()

catch_discover_tests(paptest)

# ---- End-of-file commands ----

add_folders(Test)
