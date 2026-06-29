cmake_minimum_required(VERSION 3.6)

include(FetchContent)
include(FeatureSummary)

project(pappus_proj
    LANGUAGES CXX
    )

include(CTest)

set(BUILD_TESTS_DESCRIPTION "Build unit tests.")
set(DIRECTED_ROUNDING_DESCRIPTION "Use directed rounding towards -inf or +inf for interval bounds.")

option(BUILD_TESTS ${BUILD_TESTS_DESCRIPTION} ON)
option(DIRECTED_ROUNDING ${DIRECTED_ROUNDING_DESCRIPTION} ON)

add_feature_info(BUILD_TESTS       BUILD_TESTS       ${BUILD_TESTS_DESCRIPTION})
add_feature_info(DIRECTED_ROUNDING DIRECTED_ROUNDING ${DIRECTED_ROUNDING_DESCRIPTION})

# print a status of what we found
feature_summary(WHAT ENABLED_FEATURES DESCRIPTION "Enabled features:" QUIET_ON_EMPTY)
feature_summary(WHAT DISABLED_FEATURES DESCRIPTION "Disabled features:" QUIET_ON_EMPTY)

if(BUILD_TESTS)
    enable_testing()
    FetchContent_Declare(
        aaflib
        GIT_REPOSITORY https://github.com/foolnotion/aaflib.git
        GIT_TAG master
        DOWNLOAD_DIR ${PROJECT_SOURCE_DIR}/thirdparty/aaflib
        SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/aaflib
        )
    FetchContent_MakeAvailable(aaflib)
    SET(AAF_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/aaflib)

    FetchContent_Declare(
        nanobench 
        DOWNLOAD_DIR ${PROJECT_SOURCE_DIR}/thirdparty/nanobench/
        URL https://raw.githubusercontent.com/martinus/nanobench/master/src/include/nanobench.h
        DOWNLOAD_NO_EXTRACT 1
        )
    FetchContent_MakeAvailable(nanobench)
    SET(NANOBENCH_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/nanobench)

    add_subdirectory(test)
endif()

add_subdirectory(src)
