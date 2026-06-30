cmake_minimum_required(VERSION 3.6)

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
    add_subdirectory(test)
endif()

add_subdirectory(src)
