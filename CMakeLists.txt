cmake_minimum_required(VERSION 3.6)

include(FetchContent)
include(FeatureSummary)

project(pappus_proj
    LANGUAGES CXX
    )

set(BUILD_TESTS_DESCRIPTION "Build unit tests.")
set(CRLIBM_DESCRIPTION "Use CRlibm, a correctly rounded arithmetic library")

option(BUILD_TESTS  ${BUILD_TESTS_DESCRIPTION} OFF)
option(USE_CRLIBM ${CRLIBM_DESCRIPTION}      OFF)

if(BUILD_TESTS)
    if (USE_CRLIBM)
        message(WARNING "Unit tests assume a standard libm (not CRlibm), will not build tests")
    else()
        FetchContent_Declare(
            aaflib
            GIT_REPOSITORY https://github.com/foolnotion/aaflib.git
            GIT_TAG master
            DOWNLOAD_DIR ${PROJECT_SOURCE_DIR}/thirdparty/aaflib
            SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/aaflib
            )
        FetchContent_MakeAvailable(aaflib)
        SET(AAF_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/aaflib)
        add_subdirectory(test)
    endif()
endif()

if(CMAKE_EXPORT_COMPILE_COMMANDS)
    set(CMAKE_CXX_STANDARD_INCLUDE_DIRECTORIES ${CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES})
endif()

add_subdirectory(src)
