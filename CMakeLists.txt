cmake_minimum_required(VERSION 3.6)

include(FetchContent)
include(FeatureSummary)

project(pappus_proj
    LANGUAGES CXX
    )

FetchContent_Declare(
    aaflib
    GIT_REPOSITORY https://github.com/foolnotion/aaflib.git
    GIT_TAG master
    DOWNLOAD_DIR ${PROJECT_SOURCE_DIR}/thirdparty/aaflib
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/aaflib
    )

FetchContent_MakeAvailable(aaflib)

SET(AAF_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/aaflib)

if(CMAKE_EXPORT_COMPILE_COMMANDS)
    set(CMAKE_CXX_STANDARD_INCLUDE_DIRECTORIES ${CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES})
endif()

add_subdirectory(src)
add_subdirectory(test)
