cmake_minimum_required(VERSION 3.6)

include(FetchContent)
include(FeatureSummary)

project(pappus_proj
    LANGUAGES CXX)

FetchContent_Declare(
    aaflib
    GIT_REPOSITORY git@github.com:foolnotion/aaflib.git 
    GIT_TAG master
    DOWNLOAD_DIR ${PROJECT_SOURCE_DIR}/thirdparty/aaflib
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/aaflib
    )

FetchContent_MakeAvailable(aaflib)

SET(AAF_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/aaflib)

add_subdirectory(src)
add_subdirectory(test)
