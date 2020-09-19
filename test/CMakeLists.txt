project(pappus_test)

find_package(doctest REQUIRED)

add_executable(paptest affine.cpp interval.cpp)
target_link_libraries(paptest pappus aaf.lin64 doctest::doctest)
target_compile_features(paptest PRIVATE cxx_std_17)
target_include_directories(paptest PRIVATE ${AAF_INCLUDE_DIR})

if(CMAKE_EXPORT_COMPILE_COMMANDS)
    set(CMAKE_CXX_STANDARD_INCLUDE_DIRECTORIES ${CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES})
endif()
