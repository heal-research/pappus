project(pappus_test)

find_package(doctest REQUIRED)

add_executable(paptest main.cpp)
target_link_libraries(paptest pappus aaf.lin64 doctest::doctest)
target_compile_features(paptest PRIVATE cxx_std_17)
target_include_directories(paptest PRIVATE ${AAF_INCLUDE_DIR})
