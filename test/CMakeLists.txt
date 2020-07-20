project(pappus_test)

find_package(doctest REQUIRED)

add_executable(paptest main.cpp)
target_link_libraries(paptest pappus doctest::doctest)
target_compile_features(paptest PRIVATE cxx_std_17)
