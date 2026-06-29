project(pappus_lib)

find_package(small_vector REQUIRED)
find_package(eve REQUIRED)

# header-only INTERFACE library; sources are fully templatized in headers
add_library(pappus INTERFACE)
target_include_directories(pappus INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(pappus INTERFACE gch::small_vector eve::eve)
target_compile_features(pappus INTERFACE cxx_std_23)
target_compile_definitions(pappus INTERFACE
    "$<$<BOOL:${DIRECTED_ROUNDING}>:DIRECTED_ROUNDING>"
)

if(CMAKE_EXPORT_COMPILE_COMMANDS)
    set(CMAKE_CXX_STANDARD_INCLUDE_DIRECTORIES ${CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES})
endif()
