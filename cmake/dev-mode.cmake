include(cmake/folders.cmake)

# enable implicit include directories in the compile commands
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE INTERNAL "")
if(CMAKE_EXPORT_COMPILE_COMMANDS)
    set(CMAKE_CXX_STANDARD_INCLUDE_DIRECTORIES ${CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES})
endif()

include(CTest)
if(BUILD_TESTING)
  add_subdirectory(test)
endif()

add_folders(Project)
