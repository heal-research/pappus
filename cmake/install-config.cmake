set(pappus_FOUND YES)

include(CMakeFindDependencyMacro)
find_dependency(small_vector)
find_dependency(eve)

if(pappus_FOUND)
  include("${CMAKE_CURRENT_LIST_DIR}/pappusTargets.cmake")
endif()
