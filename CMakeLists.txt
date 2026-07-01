cmake_minimum_required(VERSION 3.14)

include(cmake/prelude.cmake)

project(
    pappus
    VERSION 0.1.0
    DESCRIPTION "Interval and affine arithmetic library"
    HOMEPAGE_URL "https://github.com/foolnotion/pappus"
    LANGUAGES CXX
)

include(cmake/project-is-top-level.cmake)
include(cmake/variables.cmake)

# ---- Declare library ----

find_package(small_vector REQUIRED)
find_package(eve REQUIRED)

add_library(pappus_pappus INTERFACE)
add_library(pappus::pappus ALIAS pappus_pappus)

set_property(TARGET pappus_pappus PROPERTY EXPORT_NAME pappus)

target_include_directories(
    pappus_pappus ${warning_guard}
    INTERFACE
    "\$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
)

target_compile_features(pappus_pappus INTERFACE cxx_std_23)
target_link_libraries(pappus_pappus INTERFACE gch::small_vector eve::eve)

set(DIRECTED_ROUNDING_DESCRIPTION "Use directed rounding towards -inf or +inf for interval bounds.")
option(DIRECTED_ROUNDING ${DIRECTED_ROUNDING_DESCRIPTION} ON)
target_compile_definitions(
    pappus_pappus INTERFACE
    "$<$<BOOL:${DIRECTED_ROUNDING}>:DIRECTED_ROUNDING>"
)

# ---- Install rules ----

if(NOT CMAKE_SKIP_INSTALL_RULES)
  include(cmake/install-rules.cmake)
endif()

# ---- Developer mode ----

if(NOT pappus_DEVELOPER_MODE)
  return()
elseif(NOT PROJECT_IS_TOP_LEVEL)
  message(
      AUTHOR_WARNING
      "Developer mode is intended for developers of pappus"
  )
endif()

include(cmake/dev-mode.cmake)
