if(PROJECT_IS_TOP_LEVEL)
  set(
      CMAKE_INSTALL_INCLUDEDIR "include/pappus-${PROJECT_VERSION}"
      CACHE STRING ""
  )
  set_property(CACHE CMAKE_INSTALL_INCLUDEDIR PROPERTY TYPE PATH)
endif()

include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

# find_package(<package>) call for consumers to find this project
# should match the name of variable set in the install-config.cmake script
set(package pappus)

install(
    DIRECTORY include/
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
    COMPONENT pappus_Development
)

install(
    TARGETS pappus_pappus
    EXPORT pappusTargets
    INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
)

write_basic_package_version_file(
    "${package}ConfigVersion.cmake"
    COMPATIBILITY SameMajorVersion
    ARCH_INDEPENDENT
)

# Allow package maintainers to freely override the path for the configs
set(
    pappus_INSTALL_CMAKEDIR "${CMAKE_INSTALL_DATADIR}/${package}"
    CACHE STRING "CMake package config location relative to the install prefix"
)
set_property(CACHE pappus_INSTALL_CMAKEDIR PROPERTY TYPE PATH)
mark_as_advanced(pappus_INSTALL_CMAKEDIR)

install(
    FILES cmake/install-config.cmake
    DESTINATION "${pappus_INSTALL_CMAKEDIR}"
    RENAME "${package}Config.cmake"
    COMPONENT pappus_Development
)

install(
    FILES "${PROJECT_BINARY_DIR}/${package}ConfigVersion.cmake"
    DESTINATION "${pappus_INSTALL_CMAKEDIR}"
    COMPONENT pappus_Development
)

install(
    EXPORT pappusTargets
    NAMESPACE pappus::
    DESTINATION "${pappus_INSTALL_CMAKEDIR}"
    COMPONENT pappus_Development
)

if(PROJECT_IS_TOP_LEVEL)
  include(CPack)
endif()
