if(PROJECT_IS_TOP_LEVEL)
    set(CMAKE_INSTALL_INCLUDEDIR include/perturb CACHE PATH "")
endif()

include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

# find_package(<package>) call for consumers to find this project
set(package perturb)

install(
    DIRECTORY include/
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
    COMPONENT perturb_Development
)

install(
    TARGETS perturb
    EXPORT perturbTargets
    ARCHIVE #
    COMPONENT perturb_Development
    INCLUDES #
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
)

write_basic_package_version_file(
    "${package}ConfigVersion.cmake"
    COMPATIBILITY SameMajorVersion
)

# Allow package maintainers to freely override the path for the configs
set(
    perturb_INSTALL_CMAKEDIR "${CMAKE_INSTALL_DATADIR}/${package}"
    CACHE PATH "CMake package config location relative to the install prefix"
)
mark_as_advanced(perturb_INSTALL_CMAKEDIR)

install(
    FILES cmake/install-config.cmake
    DESTINATION "${perturb_INSTALL_CMAKEDIR}"
    RENAME "${package}Config.cmake"
    COMPONENT perturb_Development
)

install(
    FILES "${PROJECT_BINARY_DIR}/${package}ConfigVersion.cmake"
    DESTINATION "${perturb_INSTALL_CMAKEDIR}"
    COMPONENT perturb_Development
)

install(
    EXPORT perturbTargets
    NAMESPACE perturb::
    DESTINATION "${perturb_INSTALL_CMAKEDIR}"
    COMPONENT perturb_Development
)

if(PROJECT_IS_TOP_LEVEL)
    include(CPack)
endif()
