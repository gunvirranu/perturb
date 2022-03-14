cmake_minimum_required(VERSION 3.14)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
    message(FATAL_ERROR "In-source builds are not supported.")
endif()

project(test_perturb LANGUAGES CXX)

include(FetchContent)
FetchContent_Declare(
    doctest
    GIT_REPOSITORY "https://github.com/doctest/doctest"
    GIT_TAG "v2.4.8"
)
FetchContent_MakeAvailable(doctest)

add_executable(test_perturb test_perturb.cpp)
target_link_libraries(test_perturb PRIVATE perturb doctest)

target_compile_features(test_perturb PRIVATE cxx_std_11)
set_target_properties(
    perturb
    PROPERTIES
        CXX_STANDARD 11
        CXX_STANDARD_REQUIRED YES
        CXX_EXTENSIONS NO
)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    # TODO: Add -Werror once Vallado code is cleaned up
    target_compile_options(perturb PUBLIC -Wall -Wpedantic -Wextra)
elseif(MSVC)
    target_compile_options(perturb PUBLIC /W4 /WX)
    target_compile_definitions(test_perturb PUBLIC DOCTEST_CONFIG_USE_STD_HEADERS)
endif()

target_compile_definitions(
    perturb
    PUBLIC
        # PERTURB_DISABLE_IO
        PERTURB_VALLADO_SGP4_ENABLE_DEBUG=true
)

enable_testing()
include(${doctest_SOURCE_DIR}/scripts/cmake/doctest.cmake)
doctest_discover_tests(test_perturb)
