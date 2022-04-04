cmake_minimum_required(VERSION 3.14)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
    message(FATAL_ERROR "In-source builds are not supported.")
endif()

project(test_perturb LANGUAGES CXX)

# A hack to disable warnings for any targets except perturb
if(DEFINED CMAKE_CXX_CLANG_TIDY)
    set(CMAKE_CXX_CLANG_TIDY_save "${CMAKE_CXX_CLANG_TIDY}")
    set(CMAKE_CXX_CLANG_TIDY "")
endif()

include(FetchContent)
FetchContent_Declare(
    doctest
    GIT_REPOSITORY "https://github.com/doctest/doctest"
    GIT_TAG "v2.4.8"
)
FetchContent_MakeAvailable(doctest)

add_executable(test_perturb test_perturb.cpp)
target_link_libraries(test_perturb PRIVATE perturb doctest)

# Restore previous
if(DEFINED CMAKE_CXX_CLANG_TIDY_save)
    set(CMAKE_CXX_CLANG_TIDY "${CMAKE_CXX_CLANG_TIDY_save}")
endif()

target_compile_features(test_perturb PRIVATE cxx_std_11)
set_target_properties(
    perturb
    PROPERTIES
        CXX_STANDARD 11
        CXX_STANDARD_REQUIRED YES
        CXX_EXTENSIONS NO
)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    # Do nothing
elseif(MSVC)
    target_compile_definitions(test_perturb PUBLIC DOCTEST_CONFIG_USE_STD_HEADERS)
endif()

if(NOT perturb_DISABLE_IO)
    target_compile_definitions(
        perturb
        PUBLIC
            PERTURB_VALLADO_SGP4_ENABLE_DEBUG
    )
endif()

enable_testing()
include(${doctest_SOURCE_DIR}/scripts/cmake/doctest.cmake)
doctest_discover_tests(test_perturb)
