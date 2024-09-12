cmake_minimum_required(VERSION 3.14)

project(
    VeryCoolExampleProject
    DESCRIPTION "Example of using perturb from a CMake C99 application"
    LANGUAGES C
)

add_executable(VeryCoolExampleProject main.c)

# Use local perturb library
enable_language(CXX)
set(perturb_EXPORT_C_INTERFACE ON)
add_subdirectory(../../ build)

target_link_libraries(VeryCoolExampleProject perturb)
