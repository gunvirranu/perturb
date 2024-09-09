cmake_minimum_required(VERSION 3.14)

# Create your own very cool project!
project(
    VeryCoolExampleProject
    DESCRIPTION "Example of using perturb via CMake's FetchContent"
    LANGUAGES CXX
)

# A very cool executable indeed
add_executable(VeryCoolExampleProject main.cpp)

# CMake downloads and sets up perturb for you
include(FetchContent)
FetchContent_Declare(
    perturb
    GIT_REPOSITORY "https://github.com/gunvirranu/perturb"
)
FetchContent_MakeAvailable(perturb)

# Link together that bad boy as a dependency
target_link_libraries(VeryCoolExampleProject perturb)
