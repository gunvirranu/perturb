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
    GIT_TAG 5723e50a54d438876b7a02e02a326e30b51ead4c  # Release v1.0.0
)
FetchContent_MakeAvailable(perturb)

# Link together that bad boy as a dependency
target_link_libraries(VeryCoolExampleProject perturb)
