cmake_minimum_required(VERSION 3.14)

# Create your own very cool project!
project(
    VeryCoolExampleProject
    DESCRIPTION "Example of using perturb from a local source clone"
    LANGUAGES CXX
)

# A very cool executable indeed
add_executable(VeryCoolExampleProject main.cpp)

# Add perturb library root directory
# I have to use the local `build/` binary_dir cause my root one already being used
add_subdirectory(../../ build)

# Link together that bad boy as a dependency
target_link_libraries(VeryCoolExampleProject perturb)
