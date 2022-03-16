
# Examples

These examples demonstrate how easy it is to install and setup the `perturb` library with CMake. The important part for each is the `CMakeLists.txt` file for each example.  The steps for actually running these are the exact same, so they're listed below.

## Local Checkout: [`cmake-local`](cmake-local)

This example is if you prefer to manage checking out the `perturb` source repo yourself, either manually, using Git Submodules, or any other method. The only requirement is that you can point CMake to the root folder of the `perturb` library, which contains its `CMakeLists.txt`.

Since this example *is already in* in the repo, it can just point two directories up (`../../`) to the root folder. You can change this to whatever folder you clone the `perturb` repo to. For example, you might clone it in an `external/` or `third-party/` folder. 

## CMake's FetchContent: [`cmake-fetch-content`](cmake-fetch-content)

This example is if you prefer to have CMake download and checkout this repository for you during the configuration stage. This is done through a very cool feature of CMake called [`FetchContent`](https://bewagner.net/programming/2020/05/02/cmake-fetchcontent). See [this documentation page](https://cmake.org/cmake/help/latest/module/FetchContent.html) for more options, such as pinning to a specific version tag (which is recommended).

## How to Run

1. Requirements are a recent-ish version of CMake and a C++11 compiler
2. `git clone https://github.com/gunvirranu/perturb.git`
3. `cd perturb/examples/cmake-local`
3. `cmake -S . -B build`
4. `cmake --build build`
5. `./build/VeryCoolProject`
6. Be amazed at just how simple [the `CMakeLists.txt` file we used](cmake-local/CMakeLists.txt) is
