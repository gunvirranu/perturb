<div class="title-block" style="text-align: center;" align="center">

# 🌎 `perturb` 🌏

#### A modern C++11 wrapper for the SGP4 orbit propagator

[![Release][release_badge]][release]
[![Docs][docs_badge]][docs]
[![License][license_badge]][license]
[![CI][ci_badge]][ci]

</div>

Don't you just hate that awkward moment when you need to propagate Earth-centred spacecraft trajectories using the accurate-ish [SGP4][SGP4] model, but also want to avoid having to understand [source code that literally dates back to the 80s](original/SGP4.cpp)? Well you're in luck, cause `perturb` provides a simple interface so you don't have to. Throw this library a [TLE][TLE] orbital element and it'll spit out a position and velocity vector.

## Features

- [Uses Vallado's latest and greatest de-facto standard SGP4 implementation](original)
- Strictly standards compliant C++11 with no dependencies
- [Embedded support](#build-options) (i.e. no dynamic memory, no recursion, no exceptions, no virtual, no RTTI, none of that runtime funny business)
- [Built with modern CMake](#install) with a relatively standard library layout
- Has tests and CI, which is more than I can say for most of my code

## Example

See the [`examples`](examples) directory to see how you can try out this exact code snippet.

```cpp
#include <cassert>
#include <iostream>

#include <perturb/perturb.hpp>

using namespace perturb;  // I know, I know. I promise it's only for examples!

int main() {
    // Let try simulating the orbit of the International Space Station
    // Got TLE from Celestrak sometime around 2022-03-12
    std::string ISS_TLE_1 = "1 25544U 98067A   22071.78032407  .00021395  00000-0  39008-3 0  9996";
    std::string ISS_TLE_2 = "2 25544  51.6424  94.0370 0004047 256.5103  89.8846 15.49386383330227";

    // Create and initialize a satellite object from the TLE
    auto sat = Satellite::from_tle(ISS_TLE_1, ISS_TLE_2);
    assert(sat.last_error() == Sgp4Error::NONE);
    assert(sat.epoch().to_datetime().day == 12);

    // Let's see what the ISS is doing on Pi Day
    const auto t = JulianDate(DateTime { 2022, 3, 14, 1, 59, 26.535 });
    const double delta_days = t - sat.epoch();
    assert(1 < delta_days && delta_days < 3);  // It's been ~2 days since the epoch

    // Calculate the position and velocity at the chosen time
    StateVector sv;
    const auto err = sat.propagate(t, sv);
    assert(err == Sgp4Error::NONE);
    const auto &pos = sv.position, &vel = sv.velocity;

    // Conclusion: The ISS is going pretty fast (~8 km/s)
    std::cout << "Position [km]: { " << pos[0] << ", " << pos[1] << ", " << pos[2] << " }\n";
    std::cout << "Velocity [km/s]: { " << vel[0] << ", " << vel[1] << ", " << vel[2] << " }\n";
}
```

## Install

Most C++ build systems kinda suck, but CMake is the closest thing to a standard we have. This library exposes a [CMake target][modern-CMake] (named `perturb`) which hopefully should make using it very simple. To get a copy of the source code, you have a few options. In my personal order of preference:

### [CPM.cmake][CPM.cmake]

This is a simple CMake-based package manager. Just point it at this Github repo with a tagged release version and it should just work.

### [CMake's FetchContent][cmake-fetchcontent]

Point it at this Github repo with a tagged release version and you should be good to go. Here's a snippet demonstrating this:

```cmake
include(FetchContent)
# Download and setup the perturb library
FetchContent_Declare(
    perturb
    GIT_REPOSITORY "https://github.com/gunvirranu/perturb"
)
FetchContent_MakeAvailable(perturb)

# Link perturb into your project
target_link_libraries(YOUR_VERY_COOL_PROJECT_TARGET perturb)
```

Also see the [`examples`](examples) directory for a full CMake example project using this. You can make more complicated adjustments as necessary, but this demonstrates the basics.

### Git Submodules

Personally, I'm not a fan, but a lot of people seem to love these. Once you have the `perturb` repo set up as a submodule, let's say in the `external` folder, you can use it like so:

```cmake
# Add perturb library root directory
add_subdirectory(external/perturb)

# Link perturb into your project
target_link_libraries(YOUR_VERY_COOL_PROJECT_TARGET perturb)
```

### Any Other Method of Cloning

If you clone this repo manually into a specific location or have any other method of checking it out, you can then follow the exact same steps as using [Git Submodules](#git-submodules).

### Manually Downloading Source Files

You can also just download the headers from [`include/perturb`](include/perturb) and source files from [`src`](src) and then just check them straight into your project with your own build system. Not very elegant, but it gets the job done I guess.

### Other Build Systems

`perturb` _currently_ doesn't support any other build systems, primarily due to my lack of knowledge on them, but I'm open to it. If you're interested, open an issue and we can look into it.

## Basic Usage

Everything in this library is in the `perturb::` namespace. Here's a quick intro to the typical "usage flow".

I won't cover the details of [SGP4][SGP4], but in brief, it's a very popular orbit propagator for Earth-centered spacecraft. Usually, the input orbit ephemeris is through a [TLE][TLE], such as the ones provided by [Celestrek][Celestrek]. These TLE inputs can be used to construct a `perturb::Satellite` object.

A specific point in time is represented as a `perturb::JulianDate`. You can either construct one from a specific date and time via `perturb::DateTime` or offset a number of days from the `epoch()` of a satellite.

Passing in a time point to the `Satellite::propagate(...)` method of a satellite yields a `perturb::StateVector`, which contains a time-stamp, and a position and velocity vector. These vectors are just a `std::array<double, 3>`, measured in kilometres, and are represented in the [TEME][ECI-TEME] coordinate reference frame. The details of this frame can get a bit annoying, so this library does _not_ handle converting it to others. For handling Earth-centered reference frames such as TEME and transformations between them, you may be interested in the [`gelocus`][gelocus] library.

Check out [this page][perturb-docs] for some slightly more detailed documentation of the interface.

## Build Options

You *probably* don't need to configure anything.

### Disabling I/O

The TLE parsing code (which accepts a `char *` or `std::string`) may be undesired for embedded applications, either due to inefficient codegen from `sscanf` (which is used for parsing) or to avoid dynamic memory. This I/O functionality can be disabled by defining the `PERTURB_DISABLE_IO` preprocessor flag in your build system, which will strip out any mentions of I/O and strings. It is not defined by default, so the functionality is usually included.

In CMake, once you have the `perturb` target initialized, you can do:

```cmake
target_compile_definitions(perturb PRIVATE PERTURB_DISABLE_IO)
```

You could also set the `perturb_DISABLE_IO` option in CMake to `ON` before initializing the `perturb` target. This is also by default `OFF`. Setting this option will handle defining the `PERTURB_DISABLE_IO` preprocessor flag.

Do note, this will leave you with no way of parsing TLEs. You will need to pre-parse the TLE and initialize the propagator using numerical values directly. This can be done by initializing the `TwoLineElement` type however you wish and using that to construct a `Satellite` object via its constructor.

## Changelog

See [`CHANGELOG.md`](CHANGELOG.md).

## License

`perturb` is licensed under the [MIT license](LICENSE).

<!-- Links -->
[SGP4]: https://en.wikipedia.org/wiki/Simplified_perturbations_models
[TLE]: https://en.wikipedia.org/wiki/Two-line_element_set
[modern-CMake]: https://pabloariasal.github.io/2018/02/19/its-time-to-do-cmake-right
[CPM.cmake]: https://github.com/cpm-cmake/CPM.cmake
[cmake-fetchcontent]: https://bewagner.net/programming/2020/05/02/cmake-fetchcontent
[Celestrek]: https://celestrak.com
[ECI-TEME]: https://en.wikipedia.org/wiki/Earth-centered_inertial
[gelocus]: https://github.com/gunvirranu/gelocus
[perturb-docs]: https://gunvirranu.github.io/perturb

<!-- Badges -->
[release]: https://github.com/gunvirranu/perturb/releases "Release"
[release_badge]: https://img.shields.io/github/v/release/gunvirranu/perturb?color=orange&logo=github&sort=semver "Release"
[docs]: https://gunvirranu.github.io/perturb "Docs"
[docs_badge]: https://img.shields.io/badge/docs-passing-success "Docs"
[license]: #license "License"
[license_badge]: https://img.shields.io/badge/license-MIT-blue.svg "License"
[ci]: https://github.com/gunvirranu/perturb/actions "Github Actions"
[ci_badge]: https://github.com/gunvirranu/perturb/workflows/CI/badge.svg?branch=master "Github Actions"
