
# Vallado's Original SGP4 Source

This directory contains the original unmodified implementation of SGP4 provided by the great wizard of astrodynamics, David Vallado. This repo only contains the latest revision from 2020, but this and older versions can be found scattered all over the interwebs. This copy should be identical to the many others online. Note, I did convert CRLF to LF and tabs to spaces.

While I originally intended to wrap this source without _any_ modifications, I needed to make changes to disable I/O features and a few other small things. A lightly modified version is used in the library (see [`vallado_sgp.hpp`](../include/perturb/vallado_sgp4.hpp) and [`vallado_sgp4.cpp`](../src/vallado_sgp4.cpp)). You can see the exact changes with a simple `diff` or using your favourite IDE.

### Important Changes

- Put everything inside namespace `perturb::vallado_sgp4`
- Gate any I/O and string processing by the `PERTURB_DISABLE_IO` flag
  - Removes `twoline2rv` TLE parsing method as it uses `sscanf(...)`
  - Remove debug printing via `fprintf(stderr, ...)`
- Gate the verification mode by the `PERTURB_VALLADO_SGP4_ENABLE_DEBUG` flag
- Some small refactoring to fix compiler warnings and lints
- Refactor use of `strcpy` to `memcpy`
