/*
* perturb -- A modern C++11 wrapper for the SGP4 orbit propagator
* Version 1.0.0
* https://github.com/gunvirranu/perturb
*
* Licensed under the MIT License <http://opensource.org/licenses/MIT>.
* SPDX-License-Identifier: MIT
*
* Copyright (c) 2024 Gunvir Ranu
*/

//! @file
//! Exported C interface header for the C++ perturb library
//! @author Gunvir Singh Ranu
//! @version 1.0.0
//! @copyright Gunvir Singh Ranu, MIT License

#ifndef PERTURB_PERTURB_H
#define PERTURB_PERTURB_H

// Sanity checks on if you're including this header correctly
// `perturb.hpp` is the primary C++ header for this library, perturb.
// This header, `perturb.h`, is intended only for C ABI consumers.
#ifndef PERTURB_EXPORT_C_INTERFACE
#  if __cplusplus
#    error "Seems like you're using C++. Use the perturb.hpp header."
#  else
#    error "Seems like you're using C. Enable the PERTURB_EXPORT_C_INTERFACE option."
#  endif
#endif
#ifdef PERTURB_EXPORT_C_INTERFACE

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#define PERTURB_TLE_LINE_LEN 69U

// Support header inclusion from C++
#if __cplusplus
extern "C" {
#endif

typedef double perturb_real_t;

extern const size_t PERTURB_C_SATELLITE_ALLOC_SIZE;
extern const size_t PERTURB_C_TLE_ALLOC_SIZE;

struct perturb_satellite {
    void * p_satellite;
};

struct perturb_tle {
    void * p_tle;
};

struct perturb_julian_date {
    perturb_real_t jd;
    perturb_real_t jd_frac;
};

struct perturb_state_vector {
    struct perturb_julian_date epoch;
    perturb_real_t position[3];
    perturb_real_t velocity[3];
};

struct perturb_julian_date perturb_datetime_to_julian(
    uint8_t year,
    uint8_t month,
    uint8_t day,
    uint8_t hour,
    uint8_t min,
    perturb_real_t sec
);

struct perturb_julian_date perturb_add_days_to_julian(
    struct perturb_julian_date t,
    perturb_real_t days
);

struct perturb_julian_date perturb_epoch(struct perturb_satellite sat);

bool perturb_init_sat_from_tle(struct perturb_satellite sat, struct perturb_tle tle);
bool perturb_parse_tle(char * line_1, char * line_2, struct perturb_tle tle);
bool perturb_parse_tle_and_init_sat(char * line_1, char * line_2, struct perturb_satellite sat);

struct perturb_state_vector perturb_propagate_days(struct perturb_satellite sat, perturb_real_t days);
struct perturb_state_vector perturb_propagate(struct perturb_satellite sat, struct perturb_julian_date t);

#if __cplusplus
}
#endif

#endif  // PERTURB_EXPORT_C_INTERFACE
#endif  // PERTURB_PERTURB_H
