/*
 * perturb -- A modern C++11 wrapper for the SGP4 orbit propagator
 * version 1.0.0
 * https://github.com/gunvirranu/perturb
 *
 * Copyright (c) 2022 Gunvir Singh Ranu
 * SPDX-License-Identifier: MIT
 */

//! @file Internal private header for perturb's source files
//! @author Gunvir Singh Ranu
//! @version 1.0.0
//! @copyright Gunvir Singh Ranu, MIT License

#ifndef PERTURB_COMMON_PRIVATE_H
#define PERTURB_COMMON_PRIVATE_H

#define PI              3.14159265358979323846      ///< [-] Good ol' π
#define DEG_TO_RAD      (PI / 180.0)                ///< [rad / deg] Radians per degree

#define MINS_PER_DAY    (24 * 60)                   ///< [min / day] Minutes per 24 hour day
#define XP_DOT_P        (MINS_PER_DAY / (2 * PI))   ///< [min / rad] Minutes per radian of Earth rotation

#define UNUSED(x)       (void) (x)
#define ARRAY_SIZE(x)   (sizeof(x) / sizeof(x[0]))

typedef perturb_real_t  real_t;

#define FABS(x)         fabs(x)
#define FLOOR(x)        floor(x)

#endif  // PERTURB_COMMON_PRIVATE_H
