/*
 * perturb -- A modern C++11 wrapper for the SGP4 orbit propagator
 * Version 1.0.0
 * https://github.com/gunvirranu/perturb
 *
 * Licensed under the MIT License <http://opensource.org/licenses/MIT>.
 * SPDX-License-Identifier: MIT
 *
 * Copyright (c) 2022 Gunvir Singh Ranu
 */

//! @file Internal private header for perturb's source files
//! @author Gunvir Singh Ranu
//! @version 1.0.0
//! @copyright Gunvir Singh Ranu, MIT License

#ifndef PERTURB_COMMON_PRIVATE_H
#define PERTURB_COMMON_PRIVATE_H

typedef perturb_real_t real_t;

#define PI  3.14159265358979323846

#define MINS_PER_DAY  (24 * 60)

#define FABS(x)     fabs(x)
#define FLOOR(x)    floor(x)

#endif  // PERTURB_COMMON_PRIVATE_H
