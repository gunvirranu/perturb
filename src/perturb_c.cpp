/*
* perturb -- A modern C++11 wrapper for the SGP4 orbit propagator
* Version 1.0.0
* https://github.com/gunvirranu/perturb
*
* Licensed under the MIT License <http://opensource.org/licenses/MIT>.
* SPDX-License-Identifier: MIT
*
* Copyright (c) 2024 Gunvir Singh Ranu
*/

#include "perturb/perturb.h"
#include "perturb/perturb.hpp"

#ifdef PERTURB_EXPORT_C_INTERFACE

using namespace perturb;

extern "C" {

const size_t PERTURB_C_SATELLITE_ALLOC_SIZE = sizeof(Satellite);
const size_t PERTURB_C_TLE_ALLOC_SIZE = sizeof(TwoLineElement);

}

#endif  // PERTURB_EXPORT_C_INTERFACE
