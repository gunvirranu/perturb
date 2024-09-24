/*
 * perturb -- A modern C++11 wrapper for the SGP4 orbit propagator
 * Version 1.0.0
 * https://github.com/gunvirranu/perturb
 *
 * Licensed under the MIT License <http://opensource.org/licenses/MIT>.
 * SPDX-License-Identifier: MIT
 *
 * Copyright (c) 2022 Gunvir Ranu
 */

//! @file
//! Header for custom TLE (two line element) processing
//! @author Gunvir Ranu
//! @version 1.0.0
//! @copyright Gunvir Ranu, MIT License

#ifndef PERTURB_TLE_H
#define PERTURB_TLE_H

#include "perturb/perturb.h"
#include "perturb/sgp4.h"

#if __cplusplus
#  ifdef PERTURB_ENABLE_CPP_INTERFACE
namespace perturb {
namespace c_internal {
#  endif
extern "C" {
#endif

/// Possible errors when parsing a TLE.
///
/// Returned by `TwoLineElement::parse` after processing a TLE record string.
///
/// @post
/// Errors (excluding `NONE`) are guaranteed to occur in definition order.
/// Meaning that spaces are checked first, then invalid format, then invalid
/// values, and lastly checksum. This allows you to assume, for example, that
/// if there is a `CHECKSUM_MISMATCH`, then none of the previous errors
/// occurred first, and so the invalid checksum is the *only* issue.
enum perturb_tle_parse_error {
    PERTURB_TLE_PARSE_ERROR_NONE,               ///< If no issues when parsing
    PERTURB_TLE_PARSE_ERROR_SHOULD_BE_SPACE,    ///< If there is a lack of space in the TLE
    PERTURB_TLE_PARSE_ERROR_INVALID_FORMAT,     ///< If general parsing was unsuccessfully
    PERTURB_TLE_PARSE_ERROR_INVALID_VALUE,      ///< If a parsed value doesn't make sense
    PERTURB_TLE_PARSE_ERROR_CHECKSUM_MISMATCH,  ///< If the checksum doesn't match
};


/// Represents a pre-parsed TLE record.
///
/// Can be generated via `TwoLineElement::parse`, but not particularly useful
/// unless you care about the specific TLE values. If you want to parse a TLE
/// from a string and use it for SGP4 in one go, then better to use the
/// `Satellite::from_tle` methods.
///
/// The primary purpose of this type is when the `PERTURB_DISABLE_IO` flag is
/// set, as then there's no way to construct a `Satellite` from a TLE. In such
/// a case where all I/O and string processing is removed, this type still
/// allows you to construct and initialize a `Satellite` manually. However, you
/// must handle your own method of creating the `TwoLineElement`s.
struct perturb_tle {
// clang-format off
    // Line 1
    char catalog_number[6];             ///< Satellite catalog number
    char classification;                ///< Classification {U: Unclassified, C: Classified, S: Secret}
    unsigned int launch_year;           ///< International Designator - Launch year (last two digits)
    unsigned int launch_number;         ///< International Designator - Launch number of the year
    char launch_piece[4];               ///< International Designator - Piece of launch
    unsigned int epoch_year;            ///< Epoch year (last two digits)
    double epoch_day_of_year;           ///< Epoch fractional day of year
    double n_dot;                       ///< First derivative of mean motion (ballistic coefficient) [rev/day^2]
    double n_ddot;                      ///< Second derivative of mean motion [rev/day^3]
    double b_star;                      ///< B* radiation pressure coefficient [1 / (earth radii)]
    unsigned char ephemeris_type;       ///< Orbital model used to generate data (usually 0)
    unsigned int element_set_number;    ///< Element set number
    unsigned char line_1_checksum;      ///< Line 1 check-sum

    // Line 2
    double inclination;                 ///< Inclination, 0 ≤ [deg] ≤ 180
    double raan;                        ///< Right ascension of the ascending node, 0 ≤ [deg] ≤ 360
    double eccentricity;                ///< Eccentricity (0 ≤ [] ≤ 1)
    double arg_of_perigee;              ///< Argument of perigee, 0 ≤ [deg] ≤ 360
    double mean_anomaly;                ///< Mean anomaly, 0 ≤ [deg] ≤ 360
    double mean_motion;                 ///< Mean motion, 0 < [rev/day]
    unsigned long revolution_number;    ///< Revolution number at epoch, 0 ≤ [rev] ≤ 99999
    unsigned char line_2_checksum;      ///< Line 2 check-sum
// clang-format on
};

enum perturb_tle_parse_error perturb_init_sat_from_tle(
    struct perturb_tle tle, enum perturb_grav_model grav_model, struct perturb_satellite * sat
);

#ifndef PERTURB_DISABLE_IO
enum perturb_tle_parse_error perturb_parse_tle(char * line_1, char * line_2, struct perturb_tle * tle);
#endif

#ifndef PERTURB_DISABLE_IO
enum perturb_tle_parse_error perturb_parse_tle_and_init_sat(
    char * line_1, char * line_2, enum perturb_grav_model grav_model, struct perturb_satellite * sat
);
#endif

#if __cplusplus
}  // extern "C"
#  ifdef PERTURB_ENABLE_CPP_INTERFACE
}  // namespace c_internal
}  // namespace perturb
#  endif
#endif

#endif  // PERTURB_TLE_H
