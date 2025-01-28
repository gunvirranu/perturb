/*
 * perturb -- A modern C++11 wrapper for the SGP4 orbit propagator
 * version 1.0.0
 * https://github.com/gunvirranu/perturb
 *
 * Copyright (c) 2022 Gunvir Singh Ranu
 * SPDX-License-Identifier: MIT
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

#ifdef __cplusplus
#  ifdef PERTURB_ENABLE_CPP_INTERFACE
namespace perturb {
namespace c_internal {
#  endif
extern "C" {
#endif

#if !(defined(__cplusplus) && defined(PERTURB_ENABLE_CPP_INTERFACE))
/// Both lines of a TLE **must** be this length, for TLE constructors.
///
/// It is assumed that this memory can be safely accessed.
/// Lines can be longer for verification mode, but that's for internal testing
/// purposes only and doesn't pertain to general usage.
///
/// Macro definition is excluded in C++ cases to not pollute globals.
/// The C++ library redefines this as a proper type.
#  define PERTURB_TLE_LINE_LEN  69U
#endif

/// Possible errors when parsing a TLE.
///
/// Returned by `TwoLineElement::parse` after processing a TLE record string.
///
/// @post
/// Errors are guaranteed to occur in reverse definition order, in that higher
/// enum values are checked first. Meaning invalid input is checked first, then
/// spaces, and so on. This allows you to assume, for example, that if there is
/// a `CHECKSUM_MISMATCH`, none of the previous errors occured first.
enum perturb_TleParseError {
    PERTURB_TLE_PARSE_ERROR_NONE,               ///< If no issues when parsing
    PERTURB_TLE_PARSE_ERROR_CHECKSUM_MISMATCH,  ///< If the checksum doesn't match
    PERTURB_TLE_PARSE_ERROR_INVALID_VALUE,      ///< If a parsed value doesn't make sense
    PERTURB_TLE_PARSE_ERROR_INVALID_FORMAT,     ///< If general parsing was unsuccessfully
    PERTURB_TLE_PARSE_ERROR_SHOULD_BE_SPACE,    ///< If there is a lack of space in the TLE
    PERTURB_TLE_PARSE_ERROR_INVALID_INPUT,      ///< If any inputs are null pointers
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
struct perturb_TwoLineElement {
    // These are ordered same as a TLE in string format
    // clang-format off

    // Line 1 - Metadata
    char  catalog_number[6];            ///< Satellite catalog number as an unparsed string
    char  classification;               ///< Classification {U: Unclassified, C: Classified, S: Secret}

    // Line 1 - Launch
    uint8_t   launch_year;              ///< International Designator - Launch year (last two digits)
    uint16_t  launch_number;            ///< International Designator - Launch number of the year
    char      launch_piece[4];          ///< International Designator - Piece of launch

    // Line 1 - Epoch Time
    uint8_t         epoch_year;         ///< Two-digit epoch year in [1957, 2056] [year]
    perturb_real_t  epoch_day_of_year;  ///< Epoch fractional day of year in [0, 366] [day]

    // Line 1 - Trajectory
    perturb_real_t  n_dot;              ///< First derivative of mean motion (ballistic coefficient)
    perturb_real_t  n_ddot;             ///< Second derivative of mean motion [rev/day^3]
    perturb_real_t  b_star;             ///< B* radiation pressure coefficient [1 / (earth radii)]

    // Line 1 - Metadata
    uint8_t   ephemeris_type;           ///< Orbital model used to generate data (usually 0)
    uint16_t  element_set_number;       ///< Element set number
    uint8_t   line_1_checksum;          ///< Line 1 check-sum

    // Line 2 - Orbit
    perturb_real_t  inclination;        ///< Inclination, 0 ≤ [deg] ≤ 180
    perturb_real_t  raan;               ///< Right ascension of the ascending node, 0 ≤ [deg] ≤ 360
    perturb_real_t  eccentricity;       ///< Eccentricity (0 ≤ [] ≤ 1)
    perturb_real_t  arg_of_perigee;     ///< Argument of perigee, 0 ≤ [deg] ≤ 360
    perturb_real_t  mean_anomaly;       ///< Mean anomaly, 0 ≤ [deg] ≤ 360
    perturb_real_t  mean_motion;        ///< Mean motion, 0 < [rev/day]

    // Line 2 - Metadata
    uint32_t  revolution_number;        ///< Revolution number at epoch, 0 ≤ [rev] ≤ 99999
    uint8_t   line_2_checksum;          ///< Line 2 check-sum

    // clang-format on
};

enum perturb_Sgp4Error perturb_init_sat_from_tle(
    struct perturb_TwoLineElement tle,
    enum perturb_GravityModel grav_model,
    struct perturb_Satellite * sat
);

#ifndef PERTURB_DISABLE_IO
enum perturb_TleParseError perturb_parse_tle(
    const char * line_1, const char * line_2,
    struct perturb_TwoLineElement * tle
);
#endif

#ifdef __cplusplus
}  // extern "C"
#  ifdef PERTURB_ENABLE_CPP_INTERFACE
}  // namespace c_internal
}  // namespace perturb
#  endif
#endif

#endif  // PERTURB_TLE_H
