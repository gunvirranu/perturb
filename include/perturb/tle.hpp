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

#ifndef PERTURB_TLE_HPP
#define PERTURB_TLE_HPP

#include <cstddef>
#ifndef PERTURB_DISABLE_IO
#  include <string>
#endif

namespace perturb {

/// Both lines of a TLE **must** be this length, for TLE constructors.
///
/// The memory can and is accessed.
/// Lines can be longer for verification mode, but that's for internal testing
/// purposes only and doesn't pertain to general usage.
constexpr std::size_t TLE_LINE_LEN = 69;

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
enum class TLEParseError {
    NONE,               ///< If no issues when parsing
    SHOULD_BE_SPACE,    ///< If there is a lack of space in the TLE
    INVALID_FORMAT,     ///< If general parsing was unsuccessfully
    INVALID_VALUE,      ///< If a parsed value doesn't make sense
    CHECKSUM_MISMATCH,  ///< If the checksum doesn't match
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
struct TwoLineElement {
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

#ifndef PERTURB_DISABLE_IO
    /// Parse a TLE record string.
    ///
    /// You *probably* don't need this method. As I explain in the `TwoLineElement`
    /// and and `Satellite()` constructor docs, you should probably just use the
    /// `Satellite::from_tle` method directly. This method must be called on an
    /// existing `TwoLineElement` variable, so it can return the error code.
    ///
    /// This currently uses my own implementation of a parser that doesn't
    /// support every case that Vallado's impl does, so there may be the
    /// occasional false error.
    ///
    /// @post See the `perturb::TLEParseError` docs for the guaranteed error ordering.
    ///
    /// @param line_1 First line of TLE as C-string of length `perturb::TLE_LINE_LEN`
    /// @param line_2 Second line of TLE as C-string of length `perturb::TLE_LINE_LEN`
    /// @return Issues with parsing, should usually be `TLEParseError::NONE`.
    ///         The parsed values are written into the `TwoLineElement` instance.
    TLEParseError parse(const char *line_1, const char *line_2);
#endif  // PERTURB_DISABLE_IO

#ifndef PERTURB_DISABLE_IO
    /// Wrapper for `TwoLineElement::parse` that accepts C++ style strings.
    ///
    /// @param line_1 First line of TLE
    /// @param line_2 Second line of TLE
    /// @return Issues with parsing, should usually be `TLEParseError::NONE`
    TLEParseError parse(const std::string &line_1, const std::string &line_2);
#endif  // PERTURB_DISABLE_IO
};

}  // namespace perturb

#endif  // PERTURB_TLE_HPP
