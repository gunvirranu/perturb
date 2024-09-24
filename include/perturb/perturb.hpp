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

//! @file Primary C++ header for the perturb library
//! @author Gunvir Singh Ranu
//! @version 1.0.0
//! @copyright Gunvir Singh Ranu, MIT License

#ifndef PERTURB_PERTURB_HPP
#define PERTURB_PERTURB_HPP

#ifndef __cplusplus
#  error "This header is intended for C++, use perturb.h for C"
#endif

#ifndef PERTURB_ENABLE_CPP_INTERFACE
#  error "Set the PERTURB_ENABLE_CPP_INTERFACE feature flag to support the C++ interface"
#endif

#include <array>
#include <cstddef>
#ifndef PERTURB_DISABLE_IO
#  include <string>
#endif

#include "perturb/perturb.h"
#include "perturb/tle.h"
#include "perturb/sgp4.h"

/// Primary namespace for the perturb C++ library wrapper, everything is in here.
///
/// See README for a brief intro to the main library types and basic usage.
/// The C library interface is hidden in the `perturb::c_internal` namespace.
namespace perturb {

/// Alias for general floating-point type (may be 32-bit or 64-bit)
using real_t = c_internal::perturb_real_t;

/// Alias for representing 3D position and velocity vectors in [km]
using Vec3 = std::array<real_t, 3>;

/// Both lines of a TLE **must** be this length, for TLE constructors.
///
/// It is assumed that this memory can be safely accessed.
/// Lines can be longer for verification mode, but that's for internal testing
/// purposes only and doesn't pertain to general usage.
constexpr std::size_t TLE_LINE_LEN = 69;

// TODO: Add a few static asserts for enums as a sanity check

enum class TLEParseError {
    NONE,               ///< If no issues when parsing
    SHOULD_BE_SPACE,    ///< If there is a lack of space in the TLE
    INVALID_FORMAT,     ///< If general parsing was unsuccessfully
    INVALID_VALUE,      ///< If a parsed value doesn't make sense
    CHECKSUM_MISMATCH,  ///< If the checksum doesn't match
};

enum class GravModel {
    WGS72_OLD,
    WGS72,
    WGS84,
};

/// Convert C++ `GravModel` enum to internal C enum
c_internal::perturb_grav_model convert_grav_model(const GravModel model);

enum class Sgp4Error {
    NONE,
    MEAN_ELEMENTS,
    MEAN_MOTION,
    PERT_ELEMENTS,
    SEMI_LATUS_RECTUM,
    EPOCH_ELEMENTS_SUB_ORBITAL,
    DECAYED,
    INVALID_TLE,
    UNKNOWN
};

struct DateTime {
    c_internal::perturb_date_time internal;  /// Internal C data

    // TODO: Add constructor anyways
};

struct JulianDate {
    c_internal::perturb_julian_date internal;  /// Internal C data

    explicit JulianDate(c_internal::perturb_julian_date in);

    /// Construct from a Julian number of days since epoch
    explicit JulianDate(real_t jd);

    /// Construct from a Julian day number and fractional day.
    ///
    /// The "true" Julian date value is the sum of the two. The motivation for
    /// separating into two parts is to preserve floating-point precision.
    ///
    /// @param jd Larger Julian [day] value since epoch
    /// @param jd_frac Smaller fractional Julian [day] value
    explicit JulianDate(real_t jd, real_t jd_frac);

    /// Construct from a `DateTime` time point.
    ///
    /// @pre Only years from 1900 to 2100 are supported.
    /// @warning The date and time are assumed to be valid, with no checks.
    ///
    /// @param t Time point, must be from 1900 to 2100
    explicit JulianDate(DateTime t);

    /// Convert to a `DateTime` representing the same time point.
    ///
    /// @return Same time point converted to a human readable representation
    DateTime to_datetime() const;

    /// Get a normalized / canonical julian date representation
    ///
    /// Larger value is restricted to a whole number of days and the smaller
    /// value to a [0.0, 1.0) time offset. This is *not* generally needed for
    /// conversions; they will automatically normalize.
    JulianDate normalized() const;

    /// Normalize julian date representation in-place.
    ///
    /// See `JulianDate::normalized` for details.
    void normalize();

    /// Returns the difference/delta in times as a fractional number of days
    real_t operator-(const JulianDate &rhs) const;

    /// Returns another time point offset by a number of days.
    ///
    /// Offset is *only* added to `jd_frac`, so the returned value may not be
    /// "normalized". This needs to be done explicitly if desired. This could
    /// split `delta_jd` into a whole and fractional part and then add each
    /// accordingly, but I defaulted to the former method for performance. If
    /// you think this wasn't the right move, lemme know.
    ///
    /// @param delta_jd Number of fractional [day] to add
    /// @return Sum of a time point and offset, *not* normalized
    JulianDate operator+(const real_t &delta_jd) const;
    /// Add a delta number of days offset to this time point
    JulianDate &operator+=(const real_t &delta_jd);
    /// Returns another time point offset backwards by a number of days
    JulianDate operator-(const real_t &delta_jd) const;
    /// Subtracts a delta number of days offset to this time point
    JulianDate &operator-=(const real_t &delta_jd);

    /// Compare if chronologically earlier than another time point
    bool operator<(const JulianDate &rhs) const;
    /// Compare if chronologically after than another time point
    bool operator>(const JulianDate &rhs) const;
    /// Compare if earlier than or same as another time point
    bool operator<=(const JulianDate &rhs) const;
    /// Compare if after than or same as another time point
    bool operator>=(const JulianDate &rhs) const;
};

struct StateVector {
    c_internal::perturb_state_vector internal;  /// Internal C data
};

struct ClassicalOrbitalElements {
    c_internal::perturb_classical_orbital_elements internal;  /// Internal C data

    /// Construct from a `StateVector` position and velocity in TEME.
    ///
    /// @param sv A position-velocity state vector generated via SGP4
    /// TODO: param grav_model Gravity model used in SGP4 (default `GravModel::WGS72`)
    explicit ClassicalOrbitalElements(StateVector sv, GravModel grav_model = GravModel::WGS72);
};

struct TwoLineElement {
    c_internal::perturb_tle internal;  /// Internal C data

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

class Satellite {
public:
    c_internal::perturb_satellite internal;  /// Internal C data

    /// Construct from a raw SGP4 orbital record.
    ///
    /// @param sat Pre-initialized SGP4 orbital record
    explicit Satellite(c_internal::perturb_satellite sat);

    /// Construct and initialize from a pre-parsed TLE record.
    ///
    /// Does not require any string parsing, so aight for embedded. Does the
    /// same initialization steps as `perturb::sgp4::twoline2rv` and then uses
    /// `perturb::sgp4::sgp4init` to fully initialize.
    ///
    /// You probably don't want to parse a TLE string into a `TwoLineElement`
    /// and then construct a `Satellite` via this. You should instead use the
    /// `Satellite::from_tle` constructor directly.
    ///
    /// This is only useful if you don't want *any* IO (via the `PERTURB_DISABLE_IO`
    /// flag), because then this and the `Satellite(sgp4::elsetrec)` constructor are the
    /// only way of constructing a `Satellite` object. It's up to you how you wanna
    /// construct the parsed `perturb::TwoLineElement` object.
    ///
    /// @pre The TLE must be valid and contain valid values.
    ///
    /// @param tle Parsed and valid TLE
    /// TODO: param grav_model Gravity constants to use (default `GravModel::WGS72`)
    explicit Satellite(
        const TwoLineElement &tle, GravModel grav_model = GravModel::WGS72
    );

#ifndef PERTURB_DISABLE_IO
    /// Construct and initialize a `Satellite` from a TLE record.
    ///
    /// The strings are mutable because the underlying implementation in
    /// `perturb::sgp4::twoline2rv` may modify the string during parsing.
    /// Left as mutable instead of internally copying for efficiency reasons as
    /// this may be okay for the caller.
    ///
    /// @param line_1 First line of TLE as C-string of length `perturb::TLE_LINE_LEN`
    /// @param line_2 Second line of TLE as C-string of length `perturb::TLE_LINE_LEN`
    /// @param grav_model Gravity constants to use (default `GravModel::WGS72`)
    /// @return An initialized `Satellite`
    static Satellite from_tle(
        char *line_1, char *line_2, GravModel grav_model = GravModel::WGS72
    );
#endif  // PERTURB_DISABLE_IO

#ifndef PERTURB_DISABLE_IO
    /// Wrapper for `Satellite::from_tle` that accepts C++ style strings.
    ///
    /// @param line_1 First line of TLE
    /// @param line_2 Second line of TLE
    /// @param grav_model Gravity constants to use (default `GravModel::WGS72`)
    /// @return An initialized `Satellite`
    static Satellite from_tle(
        std::string &line_1, std::string &line_2, GravModel grav_model = GravModel::WGS72
    );
#endif  // PERTURB_DISABLE_IO

    /// Return the last recorded error in the internal SGP4 record
    Sgp4Error last_error() const;

    /// Return the epoch of the orbital ephemeris, likely from a TLE
    JulianDate epoch() const;

    /// Propagate the SGP4 model based on time around the epoch.
    ///
    /// @param mins_from_epoch Offset number of minutes around the epoch
    /// @param posvel Returned state vector in the TEME frame
    /// @return Issues during propagation, should usually be `Sgp4Error::NONE`
    Sgp4Error propagate_from_epoch(real_t mins_from_epoch, StateVector &sv);

    /// Propagate the SGP4 model to a specific time point.
    ///
    /// @param jd Time point in UTC or UT1
    /// @param posvel Returned state vector in the TEME frame
    /// @return Issues during propagation, should usually be `Sgp4Error::NONE`
    Sgp4Error propagate(JulianDate jd, StateVector &sv);
};

}  // namespace perturb

#endif  // PERTURB_PERTURB_HPP
