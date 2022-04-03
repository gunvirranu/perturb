/*
 * perturb -- A modern C++11 wrapper for the SGP4 orbit propagator
 * Version 0.0.0
 * https://github.com/gunvirranu/perturb
 *
 * Licensed under the MIT License <http://opensource.org/licenses/MIT>.
 * SPDX-License-Identifier: MIT
 *
 * Copyright (c) 2022 Gunvir Ranu
 */

//! @file
//! Primary header file for the perturb library
//! @author Gunvir Ranu
//! @version 0.0.0
//! @copyright Gunvir Ranu, MIT License

#ifndef PERTURB_PERTURB_HPP
#define PERTURB_PERTURB_HPP

#include "perturb/vallado_sgp4.hpp"

#include <array>
#include <cstddef>
#ifndef PERTURB_DISABLE_IO
#include <string>
#endif

/// Primary namespace for the perturb library, everything is in here.
///
/// See README for a brief intro to the main library types and basic usage.
namespace perturb {

/// Alias for representing position and velocity vectors
using Vec3 = std::array<double, 3>;

/// Both lines of a TLE **must** be this length, for TLE constructors.
///
/// The memory can and is accessed.
/// Lines can be longer for verification mode, but that's for internal testing
/// purposes only and doesn't pertain to general usage.
constexpr std::size_t TLE_LINE_LEN = 69;

/// Possible issues during SGP4 propagation or even TLE parsing.
///
/// This is important in two places:
///   1. After construction of a `Satellite`, check `Satellite::last_error`
///      for any issues with TLE parsing or SGP4 initialization.
///   2. Calling `Satellite::propagate` returns a `perturb::Sgp4Error`,
///      indicating any possible issues with propagation.
///
/// If everything is all good, the value should be `Sgp4Error::NONE`.
/// The errors `Sgp4Error::MEAN_ELEMENTS` to `SGP4Error::DECAYED` directly
/// correlate to errors in the underlying SGP4 impl, from the comments of
/// `vallado_sgp4::sgp4`. The additional `Sgp4Error::INVALID_TLE` is for issues
/// with reading the TLE strings.
enum class Sgp4Error : int {
    NONE = 0,
    MEAN_ELEMENTS,
    MEAN_MOTION,
    PERT_ELEMENTS,
    SEMI_LATUS_RECTUM,
    EPOCH_ELEMENTS_SUB_ORBITAL,
    DECAYED,
    INVALID_TLE,  // Not from base impl, added in
    UNKNOWN
};

/// Choice of gravity model / constants for the underlying SGP4 impl.
///
/// Corresponds to the `gravconsttype` type in `perturb::vallado_sgp4`.
/// Generally, WGS72 is the standard choice, despite WGS83 being the newer and
/// more accurate model. What is most important is that this is the exact same
/// as the gravity model used to generate the TLE ephemeris. This can be
/// confirmed from the source of your TLE data.
enum class GravModel {
    WGS72_OLD,
    WGS72,
    WGS84
};

/// A basic and human readable representation of a point in time.
///
/// The primary purpose of this type is to be constructed manually via
/// aggregate initialization and the converted to a `JulianDate`.
///
/// @warning
/// This type doesn't enforce a valid date and time upon construction; there
/// are no checks. Additionally, the conversion to `JulianDate` values are only
/// valid from years 1900 to 2100.
///
/// The question of what this time point represents gets complicated. In the
/// "Revisiting Spacetrack Report #3" paper from Celestrak, it is discussed
/// what this time represents. While UTC would make sense, the method
/// `vallado_sgp4::gstime_SGP4` requires UT1 time to calculate GMST. This
/// library makes the same assumption that the paper concludes:
///
/// @par
///  "The error associated with approximating UT1 with UTC is within the
///   theoretical uncertainty of the SGP4 theory itself. Except for the GMST
///   calculation, this paper and code assumes time to be realized as UTC."
struct YMDhms {
    /// Year from 1900 to 2100
    int year;
    /// Month from 1 to 12
    int month;
    /// Day from 1 to {28, 29, 30, 31} (depending on month)
    int day;
    /// Hour from 0 to 23
    int hour;
    /// Minute from 0 to 59
    int min;
    /// Fractional seconds from 0.0 to 59.999...
    double sec;
};

/// Represents a specific point in time on the Julian calendar.
///
/// Generally not constructed manually, but instead converted from a `YMDhms`
/// or from `Satellite::epoch`. Supports some basic manipulation operations.
/// For a human readable representation, can be converted back to `YMDhms`.
/// As for what time point this represents, see the comment on `YMDhms`.
///
/// Internally, this is represented as the "theoretical" sum of two double
/// precision floats. This is to preserve as much time accuracy as possible,
/// since many bits are lost due to storing the number of the day. The smaller
/// value is used to represent a more accurate time offset from that day. The
/// "normalized" value restricts the larger value to entire days and the
/// smaller to a [0.0, 1.0) time offset.
struct JulianDate {
    /// Fractional number of days since the epoch (4713 B.C.)
    double jd;
    /// Smaller fractional number of days
    double jd_frac;

    /// Construct from a Julian number of days since epoch (4713 B.C.)
    explicit JulianDate(double jd);

    /// Construct from a Julian day number and fractional day.
    ///
    /// The "true" Julian date value is the sum of the two. The motivation for
    /// separating into two parts is to preserve floating-point precision.
    ///
    /// @param jd Larger Julian day value since epoch
    /// @param jd_frac Smaller fractional Julian day value
    explicit JulianDate(double jd, double jd_frac);

    /// Construct from a `YMDhms` time point.
    ///
    /// @pre Only years from 1900 to 2100 are supported.
    /// @warning The date and time are assumed to be valid, with no checks.
    ///
    /// @param t Time point, must be from 1900 to 2100
    explicit JulianDate(YMDhms t);

    /// Convert to a `YMDhms` representing the same time point.
    ///
    /// @return Same time point converted to a human readable representation
    YMDhms to_datetime() const;

    /// Normalizes a Julian date representation to a canonical representation.
    ///
    /// Modifies the two julian day values to restric the larger value to whole
    /// days and the smaller value to a [0.0, 1.0) time offset. This is *not*
    /// needed for the conversions as they will automatically normalize, but
    /// you may want to do this for other reasons.
    void normalize();

    /// Get a normalized copy of a julian date
    JulianDate normalized() const;

    /// Returns the difference/delta in times as a fractional number of days
    double operator-(const JulianDate &rhs) const;

    /// Returns another time point offset by a number of days.
    ///
    /// Offset is *only* added to `jd_frac`, so the returned value may not be
    /// "normalized". This needs to be done explicitly if desired. This could
    /// split `delta_jd` into a whole and fractional part and then add each
    /// accordingly, but I defaulted to the former method for performance. If
    /// you think this wasn't the right move, lemme know.
    ///
    /// @param delta_jd Number of fractional days to add
    /// @return Sum of a time point and offset, *not* normalized
    JulianDate operator+(const double &delta_jd) const;

    /// Add a delta number of days offset to this time point
    JulianDate& operator+=(const double &delta_jd);

    /// Returns another time point offset backwards by a number of days
    JulianDate operator-(const double &delta_jd) const;

    /// Subtracts a delta number of days offset to this time point
    JulianDate& operator-=(const double &delta_jd);
};

/// Represents a specific orbital ephemeris for an Earth-centered trajectory.
///
/// This is the primary type in this library. Wraps the internal SGP4 record
/// type `vallado_sgp4::elsetrec`. Generally constructed via TLEs through
/// `Satellite::from_tle` constructors. Of particular importance is the
/// `Satellite::last_error` method which you can check to determine if there
/// were any issues with TLE initialization or propagation. The primary method
/// of running the SGP4 algorithm is the `Satellite::propagate` method.
class Satellite {
public:
    /// Internal SGP4 type
    vallado_sgp4::elsetrec sat_rec;

    /// Constructs from a raw SGP4 orbital record.
    ///
    /// @param sat_rec Pre-initialized SGP4 orbital record
    explicit Satellite(vallado_sgp4::elsetrec sat_rec);

#ifndef PERTURB_DISABLE_IO
    /// Construct and initialize a `Satellite` from a TLE record.
    ///
    /// The strings are mutable because the underlying implementation in
    /// `vallado_sgp4::twoline2rv` may modify the string during parsing.
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
        std::string &line_1,
        std::string &line_2,
        GravModel grav_model = GravModel::WGS72
    );
#endif  // PERTURB_ENABLE_IO

    /// Return the last recorded error in the internal SGP4 record
    Sgp4Error last_error() const;

    /// Return the epoch of the orbital ephemeris, likely from a TLE
    JulianDate epoch() const;

    /// Propagate the SGP4 model based on time around the epoch.
    ///
    /// @param mins_from_epoch Offset number of minutes around the epoch
    /// @param pos Returned position vector in km in the TEME frame
    /// @param vel Returned velocity vector in km/s in the TEME frame
    /// @return Issues during propagation, should usually be `Sgp4Error::NONE`
    Sgp4Error propagate_from_epoch(double mins_from_epoch, Vec3 &pos, Vec3 &vel);

    /// Propagate the SGP4 model to a specific time point.
    ///
    /// @param jd Time point in UTC or UT1
    /// @param pos Returned position vector in km in the TEME frame
    /// @param vel Returned velocity vector in km/s in the TEME frame
    /// @return Issues during propagation, should usually be `Sgp4Error::NONE`
    Sgp4Error propagate(JulianDate jd, Vec3 &pos, Vec3 &vel);
};
}  // namespace perturb

#endif  // PERTURB_PERTURB_HPP
