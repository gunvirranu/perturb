/**
 * perturb -- A modern C++11 wrapper for the SGP4 orbit propagator
 * Version 0.0.0
 * https://github.com/gunvirranu/perturb
 *
 * Licensed under the MIT License <http://opensource.org/licenses/MIT>.
 * SPDX-License-Identifier: MIT
 *
 * Copyright (c) 2022 Gunvir Ranu
 */

#ifndef PERTURB_PERTURB_HPP
#define PERTURB_PERTURB_HPP

// Defining preprocessor flag `PERTURB_DISABLE_IO` across the entire lib.
// Setting it removes all I/O and string related functionality.
// This may be desired if targeting for embedded. Majors impact are:
//   1. Removes TLE constructor since it relies on C-style strings and sscanf
//   2. Removes std::string helper constructor as well

#include "perturb/vallado_sgp4.hpp"

#include <array>
#include <cstddef>
#ifndef PERTURB_DISABLE_IO
#include <string>
#endif

namespace perturb {

using Vec3 = std::array<double, 3>;

// The two lines of the TLE *must* be this length as the memory is accessed.
// TLE lines are only greater for verification mode, but that's internal.
constexpr std::size_t TLE_LINE_LEN = 69;

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

enum class GravModel {
    WGS72_OLD,
    WGS72,
    WGS84
};

struct YMDhms {
    int year, month, day, hour, min;
    double sec;
};

struct JulianDate {
    double jd;

    explicit JulianDate(double jd);

    explicit JulianDate(double jd, double jd_frac);

    explicit JulianDate(YMDhms t);

    YMDhms to_datetime() const;

    double operator-(const JulianDate &rhs) const;

    JulianDate operator+(const double &delta_jd) const;

    JulianDate& operator+=(const double &delta_jd);
};

class Satellite {
public:
    vallado_sgp4::elsetrec sat_rec;

    explicit Satellite(vallado_sgp4::elsetrec sat_rec);

#ifndef PERTURB_DISABLE_IO
    static Satellite from_tle(
        char *line_1, char *line_2, GravModel grav_model = GravModel::WGS72
    );

    static Satellite from_tle(
        std::string &line_1,
        std::string &line_2,
        GravModel grav_model = GravModel::WGS72
    );
#endif  // PERTURB_ENABLE_IO

    Sgp4Error last_error() const;

    JulianDate epoch() const;

    Sgp4Error propagate_from_epoch(double mins_from_epoch, Vec3 &pos, Vec3 &vel);

    Sgp4Error propagate(JulianDate jd, Vec3 &pos, Vec3 &vel);
};
}  // namespace perturb

#endif  // PERTURB_PERTURB_HPP
