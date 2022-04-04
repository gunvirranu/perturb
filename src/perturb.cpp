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

#include "perturb/perturb.hpp"

#include <cmath>
#ifndef PERTURB_DISABLE_IO
#  include <cstring>
#endif

#include "perturb/vallado_sgp4.hpp"

namespace perturb {

static Sgp4Error convert_sgp4_error_code(const int error_code) {
    if (error_code < 0 || error_code >= static_cast<int>(Sgp4Error::UNKNOWN)) {
        return Sgp4Error::UNKNOWN;
    }
    return static_cast<Sgp4Error>(error_code);
}

static vallado_sgp4::gravconsttype convert_grav_model(const GravModel model) {
    switch (model) {
        case GravModel::WGS72_OLD: return vallado_sgp4::wgs72old;
        case GravModel::WGS72: return vallado_sgp4::wgs72;
        case GravModel::WGS84: return vallado_sgp4::wgs84;
        default: return vallado_sgp4::wgs72;
    }
}

JulianDate::JulianDate(double _jd) : jd(_jd), jd_frac(0) {}

JulianDate::JulianDate(double _jd, double _jd_frac) : jd(_jd), jd_frac(_jd_frac) {}

JulianDate::JulianDate(YMDhms t) {
    double tmp_jd, tmp_jd_frac;
    vallado_sgp4::jday_SGP4(
        t.year, t.month, t.day, t.hour, t.min, t.sec, tmp_jd, tmp_jd_frac
    );
    jd = tmp_jd;
    jd_frac = tmp_jd_frac;
}

YMDhms JulianDate::to_datetime() const {
    YMDhms t {};
    vallado_sgp4::invjday_SGP4(
        jd, jd_frac, t.year, t.month, t.day, t.hour, t.min, t.sec
    );
    return t;
}

void JulianDate::normalize() {
    // Check for fractional days included in `jd` and put them in `jd`
    const double frac_days = jd - std::floor(jd) - 0.5;
    if (std::abs(frac_days) > 1e-12) {
        jd -= frac_days;
        jd_frac += frac_days;
    }
    // Check for whole days in `jd_frac` and put them in `jd`
    if (std::abs(jd_frac) >= 1.0) {
        const double whole_days = std::floor(jd_frac);
        jd += whole_days;
        jd_frac -= whole_days;
    }
}

JulianDate JulianDate::normalized() const {
    JulianDate other = *this;
    other.normalize();
    return other;
}

double JulianDate::operator-(const JulianDate &rhs) const {
    // Grouping here is very important to preserve precision
    return (this->jd - rhs.jd) + (this->jd_frac - rhs.jd_frac);
}

JulianDate JulianDate::operator+(const double &delta_jd) const {
    // Just add entire offset to fractional value
    // Can be normalized later explicitly if needed
    return JulianDate(jd, jd_frac + delta_jd);
}

JulianDate &JulianDate::operator+=(const double &delta_jd) {
    return *this = *this + delta_jd;
}

JulianDate JulianDate::operator-(const double &delta_jd) const {
    return *this + (-delta_jd);
}

JulianDate &JulianDate::operator-=(const double &delta_jd) {
    return *this = *this - delta_jd;
}

Satellite::Satellite(const vallado_sgp4::elsetrec _sat_rec) : sat_rec(_sat_rec) {}

#ifndef PERTURB_DISABLE_IO
Satellite Satellite::from_tle(char *line_1, char *line_2, GravModel grav_model) {
    vallado_sgp4::elsetrec sat_rec {};
    const bool bad_ptrs = !line_1 || !line_2;
    if (bad_ptrs || std::strlen(line_1) < TLE_LINE_LEN
        || std::strlen(line_2) < TLE_LINE_LEN) {
        sat_rec.error = static_cast<int>(Sgp4Error::INVALID_TLE);
    } else {
        double _startmfe, _stopmfe, _deltamin;
        vallado_sgp4::twoline2rv(
            line_1, line_2, ' ', ' ', 'i', convert_grav_model(grav_model), _startmfe,
            _stopmfe, _deltamin, sat_rec
        );
    }
    return Satellite(sat_rec);
}
#endif  // PERTURB_DISABLE_IO

#ifndef PERTURB_DISABLE_IO
Satellite Satellite::from_tle(
    std::string &line_1, std::string &line_2, GravModel grav_model
) {
    if (line_1.length() < TLE_LINE_LEN || line_2.length() < TLE_LINE_LEN) {
        return from_tle(nullptr, nullptr);
    }
    return from_tle(&line_1[0], &line_2[0], grav_model);
}
#endif  // PERTURB_DISABLE_IO

Sgp4Error Satellite::last_error() const {
    return convert_sgp4_error_code(sat_rec.error);
}

JulianDate Satellite::epoch() const {
    return JulianDate(sat_rec.jdsatepoch, sat_rec.jdsatepochF);
}

Sgp4Error Satellite::propagate_from_epoch(double mins_from_epoch, Vec3 &pos, Vec3 &vel) {
    const bool is_valid =
        vallado_sgp4::sgp4(sat_rec, mins_from_epoch, pos.data(), vel.data());
    (void) is_valid;  // Unused because it is consistent with error code
    return last_error();
}

Sgp4Error Satellite::propagate(const JulianDate jd, Vec3 &pos, Vec3 &vel) {
    constexpr double MINS_PER_DAY = 24 * 60;
    const double delta_jd = jd - epoch();
    const double mins_from_epoch = delta_jd * MINS_PER_DAY;
    return propagate_from_epoch(mins_from_epoch, pos, vel);
}
}  // namespace perturb
