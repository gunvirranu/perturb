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

#include "perturb/perturb.h"

#ifdef __cplusplus
#  ifdef PERTURB_ENABLE_CPP_INTERFACE
namespace perturb {
namespace c_internal {
#  endif
extern "C" {
#endif

#include <math.h>

#include "perturb/sgp4.h"
#include "common_private.h"

struct perturb_julian_date perturb_datetime_to_julian(const struct perturb_date_time t) {
    struct perturb_julian_date jd;
    jday_SGP4(t.year, t.month, t.day, t.hour, t.min, t.sec, &jd.jd, &jd.jd_frac);
    return jd;
}

struct perturb_date_time perturb_julian_to_datetime(struct perturb_julian_date jd) {
    struct perturb_date_time t;
    invjday_SGP4(jd.jd, jd.jd_frac, &t.year, &t.month, &t.day, &t.hour, &t.min, &t.sec);
    return t;
}

struct perturb_julian_date perturb_julian_normalized(const struct perturb_julian_date t) {
    struct perturb_julian_date out = t;

    // Check for fractional days included in `jd` and put them in `jd`
    const real_t frac_days = t.jd - floor(t.jd) - 0.5;

    if (FABS(frac_days) > 1e-12) {
        out.jd -= frac_days;
        out.jd_frac += frac_days;
    }

    // Check for whole days in `jd_frac` and put them in `jd`
    if (FABS(out.jd_frac) >= 1.0) {
        const double whole_days = FLOOR(out.jd_frac);
        out.jd += whole_days;
        out.jd_frac -= whole_days;
    }
    return out;
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

ClassicalOrbitalElements::ClassicalOrbitalElements(
    StateVector sv, GravModel grav_model
) {
    double mus, _tumin, _rekm, _xke, _j2, _j3, _j4, _j3oj2;
    sgp4::getgravconst(
        convert_grav_model(grav_model), _tumin, mus, _rekm, _xke, _j2, _j3, _j4, _j3oj2
    );
    sgp4::rv2coe_SGP4(
        sv.position.data(), sv.velocity.data(), mus, semilatus_rectum, semimajor_axis,
        eccentricity, inclination, raan, arg_of_perigee, true_anomaly, mean_anomaly,
        arg_of_latitude, true_longitude, longitude_of_periapsis
    );
}

Satellite::Satellite(const sgp4::elsetrec _sat_rec) : sat_rec(_sat_rec) {}

Satellite::Satellite(const TwoLineElement &tle, GravModel grav_model) : sat_rec({}) {
    constexpr double DEG_TO_RAD = PI / 180.0;
    constexpr double XP_DOT_P = 1440.0 / (2 * PI);

    // Line 1
    std::memcpy(sat_rec.satnum, tle.catalog_number, sizeof(sat_rec.satnum));
    sat_rec.classification = tle.classification;
    // Don't bother converting to set `sat_rec.intldesg` b/c it has no effects
    sat_rec.epochyr = static_cast<int>(tle.epoch_year);
    sat_rec.epochdays = tle.epoch_day_of_year;
    sat_rec.ndot = tle.n_dot;
    sat_rec.nddot = tle.n_ddot;
    sat_rec.bstar = tle.b_star;
    sat_rec.ephtype = tle.ephemeris_type;
    sat_rec.elnum = tle.element_set_number;

    // Line 2
    sat_rec.inclo = tle.inclination;
    sat_rec.nodeo = tle.raan;
    sat_rec.ecco = tle.eccentricity;
    sat_rec.argpo = tle.arg_of_perigee;
    sat_rec.mo = tle.mean_anomaly;
    sat_rec.no_kozai = tle.mean_motion;
    sat_rec.revnum = static_cast<long>(tle.revolution_number);

    // Post-process same as how Vallado does it
    sat_rec.error = 0;
    sat_rec.no_kozai /= XP_DOT_P;  // [rad / min]
    sat_rec.ndot /= (XP_DOT_P * 1440.0);
    sat_rec.nddot /= (XP_DOT_P * 1440.0 * 1440);
    sat_rec.inclo *= DEG_TO_RAD;
    sat_rec.nodeo *= DEG_TO_RAD;
    sat_rec.argpo *= DEG_TO_RAD;
    sat_rec.mo *= DEG_TO_RAD;

    DateTime t {};
    t.year = sat_rec.epochyr + ((sat_rec.epochyr < 57) ? 2000 : 1900);
    sgp4::days2mdhms_SGP4(
        t.year, sat_rec.epochdays, t.month, t.day, t.hour, t.min, t.sec
    );
    const auto jd = JulianDate(t);
    sat_rec.jdsatepoch = jd.jd;
    sat_rec.jdsatepochF = jd.jd_frac;
    const double epoch = (jd.jd + jd.jd_frac) - 2433281.5;

    // Initialize orbit
    sgp4::sgp4init(
        convert_grav_model(grav_model), 'i', sat_rec.satnum, epoch, sat_rec.bstar,
        sat_rec.ndot, sat_rec.nddot, sat_rec.ecco, sat_rec.argpo, sat_rec.inclo,
        sat_rec.mo, sat_rec.no_kozai, sat_rec.nodeo, sat_rec
    );
}

#ifndef PERTURB_DISABLE_IO
Satellite Satellite::from_tle(char *line_1, char *line_2, GravModel grav_model) {
    sgp4::elsetrec sat_rec {};
    const bool bad_ptrs = !line_1 || !line_2;
    // FIXME: Remove `strlen` and just check last byte
    if (bad_ptrs || std::strlen(line_1) < TLE_LINE_LEN
        || std::strlen(line_2) < TLE_LINE_LEN) {
        sat_rec.error = static_cast<int>(Sgp4Error::INVALID_TLE);
    } else {
        // FIXME: Change default TLE to own parser, check downstream usage for assumptions
        double _startmfe, _stopmfe, _deltamin;
        sgp4::twoline2rv(
            line_1, line_2, ' ', ' ', 'i', convert_grav_model(grav_model), _startmfe,
            _stopmfe, _deltamin, sat_rec
        );
    }
    return Satellite(sat_rec);
}
#endif  // PERTURB_DISABLE_IO

JulianDate Satellite::epoch() const {
    return JulianDate(sat_rec.jdsatepoch, sat_rec.jdsatepochF);
}

Sgp4Error Satellite::propagate_from_epoch(double mins_from_epoch, StateVector &sv) {
    sv.epoch = epoch() + (mins_from_epoch / MINS_PER_DAY);
    const bool is_valid =
        sgp4::sgp4(sat_rec, mins_from_epoch, sv.position.data(), sv.velocity.data());
    (void) is_valid;  // Unused because it is consistent with error code
    return last_error();
}

Sgp4Error Satellite::propagate(const JulianDate jd, StateVector &sv) {
    const double delta_jd = jd - epoch();
    const double mins_from_epoch = delta_jd * MINS_PER_DAY;
    const auto err = propagate_from_epoch(mins_from_epoch, sv);
    sv.epoch = jd;  // Can save some math, ignore value from `propagate_from_epoch`
    return err;
}
}  // namespace perturb
