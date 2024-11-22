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

#ifndef __cplusplus
// There's a small chance you want the C version of perturb, but you glob compiled
// the entire src/ directory. Just ignore the .cpp files and you're golden.
#  error "Ok someone messed up haha, u may be accidentally compiling this C++ as C"
#endif

#include "perturb/perturb.hpp"

namespace perturb {

// Sanity check enums matches
#define CHECK_ENUM_MATCHES(A, B)  static_assert(static_cast<int>(A) == c_internal::PERTURB_ ## B, "Enum mismatch!")

CHECK_ENUM_MATCHES(TLEParseError::NONE, TLE_PARSE_ERROR_NONE);
CHECK_ENUM_MATCHES(TLEParseError::SHOULD_BE_SPACE, TLE_PARSE_ERROR_SHOULD_BE_SPACE);
CHECK_ENUM_MATCHES(TLEParseError::INVALID_FORMAT, TLE_PARSE_ERROR_INVALID_FORMAT);
CHECK_ENUM_MATCHES(TLEParseError::INVALID_VALUE, TLE_PARSE_ERROR_INVALID_VALUE);
CHECK_ENUM_MATCHES(TLEParseError::CHECKSUM_MISMATCH, TLE_PARSE_ERROR_CHECKSUM_MISMATCH);

CHECK_ENUM_MATCHES(GravModel::WGS72_OLD, GRAVITY_MODEL_WGS72_OLD);
CHECK_ENUM_MATCHES(GravModel::WGS72, GRAVITY_MODEL_WGS72);
CHECK_ENUM_MATCHES(GravModel::WGS84, GRAVITY_MODEL_WGS84);

CHECK_ENUM_MATCHES(Sgp4Error::NONE, SGP4_ERROR_NONE);
CHECK_ENUM_MATCHES(Sgp4Error::MEAN_ELEMENTS, SGP4_ERROR_MEAN_ELEMENTS);
CHECK_ENUM_MATCHES(Sgp4Error::MEAN_MOTION, SGP4_ERROR_MEAN_MOTION);
CHECK_ENUM_MATCHES(Sgp4Error::PERT_ELEMENTS, SGP4_ERROR_PERT_ELEMENTS);
CHECK_ENUM_MATCHES(Sgp4Error::SEMI_LATUS_RECTUM, SGP4_ERROR_SEMI_LATUS_RECTUM);
CHECK_ENUM_MATCHES(Sgp4Error::EPOCH_ELEMENTS_SUB_ORBITAL, SGP4_ERROR_EPOCH_ELEMENTS_SUB_ORBITAL);
CHECK_ENUM_MATCHES(Sgp4Error::DECAYED, SGP4_ERROR_DECAYED);
CHECK_ENUM_MATCHES(Sgp4Error::INVALID_TLE, SGP4_ERROR_INVALID_TLE);
CHECK_ENUM_MATCHES(Sgp4Error::UNKNOWN, SGP4_ERROR_UNKNOWN);

static Sgp4Error convert_sgp4_error_code(const int error_code) {
    if (error_code < 0 || error_code >= static_cast<int>(Sgp4Error::UNKNOWN)) {
        return Sgp4Error::UNKNOWN;
    }
    return static_cast<Sgp4Error>(error_code);
}

c_internal::perturb_GravityModel convert_grav_model(const GravModel model) {
    switch (model) {
        case GravModel::WGS72_OLD:  return c_internal::PERTURB_GRAVITY_MODEL_WGS72_OLD;
        case GravModel::WGS72:      return c_internal::PERTURB_GRAVITY_MODEL_WGS72;
        case GravModel::WGS84:      return c_internal::PERTURB_GRAVITY_MODEL_WGS84;
        default:                    return c_internal::PERTURB_GRAVITY_MODEL_WGS72;
    }
}

// TODO: DateTime methods

// JulianDate methods

JulianDate::JulianDate(const c_internal::perturb_JulianDate in) : internal(in) {}

JulianDate::JulianDate(const real_t jd) : internal({ jd, 0 }) {}

JulianDate::JulianDate(const real_t jd, const real_t jd_frac) : internal({ jd, jd_frac }) {}

JulianDate::JulianDate(const DateTime t) {
    this->internal = c_internal::perturb_datetime_to_julian(t.internal);
}

DateTime JulianDate::to_datetime() const {
    return DateTime { perturb_julian_to_datetime(this->internal) };
}

void JulianDate::normalize() {
    internal = c_internal::perturb_julian_normalized(internal);
}

JulianDate JulianDate::normalized() const {
    return JulianDate { c_internal::perturb_julian_normalized(internal) };
}

real_t JulianDate::operator-(const JulianDate &rhs) const {
    return c_internal::perturb_julian_subtract(this->internal, rhs.internal);
}

JulianDate JulianDate::operator+(const real_t &delta_jd) const {
    return JulianDate { c_internal::perturb_julian_add_days(internal, delta_jd) };
}

JulianDate &JulianDate::operator+=(const real_t &delta_jd) {
    return *this = *this + delta_jd;
}

JulianDate JulianDate::operator-(const real_t &delta_jd) const {
    return *this + (-delta_jd);
}

JulianDate &JulianDate::operator-=(const real_t &delta_jd) {
    return *this = *this - delta_jd;
}

bool JulianDate::operator<(const JulianDate &rhs) const {
    return (*this - rhs) < 0;
}

bool JulianDate::operator>(const JulianDate &rhs) const {
    return (*this - rhs) > 0;
}

bool JulianDate::operator<=(const JulianDate &rhs) const {
    return (*this - rhs) <= 0;
}

bool JulianDate::operator>=(const JulianDate &rhs) const {
    return (*this - rhs) >= 0;
}

// ClassicalOrbitalElements methods

ClassicalOrbitalElements::ClassicalOrbitalElements(StateVector sv, GravModel grav_model) {
    internal = c_internal::perturb_state_vector_to_orbital_elements_with_grav(
        sv.internal, convert_grav_model(grav_model)
    );
}

// TwoLineElement methods

#ifndef PERTURB_DISABLE_IO
TLEParseError TwoLineElement::parse(const char * line_1, const char * line_2) {
    // Pass parsing call directly to underlying C impl
    c_internal::perturb_Tle * tle = &this->internal;
    const c_internal::perturb_TleParseError err = c_internal::perturb_parse_tle(line_1, line_2,  tle);
    return static_cast<TLEParseError>(err);
}
#endif  // PERTURB_DISABLE_IO

#ifndef PERTURB_DISABLE_IO
TLEParseError TwoLineElement::parse(
    const std::string &line_1, const std::string &line_2
) {
    // Delegate to C-sting parsing, which assumes a minimum length
    if (line_1.length() < TLE_LINE_LEN || line_2.length() < TLE_LINE_LEN) {
        return TLEParseError::INVALID_FORMAT;
    }
    return this->parse(line_1.c_str(), line_2.c_str());
}
#endif  // PERTURB_DISABLE_IO

// Satellite methods

Satellite::Satellite(c_internal::perturb_Satellite sat) : internal(sat) {}

Satellite::Satellite(const TwoLineElement &tle, const GravModel grav_model) {
    const auto _err = c_internal::perturb_init_sat_from_tle(
        tle.internal, convert_grav_model(grav_model), &internal
    );
    static_cast<void>(_err);  // Unused, accessible via `last_error`
}

#ifndef PERTURB_DISABLE_IO
Satellite Satellite::from_tle(char * line_1, char * line_2, GravModel grav_model) {
    c_internal::perturb_Satellite sat {};
    c_internal::perturb_parse_tle_and_init_sat(line_1, line_2, convert_grav_model(grav_model), &sat);
    return Satellite(sat);
}
#endif  // PERTURB_DISABLE_IO

#ifndef PERTURB_DISABLE_IO
Satellite Satellite::from_tle(
    std::string &line_1, std::string &line_2, GravModel grav_model
) {
    if (line_1.length() < TLE_LINE_LEN || line_2.length() < TLE_LINE_LEN) {
        return from_tle(nullptr, nullptr);
    }
    // FIXME: Find a way to remove usage of &str[0]
    return from_tle(&line_1[0], &line_2[0], grav_model);
}
#endif  // PERTURB_DISABLE_IO

Sgp4Error Satellite::last_error() const {
    return convert_sgp4_error_code(internal.error);
}

JulianDate Satellite::epoch() const {
    return JulianDate(c_internal::perturb_epoch(internal));
}

Sgp4Error Satellite::propagate_from_epoch(const real_t mins_from_epoch, StateVector &sv) {
    const real_t days_from_epoch = (mins_from_epoch / (24 * 60));
    const auto err = c_internal::perturb_propagate_days_from_epoch(internal, days_from_epoch, &sv.internal);
    return convert_sgp4_error_code(err);
}

Sgp4Error Satellite::propagate(const JulianDate jd, StateVector &sv) {
    const auto err = c_internal::perturb_propagate(internal, jd.internal, &sv.internal);
    return convert_sgp4_error_code(err);
}

}  // namespace perturb
