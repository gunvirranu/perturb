#include "perturb/perturb.hpp"

#include <cstring>

#include "perturb/vallado_sgp4.hpp"

namespace perturb {

static Sgp4Error convert_sgp4_error_code(const int error_code) {
    if (error_code < 0 || error_code >= static_cast<int>(Sgp4Error::UNKNOWN)) {
        return Sgp4Error::UNKNOWN;
    }
    return static_cast<Sgp4Error>(error_code);
}

JulianDate::JulianDate(double jd) : jd(jd) {}

JulianDate::JulianDate(double jd, double jd_frac) : jd(jd + jd_frac) {}

JulianDate::JulianDate(YMDhms t) {
    double tmp_jd, tmp_jd_frac;
    vallado_sgp4::jday_SGP4(t.year, t.month, t.day, t.hour, t.min, t.sec, tmp_jd, tmp_jd_frac);
    jd = tmp_jd + tmp_jd_frac;
}

YMDhms JulianDate::to_datetime() const {
    YMDhms t {};
    vallado_sgp4::invjday_SGP4(jd, 0.0, t.year, t.month, t.day, t.hour, t.min, t.sec);
    return t;
}

double JulianDate::operator-(const JulianDate &rhs) const {
    return this->jd - rhs.jd;
}

JulianDate JulianDate::operator+(const double &delta_jd) const {
    return JulianDate(jd, 0.0 + delta_jd);
}

JulianDate &JulianDate::operator+=(const double &delta_jd) {
    return *this = *this + delta_jd;
}

Satellite::Satellite(const vallado_sgp4::elsetrec sat_rec) : sat_rec(sat_rec) {}

Satellite Satellite::from_tle(char *line_1, char *line_2) {
    double _startmfe, _stopmfe, _deltamin;
    vallado_sgp4::elsetrec sat_rec {};
    const bool bad_ptrs = !line_1 || !line_2;
    if (bad_ptrs || std::strlen(line_1) < TLE_LINE_LEN || std::strlen(line_2) < TLE_LINE_LEN) {
        sat_rec.error = static_cast<int>(Sgp4Error::INVALID_TLE);
    } else {
        vallado_sgp4::twoline2rv(
            line_1, line_2, ' ', ' ', 'i',
            vallado_sgp4::wgs72, _startmfe, _stopmfe, _deltamin, sat_rec
        );
    }
    return Satellite(sat_rec);
}

Satellite Satellite::from_tle(std::string &line_1, std::string &line_2) {
    if (line_1.length() < TLE_LINE_LEN || line_2.length() < TLE_LINE_LEN) {
        vallado_sgp4::elsetrec sat_rec;
        sat_rec.error = static_cast<int>(Sgp4Error::INVALID_TLE);
        return Satellite(sat_rec);
    }
    return from_tle(&line_1[0], &line_2[0]);
}

Sgp4Error Satellite::last_error() const {
    return convert_sgp4_error_code(sat_rec.error);
}

JulianDate Satellite::epoch() const {
    return JulianDate(sat_rec.jdsatepoch, sat_rec.jdsatepochF);
}

Sgp4Error Satellite::propagate_from_epoch(double mins_from_epoch, Vec3 &pos, Vec3 &vel) {
    const bool is_valid = vallado_sgp4::sgp4(sat_rec, mins_from_epoch, pos.data(), vel.data());
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
