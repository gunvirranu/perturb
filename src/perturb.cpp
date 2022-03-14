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

Satellite::Satellite(const vallado_sgp4::elsetrec sat_rec) : sat_rec(sat_rec) {}

Satellite Satellite::from_tle(
    const std::array<char, TLE_LINE_LEN> &line_1,
    const std::array<char, TLE_LINE_LEN> &line_2
) {
    char line_buf_1[130], line_buf_2[130];
    std::memcpy(line_buf_1, line_1.data(), TLE_LINE_LEN);
    std::memcpy(line_buf_2, line_2.data(), TLE_LINE_LEN);
    double _startmfe, _stopmfe, _deltamin;
    vallado_sgp4::elsetrec sat_rec;
    vallado_sgp4::twoline2rv(
        line_buf_1, line_buf_2, ' ', ' ', 'i',
        vallado_sgp4::wgs72, _startmfe, _stopmfe, _deltamin, sat_rec
    );
    return Satellite(sat_rec);
}

Sgp4Error Satellite::last_error() const {
    return convert_sgp4_error_code(sat_rec.error);
}


JulianDate Satellite::epoch() const {
    return JulianDate(sat_rec.jdsatepoch, sat_rec.jdsatepochF);
}

Sgp4Error Satellite::propogate_from_epoch(double mins_from_epoch, Vec3 &pos, Vec3 &vel) {
    const bool is_valid = vallado_sgp4::sgp4(sat_rec, mins_from_epoch, pos.data(), vel.data());
    (void) is_valid;  // Unused because it is consistent with error code
    return last_error();
}

Sgp4Error Satellite::propogate(const JulianDate jd, Vec3 &pos, Vec3 &vel) {
    constexpr double MINS_PER_DAY = 24 * 60;
    const double delta_jd = jd - epoch();
    const double mins_from_epoch = delta_jd * MINS_PER_DAY;
    return propogate_from_epoch(mins_from_epoch, pos, vel);
}
}  // namespace perturb
