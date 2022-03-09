#include "perturb/perturb.hpp"

#include "perturb/vallado_sgp4.hpp"

namespace perturb {

Sgp4Error convert_sgp4_error_code(const int error_code) {
    if (error_code < 0 || error_code >= static_cast<int>(Sgp4Error::UNKNOWN)) {
        return Sgp4Error::UNKNOWN;
    }
    return static_cast<Sgp4Error>(error_code);
}

JulianDate JulianDate::from_datetime(
    const int year, const int month, const int day,
    const int hour, const int min, const double sec
) {
    double tmp_jd, tmp_jd_frac;
    vallado_sgp4::jday_SGP4(year, month, day, hour, min, sec, tmp_jd, tmp_jd_frac);
    double jd = tmp_jd + tmp_jd_frac;
    return JulianDate { jd };
}

void JulianDate::to_datetime(
    int &year, int &month, int &day, int &hour, int &min, double &sec
) const {
    vallado_sgp4::invjday_SGP4(jd, 0.0, year, month, day, hour, min, sec);
}

}  // namespace perturb
