#ifndef PERTURB_PERTURB_HPP
#define PERTURB_PERTURB_HPP

#include "perturb/vallado_sgp4.hpp"

#include <array>
#include <cstddef>
#include <string>

namespace perturb {

using Vec3 = std::array<double, 3>;

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

    static Satellite from_tle(char *line_1, char *line_2);

    static Satellite from_tle(std::string &line_1, std::string &line_2);

    Sgp4Error last_error() const;

    JulianDate epoch() const;

    Sgp4Error propogate_from_epoch(double mins_from_epoch, Vec3 &pos, Vec3 &vel);

    Sgp4Error propogate(JulianDate jd, Vec3 &pos, Vec3 &vel);
};
}  // namespace perturb

#endif  // PERTURB_PERTURB_HPP
