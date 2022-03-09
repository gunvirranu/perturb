#ifndef PERTURB_PERTURB_HPP
#define PERTURB_PERTURB_HPP

#include "perturb/vallado_sgp4.hpp"

#include <array>
#include <cstddef>

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
    UNKNOWN
};

struct JulianDate {
    double jd;

    static JulianDate from_datetime(int year, int month, int day, int hour, int min, double sec);

    void to_datetime(int &year, int &month, int &day, int &hour, int &min, double &sec) const;
};

class Satellite {
public:
    vallado_sgp4::elsetrec sat_rec;

    explicit Satellite(vallado_sgp4::elsetrec sat_rec);

    static Satellite from_tle(
        const std::array<char, TLE_LINE_LEN> &line_1,
        const std::array<char, TLE_LINE_LEN> &line_2
    );

    Sgp4Error propogate_from_epoch(double mins_from_epoch, Vec3 &pos, Vec3 &vel);

    Sgp4Error propogate(JulianDate jd, Vec3 &pos, Vec3 &vel);
};
}  // namespace perturb

#endif  // PERTURB_PERTURB_HPP
