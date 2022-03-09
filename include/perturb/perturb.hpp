#ifndef PERTURB_PERTURB_HPP
#define PERTURB_PERTURB_HPP

namespace perturb {

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

}  // namespace perturb

#endif  // PERTURB_PERTURB_HPP
