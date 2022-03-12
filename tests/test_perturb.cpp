#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

#include "perturb/perturb.hpp"

using namespace::perturb;

TEST_CASE(
    "test_julian_date"
    * doctest::description("Check julian date conversions from 1901 to 2099")
) {
    const auto JD_START = JulianDate(2415750.5);
    const auto JD_END = JulianDate(2488068.5);
    constexpr int N_CHECKS = 123479;  // Shouldn't evenly divide range
    const double DELTA_JD = (JD_END - JD_START) / N_CHECKS;

    for (int i = 0; i < N_CHECKS; ++i) {
        const auto jd = JulianDate(JD_START.jd + i * DELTA_JD);
        const YMDhms t = jd.to_datetime();
        CHECK_MESSAGE(
            jd.jd == JulianDate(t).jd,
            t.year, "-", t.month, "-", t.day, " ", t.hour, ":", t.min, ":", t.sec
        );
    }
}
