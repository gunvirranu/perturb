#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

#include "perturb/perturb.hpp"

using namespace::perturb;

TEST_CASE(
    "test_julian_date"
    * doctest::description("Check julian date conversions from 1901 to 2099")
) {
    constexpr double JD_START = 2415750.5;
    constexpr double JD_END = 2488068.5;
    constexpr std::size_t N_CHECKS = 123479;  // Shouldn't evenly divide range
    constexpr double DELTA_JD = (JD_END - JD_START) / N_CHECKS;

    for (std::size_t i = 0; i < N_CHECKS; ++i) {
        const auto jd = JulianDate(JD_START + i * DELTA_JD);
        const YMDhms t = jd.to_datetime();
        CHECK_MESSAGE(
            jd.jd == JulianDate(t).jd,
            t.year, "-", t.month, "-", t.day, " ", t.hour, ":", t.min, ":", t.sec
        );
    }
}
