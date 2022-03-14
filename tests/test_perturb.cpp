#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

#include <array>
#include <cstring>
#include <cmath>

#include "perturb/perturb.hpp"

using namespace perturb;

std::array<char, TLE_LINE_LEN> str_to_arr(const char *str) {
    std::array<char, TLE_LINE_LEN> line;
    std::memcpy(line.data(), str, TLE_LINE_LEN);
    return line;
}

double norm(const Vec3 &v) {
    return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

TEST_CASE("test_julian_date_type") {
    const auto t = YMDhms { 2022, 3, 14, 0, 31, 19.3 };
    const auto jd = JulianDate(t);
    auto t2 = t;
    t2.day = 17;
    t2.hour = 15;
    t2.min = 45;
    const auto jd2 = JulianDate(t2);
    const double dt = jd2 - jd;
    const double EXPECTED_DT = (t2.day - t.day) + ((t2.hour - t.hour) + (t2.min - t.min) / 60.0) / 24.0;
    CHECK(dt == doctest::Approx(EXPECTED_DT).epsilon(1e-8));

    const auto jd3 = jd + dt;
    CHECK(jd2.jd == jd3.jd);

    auto jd4 = jd;
    jd4 += dt;
    CHECK(jd3.jd == jd4.jd);
}

TEST_CASE(
    "test_julian_date_conversions"
    * doctest::description("Check julian date conversions from 1901 to 2099")
) {
    const auto JD_START = JulianDate(2415750.5);    // Around start of 1901
    const auto JD_END = JulianDate(2488068.5);      // Around end of 2099
    constexpr int N_CHECKS = 123479;  // Shouldn't evenly divide range
    const double DELTA_JD = (JD_END - JD_START) / N_CHECKS;

    for (int i = 0; i < N_CHECKS; ++i) {
        const auto jd = JD_START + i * DELTA_JD;
        const YMDhms t = jd.to_datetime();
        CHECK_MESSAGE(
            jd.jd == JulianDate(t).jd,
            t.year, "-", t.month, "-", t.day, " ", t.hour, ":", t.min, ":", t.sec
        );
    }
}

TEST_CASE("test_sgp4_iss_tle") {
    // Pulled sometime around 2022-03-12
    const auto ISS_TLE_1 = str_to_arr("1 25544U 98067A   22071.78032407  .00021395  00000-0  39008-3 0  9996");
    const auto ISS_TLE_2 = str_to_arr("2 25544  51.6424  94.0370 0004047 256.5103  89.8846 15.49386383330227");
    auto sat = Satellite::from_tle(ISS_TLE_1, ISS_TLE_2);
    REQUIRE(sat.last_error() == Sgp4Error::NONE);

    // Check that the epoch is correct based of manual calculations
    SUBCASE("test_epoch") {
        const YMDhms epoch = sat.epoch().to_datetime();
        CHECK(epoch.year == 2022);
        CHECK(epoch.month == 3);
        CHECK(epoch.day == 12);
        CHECK(epoch.hour == 18);
        CHECK(epoch.min == 43);
        CHECK(epoch.sec == doctest::Approx(40).epsilon(1e-5));
    }

    // Check height above Earth and speed roughly stay the same over a week.
    // The ISS is in a pretty circular orbit, but still need a ~5 % tolerance.
    SUBCASE("test_propagation_distance_speed") {
        // Propagate every minute for the next week
        constexpr double AVG_EARTH_RADIUS = 6371;       // km
        constexpr double AVG_ISS_HEIGHT = 410;          // km
        constexpr double AVG_ISS_SPEED = 8;             // km / s
        constexpr double CHECK_EVERY_MINS = 1;          // Every 1 minute
        constexpr double CHECK_FOR_MINS = 7 * 24 * 60;  // For 1 week
        double mins = 0;
        while (mins < CHECK_FOR_MINS) {
            Vec3 pos, vel;
            const auto err = sat.propogate_from_epoch(mins, pos, vel);
            CHECK(err == Sgp4Error::NONE);
            const double dist = norm(pos) - AVG_EARTH_RADIUS;
            const double speed = norm(vel);
            CHECK(dist == doctest::Approx(AVG_ISS_HEIGHT).epsilon(0.05));
            CHECK(speed == doctest::Approx(AVG_ISS_SPEED).epsilon(0.05));
            mins += CHECK_EVERY_MINS;
        }
    }

    // Check position and velocity vectors are roughly the same after entire orbits
    SUBCASE("test_whole_orbit") {
        constexpr double AVG_ISS_ORBITAL = 92.8;    // minutes
        constexpr double EPS = 0.1;                 // Acceptable relative delta
        constexpr int CHECK_N_ORBITS = 1000;        // Number of consecutive orbits
        Vec3 pos_1, vel_1, pos_2, vel_2;
        for (int i = 0; i < CHECK_N_ORBITS; ++i) {
            const double t = i * AVG_ISS_ORBITAL;
            sat.propogate_from_epoch(t, pos_1, vel_1);
            sat.propogate_from_epoch(t + AVG_ISS_ORBITAL, pos_2, vel_2);
            CHECK(vel_1[0] == doctest::Approx(vel_2[0]).epsilon(EPS));
            CHECK(vel_1[1] == doctest::Approx(vel_2[1]).epsilon(EPS));
            CHECK(vel_1[2] == doctest::Approx(vel_2[2]).epsilon(EPS));
            CHECK(vel_1[0] == doctest::Approx(vel_2[0]).epsilon(EPS));
            CHECK(vel_1[1] == doctest::Approx(vel_2[1]).epsilon(EPS));
            CHECK(vel_1[2] == doctest::Approx(vel_2[2]).epsilon(EPS));
        }
    }

    // Check position and velocity is roughly inverse after half orbits
    SUBCASE("test_half_orbit") {
        constexpr double AVG_ISS_ORBITAL = 92.8;    // minutes
        constexpr double EPS = 0.05;                // Acceptable relative delta
        constexpr int CHECK_N_ORBITS = 1000;           // Number of consecutive orbits
        Vec3 pos_1, vel_1, pos_2, vel_2;
        for (int i = 0; i < CHECK_N_ORBITS; ++i) {
            const auto t = i * AVG_ISS_ORBITAL;
            sat.propogate_from_epoch(t, pos_1, vel_1);
            sat.propogate_from_epoch(t + AVG_ISS_ORBITAL / 2.0, pos_2, vel_2);
            CHECK(vel_1[0] == doctest::Approx(-vel_2[0]).epsilon(EPS));
            CHECK(vel_1[1] == doctest::Approx(-vel_2[1]).epsilon(EPS));
            CHECK(vel_1[2] == doctest::Approx(-vel_2[2]).epsilon(EPS));
            CHECK(vel_1[0] == doctest::Approx(-vel_2[0]).epsilon(EPS));
            CHECK(vel_1[1] == doctest::Approx(-vel_2[1]).epsilon(EPS));
            CHECK(vel_1[2] == doctest::Approx(-vel_2[2]).epsilon(EPS));
        }
    }
}
