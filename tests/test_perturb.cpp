#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

#include <array>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <fstream>

#include "perturb/perturb.hpp"

using namespace perturb;

#ifndef PERTURB_DISABLE_IO
Satellite sat_from_verif_tle(
    std::string &line_1, std::string &line_2,
    double &startmfe, double &stopmfe, double &deltamin
) {
    constexpr char OPS_MODE = 'a';
    REQUIRE(line_1.length() >= TLE_LINE_LEN);
    REQUIRE(line_2.length() >= TLE_LINE_LEN);
    vallado_sgp4::elsetrec sat_rec;
    vallado_sgp4::twoline2rv(
        &line_1[0], &line_2[0], 'v', 'e', OPS_MODE,
        vallado_sgp4::wgs72, startmfe, stopmfe, deltamin, sat_rec
    );
    return Satellite(sat_rec);
}
#endif  // PERTURB_DISABLE_IO

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

#ifndef PERTURB_DISABLE_IO
TEST_CASE("test_sgp4_iss_tle") {
    // Pulled sometime around 2022-03-12
    std::string ISS_TLE_1("1 25544U 98067A   22071.78032407  .00021395  00000-0  39008-3 0  9996");
    std::string ISS_TLE_2("2 25544  51.6424  94.0370 0004047 256.5103  89.8846 15.49386383330227");
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
            const auto err = sat.propagate_from_epoch(mins, pos, vel);
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
            sat.propagate_from_epoch(t, pos_1, vel_1);
            sat.propagate_from_epoch(t + AVG_ISS_ORBITAL, pos_2, vel_2);
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
        constexpr int CHECK_N_ORBITS = 1000;        // Number of consecutive orbits
        Vec3 pos_1, vel_1, pos_2, vel_2;
        for (int i = 0; i < CHECK_N_ORBITS; ++i) {
            const auto t = i * AVG_ISS_ORBITAL;
            sat.propagate_from_epoch(t, pos_1, vel_1);
            sat.propagate_from_epoch(t + AVG_ISS_ORBITAL / 2.0, pos_2, vel_2);
            CHECK(vel_1[0] == doctest::Approx(-vel_2[0]).epsilon(EPS));
            CHECK(vel_1[1] == doctest::Approx(-vel_2[1]).epsilon(EPS));
            CHECK(vel_1[2] == doctest::Approx(-vel_2[2]).epsilon(EPS));
            CHECK(vel_1[0] == doctest::Approx(-vel_2[0]).epsilon(EPS));
            CHECK(vel_1[1] == doctest::Approx(-vel_2[1]).epsilon(EPS));
            CHECK(vel_1[2] == doctest::Approx(-vel_2[2]).epsilon(EPS));
        }
    }
}
#endif  // PERTURB_DISABLE_IO

// Can't run verification mode without debug mode
#ifdef PERTURB_VALLADO_SGP4_ENABLE_DEBUG
TEST_CASE(
    "test_sgp4_verification_mode"
    * doctest::description("Run verification mode based off Vallado's test code")
) {
    constexpr double PI = 3.14159265358979323846;
    constexpr double RAD_TO_DEG = 180 / PI;
    constexpr double DEG_TO_RAD = PI / 180;

    std::ifstream in_file("SGP4-VER.TLE");
    FILE *out_file = std::fopen("generated-tcppver.out", "w");
    REQUIRE(out_file != nullptr);

    Vec3 pos, vel;
    std::string line_1, line_2;
    while (std::getline(in_file, line_1)) {
        if (line_1[0] == '#') {
            continue;
        }
        REQUIRE(std::getline(in_file, line_2));

        double startmfe, stopmfe, deltamin;
        auto sat = sat_from_verif_tle(line_1, line_2, startmfe, stopmfe, deltamin);

        sat.propagate_from_epoch(0.0, pos, vel);  // Initialize maybe??
        std::fprintf(out_file, "%s xx\n", sat.sat_rec.satnum);
        std::fprintf(
            out_file,
            " %16.8f %16.8f %16.8f %16.8f %12.9f %12.9f %12.9f\n",
            sat.sat_rec.t, pos[0], pos[1], pos[2], vel[0], vel[1], vel[2]
        );

        double tsince = startmfe;
        if (std::fabs(tsince) > 1e-8) {
            tsince -= deltamin;
        }

        while ((tsince < stopmfe) && sat.last_error() == Sgp4Error::NONE) {
            tsince = std::min(tsince + deltamin, stopmfe);
            if (sat.propagate_from_epoch(tsince, pos, vel) != Sgp4Error::NONE) {
                continue;
            }

            // FIXME: Change to this when fractional JD added to `JulianDate`
            // const auto jd = sat.epoch() + tsince / 1440.0;
            // const auto ymdhms = jd.to_datetime();

            int year, mon, day, hr, min;
            double sec;
            double jd = sat.sat_rec.jdsatepoch;
            double jd_frac = sat.sat_rec.jdsatepochF + tsince / 1440.0;
            if (jd_frac < 0) {
                jd -= 1;
                jd_frac += 1;
            }
            vallado_sgp4::invjday_SGP4(jd, jd_frac, year, mon, day, hr, min, sec);
            const YMDhms ymdhms { year, mon, day, hr, min, sec };

            std::fprintf(
                out_file,
                " %16.8f %16.8f %16.8f %16.8f %12.9f %12.9f %12.9f",
                tsince, pos[0], pos[1], pos[2], vel[0], vel[1], vel[2]
            );

            double p, a, ecc, incl, node, argp, nu, m, arglat, truelon, lonper;
            vallado_sgp4::rv2coe_SGP4(
                pos.data(), vel.data(), sat.sat_rec.mus,
                p, a, ecc, incl, node, argp, nu, m, arglat, truelon, lonper
            );
            std::fprintf(
                out_file,
                " %14.6f %8.6f %10.5f %10.5f %10.5f %10.5f %10.5f %5i%3i%3i %2i:%2i:%9.6f\n",
                a, ecc, incl * RAD_TO_DEG, node * RAD_TO_DEG, argp * RAD_TO_DEG, nu * RAD_TO_DEG, m * RAD_TO_DEG,
                ymdhms.year, ymdhms.month, ymdhms.day, ymdhms.hour, ymdhms.min, ymdhms.sec
            );
        }
    }

    // Useless scope so it can be visually folded
    {
        vallado_sgp4::elsetrec sat_rec {};
        std::strcpy(sat_rec.satnum, "8195");
        sat_rec.jdsatepoch = 2453911.0;
        sat_rec.jdsatepochF = 0.8321544402;
        sat_rec.no_kozai = 2.00491383;
        sat_rec.ecco = 0.6877146;
        sat_rec.inclo = 64.1586;
        sat_rec.nodeo = 279.0717;
        sat_rec.argpo = 264.7651;
        sat_rec.mo = 20.2257;
        sat_rec.nddot = 0.00000e0;
        sat_rec.bstar = 0.11873e-3;
        sat_rec.ndot = 0.00000099;
        sat_rec.elnum = 813;
        sat_rec.revnum = 22565;
        sat_rec.classification = 'U';
        sat_rec.ephtype = 0;
        std::strcpy(sat_rec.intldesg, "          ");

        constexpr double XPDOTP = 1440.0 / (2 * PI);
        sat_rec.no_kozai = sat_rec.no_kozai / XPDOTP;
        sat_rec.ndot = sat_rec.ndot / (XPDOTP * 1440.0);
        sat_rec.nddot = sat_rec.nddot / (XPDOTP * 1440.0 * 1440);
        sat_rec.inclo = sat_rec.inclo * DEG_TO_RAD;
        sat_rec.nodeo = sat_rec.nodeo * DEG_TO_RAD;
        sat_rec.argpo = sat_rec.argpo * DEG_TO_RAD;
        sat_rec.mo = sat_rec.mo * DEG_TO_RAD;

        double startmfe = 0;
        double stopmfe = 2880;
        double deltamin = 120;
        vallado_sgp4::sgp4init(
            vallado_sgp4::wgs72, 'a', sat_rec.satnum, sat_rec.jdsatepoch + sat_rec.jdsatepochF - 2433281.5,
            sat_rec.bstar, sat_rec.ndot, sat_rec.nddot, sat_rec.ecco, sat_rec.argpo, sat_rec.inclo,
            sat_rec.mo, sat_rec.no_kozai, sat_rec.nodeo, sat_rec
        );
        auto sat = Satellite(sat_rec);

        double tsince = startmfe;
        while ((tsince < stopmfe) && sat.last_error() == Sgp4Error::NONE) {
            tsince = std::min(tsince + deltamin, stopmfe);
            CHECK(sat.propagate_from_epoch(tsince, pos, vel) == Sgp4Error::NONE);
            std::fprintf(
                out_file,
                " %16.8f %16.8f %16.8f %16.8f %12.9f %12.9f %12.9f",
                tsince, pos[0], pos[1], pos[2], vel[0], vel[1], vel[2]
            );
        }
    }
    std::fclose(out_file);
}
#endif  // PERTURB_VALLADO_SGP4_ENABLE_DEBUG
