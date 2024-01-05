#define DOCTEST_CONFIG_TREAT_CHAR_STAR_AS_STRING
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include <array>
#include <cmath>
#include <fstream>

#include "perturb/perturb.hpp"
#include "perturb/tle.hpp"

using namespace perturb;

using doctest::Approx;

#define CHECK_VEC(a, b, eps, scl)                            \
    CHECK((a)[0] == Approx((b)[0]).scale(scl).epsilon(eps)); \
    CHECK((a)[1] == Approx((b)[1]).scale(scl).epsilon(eps)); \
    CHECK((a)[2] == Approx((b)[2]).scale(scl).epsilon(eps))

#ifndef PERTURB_DISABLE_IO
Satellite sat_from_verif_tle(
    std::string &line_1, std::string &line_2, double &startmfe, double &stopmfe,
    double &deltamin
) {
    // Specific to verification TLEs
    constexpr char RUN_TYPE = 'v';
    constexpr char INPUT_TYPE = 'e';
    constexpr char OPS_MODE = 'a';

    // Check received string buffers are of appropriate length
    REQUIRE(line_1.length() >= TLE_LINE_LEN);
    REQUIRE(line_2.length() >= TLE_LINE_LEN);

    // Initialize empty `sat_rec` and let `twoline2rv` fill out
    sgp4::elsetrec sat_rec;
    sgp4::twoline2rv(
        &line_1[0], &line_2[0], RUN_TYPE, INPUT_TYPE, OPS_MODE, sgp4::wgs72, startmfe,
        stopmfe, deltamin, sat_rec
    );

    // Construct and return `Satellite` using pre-parsed `sat_rec`
    return Satellite(sat_rec);
}
#endif  // PERTURB_DISABLE_IO

double norm(const Vec3 &v) {
    return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

TEST_CASE("test_julian_date_type") {
    constexpr double EPS = 1e-10;

    const auto t = DateTime { 2022, 3, 14, 0, 31, 19.3 };
    const auto jd = JulianDate(t);

    // Check subtraction of two JDs
    auto t2 = t;
    t2.day = 17;
    t2.hour = 15;
    t2.min = 45;
    const auto jd2 = JulianDate(t2);
    const double dt = jd2 - jd;
    const double expected_dt =
        (t2.day - t.day) + ((t2.hour - t.hour) + (t2.min - t.min) / 60.0) / 24.0;
    CHECK(dt == Approx(expected_dt).epsilon(EPS));

    // Check that normalization works
    const auto jd_unnorm = jd + dt;
    const auto jd_norm = jd_unnorm.normalized();
    CHECK(jd_norm.jd - 0.5 == std::floor(jd_norm.jd));
    CHECK(0 <= jd_norm.jd_frac);
    CHECK(jd_norm.jd_frac < 1);

    // Check addition of JD and offset
    const auto jd3 = (jd + dt).normalized();
    CHECK(jd2.jd == jd3.jd);
    CHECK(jd2.jd_frac == Approx(jd3.jd_frac).epsilon(EPS));

    // Check addition assignment
    auto jd4 = jd;
    jd4 += dt;
    jd4.normalize();
    CHECK(jd3.jd == jd4.jd);
    CHECK(jd3.jd_frac == jd4.jd_frac);

    // Check subtraction of JD and offset
    const auto jd5 = (jd3 - dt).normalized();
    CHECK(jd5.jd == jd.jd);
    CHECK(jd5.jd_frac == Approx(jd.jd_frac).epsilon(EPS));

    // Check ordering overloads
    CHECK(jd < jd2);
    CHECK_FALSE(jd2 < jd);
    CHECK(jd2 > jd);
    CHECK_FALSE(jd > jd2);
    CHECK(jd <= jd2);
    CHECK_FALSE(jd2 <= jd);
    CHECK(jd <= jd);
    CHECK(jd2 <= jd2);
    CHECK(jd2 >= jd);
    CHECK_FALSE(jd >= jd2);
    CHECK(jd >= jd);
    CHECK(jd2 >= jd2);
}

TEST_CASE(
    "test_julian_date_conversions"
    * doctest::description("Check julian date conversions from 1901 to 2099")
) {
    const auto JD_START = JulianDate(2415750.5);  // Around start of 1901
    const auto JD_END = JulianDate(2488068.5);    // Around end of 2099

    constexpr int N_CHECKS = 123479;  // Shouldn't evenly divide range
    const double DELTA_JD = (JD_END - JD_START) / N_CHECKS;

    for (int i = 0; i < N_CHECKS; ++i) {
        const auto jd = (JD_START + i * DELTA_JD).normalized();
        const DateTime t = jd.to_datetime();
        const auto jd_conv = JulianDate(t);
        CHECK_MESSAGE(jd.jd == jd_conv.jd, t.year, "-", t.month, "-", t.day);
        CHECK_MESSAGE(
            jd.jd_frac == Approx(jd_conv.jd_frac).epsilon(1e-12), t.hour, ":", t.min,
            ":", t.sec
        );
    }
}

#ifndef PERTURB_DISABLE_IO
TEST_CASE("test_tle_parse") {
    const char *TLE_1 =
        "1 25544U 98067A   22071.78032407  .00021395  00000-0  39008-3 0  9996";
    const char *TLE_2 =
        "2 25544  51.6424  94.0370 0004047 256.5103  89.8846 15.49386383330227";

    SUBCASE("test_standard") {
        TwoLineElement tle {};
        const auto err1 = tle.parse(TLE_1, TLE_2);
        CHECK(err1 == TLEParseError::NONE);

        CHECK(tle.catalog_number == "25544");
        CHECK(tle.classification == 'U');
        CHECK(tle.launch_year == 98U);
        CHECK(tle.launch_number == 67U);
        CHECK(tle.launch_piece == "A");
        CHECK(tle.epoch_year == 22U);
        CHECK(tle.epoch_day_of_year == 71.78032407);
        CHECK(tle.n_dot == 0.00021395);
        CHECK(tle.n_ddot == 0.0e0);
        CHECK(tle.b_star == 0.39008e-3);
        CHECK(tle.ephemeris_type == 0U);
        CHECK(tle.element_set_number == 999U);
        CHECK(tle.line_1_checksum == 6U);

        CHECK(tle.inclination == 51.6424);
        CHECK(tle.raan == 94.0370);
        CHECK(tle.eccentricity == 0.0004047);
        CHECK(tle.arg_of_perigee == 256.5103);
        CHECK(tle.mean_anomaly == 89.8846);
        CHECK(tle.mean_motion == 15.49386383);
        CHECK(tle.revolution_number == 33022UL);
        CHECK(tle.line_2_checksum == 7U);

        const auto err2 = tle.parse(
            "1 25544U 98067 BA 22071.78032407  .00021395 .00000-0 .39008-3 0 39999",
            "2 25544  51.6424  94.0370 0004047 256.5103  89.8846  5.49386383 30223"
        );
        CHECK(err2 == TLEParseError::NONE);
        CHECK(tle.launch_piece == "BA");
        CHECK(tle.n_ddot == 0.0e0);
        CHECK(tle.b_star == 0.39008e-3);
        CHECK(tle.mean_motion == 5.49386383);
        CHECK(tle.revolution_number == 3022UL);
        CHECK(tle.line_2_checksum == 3U);
    }

    SUBCASE("test_line_1_errors") {
        TwoLineElement tle {};
        const auto err1 = tle.parse(
            "1 25544U*98067A   22071.78032407  .00021395  00000-0  39008-3 0  9996",
            TLE_2
        );
        CHECK(err1 == TLEParseError::SHOULD_BE_SPACE);

        const auto err2 = tle.parse(
            "1 25544U 98067A   22071.78032407  .00021395  00000-0  39008-3 0  9990",
            TLE_2
        );
        CHECK(err2 == TLEParseError::CHECKSUM_MISMATCH);

        const auto err3 = tle.parse(
            "1 25544U 98067A   22071.78*32407  .00021395  00000-0  39008-3 0  9996",
            TLE_2
        );
        CHECK(err3 == TLEParseError::INVALID_FORMAT);
    }

    SUBCASE("test_line_2_errors") {
        TwoLineElement tle {};
        const auto err1 = tle.parse(
            TLE_1,
            "2 25544  51.6424* 94.0370 0004047 256.5103  89.8846 15.49386383330227"
        );
        CHECK(err1 == TLEParseError::SHOULD_BE_SPACE);

        const auto err2 = tle.parse(
            TLE_1,
            "2 25544  51.6424  94.0370 0004047 256.5103  89.8846 15.49386383330220"
        );
        CHECK(err2 == TLEParseError::CHECKSUM_MISMATCH);

        const auto err3 = tle.parse(
            TLE_1,
            "2 25544  51.6424  94.0370 0004047 256.5103  89.8846 15.493*6383330227"
        );
        CHECK(err3 == TLEParseError::INVALID_FORMAT);
    }
}
#endif  // PERTURB_DISABLE_IO

#ifndef PERTURB_DISABLE_IO
TEST_CASE("test_sgp4_iss_tle") {
    // Pulled sometime around 2022-03-12
    std::string ISS_TLE_1(
        "1 25544U 98067A   22071.78032407  .00021395  00000-0  39008-3 0  9996"
    );
    std::string ISS_TLE_2(
        "2 25544  51.6424  94.0370 0004047 256.5103  89.8846 15.49386383330227"
    );

    auto sat = Satellite::from_tle(ISS_TLE_1, ISS_TLE_2);
    REQUIRE(sat.last_error() == Sgp4Error::NONE);

    // Check that the epoch is correct based off manual calculations
    SUBCASE("test_epoch") {
        const DateTime epoch = sat.epoch().to_datetime();

        CHECK(epoch.year == 2022);
        CHECK(epoch.month == 3);
        CHECK(epoch.day == 12);
        CHECK(epoch.hour == 18);
        CHECK(epoch.min == 43);
        CHECK(epoch.sec == Approx(40).epsilon(1e-5));
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
            StateVector sv;
            const auto err = sat.propagate_from_epoch(mins, sv);
            CHECK(err == Sgp4Error::NONE);

            const double dist = norm(sv.position) - AVG_EARTH_RADIUS;
            CHECK(dist == Approx(AVG_ISS_HEIGHT).epsilon(0.05));

            const double speed = norm(sv.velocity);
            CHECK(speed == Approx(AVG_ISS_SPEED).epsilon(0.05));

            mins += CHECK_EVERY_MINS;
        }
    }

    // Check position and velocity vectors are roughly the same after entire orbits
    SUBCASE("test_whole_orbit") {
        constexpr double AVG_ISS_ORBITAL = 92.8;  // minutes
        constexpr int CHECK_N_ORBITS = 1000;      // Number of consecutive orbits

        for (int i = 0; i < CHECK_N_ORBITS; ++i) {
            const double t = i * AVG_ISS_ORBITAL;

            StateVector sv_1, sv_2;
            sat.propagate_from_epoch(t, sv_1);
            sat.propagate_from_epoch(t + AVG_ISS_ORBITAL, sv_2);

            const auto &pos_1 = sv_1.position, &pos_2 = sv_2.position;
            const auto &vel_1 = sv_1.velocity, &vel_2 = sv_2.velocity;
            CHECK_VEC(pos_1, pos_2, 0.05, 1000);
            CHECK_VEC(vel_1, vel_2, 0.05, 5);
        }
    }

    // Check position and velocity is roughly inverse after half orbits
    SUBCASE("test_half_orbit") {
        constexpr double AVG_ISS_ORBITAL = 92.8;  // minutes
        constexpr int CHECK_N_ORBITS = 1000;      // Number of consecutive orbits

        for (int i = 0; i < CHECK_N_ORBITS; ++i) {
            const auto t = i * AVG_ISS_ORBITAL;

            StateVector sv_1, sv_2;
            sat.propagate_from_epoch(t, sv_1);
            sat.propagate_from_epoch(t + AVG_ISS_ORBITAL / 2.0, sv_2);

            // Flip `sv_2` so it can be compared against `sv_1`
            for (std::size_t j = 0; j < 3; ++j) {
                sv_2.position[j] *= -1;
                sv_2.velocity[j] *= -1;
            }

            const auto &pos_1 = sv_1.position, &pos_2 = sv_2.position;
            const auto &vel_1 = sv_1.velocity, &vel_2 = sv_2.velocity;
            CHECK_VEC(pos_1, pos_2, 0.05, 1000);
            CHECK_VEC(vel_1, vel_2, 0.05, 5);
        }
    }
}
#endif  // PERTURB_DISABLE_IO

// Can't run verification mode without debug mode
#ifdef PERTURB_SGP4_ENABLE_DEBUG
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

    std::string line_1, line_2;
    while (std::getline(in_file, line_1)) {
        if (line_1[0] == '#') {
            continue;
        }
        REQUIRE(std::getline(in_file, line_2));

        double startmfe, stopmfe, deltamin;
        auto sat = sat_from_verif_tle(line_1, line_2, startmfe, stopmfe, deltamin);

        StateVector sv;
        sat.propagate_from_epoch(0.0, sv);  // Initialize maybe??
        auto pos = sv.position;
        auto vel = sv.velocity;

        std::fprintf(out_file, "%s xx\n", sat.sat_rec.satnum);
        std::fprintf(
            out_file, " %16.8f %16.8f %16.8f %16.8f %12.9f %12.9f %12.9f\n",
            sat.sat_rec.t, pos[0], pos[1], pos[2], vel[0], vel[1], vel[2]
        );

        double tsince = startmfe;
        if (std::fabs(tsince) > 1e-8) {
            tsince -= deltamin;
        }

        while ((tsince < stopmfe) && sat.last_error() == Sgp4Error::NONE) {
            tsince = std::min(tsince + deltamin, stopmfe);

            if (sat.propagate_from_epoch(tsince, sv) != Sgp4Error::NONE) {
                continue;
            }
            pos = sv.position;
            vel = sv.velocity;

            const auto jd = (sat.epoch() + tsince / 1440.0).normalized();
            const auto ymdhms = jd.to_datetime();

            std::fprintf(
                out_file, " %16.8f %16.8f %16.8f %16.8f %12.9f %12.9f %12.9f", tsince,
                pos[0], pos[1], pos[2], vel[0], vel[1], vel[2]
            );

            double p, a, ecc, incl, node, argp, nu, m, arglat, truelon, lonper;
            sgp4::rv2coe_SGP4(
                pos.data(), vel.data(), sat.sat_rec.mus, p, a, ecc, incl, node, argp, nu,
                m, arglat, truelon, lonper
            );

            std::fprintf(
                out_file,
                " %14.6f %8.6f %10.5f %10.5f %10.5f %10.5f %10.5f %5i%3i%3i %2i:%2i:%9.6f\n",
                a, ecc, incl * RAD_TO_DEG, node * RAD_TO_DEG, argp * RAD_TO_DEG,
                nu * RAD_TO_DEG, m * RAD_TO_DEG, ymdhms.year, ymdhms.month, ymdhms.day,
                ymdhms.hour, ymdhms.min, ymdhms.sec
            );
        }
    }

    // Useless scope so it can be visually folded
    {
        sgp4::elsetrec sat_rec {};
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
        sgp4::sgp4init(
            sgp4::wgs72, 'a', sat_rec.satnum,
            sat_rec.jdsatepoch + sat_rec.jdsatepochF - 2433281.5, sat_rec.bstar,
            sat_rec.ndot, sat_rec.nddot, sat_rec.ecco, sat_rec.argpo, sat_rec.inclo,
            sat_rec.mo, sat_rec.no_kozai, sat_rec.nodeo, sat_rec
        );
        auto sat = Satellite(sat_rec);

        double tsince = startmfe;
        while ((tsince < stopmfe) && sat.last_error() == Sgp4Error::NONE) {
            tsince = std::min(tsince + deltamin, stopmfe);

            StateVector sv;
            CHECK(sat.propagate_from_epoch(tsince, sv) == Sgp4Error::NONE);

            const auto &pos = sv.position;
            const auto &vel = sv.velocity;
            std::fprintf(
                out_file, " %16.8f %16.8f %16.8f %16.8f %12.9f %12.9f %12.9f", tsince,
                pos[0], pos[1], pos[2], vel[0], vel[1], vel[2]
            );
        }
    }
    std::fclose(out_file);
}
#endif  // PERTURB_SGP4_ENABLE_DEBUG

#ifdef PERTURB_SGP4_ENABLE_DEBUG
#  define CHECK_AB_MEMBER(x)          CHECK(a.x == b.x)
#  define CHECK_AB_MEMBER_EPS(x, eps) CHECK(a.x == Approx(b.x).epsilon(eps))

TEST_CASE(
    "test_tle_parser_with_verif_mode"
    * doctest::description("Compare TLE parser against Vallado on all verif TLEs")
) {
    std::ifstream in_file("SGP4-VER.TLE");
    std::string line_1, line_2;

    while (std::getline(in_file, line_1)) {
        if (line_1[0] == '#') {
            continue;
        }
        REQUIRE(std::getline(in_file, line_2));

        line_1.resize(TLE_LINE_LEN);
        line_2.resize(TLE_LINE_LEN);
        INFO("TLE: ", line_1, "\n       ", line_2);

        // Parse into `TLE` type and construct `Satellite` as `sat_tle`
        TwoLineElement tle {};
        const auto err = tle.parse(line_1, line_2);

        // Ignore `NONE`, `CHECKSUM_MISMATCH`, and `INVALID_VALUE`
        CHECK(err != TLEParseError::SHOULD_BE_SPACE);
        CHECK(err != TLEParseError::INVALID_FORMAT);
        // `TLE::parse` can't handle some odd parsing cases
        if (err == TLEParseError::INVALID_VALUE) {
            continue;
        }
        auto sat_tle = Satellite(tle);

        // Parse and construct `sat_orig` using Vallado's impl
        auto sat_orig = Satellite::from_tle(line_1, line_2);
        CHECK(sat_orig.last_error() != Sgp4Error::INVALID_TLE);

        // Correct some unimportant differences
        sat_tle.sat_rec.elnum = (sat_tle.sat_rec.elnum * 10 + tle.line_1_checksum);
        sat_tle.sat_rec.revnum = (sat_tle.sat_rec.revnum * 10 + tle.line_2_checksum);

        // Check that a subset of the member vars match. Others are based off these.
        const auto &a = sat_tle.sat_rec, &b = sat_orig.sat_rec;
        // Line 1
        CHECK_AB_MEMBER(satnum);
        CHECK_AB_MEMBER(classification);
        // `Satellite(const TLE &)` doesn't set `sat_rec.intldesg`, so ignore
        // CHECK_SAT_MEMBER(intldesg);
        CHECK_AB_MEMBER(epochyr);
        CHECK_AB_MEMBER(epochdays);
        CHECK_AB_MEMBER(epochdays);
        CHECK_AB_MEMBER(ndot);
        CHECK_AB_MEMBER_EPS(nddot, 1e-16);
        CHECK_AB_MEMBER_EPS(bstar, 1e-16);
        CHECK_AB_MEMBER(ephtype);
        CHECK_AB_MEMBER(elnum);
        // Line 2
        CHECK_AB_MEMBER(inclo);
        CHECK_AB_MEMBER(nodeo);
        CHECK_AB_MEMBER(ecco);
        CHECK_AB_MEMBER(argpo);
        CHECK_AB_MEMBER(mo);
        CHECK_AB_MEMBER(no_kozai);
        CHECK_AB_MEMBER(revnum);

        // Try propagating to check that output predictions match
        for (const double mins : { 0.0, 0.5, 5.0, 30.0, 1440.0, 20000.0 }) {
            CAPTURE(mins);

            StateVector sv_a {}, sv_b {};
            (void) sat_tle.propagate_from_epoch(mins, sv_a);
            (void) sat_orig.propagate_from_epoch(mins, sv_b);

            CHECK(sv_a.epoch.jd == sv_b.epoch.jd);
            CHECK(sv_a.epoch.jd_frac == sv_b.epoch.jd_frac);
            CHECK_VEC(sv_a.position, sv_b.position, 1e-16, 1);
            CHECK_VEC(sv_a.velocity, sv_b.velocity, 1e-14, 1);
        }
    }
}
#endif  // PERTURB_SGP4_ENABLE_DEBUG
