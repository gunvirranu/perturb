/*
 * perturb -- A modern C++11 wrapper for the SGP4 orbit propagator
 * Version 1.0.0
 * https://github.com/gunvirranu/perturb
 *
 * Licensed under the MIT License <http://opensource.org/licenses/MIT>.
 * SPDX-License-Identifier: MIT
 *
 * Copyright (c) 2022 Gunvir Ranu
 */

#include "perturb/tle.hpp"

#include <array>
#ifndef PERTURB_DISABLE_IO
#  include <cctype>
#  include <cmath>
#  include <cstdio>
#  include <cstring>
#endif

namespace perturb {

#ifndef PERTURB_DISABLE_IO
static unsigned int calc_tle_line_checksum(const char *line) {
    unsigned int checksum = 0U;
    for (std::size_t i = 0; i < (TLE_LINE_LEN - 1); ++i) {
        if (std::isdigit(line[i])) {
            checksum += static_cast<unsigned int>(line[i] - '0');
        }
        if (line[i] == '-') {
            checksum += 1U;
        }
    }
    return (checksum % 10U);
}
#endif  // PERTURB_DISABLE_IO

#ifndef PERTURB_DISABLE_IO
// FIXME: Use a more robust parsing method. I wish string_view existed :(
TLEParseError TwoLineElement::parse(const char *line_1, const char *line_2) {
    // Make sure there are spaces in the right places
    constexpr std::array<int, 8> LINE_1_SPACES = { 2, 9, 18, 33, 44, 53, 62, 64 };
    for (const int i : LINE_1_SPACES) {
        if (line_1[i - 1] != ' ') {
            return TLEParseError::SHOULD_BE_SPACE;
        }
    }
    constexpr std::array<int, 7> LINE_2_SPACES = { 2, 8, 17, 26, 34, 43, 52 };
    for (const int i : LINE_2_SPACES) {
        if (line_2[i - 1] != ' ') {
            return TLEParseError::SHOULD_BE_SPACE;
        }
    }

    // Parse format

    // Line 1
    constexpr auto LINE_1_FMT_STR =
        "%1hhu %5s %1c %2u %3u %3s %2u %12lf %10lf %6lf %2d %6lf %2d %1hhu %4u %n %1hhu";
    unsigned char line1_num;
    int n_ddot_exp, b_star_exp, l1_pre_checksum;
    int l1_scanned = std::sscanf(  // Brooo C++ sucks at string processing :(
        line_1, LINE_1_FMT_STR, &line1_num, this->catalog_number, &this->classification,
        &this->launch_year, &this->launch_number, this->launch_piece, &this->epoch_year,
        &this->epoch_day_of_year, &this->n_dot, &this->n_ddot, &n_ddot_exp,
        &this->b_star, &b_star_exp, &this->ephemeris_type, &this->element_set_number,
        &l1_pre_checksum, &this->line_1_checksum
    );

    // Line 2
    constexpr auto LINE_2_FMT_STR_NO_SPACE =
        "%1hhu %5s %8lf %8lf %7lu %8lf %8lf %11lf %5lu %n %1hhu";
    constexpr auto LINE_2_FMT_STR_WT_SPACE =
        "%1hhu %5s %8lf %8lf %7lu %8lf %8lf %10lf %5lu %n %1hhu";
    unsigned char line2_num;
    char catlog_num_line2[6];
    unsigned long eccentricity_int;
    int l2_scanned, l2_pre_checksum;
    // :( Yeah, I know, this is annoying...
    if (line_2[52] != ' ') {
        l2_scanned = std::sscanf(
            line_2, LINE_2_FMT_STR_NO_SPACE, &line2_num, catlog_num_line2,
            &this->inclination, &this->raan, &eccentricity_int, &this->arg_of_perigee,
            &this->mean_anomaly, &this->mean_motion, &this->revolution_number,
            &l2_pre_checksum, &this->line_2_checksum
        );
    } else {
        l2_scanned = std::sscanf(
            line_2, LINE_2_FMT_STR_WT_SPACE, &line2_num, catlog_num_line2,
            &this->inclination, &this->raan, &eccentricity_int, &this->arg_of_perigee,
            &this->mean_anomaly, &this->mean_motion, &this->revolution_number,
            &l2_pre_checksum, &this->line_2_checksum
        );
    }

    this->catalog_number[5] = '\0';
    this->launch_piece[3] = '\0';
    catlog_num_line2[5] = '\0';

    // Handle some annoying special cases
    const bool fix_l1_checksum = (l1_scanned == 15)
        && (l1_pre_checksum >= static_cast<int>(TLE_LINE_LEN))
        && (line_1[TLE_LINE_LEN - 5] == ' ');
    if (fix_l1_checksum) {
        // Element set number often doesn't have leading zeroes,
        // so `sscanf` over-consumes and eats up the final checksum.
        this->element_set_number /= 10U;  // Chop of least-significant digit
        this->line_1_checksum =
            static_cast<unsigned char>(line_1[TLE_LINE_LEN - 1] - '0');
        l1_scanned += 1;
    }
    const bool fix_l2_checksum = (l2_scanned == 9)
        && (l2_pre_checksum >= static_cast<int>(TLE_LINE_LEN))
        && (line_2[TLE_LINE_LEN - 6] == ' ');
    if (fix_l2_checksum) {
        // If revolution number doens't have leading zero, so same issue as above
        this->revolution_number /= 10U;
        this->line_2_checksum =
            static_cast<unsigned char>(line_2[TLE_LINE_LEN - 1] - '0');
        l2_scanned += 1;
    }
    if (line_1[44] != '.') {
        n_ddot_exp -= 5;
    }
    if (line_1[53] != '.') {
        b_star_exp -= 5;
    }

    // Check that the correct number of values were parsed
    if (l1_scanned != 16 || l2_scanned != 10) {
        return TLEParseError::INVALID_FORMAT;
    }

    // Post-process
    this->n_ddot *= std::pow(10.0, n_ddot_exp);
    this->b_star *= std::pow(10.0, b_star_exp);
    this->eccentricity = static_cast<double>(eccentricity_int) / 1.0e7;

    // Check that valid values were parsed
    bool valid_vals = true;
    // Line 1
    valid_vals &= (line1_num == 1);
    const char clsf = this->classification;
    valid_vals &= (clsf == 'U') || (clsf == 'C') || (clsf == 'S');
    valid_vals &= (this->launch_year < 100U) && (this->epoch_year < 100U);
    valid_vals &= (1.0 <= this->epoch_day_of_year) && (this->epoch_day_of_year <= 366.0);
    valid_vals &= (-15 < n_ddot_exp) && (n_ddot_exp < 10);
    valid_vals &= (-15 < b_star_exp) && (b_star_exp < 10);
    valid_vals &= (this->ephemeris_type == 0U);
    valid_vals &= (this->element_set_number < 10000U);
    // Line 2
    valid_vals &= (line2_num == 2);
    valid_vals &= (std::strcmp(this->catalog_number, catlog_num_line2) == 0);
    valid_vals &= (0.0 <= this->inclination) && (this->inclination <= 180.0);
    valid_vals &= (0.0 <= this->raan) && (this->raan <= 360.0);
    valid_vals &= (0.0 <= this->arg_of_perigee) && (this->arg_of_perigee <= 360.0);
    valid_vals &= (0.0 <= this->mean_anomaly) && (this->mean_anomaly <= 360.0);
    if (!valid_vals) {
        return TLEParseError::INVALID_VALUE;
    }

    // Calculate and compare checksums
    bool checksum_matches = true;
    checksum_matches &= (calc_tle_line_checksum(line_1) == this->line_1_checksum);
    checksum_matches &= (calc_tle_line_checksum(line_2) == this->line_2_checksum);
    if (!checksum_matches) {
        return TLEParseError::CHECKSUM_MISMATCH;
    }

    return TLEParseError::NONE;
}
#endif  // PERTURB_DISABLE_IO

#ifndef PERTURB_DISABLE_IO
TLEParseError TwoLineElement::parse(
    const std::string &line_1, const std::string &line_2
) {
    if (line_1.length() < TLE_LINE_LEN || line_2.length() < TLE_LINE_LEN) {
        return TLEParseError::INVALID_FORMAT;
    }
    return this->parse(line_1.c_str(), line_2.c_str());
}
#endif  // PERTURB_DISABLE_IO

}  // namespace perturb
