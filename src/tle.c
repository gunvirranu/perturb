/*
 * perturb -- A modern C++11 wrapper for the SGP4 orbit propagator
 * version 1.0.0
 * https://github.com/gunvirranu/perturb
 *
 * Copyright (c) 2022 Gunvir Singh Ranu
 * SPDX-License-Identifier: MIT
 */

#include "perturb/tle.h"

#include <stddef.h>
#ifndef PERTURB_DISABLE_IO
#  include <stdio.h>
#  include <inttypes.h>
#endif

#include "common_private.h"

#ifdef __cplusplus
#  error "I kindly request you compile this as C instead of C++"
#endif

static bool check_for_missing_spaces(
    const char * line,
    const size_t * spaces,
    const size_t n_spaces
) {
    for (size_t i = 0U; i < n_spaces; ++i) {
        const size_t idx_space = spaces[i];

        if (line[idx_space] != ' ') {
            return true;
        }
    }
    return false;
}

enum perturb_Sgp4Error perturb_init_sat_from_tle(
    const struct perturb_TwoLineElement tle,
    const enum perturb_GravityModel grav_model,
    struct perturb_Satellite * sat
) {
    if (sat == NULL) {
        return PERTURB_SGP4_ERROR_INVALID_INPUT;
    }

    // FIXME: impl
    return PERTURB_SGP4_ERROR_INVALID_INPUT;
}

#ifndef PERTURB_DISABLE_IO
enum perturb_TleParseError perturb_parse_tle(
    const char * line_1, const char * line_2,
    struct perturb_TwoLineElement * tle
) {
    const bool bad_ptrs = (line_1 == NULL) || (line_2 == NULL) || (tle == NULL);
    if (bad_ptrs) {
        return PERTURB_TLE_PARSE_ERROR_INVALID_INPUT;
    }

    // Make sure there are spaces in the right places
    const size_t LINE_1_SPACES[] = { 1, 8, 17, 32, 43, 52, 61, 64 };
    const size_t LINE_2_SPACES[] = { 1, 7, 16, 25, 33, 42, 51 };

    const bool is_space_missing = (
        check_for_missing_spaces(line_1, LINE_1_SPACES, ARRAY_SIZE(LINE_1_SPACES)) ||
        check_for_missing_spaces(line_2, LINE_2_SPACES, ARRAY_SIZE(LINE_2_SPACES))
    );
    if (is_space_missing) {
        return PERTURB_TLE_PARSE_ERROR_SHOULD_BE_SPACE;
    }

    // Parse format - Line 1
    // TODO: Use a more robust parsing method. I wish from_chars existed :(
    const char LINE_1_FMT_STR[] = (
        "%1c %5s %1c "
        "%2" SCNu8 " %3" SCNu16 " %3s "
        "%2" SCNu8 " %12lf "
        "%10lf %6lf %2" SCNd8 " "
        "%6lf %2" SCNd8 " "
        "%1" SCNu8 " %4" SCNu16 " "
        "%n %1" SCNu8
    );

    // FIXME: swap these to use normal int types and then convert to fixed-size for saving
    char line1_num = '0';
    int8_t n_ddot_exp = 0;
    int8_t b_star_exp = 0;
    unsigned int line1_pre_checksum = 0U;

    int line1_scanned = sscanf(  // bruh C and C++ both suck at string processing :(
        line_1, LINE_1_FMT_STR,
        &line1_num, tle->catalog_number, &tle->classification,
        &tle->launch_year, &tle->launch_number, tle->launch_piece,
        &tle->epoch_year, &tle->epoch_day_of_year,
        &tle->n_dot, &tle->n_ddot, &n_ddot_exp,
        &tle->b_star, &b_star_exp,
        &tle->ephemeris_type, &tle->element_set_number,
        &line1_pre_checksum, &tle->line_1_checksum
    );

    (void) line1_scanned;

    return PERTURB_TLE_PARSE_ERROR_INVALID_INPUT;
}
#endif

#ifndef PERTURB_DISABLE_IO
enum perturb_TleParseError perturb_parse_tle_and_init_sat(
    char * line_1, char * line_2,
    enum perturb_GravityModel grav_model,
    struct perturb_Satellite * sat
) {
    // FIXME: impl
}
#endif

// FIXME: clean up
#if 1
static unsigned int calc_tle_line_checksum(const char *line) {
    unsigned int checksum = 0U;
    for (size_t i = 0; i < (PERTURB_TLE_LINE_LEN - 1); ++i) {
        if (std::isdigit(line[i])) {
            checksum += static_cast<unsigned int>(line[i] - '0');
        }
        if (line[i] == '-') {
            checksum += 1U;
        }
    }
    return (checksum % 10U);
}

#ifndef PERTURB_DISABLE_IO
TLEParseError TwoLineElement::parse(const char *line_1, const char *line_2) {
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
    // FIXME: Change to `strncmp`
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
#endif