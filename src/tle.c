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
#  include <ctype.h>
#  include <inttypes.h>
#  include <math.h>
#  include <stdio.h>
#  include <string.h>
#endif

#include "common_private.h"

#ifdef __cplusplus
#  error "I kindly request you compile this as C instead of C++"
#endif

static const real_t JULIAN_DATE_1950 = 2433281.5;  ///< [julian date] Jan 1, 1950  00:00:00.000

#ifndef PERTURB_DISABLE_IO
/// Check a string has spaces in specific indices
///
/// @return True if a space is missing (i.e. index occupied by non-space), False if all good.
static bool check_for_missing_spaces(
    const char * line,      ///< Some string, no length checking!
    const size_t * spaces,  ///< Array of specific indices in line to check (assumes valid)
    const size_t n_spaces   ///< Length of `spaces` array
) {
    for (size_t i = 0U; i < n_spaces; ++i)
    {
        // SAFETY: Assume `n_spaces` is accurate and valid length of `spaces`
        const size_t idx_space = spaces[i];

        // SAFETY: Assume indices from `spaces` are valid. Should be fixed at compile-time.
        if (line[idx_space] != ' ')
        {
            return true;
        }
    }
    return false;
}
#endif  // PERTURB_DISABLE_IO

#ifndef PERTURB_DISABLE_IO
/// Compute [0, 10] integer checksum of a TLE line based on their standard rules
static uint8_t calc_tle_line_checksum(const char * line)
{
    unsigned int checksum = 0U;
    for (size_t i = 0U; i < (PERTURB_TLE_LINE_LEN - 1U); ++i)
    {
        // SAFETY: Assumes `line` is valid buffer of at least length `PERTURB_TLE_LINE_LEN`
        if (isdigit(line[i]))
        {
            checksum += (unsigned int) (line[i] - '0');
        }
        else if (line[i] == '-')
        {
            checksum += 1U;
        }
        else
        {
            // Nothing
        }
    }
    return (uint8_t) (checksum % 10U);
}
#endif  // PERTURB_DISABLE_IO

enum perturb_Sgp4Error perturb_init_sat_from_tle(
    const struct perturb_TwoLineElement tle,
    const enum perturb_GravityModel grav_model,
    struct perturb_Satellite * sat
) {
    if (sat == NULL)
    {
        return PERTURB_SGP4_ERROR_INVALID_INPUT;
    }

    // Line 1
    memcpy(sat->satnum, tle.catalog_number, sizeof(sat->satnum));
    sat->classification = tle.classification;
    UNUSED(sat->intldesg);  // Don't bother converting to set `sat.intldesg` b/c it has no effects
    sat->epochyr = tle.epoch_year;
    sat->epochdays = tle.epoch_day_of_year;
    sat->ndot = tle.n_dot;
    sat->nddot = tle.n_ddot;
    sat->bstar = tle.b_star;
    sat->ephtype = tle.ephemeris_type;
    sat->elnum = tle.element_set_number;

    // Line 2
    sat->inclo = tle.inclination;
    sat->nodeo = tle.raan;
    sat->ecco = tle.eccentricity;
    sat->argpo = tle.arg_of_perigee;
    sat->mo = tle.mean_anomaly;
    sat->no_kozai = tle.mean_motion;
    sat->revnum = tle.revolution_number;

    // Post-process same as how Vallado does it
    sat->error = 0;
    sat->no_kozai /= XP_DOT_P;  // [rad/min]
    sat->ndot  /= (XP_DOT_P * MINS_PER_DAY);
    sat->nddot /= ((XP_DOT_P * MINS_PER_DAY) * MINS_PER_DAY);
    sat->inclo *= DEG_TO_RAD;
    sat->nodeo *= DEG_TO_RAD;
    sat->argpo *= DEG_TO_RAD;
    sat->mo    *= DEG_TO_RAD;

    // Convert (year & fractional day of year) to (date & time)
    struct perturb_DateTime t = { 0 };
    t.year = (uint16_t) (sat->epochyr + ((sat->epochyr < 57) ? 2000 : 1900));
    days2mdhms_SGP4(t.year, sat->epochdays, &t.month, &t.day, &t.hour, &t.min, &t.sec);

    // Convert (date & time) to (julian date)
    const struct perturb_JulianDate jd = perturb_datetime_to_julian(t);
    // TODO: Make Satellite just store `JulianDate` instead of seperate `jdsatepoch`
    sat->jdsatepoch = jd.jd;
    sat->jdsatepochF = jd.jd_frac;
    const real_t days_since_1950 = (jd.jd - JULIAN_DATE_1950) + jd.jd_frac;

    // Initialize orbit
    const bool ret = sgp4init(
        grav_model, 'i', sat->satnum, days_since_1950, sat->bstar,
        sat->ndot, sat->nddot, sat->ecco, sat->argpo, sat->inclo,
        sat->mo, sat->no_kozai, sat->nodeo, sat
    );
    UNUSED(ret);  // Not used, error code is saved internally

    return (enum perturb_Sgp4Error) sat->error;
}

#ifndef PERTURB_DISABLE_IO
enum perturb_TleParseError perturb_parse_tle(
    const char * line_1, const char * line_2,
    struct perturb_TwoLineElement * tle
){
    const bool bad_ptrs = (line_1 == NULL) || (line_2 == NULL) || (tle == NULL);
    if (bad_ptrs)
    {
        return PERTURB_TLE_PARSE_ERROR_INVALID_INPUT;
    }

    // Make sure there are spaces in the right places
    const size_t LINE_1_SPACES[] = { 1, 8, 17, 32, 43, 52, 61, 64 };
    const size_t LINE_2_SPACES[] = { 1, 7, 16, 25, 33, 42, 51 };

    const bool is_space_missing = (
        check_for_missing_spaces(line_1, LINE_1_SPACES, ARRAY_SIZE(LINE_1_SPACES)) ||
        check_for_missing_spaces(line_2, LINE_2_SPACES, ARRAY_SIZE(LINE_2_SPACES))
    );
    if (is_space_missing)
    {
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

    // Parse format - Line 2
    const char LINE_2_FMT_STR_NO_SPACE[] = (
        "%1hhu %5s "
        "%8lf %8lf %7lu %8lf "
        "%8lf %11lf %5" SCNu32 " "
        "%n %1hhu"
    );
    const char LINE_2_FMT_STR_WT_SPACE[] = (
        "%1hhu %5s "
        "%8lf %8lf %7lu %8lf "
        "%8lf %10lf %5" SCNu32 " "  // Only diff from `LINE_2_FMT_STR_NO_SPACE`
        "%n %1hhu"
    );

    char line2_num = '0';
    char catlog_num_line2[6] = "-----";
    unsigned long eccentricity_int = 0U;
    unsigned int line2_pre_checksum = 0U;
    int line2_scanned = 0;

    // :( yeah, I know, this is annoying...
    if (line_2[52] != ' ')
    {
        line2_scanned = sscanf(
            line_2, LINE_2_FMT_STR_NO_SPACE,
            &line2_num, catlog_num_line2,
            &tle->inclination, &tle->raan, &eccentricity_int, &tle->arg_of_perigee,
            &tle->mean_anomaly, &tle->mean_motion, &tle->revolution_number,
            &line2_pre_checksum, &tle->line_2_checksum
        );
    }
    else
    {
        line2_scanned = sscanf(
            line_2, LINE_2_FMT_STR_WT_SPACE,
            &line2_num, catlog_num_line2,
            &tle->inclination, &tle->raan, &eccentricity_int, &tle->arg_of_perigee,
            &tle->mean_anomaly, &tle->mean_motion, &tle->revolution_number,
            &line2_pre_checksum, &tle->line_2_checksum
        );
    }

    // Ensure all strings are null terminated
    tle->catalog_number[5] = '\0';
    tle->launch_piece[3] = '\0';
    catlog_num_line2[5] = '\0';

    // Handle some annoying special cases
    const bool fix_l1_checksum = (
        (line1_scanned == 15) &&
        (line1_pre_checksum >= PERTURB_TLE_LINE_LEN) &&
        (line_1[PERTURB_TLE_LINE_LEN - 5] == ' ')
    );

    // Element set number often doesn't have leading zeroes,
    // so `sscanf` over-consumes and eats up the final checksum.
    if (fix_l1_checksum)
    {
        tle->element_set_number /= 10U;  // Chop of least-significant digit
        tle->line_1_checksum = (uint8_t) (line_1[PERTURB_TLE_LINE_LEN - 1] - '0');
        line1_scanned += 1;
    }

    const bool fix_l2_checksum = (
        (line2_scanned == 9) &&
        (line2_pre_checksum >= PERTURB_TLE_LINE_LEN) &&
        (line_2[PERTURB_TLE_LINE_LEN - 6] == ' ')
    );

    // If revolution number doens't have leading zero, so same issue as above
    if (fix_l2_checksum)
    {
        tle->revolution_number /= 10U;
        tle->line_2_checksum = (uint8_t) (line_2[PERTURB_TLE_LINE_LEN - 1] - '0');
        line2_scanned += 1;
    }

    // Check that the correct number of values were parsed
    if (line1_scanned != 16 || line2_scanned != 10)
    {
        return PERTURB_TLE_PARSE_ERROR_INVALID_FORMAT;
    }

    // Handle implicit decimal for exponents
    if (line_1[44] != '.')
    {
        n_ddot_exp -= 5;
    }
    if (line_1[53] != '.')
    {
        b_star_exp -= 5;
    }

    // Check that valid values were parsed
    bool valid_vals = true;

    // Line 1
    valid_vals &= (line1_num == 1);
    const char clsf = tle->classification;
    valid_vals &= (clsf == 'U') || (clsf == 'C') || (clsf == 'S');
    valid_vals &= (tle->launch_year < 100U) && (tle->epoch_year < 100U);
    valid_vals &= (1.0 <= tle->epoch_day_of_year) && (tle->epoch_day_of_year <= 366.0);
    valid_vals &= (-15 < n_ddot_exp) && (n_ddot_exp < 10);
    valid_vals &= (-15 < b_star_exp) && (b_star_exp < 10);
    valid_vals &= (tle->ephemeris_type == 0U);
    valid_vals &= (tle->element_set_number < 10000U);

    // Line 2
    valid_vals &= (line2_num == 2);
    valid_vals &= (strcmp(tle->catalog_number, catlog_num_line2) == 0);  // TODO: Change to `strncmp`
    valid_vals &= (0.0 <= tle->inclination) && (tle->inclination <= 180.0);
    valid_vals &= (0.0 <= tle->raan) && (tle->raan <= 360.0);
    valid_vals &= (0.0 <= tle->arg_of_perigee) && (tle->arg_of_perigee <= 360.0);
    valid_vals &= (0.0 <= tle->mean_anomaly) && (tle->mean_anomaly <= 360.0);

    if (!valid_vals)
    {
        return PERTURB_TLE_PARSE_ERROR_INVALID_VALUE;
    }

    // Post-process
    tle->n_ddot *= pow(10.0, n_ddot_exp);
    tle->b_star *= pow(10.0, b_star_exp);
    tle->eccentricity = ((real_t) eccentricity_int) / ((real_t) 1.0e7);

    // Calculate and compare checksums
    const bool checksum_matches = (
        (calc_tle_line_checksum(line_1) == tle->line_1_checksum) &&
        (calc_tle_line_checksum(line_2) == tle->line_2_checksum)
    );
    if (!checksum_matches)
    {
        return PERTURB_TLE_PARSE_ERROR_CHECKSUM_MISMATCH;
    }

    return PERTURB_TLE_PARSE_ERROR_NONE;
}
#endif  // PERTURB_DISABLE_IO
