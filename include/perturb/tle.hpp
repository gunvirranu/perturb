//! @file
//! Header for custom TLE (two line element) processing
//! @author Gunvir Ranu
//! @version 0.0.0
//! @copyright Gunvir Ranu, MIT License

#ifndef PERTURB_TLE_HPP
#define PERTURB_TLE_HPP

#include <cstddef>
#ifndef PERTURB_DISABLE_IO
#  include <string>
#endif

namespace perturb {

/// Both lines of a TLE **must** be this length, for TLE constructors.
///
/// The memory can and is accessed.
/// Lines can be longer for verification mode, but that's for internal testing
/// purposes only and doesn't pertain to general usage.
constexpr std::size_t TLE_LINE_LEN = 69;

enum class TLEParseError {
    NONE,
    SHOULD_BE_SPACE,
    INVALID_FORMAT,
    INVALID_VALUE,
    CHECKSUM_MISMATCH,
};

struct TwoLineElement {
    // Line 1
    char catalog_number[6];
    char classification;
    unsigned int launch_year;
    unsigned int launch_number;
    char launch_piece[4];
    unsigned int epoch_year;
    double epoch_day_of_year;
    double n_dot;
    double n_ddot;
    double b_star;
    unsigned char ephemeris_type;
    unsigned int element_set_number;
    unsigned char line_1_checksum;

    // Line 2
    double inclination;
    double raan;
    double eccentricity;
    double arg_of_perigee;
    double mean_anomaly;
    double mean_motion;
    unsigned long revolution_number;
    unsigned char line_2_checksum;

#ifndef PERTURB_DISABLE_IO
    TLEParseError parse(const char *line_1, const char *line_2);
#endif  // PERTURB_DISABLE_IO

#ifndef PERTURB_DISABLE_IO
    TLEParseError parse(const std::string &line_1, const std::string &line_2);
#endif  // PERTURB_DISABLE_IO
};

}  // namespace perturb

#endif  // PERTURB_TLE_HPP
