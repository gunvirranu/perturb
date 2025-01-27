/*
 * perturb -- A modern C++11 wrapper for the SGP4 orbit propagator
 * version 1.0.0
 * https://github.com/gunvirranu/perturb
 *
 * Copyright (c) 2022 Gunvir Singh Ranu
 * SPDX-License-Identifier: MIT
 */

#include "perturb/perturb.h"

#include <math.h>

#include "perturb/sgp4.h"
#include "common_private.h"

#ifdef __cplusplus
// You might be glob compiling this the src/ directory as C++.
// Build the .c files as C and .cpp as C++, should be good.
#  error "Hah someone messed up, why r u compiling this C as C++"
#endif

// FIXME: Reformat files for C and C++

struct perturb_JulianDate perturb_datetime_to_julian(const struct perturb_DateTime t)
{
    struct perturb_JulianDate jd;
    jday_SGP4(t.year, t.month, t.day, t.hour, t.min, t.sec, &jd.jd, &jd.jd_frac);
    return jd;
}

struct perturb_DateTime perturb_julian_to_datetime(const struct perturb_JulianDate jd)
{
    struct perturb_DateTime t;
    invjday_SGP4(jd.jd, jd.jd_frac, &t.year, &t.month, &t.day, &t.hour, &t.min, &t.sec);
    return t;
}

struct perturb_JulianDate perturb_julian_normalized(const struct perturb_JulianDate t)
{
    struct perturb_JulianDate out = t;

    // Check for fractional days included in `jd` and put them in `jd`
    const real_t frac_days = t.jd - floor(t.jd) - 0.5;

    if (FABS(frac_days) > 1e-12)
    {
        out.jd -= frac_days;
        out.jd_frac += frac_days;
    }

    // Check for whole days in `jd_frac` and put them in `jd`
    if (FABS(out.jd_frac) >= 1.0)
    {
        const real_t whole_days = FLOOR(out.jd_frac);
        out.jd += whole_days;
        out.jd_frac -= whole_days;
    }
    return out;
}

struct perturb_JulianDate perturb_julian_add_days(
    const struct perturb_JulianDate t, const perturb_real_t days
) {
    struct perturb_JulianDate t_new = t;
    // Just add entire offset to fractional value
    // Can be normalized later explicitly if needed
    t_new.jd_frac += days;
    return t_new;
}

/// lhs - rhs
perturb_real_t perturb_julian_subtract(
    const struct perturb_JulianDate lhs, const struct perturb_JulianDate rhs
) {
    // Grouping here is important to preserve precision
    return (lhs.jd - rhs.jd) + (lhs.jd_frac - rhs.jd_frac);
}

struct perturb_OrbitalElements perturb_state_vector_to_orbital_elements(const struct perturb_StateVector sv)
{
    // FIXME: Explain why this default
    const enum perturb_GravityModel grav_model = PERTURB_GRAVITY_MODEL_WGS72;
    return perturb_state_vector_to_orbital_elements_with_grav(sv, grav_model);
}

struct perturb_OrbitalElements perturb_state_vector_to_orbital_elements_with_grav(
    struct perturb_StateVector sv, enum perturb_GravityModel grav_model
) {
    // Fetch constants that depend on gravity model
    real_t mus, _tumin, _rekm, _xke, _j2, _j3, _j4, _j3oj2;
    getgravconst(grav_model, &_tumin, &mus, &_rekm, &_xke, &_j2, &_j3, &_j4, &_j3oj2);

    struct perturb_OrbitalElements elems = { 0 };
    rv2coe_SGP4(
        sv.position, sv.velocity, mus, &elems.semilatus_rectum, &elems.semimajor_axis,
        &elems.eccentricity, &elems.inclination, &elems.raan, &elems.arg_of_perigee, &elems.true_anomaly, &elems.mean_anomaly,
        &elems.arg_of_latitude, &elems.true_longitude, &elems.longitude_of_periapsis
    );
    return elems;
}

#ifdef __cplusplus
}  // extern "C"
#  ifdef PERTURB_ENABLE_CPP_INTERFACE
}  // namespace c_internal
}  // namespace perturb
#  endif
#endif
