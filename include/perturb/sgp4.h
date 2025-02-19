/*
 * perturb -- A modern C++11 wrapper for the SGP4 orbit propagator
 * version 1.0.0
 * https://github.com/gunvirranu/perturb
 *
 * Copyright (c) 2024 Gunvir Singh Ranu
 * SPDX-License-Identifier: MIT
 */

//! @file Header for the internal modified SGP4 impl from Vallado
//! @author David Vallado
//! @author Gunvir Singh Ranu
//! @version SGP4 Version 2020-07-13
//! @date 2024-09-12

/**    ----------------------------------------------------------------
 *
 *                                 SGP4.h
 *
 *    this file contains the sgp4 procedures for analytical propagation
 *    of a satellite. the code was originally released in the 1980 and 1986
 *    spacetrack papers. a detailed discussion of the theory and history
 *    may be found in the 2006 aiaa paper by vallado, crawford, hujsak,
 *    and kelso.
 *
 *    current :
 *              12 mar 20  david vallado
 *                           chg satnum to string for alpha 5 or 9-digit
 *    changes :
 *               7 dec 15  david vallado
 *                           fix jd, jdfrac
 *               3 nov 14  david vallado
 *                           update to msvs2013 c++
 *              30 Dec 11  david vallado
 *                           consolidate updated code version
 *              30 Aug 10  david vallado
 *                           delete unused variables in initl
 *                           replace pow inetger 2, 3 with multiplies for speed
 *               3 Nov 08  david vallado
 *                           put returns in for error codes
 *              29 sep 08  david vallado
 *                           fix atime for faster operation in dspace
 *                           add operationmode for afspc (a) or improved (i)
 *                           performance mode
 *              20 apr 07  david vallado
 *                           misc fixes for constants
 *              11 aug 06  david vallado
 *                           chg lyddane choice back to strn3, constants, misc doc
 *              15 dec 05  david vallado
 *                           misc fixes
 *              26 jul 05  david vallado
 *                           fixes for paper
 *                           note that each fix is preceded by a
 *                           comment with "sgp4fix" and an explanation of
 *                           what was changed
 *              10 aug 04  david vallado
 *                           2nd printing baseline working
 *              14 may 01  david vallado
 *                           2nd edition baseline
 *                     80  norad
 *                           original baseline
 *       ----------------------------------------------------------------      */

#ifndef PERTURB_SGP4_H
#define PERTURB_SGP4_H

#include <stdbool.h>
#include <stdint.h>

#include "perturb/perturb.h"

// Define `PERTURB_SGP4_ENABLE_DEBUG` to enable Vallado's verification mode.
// Generally, this is unwanted. Only enabled for tests. Requires IO.
#if (defined(PERTURB_SGP4_ENABLE_DEBUG) && defined(PERTURB_DISABLE_IO))
#  error "Cannot enable SGP4 debug without I/O functionality"
#endif

#ifdef __cplusplus
#  ifdef PERTURB_ENABLE_CPP_INTERFACE
namespace perturb {
namespace c_internal {
#  endif
extern "C" {
#endif

/// Possible issues during SGP4 propagation or even TLE parsing.
///
/// This is important in two places:
///   1. After construction of a `Satellite`, check `Satellite::last_error`
///      for any issues with TLE parsing or SGP4 initialization.
///   2. Calling `Satellite::propagate` returns a `perturb::Sgp4Error`,
///      indicating any possible issues with propagation.
///
/// If everything is all good, the value should be `Sgp4Error::NONE`.
/// The errors `Sgp4Error::MEAN_ELEMENTS` to `SGP4Error::DECAYED` directly
/// correlate to errors in the underlying SGP4 impl, from the comments of
/// `perturb::sgp4::sgp4`. The additional `Sgp4Error::INVALID_TLE` is for issues
/// with reading the TLE strings.
enum perturb_Sgp4Error {
    PERTURB_SGP4_ERROR_NONE,
    PERTURB_SGP4_ERROR_MEAN_ELEMENTS,
    PERTURB_SGP4_ERROR_MEAN_MOTION,
    PERTURB_SGP4_ERROR_PERT_ELEMENTS,
    PERTURB_SGP4_ERROR_SEMI_LATUS_RECTUM,
    PERTURB_SGP4_ERROR_EPOCH_ELEMENTS_SUB_ORBITAL,
    PERTURB_SGP4_ERROR_DECAYED,
    PERTURB_SGP4_ERROR_INVALID_INPUT,  ///< Not from base impl, added in
    PERTURB_SGP4_ERROR_UNKNOWN
};

/// Represents a specific orbital ephemeris for an Earth-centered trajectory.
///
/// This is the primary type in this library. Wraps the internal SGP4 record
/// type `perturb::sgp4::elsetrec`. Generally constructed via TLEs through
/// `Satellite::from_tle` constructors. Of particular importance is the
/// `Satellite::last_error` method which you can check to determine if there
/// were any issues with TLE initialization or propagation. The primary method
/// of running the SGP4 algorithm is the `Satellite::propagate` method.
struct perturb_Satellite {
    char satnum[6];

    int epochyr;
    int epochtynumrev;  ///< ???

    int error;  ///< Latched error code from previous function call

    bool operationmode;  ///< True if improved (i), False if AFSPC (a)
    bool init;           ///< True if initialised
    char method;         ///< SGP4 algorith method: deep space (d) or not (n)

    /* Near Earth */
    bool isimp;
    perturb_real_t aycof, con41, cc1, cc4, cc5, d2, d3, d4, delmo, eta, argpdot, omgcof,
        sinmao, t, t2cof, t3cof, t4cof, t5cof, x1mth2, x7thm1, mdot, nodedot, xlcof,
        xmcof, nodecf;

    /* Deep Space */
    int irez;  ///< Flag for resonance (0 = none, 1 = one day, 2 = half day)
    perturb_real_t d2201, d2211, d3210, d3222, d4410, d4422, d5220, d5232, d5421, d5433,
        dedt, del1, del2, del3, didt, dmdt, dnodt, domdt, e3, ee2, peo, pgho, pho, pinco,
        plo, se2, se3, sgh2, sgh3, sgh4, sh2, sh3, si2, si3, sl2, sl3, sl4, gsto, xfact,
        xgh2, xgh3, xgh4, xh2, xh3, xi2, xi3, xl2, xl3, xl4, xlamo, zmol, zmos, atime,
        xli, xni;

    perturb_real_t a, altp, alta, epochdays, jdsatepoch, jdsatepochF, nddot, ndot, bstar,
        rcse, inclo, nodeo, ecco, argpo, mo, no_kozai;

    // sgp4fix add new variables from tle
    char classification, intldesg[11];
    int ephtype;

    long elnum, revnum;

    // sgp4fix add unkozai'd variable
    perturb_real_t no_unkozai;  ///< Mean motion of satellite

    // sgp4fix add singly averaged variables
    perturb_real_t am, em, im, Om, om, mm, nm;

    // sgp4fix add constant parameters to eliminate mutliple calls during execution
    perturb_real_t tumin, mus, radiusearthkm, xke, j2, j3, j4, j3oj2;
};

struct perturb_JulianDate perturb_epoch(const struct perturb_Satellite * sat);

enum perturb_Sgp4Error perturb_propagate(
    struct perturb_Satellite * sat, struct perturb_JulianDate t, struct perturb_StateVector * sv
);
enum perturb_Sgp4Error perturb_propagate_days_from_epoch(
    struct perturb_Satellite * sat, perturb_real_t days_since_epoch, struct perturb_StateVector * sv
);

bool sgp4init(
    enum perturb_GravityModel whichconst, char opsmode, const char satn[5], const perturb_real_t epoch,
    const perturb_real_t xbstar, const perturb_real_t xndot, const perturb_real_t xnddot, const perturb_real_t xecco,
    const perturb_real_t xargpo, const perturb_real_t xinclo, const perturb_real_t xmo, const perturb_real_t xno,
    const perturb_real_t xnodeo, struct perturb_Satellite * satrec
);

bool sgp4(
    // no longer need gravconsttype whichconst, all data contained in satrec
    struct perturb_Satellite * satrec, perturb_real_t tsince, perturb_real_t r[3], perturb_real_t v[3]
);

void getgravconst(
    enum perturb_GravityModel whichconst, perturb_real_t * tumin, perturb_real_t * mus, perturb_real_t * radiusearthkm,
    perturb_real_t * xke, perturb_real_t * j2, perturb_real_t * j3, perturb_real_t * j4, perturb_real_t * j3oj2
);

#ifndef PERTURB_DISABLE_IO
// older sgp4io methods
void twoline2rv(
    char longstr1[130], char longstr2[130], char typerun, char typeinput, char opsmode,
    enum perturb_GravityModel whichconst, perturb_real_t * startmfe, perturb_real_t * stopmfe, perturb_real_t * deltamin,
    struct perturb_Satellite * satrec
);
#endif  // PERTURB_DISABLE_IO

perturb_real_t gstime_SGP4(perturb_real_t jdut1);

perturb_real_t sgn_SGP4(perturb_real_t x);

perturb_real_t mag_SGP4(perturb_real_t x[3]);

void cross_SGP4(perturb_real_t vec1[3], perturb_real_t vec2[3], perturb_real_t outvec[3]);

perturb_real_t dot_SGP4(perturb_real_t x[3], perturb_real_t y[3]);

perturb_real_t angle_SGP4(perturb_real_t vec1[3], perturb_real_t vec2[3]);

void newtonnu_SGP4(perturb_real_t ecc, perturb_real_t nu, perturb_real_t * e0, perturb_real_t * m);

perturb_real_t asinh_SGP4(perturb_real_t xval);

void rv2coe_SGP4(
    const perturb_real_t r[3], const perturb_real_t v[3], perturb_real_t mus, perturb_real_t * p, perturb_real_t * a, perturb_real_t * ecc,
    perturb_real_t * incl, perturb_real_t * omega, perturb_real_t * argp, perturb_real_t * nu, perturb_real_t * m, perturb_real_t * arglat,
    perturb_real_t * truelon, perturb_real_t * lonper
);

void jday_SGP4(
    uint16_t year, uint8_t mon, uint8_t day, uint8_t hr, uint8_t minute, perturb_real_t sec,
    perturb_real_t * jd, perturb_real_t * jdFrac
);

void days2mdhms_SGP4(
    uint16_t year, perturb_real_t days, uint8_t * mon, uint8_t * day, uint8_t * hr, uint8_t * minute, perturb_real_t * sec
);

void invjday_SGP4(
    perturb_real_t jd, perturb_real_t jdFrac, uint16_t * year, uint8_t * mon, uint8_t * day, uint8_t * hr, uint8_t * minute,
    perturb_real_t * sec
);

#ifdef __cplusplus
}  // extern "C"
#  ifdef PERTURB_ENABLE_CPP_INTERFACE
}  // namespace c_internal
}  // namespace perturb
#  endif
#endif

#endif  // PERTURB_SGP4_H
