/*
 * Implementation of the ideal gas equation of state (EOS) for EOSlib.
 *
 * Author: Christian Reinhardt
 * Date: 15.06.2022
 * Modified: 26.10.2026 (1.1.0)
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "igeos.h"

IGEOSMAT *igeosInitMat(int iMat, double dConstGamma, double dMeanMolMass, double dKpcUnit, double dMsolUnit) {
    /* Gasoline/pkdgrav3 unit constants. Must match the hydro code and the other EOS libraries exactly. */
    const double KBOLTZ = 1.38e-16;      /* Boltzmann constant in cgs */
    const double MHYDR = 1.67e-24;       /* mass of hydrogen atom in grams */
    const double MSOLG = 1.99e33;        /* solar mass in grams */
    const double GCGS = 6.67e-8;         /* G in cgs */
    const double KPCCM = 3.085678e21;    /* kiloparsec in centimeters */
    IGEOSMAT *Mat;

    Mat = (IGEOSMAT *) calloc(1, sizeof(IGEOSMAT));
    assert(Mat != NULL);

    /*
     * Unit system (the same rule in EOSlib and every sublibrary): both units > 0 and finite are code
     * units, both <= 0 is cgs; anything else (exactly one of them <= 0, NaN, inf) is a caller error.
     */
    if (!(((dKpcUnit > 0.0) && (dMsolUnit > 0.0) && isfinite(dKpcUnit) && isfinite(dMsolUnit)) || ((dKpcUnit <= 0.0) && (dMsolUnit <= 0.0)))) {
        fprintf(stderr, "igeosInitMat: invalid unit system dKpcUnit= %g dMsolUnit= %g (give both > 0 for code units or both <= 0 for cgs).\n", dKpcUnit, dMsolUnit);
        assert(0);
        exit(1);
    }

    Mat->iMat = iMat;
    assert(Mat->iMat == 0);

    /* These two parameters define the unit system. */
    Mat->dKpcUnit = dKpcUnit;
    Mat->dMsolUnit = dMsolUnit;

    if (dKpcUnit <= 0.0 && dMsolUnit <= 0.0)
    {
        /* Units are not converted, so the code units are cgs. */
        Mat->dErgPerGmUnit = 1.0;
        Mat->dGmPerCcUnit = 1.0;
        Mat->dSecUnit = 1.0;
    } else {
        /* Conversion factors from code units to cgs, assuming G == 1 in code units. */
        /* code energy per unit mass --> erg per g */
        Mat->dErgPerGmUnit = GCGS*Mat->dMsolUnit*MSOLG/(Mat->dKpcUnit*KPCCM);
        /* code density --> g per cc */
        Mat->dGmPerCcUnit = (Mat->dMsolUnit*MSOLG)/pow(Mat->dKpcUnit*KPCCM,3.0);
        /* code time --> seconds */
        Mat->dSecUnit = sqrt(1/(Mat->dGmPerCcUnit*GCGS));
    }

    /* k_B/m_H in code units (specific energy per Kelvin); same convention as pkdgrav3. */
    Mat->dGasConst = (KBOLTZ/MHYDR)/Mat->dErgPerGmUnit;

    Mat->dConstGamma = dConstGamma;
    Mat->dMeanMolMass = dMeanMolMass;

    /* Specific heat at constant volume, cv = k_B/((gamma-1) mu m_H), in code units. */
    Mat->cv = Mat->dGasConst/((Mat->dConstGamma-1.0)*Mat->dMeanMolMass);

    /* Reference density, hard coded and used as a default value in ballic (cgs -> code units). */
    Mat->rho0 = 0.001/Mat->dGmPerCcUnit;

    return Mat;
}

/* Free memory. Returns IGEOS_SUCCESS. */
int igeosFinalizeMat(IGEOSMAT *Mat) {
    free(Mat);
    return IGEOS_SUCCESS;
}

/* Print material data. */
void igeosPrintMat(IGEOSMAT *Mat, FILE *fp) {
    assert(Mat != NULL);

    fprintf(fp,"# Material: %i (IDEALGAS)\n", Mat->iMat);

    fprintf(fp,"# dConstGamma:  %g\n", Mat->dConstGamma);
    fprintf(fp,"# dMeanMolMass: %g\n", Mat->dMeanMolMass);
    fprintf(fp,"# rho0:         %g\n", Mat->rho0);
    fprintf(fp,"# cv:           %g\n", Mat->cv);
    fprintf(fp,"# Failed lookups and fallbacks (nLookupFail): %ld\n", Mat->nLookupFail);
}

/*
 * The cold curve of the ideal gas is u = 0 (T_min = 0). A state below it (u < 0, or u NaN) falls
 * back to (rho, T = 0) like a state below the cold curve of a tabulated EOS (ANEOSTofRhoU,
 * tillTempRhoU): T = 0, P = 0, c = 0, counted once in nLookupFail. Up to 1.1.0 u < 0 gave T < 0,
 * P < 0 and c = NaN (a NaN time step in a hydro code). Returns the u to evaluate with.
 */
static double igeosColdClamp(IGEOSMAT *Mat, const char *chFunc, double rho, double u) {
    if (u >= 0.0) return u;
    Mat->nLookupFail++;
    IGEOS_WARN("%s: u= %g below the cold curve u = 0 at rho= %g (iMat= %i), using T = 0.\n", chFunc, u, rho, Mat->iMat);
    return 0.0;
}

/* Calculate P(rho, u). P(rho, 0) = 0 below the cold curve (counted). */
double igeosPofRhoU(IGEOSMAT *Mat, double rho, double u) {
    u = igeosColdClamp(Mat, __func__, rho, u);
    return (Mat->dConstGamma-1.0)*rho*u;
}

/* Calculate P(rho, T). IGEOS_FAIL_VALUE (counted) for T < 0 or NaN, like T below a table. */
double igeosPofRhoT(IGEOSMAT *Mat, double rho, double T) {
    if (!(T >= 0.0)) {
        Mat->nLookupFail++;
        IGEOS_WARN("igeosPofRhoT: T= %g < 0 (iMat= %i).\n", T, Mat->iMat);
        return IGEOS_FAIL_VALUE;
    }
    return (Mat->dConstGamma-1.0)*rho*Mat->cv*T;
}

/* Calculate cs(rho, u). 0 below the cold curve (counted). */
double igeosCofRhoU(IGEOSMAT *Mat, double rho, double u) {
    u = igeosColdClamp(Mat, __func__, rho, u);
    return sqrt(Mat->dConstGamma*(Mat->dConstGamma-1.0)*u);
}

/* Calculate P(rho, u) and cs(rho, u) (a state below the cold curve counted once). */
double igeosPCofRhoU(IGEOSMAT *Mat, double rho, double u, double *c) {
    u = igeosColdClamp(Mat, __func__, rho, u);
    *c = sqrt(Mat->dConstGamma*(Mat->dConstGamma-1.0)*u);
    return (Mat->dConstGamma-1.0)*rho*u;
}

/* Calculate P(rho, u), cs(rho, u) and T(rho, u) (a state below the cold curve counted once). */
double igeosPCTofRhoU(IGEOSMAT *Mat, double rho, double u, double *c, double *T) {
    u = igeosColdClamp(Mat, __func__, rho, u);
    *c = sqrt(Mat->dConstGamma*(Mat->dConstGamma-1.0)*u);
    *T = u/Mat->cv;
    return (Mat->dConstGamma-1.0)*rho*u;
}

/*
 * Calculate u2(rho1, u1, rho2) so that (rho1, u1) and (rho2, u2) are on the same isentrope.
 * u1 below the cold curve: the cold curve at rho2, u2 = 0 (counted, like ANEOSisentropicU).
 * IGEOS_FAIL_VALUE (counted) if rho1 or rho2 is not positive and finite.
 */
double igeosIsentropicU(IGEOSMAT *Mat, double rho1, double u1, double rho2) {
    if (!(rho1 > 0.0) || !(rho2 > 0.0) || !isfinite(rho1) || !isfinite(rho2)) {
        Mat->nLookupFail++;
        IGEOS_WARN("igeosIsentropicU: invalid density rho1= %g rho2= %g (iMat= %i).\n", rho1, rho2, Mat->iMat);
        return IGEOS_FAIL_VALUE;
    }
    u1 = igeosColdClamp(Mat, __func__, rho1, u1);
    return u1*pow(rho2/rho1, Mat->dConstGamma-1.0);
}

/* Calculate T(rho, u). Never fails: T = 0 below the cold curve (counted). */
double igeosTofRhoU(IGEOSMAT *Mat, double rho, double u) {
    u = igeosColdClamp(Mat, __func__, rho, u);
    return u/Mat->cv;
}

/* Calculate u(rho, T). IGEOS_FAIL_VALUE (counted) for T < 0 or NaN, like T below a table. */
double igeosUofRhoT(IGEOSMAT *Mat, double rho, double T) {
    if (!(T >= 0.0)) {
        Mat->nLookupFail++;
        IGEOS_WARN("igeosUofRhoT: T= %g < 0 (iMat= %i).\n", T, Mat->iMat);
        return IGEOS_FAIL_VALUE;
    }
    return Mat->cv*T;
}

/*
 * Calculate rho(P, T). IGEOS_FAIL_VALUE (counted) for T <= 0 or P < 0 (or NaN): up to 1.1.0
 * T = 0 gave inf and P < 0 a negative density, both uncounted.
 */
double igeosRhoofPT(IGEOSMAT *Mat, double P, double T) {
    if (!(T > 0.0) || !(P >= 0.0) || !isfinite(P) || !isfinite(T)) {
        Mat->nLookupFail++;
        IGEOS_WARN("igeosRhoofPT: no density for P= %g T= %g (iMat= %i).\n", P, T, Mat->iMat);
        return IGEOS_FAIL_VALUE;
    }
    return P / ((Mat->dConstGamma - 1.0) * Mat->cv * T);
}

/*
 * Calculate rho(u, T). Undefined for an ideal gas with constant cv: u = cv*T
 * does not depend on rho. Calling it is a programming error and therefore
 * fatal (EOSlib convention for unimplemented functions).
 */
double igeosRhoofUT(IGEOSMAT *Mat, double U, double T) {
    fprintf(stderr, "igeosRhoofUT: rho(u, T) is undefined for an ideal gas (iMat= %i, u= %g, T= %g).\n", Mat->iMat, U, T);
    assert(0);
    exit(1);
    return IGEOS_FAIL_VALUE;   /* not reached */
}

/* Calculate dP/drho(rho, u), at (rho, T = 0) below the cold curve (counted). */
double igeosdPdRho(IGEOSMAT *Mat, double rho, double u) {
    u = igeosColdClamp(Mat, __func__, rho, u);
    return (Mat->dConstGamma-1.0)*u;
}

/* Calculate dP/du(rho, u), at (rho, T = 0) below the cold curve (counted, as for a table). */
double igeosdPdU(IGEOSMAT *Mat, double rho, double u) {
    (void) igeosColdClamp(Mat, __func__, rho, u);
    return (Mat->dConstGamma-1.0)*rho;
}

/*
 * Calculate the specific entropy s(rho, T) of an ideal gas with constant cv:
 *
 *   s = cv*ln(T) - (gamma-1)*cv*ln(rho) + const
 *
 * where (gamma-1)*cv = k_B/(mu*m_H). The result is returned in cgs
 * (erg/g/K), consistent with the tabulated EOS libraries, and the additive
 * constant is set to zero. Only entropy differences are meaningful.
 * Returns IGEOS_FAIL_VALUE (counted in nLookupFail) for rho <= 0 or T <= 0.
 */
double igeosSofRhoT(IGEOSMAT *Mat, double rho, double T) {
    double cv_cgs = Mat->cv * Mat->dErgPerGmUnit;
    double rho_cgs = rho * Mat->dGmPerCcUnit;

    if (!(rho_cgs > 0.0) || !(T > 0.0)) {
        Mat->nLookupFail++;
        IGEOS_WARN("igeosSofRhoT: entropy undefined for rho= %g T= %g (iMat= %i).\n", rho, T, Mat->iMat);
        return IGEOS_FAIL_VALUE;
    }

    return cv_cgs * (log(T) - (Mat->dConstGamma - 1.0) * log(rho_cgs));
}

/*
 * Calculate the specific entropy s(rho, u). IGEOS_FAIL_VALUE for u <= 0 (s(rho, T = 0) is
 * undefined) or rho <= 0, counted once (not as a T fallback as well).
 */
double igeosSofRhoU(IGEOSMAT *Mat, double rho, double u) {
    return igeosSofRhoT(Mat, rho, u/Mat->cv);
}