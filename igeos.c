/*
 * Implementation of the ideal gas equation of state (EOS) for EOSlib.
 *
 * Author: Christian Reinhardt
 * Date: 15.06.2022
 * Modified: 28.10.2026 (1.1.0)
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

    /*
     * Unit system (the same rule in EOSlib and every sublibrary): both units > 0 and finite are code
     * units, both <= 0 is cgs; anything else (exactly one of them <= 0, NaN, inf) is a caller error.
     */
    if (!(((dKpcUnit > 0.0) && (dMsolUnit > 0.0) && isfinite(dKpcUnit) && isfinite(dMsolUnit)) || ((dKpcUnit <= 0.0) && (dMsolUnit <= 0.0)))) {
        fprintf(stderr, "igeosInitMat: invalid unit system dKpcUnit= %g dMsolUnit= %g (give both > 0 for code units or both <= 0 for cgs).\n", dKpcUnit, dMsolUnit);
        assert(0);
        exit(1);
    }

    /*
     * The material parameters (fatal, like a defective table): cv = k_B/((gamma-1) mu m_H) needs
     * gamma > 1 and mu > 0, both finite. Up to 1.1.0 gamma <= 1 gave cv = inf or cv < 0, mu <= 0 or
     * NaN cv = inf or NaN, and every lookup after that silently returned garbage.
     */
    if (!(dConstGamma > 1.0) || !isfinite(dConstGamma) || !(dMeanMolMass > 0.0) || !isfinite(dMeanMolMass)) {
        fprintf(stderr, "igeosInitMat: invalid ideal gas parameters dConstGamma= %g dMeanMolMass= %g (need gamma > 1 and mu > 0, both finite).\n",
                dConstGamma, dMeanMolMass);
        assert(0);
        exit(1);
    }

    Mat = (IGEOSMAT *) calloc(1, sizeof(IGEOSMAT));
    assert(Mat != NULL);

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

    /* Dimensional values in cgs with the unit, like ANEOSPrintMat, scvheosPrintMat and reos3PrintMat. */
    fprintf(fp,"# Material: %i (IDEALGAS)\n", Mat->iMat);
    fprintf(fp,"# dKpcUnit: %g dMsolUnit: %g\n", Mat->dKpcUnit, Mat->dMsolUnit);
    fprintf(fp,"# dConstGamma:  %g\n", Mat->dConstGamma);
    fprintf(fp,"# dMeanMolMass: %g (in m_H)\n", Mat->dMeanMolMass);
    fprintf(fp,"# Reference density rho0 (cgs units): %g g/cc\n", Mat->rho0*Mat->dGmPerCcUnit);
    fprintf(fp,"# cv (cgs units): %g erg/g/K\n", Mat->cv*Mat->dErgPerGmUnit);
    fprintf(fp,"# Failed lookups and fallbacks (nLookupFail): %ld\n", Mat->nLookupFail);
}

/*
 * A valid density: positive and finite. Like a density outside the table of a tabulated EOS, an
 * invalid rho is a failure of every function that takes one (counted once, by the caller), and
 * T(rho, u) falls back to T_min = 0. Up to 1.1.0 rho <= 0 gave a negative or zero P and c, and
 * EOSdUdRho a NaN at rho = 0, uncounted.
 */
static inline int igeosValidRho(double rho) {
    return (rho > 0.0) && isfinite(rho);
}

/*
 * The state check of the (rho, u) functions: returns 0 and the u to evaluate with, or -1 for an
 * invalid rho (counted). The cold curve of the ideal gas is u = 0 (T_min = 0). A state below it
 * (u < 0, or u NaN) falls back to (rho, T = 0) like a state below the cold curve of a tabulated EOS
 * (ANEOSTofRhoU, tillTempRhoU): T = 0, P = 0, c = 0, counted once in nLookupFail. Up to 1.1.0
 * u < 0 gave T < 0, P < 0 and c = NaN (a NaN time step in a hydro code). An invalid rho is counted
 * once, not as a cold-curve fallback as well.
 */
static int igeosCheckState(IGEOSMAT *Mat, const char *chFunc, double rho, double *pu) {
    if (!igeosValidRho(rho)) {
        Mat->nLookupFail++;
        IGEOS_WARN("%s: invalid density rho= %g (u= %g, iMat= %i).\n", chFunc, rho, *pu, Mat->iMat);
        return -1;
    }
    if (*pu >= 0.0) return 0;
    Mat->nLookupFail++;
    IGEOS_WARN("%s: u= %g below the cold curve u = 0 at rho= %g (iMat= %i), using T = 0.\n", chFunc, *pu, rho, Mat->iMat);
    *pu = 0.0;
    return 0;
}

/*
 * The state check of the (rho, T) functions: 0, or -1 (counted) for an invalid rho or T < 0
 * (below T_min = 0) or a NaN T, like (rho, T) outside a table.
 */
static int igeosCheckRhoT(IGEOSMAT *Mat, const char *chFunc, double rho, double T) {
    if (igeosValidRho(rho) && T >= 0.0) return 0;
    Mat->nLookupFail++;
    IGEOS_WARN("%s: invalid state rho= %g T= %g (iMat= %i).\n", chFunc, rho, T, Mat->iMat);
    return -1;
}

/* Calculate P(rho, u). P(rho, 0) = 0 below the cold curve (counted); IGEOS_FAIL_VALUE for an invalid rho. */
double igeosPofRhoU(IGEOSMAT *Mat, double rho, double u) {
    if (igeosCheckState(Mat, __func__, rho, &u)) return IGEOS_FAIL_VALUE;
    return (Mat->dConstGamma-1.0)*rho*u;
}

/* Calculate P(rho, T). IGEOS_FAIL_VALUE (counted) for T < 0, NaN or an invalid rho, like (rho, T) outside a table. */
double igeosPofRhoT(IGEOSMAT *Mat, double rho, double T) {
    if (igeosCheckRhoT(Mat, __func__, rho, T)) return IGEOS_FAIL_VALUE;
    return (Mat->dConstGamma-1.0)*rho*Mat->cv*T;
}

/* Calculate cs(rho, u). 0 below the cold curve (counted); IGEOS_FAIL_VALUE for an invalid rho. */
double igeosCofRhoU(IGEOSMAT *Mat, double rho, double u) {
    if (igeosCheckState(Mat, __func__, rho, &u)) return IGEOS_FAIL_VALUE;
    return sqrt(Mat->dConstGamma*(Mat->dConstGamma-1.0)*u);
}

/* Calculate P(rho, u) and cs(rho, u) (a bad state counted once); both IGEOS_FAIL_VALUE for an invalid rho. */
double igeosPCofRhoU(IGEOSMAT *Mat, double rho, double u, double *c) {
    if (igeosCheckState(Mat, __func__, rho, &u)) {
        *c = IGEOS_FAIL_VALUE;
        return IGEOS_FAIL_VALUE;
    }
    *c = sqrt(Mat->dConstGamma*(Mat->dConstGamma-1.0)*u);
    return (Mat->dConstGamma-1.0)*rho*u;
}

/*
 * Calculate P(rho, u), cs(rho, u) and T(rho, u) (a bad state counted once). For an invalid rho P
 * and c are IGEOS_FAIL_VALUE and T = 0 (T(rho, u) never fails), as in the tabulated libraries.
 */
double igeosPCTofRhoU(IGEOSMAT *Mat, double rho, double u, double *c, double *T) {
    if (igeosCheckState(Mat, __func__, rho, &u)) {
        *c = IGEOS_FAIL_VALUE;
        *T = 0.0;
        return IGEOS_FAIL_VALUE;
    }
    *c = sqrt(Mat->dConstGamma*(Mat->dConstGamma-1.0)*u);
    *T = u/Mat->cv;
    return (Mat->dConstGamma-1.0)*rho*u;
}

/*
 * Calculate u2(rho1, u1, rho2) so that (rho1, u1) and (rho2, u2) are on the same isentrope.
 * u1 below the cold curve: the cold curve at rho2, u2 = 0 (counted, like ANEOSisentropicU).
 * IGEOS_FAIL_VALUE (counted once) if rho1 or rho2 is not positive and finite.
 */
double igeosIsentropicU(IGEOSMAT *Mat, double rho1, double u1, double rho2) {
    if (!igeosValidRho(rho2)) {
        Mat->nLookupFail++;
        IGEOS_WARN("igeosIsentropicU: invalid density rho1= %g rho2= %g (iMat= %i).\n", rho1, rho2, Mat->iMat);
        return IGEOS_FAIL_VALUE;
    }
    if (igeosCheckState(Mat, __func__, rho1, &u1)) return IGEOS_FAIL_VALUE;
    return u1*pow(rho2/rho1, Mat->dConstGamma-1.0);
}

/* Calculate T(rho, u). Never fails: T = 0 below the cold curve or for an invalid rho (counted). */
double igeosTofRhoU(IGEOSMAT *Mat, double rho, double u) {
    if (igeosCheckState(Mat, __func__, rho, &u)) return 0.0;
    return u/Mat->cv;
}

/* Calculate u(rho, T). IGEOS_FAIL_VALUE (counted) for T < 0, NaN or an invalid rho, like (rho, T) outside a table. */
double igeosUofRhoT(IGEOSMAT *Mat, double rho, double T) {
    if (igeosCheckRhoT(Mat, __func__, rho, T)) return IGEOS_FAIL_VALUE;
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

/* Calculate dP/drho(rho, u), at (rho, T = 0) below the cold curve (counted); IGEOS_FAIL_VALUE for an invalid rho. */
double igeosdPdRho(IGEOSMAT *Mat, double rho, double u) {
    if (igeosCheckState(Mat, __func__, rho, &u)) return IGEOS_FAIL_VALUE;
    return (Mat->dConstGamma-1.0)*u;
}

/* Calculate dP/du(rho, u), at (rho, T = 0) below the cold curve (counted, as for a table); IGEOS_FAIL_VALUE for an invalid rho. */
double igeosdPdU(IGEOSMAT *Mat, double rho, double u) {
    if (igeosCheckState(Mat, __func__, rho, &u)) return IGEOS_FAIL_VALUE;
    return (Mat->dConstGamma-1.0)*rho;
}

/*
 * The (rho, T) derivatives of P = (gamma-1) rho cv T: dP/drho|T = (gamma-1) cv T and
 * dP/dT|rho = (gamma-1) rho cv. Same domain as igeosPofRhoT: IGEOS_FAIL_VALUE (counted) for
 * T < 0, NaN or an invalid rho, T = 0 (the cold curve) is valid. Up to 1.1.0 EOSdPdRhoofRhoT and
 * EOSdPdT were fatal (not implemented) for the ideal gas.
 */
double igeosdPdRhoofRhoT(IGEOSMAT *Mat, double rho, double T) {
    if (igeosCheckRhoT(Mat, __func__, rho, T)) return IGEOS_FAIL_VALUE;
    return (Mat->dConstGamma-1.0)*Mat->cv*T;
}

double igeosdPdTofRhoT(IGEOSMAT *Mat, double rho, double T) {
    if (igeosCheckRhoT(Mat, __func__, rho, T)) return IGEOS_FAIL_VALUE;
    return (Mat->dConstGamma-1.0)*rho*Mat->cv;
}

/*
 * The cold curve u(rho, T = 0) = 0 for a valid rho; IGEOS_FAIL_VALUE (counted) for an invalid one,
 * like u_cold outside a table (ANEOSUCold, reos3UCold, scvheosUCold). EOSUCold returned 0 for any
 * rho up to 1.1.0, so EOSisbelowColdCurve answered for rho <= 0 or NaN instead of failing.
 */
double igeosUCold(IGEOSMAT *Mat, double rho) {
    if (!igeosValidRho(rho)) {
        Mat->nLookupFail++;
        IGEOS_WARN("igeosUCold: invalid density rho= %g (iMat= %i).\n", rho, Mat->iMat);
        return IGEOS_FAIL_VALUE;
    }
    return 0.0;
}

/*
 * Calculate the specific entropy s(rho, T) of an ideal gas with constant cv:
 *
 *   s = cv*ln(T) - (gamma-1)*cv*ln(rho) + const
 *
 * where (gamma-1)*cv = k_B/(mu*m_H). The result is returned in cgs
 * (erg/g/K), consistent with the tabulated EOS libraries, and the additive
 * constant is set to zero. Only entropy differences are meaningful.
 * Returns IGEOS_FAIL_VALUE (counted in nLookupFail) for an invalid rho or T <= 0.
 */
double igeosSofRhoT(IGEOSMAT *Mat, double rho, double T) {
    double cv_cgs = Mat->cv * Mat->dErgPerGmUnit;
    double rho_cgs = rho * Mat->dGmPerCcUnit;

    if (!igeosValidRho(rho) || !(T > 0.0)) {
        Mat->nLookupFail++;
        IGEOS_WARN("igeosSofRhoT: entropy undefined for rho= %g T= %g (iMat= %i).\n", rho, T, Mat->iMat);
        return IGEOS_FAIL_VALUE;
    }

    return cv_cgs * (log(T) - (Mat->dConstGamma - 1.0) * log(rho_cgs));
}

/*
 * Calculate the specific entropy s(rho, u). IGEOS_FAIL_VALUE for u <= 0 (s(rho, T = 0) is
 * undefined) or an invalid rho, counted once (not as a T fallback as well).
 */
double igeosSofRhoU(IGEOSMAT *Mat, double rho, double u) {
    return igeosSofRhoT(Mat, rho, u/Mat->cv);
}