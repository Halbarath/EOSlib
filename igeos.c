/*
 * Implementation of the ideal gas equation of state (EOS) for EOSlib.
 *
 * Author: Christian Reinhardt
 * Date: 15.06.2022
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
}

/* Calculate P(rho, u). */
double igeosPofRhoU(IGEOSMAT *Mat, double rho, double u) {
    return (Mat->dConstGamma-1.0)*rho*u;
}

/* Calculate P(rho, T). */
double igeosPofRhoT(IGEOSMAT *Mat, double rho, double T) {
    return (Mat->dConstGamma-1.0)*rho*Mat->cv*T;
}

/* Calculate cs(rho, u). */
double igeosCofRhoU(IGEOSMAT *Mat, double rho, double u) {
    return sqrt(Mat->dConstGamma*(Mat->dConstGamma-1.0)*u);
}

/* Calculate P(rho, u) and cs(rho, u). */
double igeosPCofRhoU(IGEOSMAT *Mat, double rho, double u, double *c) {
   *c = igeosCofRhoU(Mat, rho, u);
   return igeosPofRhoU(Mat, rho, u);
}

/* Calculate u2(rho1, u1, rho2) so that (rho1, u1) and (rho2, u2) are on the same isentrope. */
double igeosIsentropicU(IGEOSMAT *Mat, double rho1, double u1, double rho2) {
    return u1*pow(rho2/rho1, Mat->dConstGamma-1.0);
}

/* Calculate T(rho, u). */
double igeosTofRhoU(IGEOSMAT *Mat, double rho, double u) {
    return u/Mat->cv;
}

/* Calculate u(rho, T). */
double igeosUofRhoT(IGEOSMAT *Mat, double rho, double T) {
    return Mat->cv*T;
}

/* Calculate rho(P, T). */
double igeosRhoofPT(IGEOSMAT *Mat, double P, double T) {
    return P / ((Mat->dConstGamma - 1.0) * Mat->cv * T);
}

/* Calculate rho(U, T). */
double igeosRhoofUT(IGEOSMAT *Mat, double U, double T) {
    // This is undefined, as cv is constant (both u(rho,T) and T(rho,u) do not use rho)
    assert(0);
    return 0;
}

/* Calculate dP/drho(rho, u). */
double igeosdPdRho(IGEOSMAT *Mat, double rho, double u) {
    return (Mat->dConstGamma-1.0)*u;
}

/* Calculate dP/du(rho, u). */
double igeosdPdU(IGEOSMAT *Mat, double rho, double u) {
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
 */
double igeosSofRhoT(IGEOSMAT *Mat, double rho, double T) {
    double cv_cgs = Mat->cv * Mat->dErgPerGmUnit;
    double rho_cgs = rho * Mat->dGmPerCcUnit;

    if (!(rho_cgs > 0.0) || !(T > 0.0)) return -1e50;

    return cv_cgs * (log(T) - (Mat->dConstGamma - 1.0) * log(rho_cgs));
}

/* Calculate the specific entropy s(rho, u). */
double igeosSofRhoU(IGEOSMAT *Mat, double rho, double u) {
    return igeosSofRhoT(Mat, rho, igeosTofRhoU(Mat, rho, u));
}
