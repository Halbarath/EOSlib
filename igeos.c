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
    const double KBOLTZ = 1.38e-16;      /* bolzman constant in cgs */
    const double MHYDR = 1.67e-24;       /* mass of hydrogen atom in grams */
    const double MSOLG = 1.99e33;        /* solar mass in grams */
    const double GCGS = 6.67e-8;         /* G in cgs */
    const double KPCCM = 3.085678e21;    /* kiloparsec in centimeters */
    const double NA = 6.022e23;          /* Avogadro's number */
    IGEOSMAT *Mat;

    Mat = (IGEOSMAT *) calloc(1, sizeof(IGEOSMAT));
    assert(Mat != NULL);

    Mat->iMat = iMat;
    assert(Mat->iMat == 0);

	/* This two parameters define the unit system. */
    Mat->dKpcUnit = dKpcUnit;
    Mat->dMsolUnit = dMsolUnit;

	if (dKpcUnit <= 0.0 && dMsolUnit <= 0.0)
	{
		/* In this case units are not converted, so the code units are cgs. */
		Mat->dGasConst = KBOLTZ;   // dGasConst is NOT KBOLZ !!!!!!!
		Mat->dErgPerGmUnit = 1.0;
		Mat->dGmPerCcUnit = 1.0;
		Mat->dSecUnit = 1.0;
	} else {
		/* Convert kboltz/mhydrogen to system units, assuming that G == 1. */
		Mat->dGasConst = Mat->dKpcUnit*KPCCM*KBOLTZ
		/MHYDR/GCGS/Mat->dMsolUnit/MSOLG;
		/* code energy per unit mass --> erg per g */
		Mat->dErgPerGmUnit = GCGS*Mat->dMsolUnit*MSOLG/(Mat->dKpcUnit*KPCCM);
		/* code density --> g per cc */
		Mat->dGmPerCcUnit = (Mat->dMsolUnit*MSOLG)/pow(Mat->dKpcUnit*KPCCM,3.0);
		/* code time --> seconds */
		Mat->dSecUnit = sqrt(1/(Mat->dGmPerCcUnit*GCGS));
	}
        
    Mat->dConstGamma = dConstGamma;
    Mat->dMeanMolMass = dMeanMolMass;

    /* cv = kb/mp */
    Mat->cv = KBOLTZ/((Mat->dConstGamma-1.0)*MHYDR*Mat->dMeanMolMass);

    /* This value is hard coded and used as a default value in ballic. */
    Mat->rho0 = 0.001;

    /* Convert to code units. */
    Mat->rho0 /= Mat->dGmPerCcUnit;
    Mat->cv /= Mat->dErgPerGmUnit;

    return Mat;
}

/* Free memory. */
int igeosFinalizeMat(IGEOSMAT *Mat) {
    free(Mat);
    return 0;
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
    assert(0);
    return 0;
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
    assert(0);
    return 0;
}

/* Calculate rho(U, T). */
double igeosRhoofUT(IGEOSMAT *Mat, double P, double T) {
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


