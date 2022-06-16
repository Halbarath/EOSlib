/*
 * Test the implementation of the ideal gas equation of state (EOS) for EOSlib.
 *
 * Author: Christian Reinhardt
 * Date: 16.06.2022
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "EOSlib.h"
#include "igeos.h"

int main(int argc, char *argv[]) {
	EOSMATERIAL *Mat1;
    IGEOSMAT *Mat2;
    /* Ideal gas with gamma = 5/3 and mu = 1.0. */
    double dConstGamma = 5.0/3.0;
    double dMeanMolMass = 1.0;
    /* L_unit = 1 RE, v_unit = 1 km/s. */
    double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
    double *rhoAxis;
    double *uAxis;
    double rho_min = 1e-3;
    double rho_max = 1e2;
    double u_min = 1e-3;
    double u_max = 1e2;
    int nRho = 1000;
    int nU = 1000;
    double rho;
    double u;
    double eps = 1e-15;
    struct igeosParam param;

    param.dConstGamma = dConstGamma;
    param.dMeanMolMass = dMeanMolMass;

    /* Initialize the ideal gas implemented in the Tillotson EOS library. */
	Mat1 = EOSinitMaterial(0, dKpcUnit, dMsolUnit, &param);

    Mat2 = igeosInitMat(0, dConstGamma, dMeanMolMass, dKpcUnit, dMsolUnit);

    rhoAxis = (double *) calloc(nRho, sizeof(double));
    uAxis = (double *) calloc(nU, sizeof(double));

    for (int i=0; i<nRho; i++) {
        rhoAxis[i] = rho_min + i*(rho_max-rho_min)/(nRho-1);
    }

    for (int j=0; j<nU; j++) {
        uAxis[j] = u_min + j*(u_max-u_min)/(nU-1);
    }
 
    fprintf(stderr, "Compare EOS...\n");

    for (int i=0; i<nRho; i++) {
        for (int j=0; j<nU; j++) {
            rho = rhoAxis[i];
            u = uAxis[j];

            //fprintf(stderr, "i= %i j= %i rho= %g u= %g\n", i, j, rho, u);

            /* P(rho, u). */
            assert((igeosPofRhoU(Mat2, rho, u)-EOSPofRhoU(Mat1, rho, u))/igeosPofRhoU(Mat2, rho, u) <= eps);

            /* cs(rho, u). */
            assert((igeosCofRhoU(Mat2, rho, u)-EOSCofRhoU(Mat1, rho, u))/igeosCofRhoU(Mat2, rho, u) <= eps);
            /* T(rho, u). */
            assert((igeosTofRhoU(Mat2, rho, u)-EOSTofRhoU(Mat1, rho, u))/igeosTofRhoU(Mat2, rho, u) <= eps);

            /* Calculate dP/drho(rho, u). */
            assert((igeosdPdRho(Mat2, rho, u)-EOSdPdRho(Mat1, rho, u))/igeosdPdRho(Mat2, rho, u) <= eps);

            /* Calculate dP/du(rho, u). */
            assert((igeosdPdU(Mat2, rho, u)-EOSdPdU(Mat1, rho, u))/igeosdPdU(Mat2, rho, u) <= eps);
        }
    }

    fprintf(stderr, "Done.\n");

    /* Test isentropic evolution. */
    double u1 = EOSIsentropic(Mat1, rhoAxis[0], uAxis[0], rhoAxis[nRho-1]);
    double u2 = igeosIsentropicU(Mat2, rhoAxis[0], uAxis[0], rhoAxis[nRho-1]);
    
    assert((u2-u1)/u2 < eps);


    EOSfinalizeMaterial(Mat1);
    igeosFinalizeMat(Mat2);
}


