/*
 * Equation Of State library EOSlib
 *
 * This is a wrapper to be used in Gasoline and ballic
 * It wraps around the Tillotson material library and the ANEOSmaterial library
 * but can be extended with additional material libraries
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "EOSlib.h"

/*
 * Initialization of the material structures
 * iMat: Material number
 * dKpcUnit, dMsolUnit: Unit conversion factors
 * additional_data: void pointer that can be used to give additional options:
 *   ideal gas: not used
 *   Tillotson: if not NULL pointer, lookup table is initialized
 *   ANEOS: not used
 */
EOSMATERIAL *EOSinitMaterial(int iMat, double dKpcUnit, double dMsolUnit, const void * additional_data)
{
	EOSMATERIAL *material;
	material = (EOSMATERIAL *) calloc(1, sizeof(EOSMATERIAL));
	material->iMat = iMat;
	material->canDoIsentropic = 0;

	if (iMat == iMATIDEALGAS)
	{
		// not implemented
	} else if (iMat>=iMATTILLOTSONMIN && iMat<=iMATTILLOTSONMAX)
	{
		/* Check if the Tillotson library has the right version. */
		if (TILL_VERSION_MAJOR != 3 || TILL_VERSION_MINOR < 4) {
			fprintf(stderr, "EOSinitMaterial: Tillotson library has the wrong version (%s)\n", TILL_VERSION_TEXT);
			exit(1);
		}
		material->matType = EOSTILLOTSON;
		material->tillmaterial = tillInitMaterial(iMat, dKpcUnit, dMsolUnit);
		if (additional_data != NULL)
		{
			tillInitLookup(material->tillmaterial, 1000, 1000, 1e-4, 200.0, 1200.0);
			material->canDoIsentropic = 1;
		}
		material->rho0 = material->tillmaterial->rho0;
	} else if (iMat>=iMATANEOSMIN && iMat <=iMATANEOSMAX)
	{
		/* Check if the ANEOS library has the right version. */
		if (ANEOS_VERSION_MAJOR != 1) {
			fprintf(stderr, "EOSinitMaterial: ANEOS library has the wrong version (%s)\n", ANEOS_VERSION_TEXT);
			exit(1);
		}
		material->matType = EOSANEOS;
		material->ANEOSmaterial = ANEOSinitMaterial(iMat, dKpcUnit, dMsolUnit);
		material->rho0 = ANEOSgetRho0(material->ANEOSmaterial);
		material->canDoIsentropic = 1;
	}
	
	return material;
}

/*
 * Finalization of the material structures
 */
void EOSfinalizeMaterial(EOSMATERIAL *material)
{
	switch(material->matType)
	{
		case EOSIDEALGAS:
		// not implemented
			break;
		case EOSTILLOTSON:
			tillFinalizeMaterial(material->tillmaterial);
			break;
		case EOSANEOS:
			ANEOSfinalizeMaterial(material->ANEOSmaterial);
			break;
		default:
		assert(0);
	}

	free(material);
}

/*
 * Calculates the pressure P(rho,u) for a material
 */
double EOSPofRhoU(EOSMATERIAL *material, double rho, double u)
{
	double P = 0;
	
	switch(material->matType)
	{
		case EOSIDEALGAS:
		// not implemented
			break;
		case EOSTILLOTSON:
			P = tillPressure(material->tillmaterial, rho, u);
			break;
		case EOSANEOS:
			P = ANEOSPofRhoU(material->ANEOSmaterial, rho, u);
			break;
		default:
		assert(0);
	}
	
	return P;
}

/*
 * Calculates the sound speed c(rho,u) for a material
 */
double EOSCofRhoU(EOSMATERIAL *material, double rho, double u)
{
	double c = 0;
	double P = 0;
	
	switch(material->matType)
	{
		case EOSIDEALGAS:
		// not implemented
			break;
		case EOSTILLOTSON:
			P = tillPressureSound(material->tillmaterial, rho, u, &c);
			break;
		case EOSANEOS:
			c = ANEOSCofRhoU(material->ANEOSmaterial, rho, u);
			break;
		default:
		assert(0);
	}
	
	return c;
}

/*
 * Calculates the pressure P(rho,u) and the sound speed c(rho,u) for a material
 * This should be used wherever both the pressure and the sound speed are needed for the same rho and u
 * as it is faster for some material models (esp. ANEOS)
 */
double EOSPCofRhoU(EOSMATERIAL *material, double rho, double u, double *c)
{
	double P = 0;
	double T = 0;
	
	switch(material->matType)
	{
		case EOSIDEALGAS:
		// not implemented
			break;
		case EOSTILLOTSON:
			P = tillPressureSound(material->tillmaterial, rho, u, c);
			break;
		case EOSANEOS:
			T = ANEOSTofRhoU(material->ANEOSmaterial, rho, u);
			P = ANEOSPofRhoT(material->ANEOSmaterial, rho, T);
			*c = ANEOSCofRhoT(material->ANEOSmaterial, rho, T);
			break;
		default:
		assert(0);
	}
	
	return P;
}

/*
 * Calculates the internal energy u2(rho1, u1, rho2) after an isentropic evolution for a material
 */
double EOSIsentropic(EOSMATERIAL *material, double rho1, double u1, double rho2)
{
	double u2;
	switch(material->matType)
	{
		case EOSIDEALGAS:
		// not implemented
			break;
		case EOSTILLOTSON:
			u2 = tillLookupU(material->tillmaterial, rho1, u1, rho2, 0); // last argument is actually a particle number
			break;
		case EOSANEOS:
			u2 = ANEOSisentropicU(material->ANEOSmaterial, rho1, u1, rho2);
			break;
		default:
		assert(0);
	}

	return u2;
}

/*
 * Calculates the temperature T(rho,u) for a material
 */
double EOSTofRhoU(EOSMATERIAL *material, double rho, double u)
{
	double T;
	switch(material->matType)
	{
		case EOSIDEALGAS:
		// not implemented
			break;
		case EOSTILLOTSON:
			T = tillTempRhoU(material->tillmaterial, rho, u);
			break;
		case EOSANEOS:
			T = ANEOSTofRhoU(material->ANEOSmaterial, rho, u);
			break;
		default:
		assert(0);
	}

	return T;
}

/*
 * Calculates the internal energy u(rho,T) for a material
 */
double EOSUofRhoT(EOSMATERIAL *material, double rho, double T)
{
	double u;
	switch(material->matType)
	{
		case EOSIDEALGAS:
		// not implemented
			break;
		case EOSTILLOTSON:
			u = tillURhoTemp(material->tillmaterial, rho, T);
			break;
		case EOSANEOS:
			u = ANEOSUofRhoT(material->ANEOSmaterial, rho, T);
			break;
		default:
		assert(0);
	}

	return u;
}

/*
 * Calculates the density rho(p,T) for a material
 */
double EOSRhoofPT(EOSMATERIAL *material, double p, double T)
{
	double rho;
	switch(material->matType)
	{
		case EOSIDEALGAS:
		// not implemented
			break;
		case EOSTILLOTSON:
			rho = tillRhoPTemp(material->tillmaterial, p, T);
			break;
		case EOSANEOS:
			rho = ANEOSRhoofPT(material->ANEOSmaterial, p, T);
			break;
		default:
		assert(0);
	}

	return rho;
}

/*
 * Calculates the density rho(u,T) for a material
 */
double EOSRhoofUT(EOSMATERIAL *material, double u, double T)
{
	double rho = 0;
	switch(material->matType)
	{
		case EOSIDEALGAS:
		// not implemented
			break;
		case EOSTILLOTSON:
			// not implemented
			break;
		case EOSANEOS:
			rho = ANEOSRhoofUT(material->ANEOSmaterial, u, T);
			break;
		default:
		assert(0);
	}

	return rho;
}

/*
 * Tests if the given combination is below the cold curve, returns 1 if so and 0 if not
 */
int EOSisbelowColdCurve(EOSMATERIAL *material, double rho, double u)
{
	double ucold = EOSUofRhoT(material, rho, 1e-4);
	return (u < ucold);
}

/*
 * Calculates the derivative dPdrho(rho,u) for a material
 */
double EOSdPdRho(EOSMATERIAL *material, double rho, double u)
{
	double dPdRho;
	switch(material->matType)
	{
		case EOSIDEALGAS:
		// not implemented
			break;
		case EOSTILLOTSON:
			dPdRho = tilldPdrho(material->tillmaterial, rho, u);
			break;
		case EOSANEOS:
			dPdRho = ANEOSdPdRhoofRhoU(material->ANEOSmaterial, rho, u);
			break;
		default:
		assert(0);
	}

	return dPdRho;
}

/*
 * Calculates the derivative dPdu(rho,u) for a material
 */
double EOSdPdU(EOSMATERIAL *material, double rho, double u)
{
	double dPdU;
	switch(material->matType)
	{
		case EOSIDEALGAS:
		// not implemented
			break;
		case EOSTILLOTSON:
			dPdU = tilldPdu(material->tillmaterial, rho, u);
			break;
		case EOSANEOS:
			dPdU = ANEOSdPdUofRhoU(material->ANEOSmaterial, rho, u);
			break;
		default:
		assert(0);
	}

	return dPdU;
}

/*
 * Calculates the derivative dudrho(rho,u) for a material
 */
double EOSdUdRho(EOSMATERIAL *material, double rho, double u)
{
	double dUdRho;
	switch(material->matType)
	{
		case EOSIDEALGAS:
		// not implemented
			break;
		case EOSTILLOTSON:
			dUdRho = tilldudrho(material->tillmaterial, rho, u);
			break;
		case EOSANEOS:
			dUdRho = ANEOSdUdRhoofRhoU(material->ANEOSmaterial, rho, u);
			break;
		default:
		assert(0);
	}

	return dUdRho;
}

/*
 * Given rho1, u1 (material 1) solve for rho2, u2 (material 2) at the interface between two material.
 *
 * The b.c. are:
 *
 * P1(rho1,u1)=P2(rho2,u2) and T1(rho1,u1)=T2(rho2,u2)
 *
 * Since P = P(rho, T) and T1=T2 we solve for P2(rho2)-P1=0.
 *
 * Returns 0 if successful or -1 if not.
 */
int EOSSolveBC(EOSMATERIAL *material1, EOSMATERIAL *material2, double rho1, double u1,
 double *prho2, double *pu2)
{
	double P, T;
    double a, ua, Pa, b, ub, Pb, c, uc, Pc;
    int iRet;

    iRet = -1;
    Pc = 0.0;

    /* Calculate P and T in material 1. */
	P = EOSPofRhoU(material1, rho1, u1);
	T = EOSTofRhoU(material1, rho1, u1);

    /*
     * We use rho1 as an upper limit for rho2 assuming that the denser component is in the inner shell.
     */
    a = 90; // hard coded maximum
	ua = EOSUofRhoT(material2, a, T);
	Pa = EOSPofRhoU(material2, a, ua);

	b = 1e-8; // hard coded minimum
	ub = EOSUofRhoT(material2, b, T);
	Pb = EOSPofRhoU(material2, b, ub);

    /*
     * Assert that the root is bracketed by (a, b).
     */
    if (Pa < P || Pb > P)
    {
        return iRet;
    }

    /*
     * Root bracketed by (a,b).
     */
    while (2*(Pa-Pb)/(Pa+Pb) > 1e-10) {
        c = 0.5*(a + b);
		uc = EOSUofRhoT(material2, c, T);
		Pc = EOSPofRhoU(material2, c, uc);

        if (Pc < P) {
            b = c;
            Pb = Pc;
        }
        else {
            a = c;
            Pa = Pc;
        }
    }

    /*
     * Return values.
     */
    *prho2 = c;
    *pu2 = uc;

    iRet = 0;

    return iRet;
}

/*
 * Calculate the coefficient
 *
 * f_ij := rho_mat1(P, T)/rho_mat2(P, T)
 *
 * needed to correct the density at a material interface.
 */
double EOSWoolfsonCoeff(EOSMATERIAL *material1, EOSMATERIAL *material2, double P, double T)
{
    double rho1;
    double rho2;

    // Check if there is indeed a material interface.
    if (material1->iMat == material2->iMat)
	{
        return 1.0;
	}
    /*
     * In the low density region a density correction can be problematic (dPdrho not monotonic,
     * interpretation of mixed phases in the Tillotson EOS difficult), so the correction factor
     * is one in this case.
     */

    if (P < WOOLFSON_MIN_PRESSURE) {
        fprintf(stderr, "Warning: Pressure is close to zero so no density correction is done.\n");
        return 1.0;
    }

	rho1 = EOSRhoofPT(material1, P, T);
	rho2 = EOSRhoofPT(material2, P, T);


    // If the density is unphysical return 1 so the density is not corrected.
    if ((rho1 <= 0.0) || (rho2 <= 0.0))
    {
        fprintf(stderr, "CalcWoolfsonCoeff: unphysical density (rho1= %g, rho2= %g).\n", rho1, rho2);
        //        return 1.0;
        return -1.0;
    }

    return (rho1/rho2);
}
