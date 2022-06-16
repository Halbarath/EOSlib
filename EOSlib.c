/*  This file is part of EOSlib.
 *  Copyright (c) 2020-2021 Thomas Meier & Christian Reinhardt
 *
 *  EOSlib is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  EOSlib is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with EOSlib.  If not, see <http://www.gnu.org/licenses/>.
 */

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
#include <string.h>
#include <string.h>
#include "EOSlib.h"

/*
 * Initialization of the material structures
 * iMat: Material number
 * dKpcUnit, dMsolUnit: Unit conversion factors
 * additional_data: void pointer that can be used to give additional options, not used at the moment
 */
EOSMATERIAL *EOSinitMaterial(int iMat, double dKpcUnit, double dMsolUnit, const void * additional_data)
{
    EOSMATERIAL *material;
    material = (EOSMATERIAL *) calloc(1, sizeof(EOSMATERIAL));
    material->iMat = iMat;
    material->bEntropyTableInit = EOS_FALSE;

    if (iMat == MAT_IDEALGAS)
    {
        /* Check if the Tillotson library has the right version. */
        if (IGEOS_VERSION_MAJOR != 1) {
            fprintf(stderr, "EOSinitMaterial: Ideal gas EOS library has the wrong version (%s)\n", IGEOS_VERSION_TEXT);
            exit(1);
        }

        assert(additional_data != NULL);

        /* Pass additional parameter to the ideal gas EOS library. */
        struct igeosParam *Param;
        Param = (struct igeosParam *) additional_data;

        material->matType = EOSIDEALGAS;
        material->igeosmaterial = igeosInitMat(iMat, Param->dConstGamma, Param->dMeanMolMass, dKpcUnit, dMsolUnit);
        material->rho0 = material->igeosmaterial->rho0;
        material->minSoundSpeed = 0;
        /* No entropy look up table required. */
        material->bEntropyTableInit = EOS_TRUE;
        strcpy(material->MatString, "IDEAL GAS");
    } else if (iMat>=MAT_TILLOTSON_MIN && iMat<=MAT_TILLOTSON_MAX)
    {
        /* Check if the Tillotson library has the right version. */
        if (TILL_VERSION_MAJOR != 3 || TILL_VERSION_MINOR < 4) {
            fprintf(stderr, "EOSinitMaterial: Tillotson library has the wrong version (%s)\n", TILL_VERSION_TEXT);
            exit(1);
        }
        material->matType = EOSTILLOTSON;
        material->tillmaterial = tillInitMaterial(iMat, dKpcUnit, dMsolUnit);
        material->rho0 = material->tillmaterial->rho0;
        material->minSoundSpeed = sqrt(material->tillmaterial->A/material->tillmaterial->rho0);
        tilliMatString(material->tillmaterial, material->MatString);
    } else if (iMat>=MAT_ANEOS_MIN && iMat <=MAT_ANEOS_MAX)
    {
        /* Check if the ANEOS library has the right version. */
        if (ANEOS_VERSION_MAJOR != 1) {
            fprintf(stderr, "EOSinitMaterial: ANEOS library has the wrong version (%s)\n", ANEOS_VERSION_TEXT);
            exit(1);
        }
        material->matType = EOSANEOS;
        material->ANEOSmaterial = ANEOSinitMaterial(iMat, dKpcUnit, dMsolUnit);
        material->rho0 = ANEOSgetRho0(material->ANEOSmaterial);
        // The entropy look up table is initialized in ANEOSinitMaterial
        material->bEntropyTableInit = EOS_TRUE;
        material->minSoundSpeed = ANEOSCofRhoT(material->ANEOSmaterial, material->rho0, material->ANEOSmaterial->TAxis[0]);
        strcpy(material->MatString, material->ANEOSmaterial->matName);
    } else {
        fprintf(stderr, "EOSinitMaterial: iMat %i does not exist.\n",iMat);
    }

    return material;
}

/*
 * Initialize the lookup tables needed for the isentropoic evolution
 */
void EOSinitIsentropicLookup(EOSMATERIAL *material, const void * additional_data)
{
    switch(material->matType)
    {
        case EOSIDEALGAS:
            if (material->bEntropyTableInit != EOS_TRUE) material->bEntropyTableInit = EOS_TRUE;
            break;
        case EOSTILLOTSON:
            tillInitLookup(material->tillmaterial, 1000, 1000, 1e-4, 200.0, 1200.0);
            material->bEntropyTableInit = EOS_TRUE;
            break;
        case EOSANEOS:
            // nothing to do
            if (material->bEntropyTableInit != EOS_TRUE) material->bEntropyTableInit = EOS_TRUE;
            break;
        default:
            assert(0);
    }
}

/*
 * Finalization of the material structures
 */
void EOSfinalizeMaterial(EOSMATERIAL *material)
{
    switch(material->matType)
    {
        case EOSIDEALGAS:
            igeosFinalizeMat(material->igeosmaterial);
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
 * Print material data, e.g., EOS coefficients or the size of the EOS table.
 */
void EOSPrintMat(EOSMATERIAL *material, FILE *fp)
{
    switch(material->matType)
    {
        case EOSIDEALGAS:
            igeosPrintMat(material->igeosmaterial, fp);
            break;
        case EOSTILLOTSON:
            tillPrintMat(material->tillmaterial, fp);
            break;
        case EOSANEOS:
            ANEOSPrintMat(material->ANEOSmaterial, fp);
            break;
        default:
            assert(0);
    }
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
            P = igeosPofRhoU(material->igeosmaterial, rho, u);
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
 * Calculates the pressure P(rho,T) for a material
 */
double EOSPofRhoT(EOSMATERIAL *material, double rho, double T)
{
    double P = 0;

    switch(material->matType)
    {
        case EOSIDEALGAS:
            // not implemented
            assert(0);
            break;
        case EOSTILLOTSON:
            // not implemented
            assert(0);
            break;
        case EOSANEOS:
            P = ANEOSPofRhoT(material->ANEOSmaterial, rho, T);
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

    switch(material->matType)
    {
        case EOSIDEALGAS:
            c = igeosCofRhoU(material->igeosmaterial, rho, u);
            break;
        case EOSTILLOTSON:
            tillPressureSound(material->tillmaterial, rho, u, &c);
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
            P = igeosPCofRhoU(material->igeosmaterial, rho, u, c);
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
    double u2 = 0;
    switch(material->matType)
    {
        case EOSIDEALGAS:
            u2 = igeosIsentropicU(material->igeosmaterial, rho1, u1, rho2);
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
    double T = 0;
    switch(material->matType)
    {
        case EOSIDEALGAS:
            T = igeosTofRhoU(material->igeosmaterial, rho, u);
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
    double u = 0;
    switch(material->matType)
    {
        case EOSIDEALGAS:
            u = igeosUofRhoT(material->igeosmaterial, rho, T);
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
    double rho = 0;
    switch(material->matType)
    {
        case EOSIDEALGAS:
            rho = igeosRhoofPT(material->igeosmaterial, p, T);
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
            assert(0);
            break;
        case EOSTILLOTSON:
            // not implemented
            assert(0);
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
 * Tests if the given combination is below the cold curve, returns EOS_TRUE if so and EOS_FALSE if not
 */
int EOSisbelowColdCurve(EOSMATERIAL *material, double rho, double u)
{
    double ucold = EOSUCold(material, rho);
    return (u < ucold);
}

/*
 * Tests if the given combination is inside the lookup tables
 * returns EOS_SUCCESS if so, EOS_FAIL if not
 */
int EOSIsInTable(EOSMATERIAL *material, double rho, double u)
{
    int iret = EOS_FAIL;
    switch(material->matType)
    {
        case EOSIDEALGAS:
            iret = EOS_SUCCESS;
            break;
        case EOSTILLOTSON:
            iret = tillIsInTable(material->tillmaterial, rho, u);
            if (iret == TILL_LOOKUP_SUCCESS) return EOS_SUCCESS;
            if (iret == TILL_LOOKUP_OUTSIDE_RHOMIN) return EOS_OUTSIDE_RHOMIN;
            if (iret == TILL_LOOKUP_OUTSIDE_RHOMAX) return EOS_OUTSIDE_RHOMAX;
            if (iret == TILL_LOOKUP_OUTSIDE_VMIN) return EOS_OUTSIDE_VMIN;
            /* Allow particles to have v > v_max. */
            if (iret == TILL_LOOKUP_OUTSIDE_VMAX) return EOS_SUCCESS;
            break;
        case EOSANEOS:
            if (rho < material->ANEOSmaterial->rhoAxis[0]/material->ANEOSmaterial->CodeUnitstoCGSforRho)
            {
                return EOS_OUTSIDE_RHOMIN;
            }
            if (rho >= material->ANEOSmaterial->rhoAxis[material->ANEOSmaterial->nRho-1]/material->ANEOSmaterial->CodeUnitstoCGSforRho)
            {
                return EOS_OUTSIDE_RHOMAX;
            }
            if (u < ANEOSUofRhoT(material->ANEOSmaterial, rho, material->ANEOSmaterial->TAxis[0]))
            {
                return EOS_OUTSIDE_VMIN;
            }
            if (u > ANEOSUofRhoT(material->ANEOSmaterial, rho, material->ANEOSmaterial->TAxis[material->ANEOSmaterial->nT-1]*0.9999))
            {
                return EOS_OUTSIDE_VMAX;
            }
            return EOS_SUCCESS;
            break;
        default:
            assert(0);
    }

    return iret;
}


/*
 * Calculates the derivative dPdrho(rho,u) for a material
 */
double EOSdPdRho(EOSMATERIAL *material, double rho, double u)
{
    double dPdRho = 0;
    switch(material->matType)
    {
        case EOSIDEALGAS:
            dPdRho = igeosdPdRho(material->igeosmaterial, rho, u);
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
    double dPdU = 0;
    switch(material->matType)
    {
        case EOSIDEALGAS:
            dPdU = igeosdPdU(material->igeosmaterial, rho, u);
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
    double dUdRho = 0;
    switch(material->matType)
    {
        case EOSIDEALGAS:
            dUdRho = igeosPofRhoU(material->igeosmaterial, rho, u) / (rho * rho);
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
 * Calculates the derivative dPdrho(rho,T) for a material
 */
double EOSdPdRhoatT(EOSMATERIAL *material, double rho, double T)
{
    double dPdRho = 0;
    switch(material->matType)
    {
        case EOSIDEALGAS:
            // not implemented
            assert(0);
            break;
        case EOSTILLOTSON:
            // not implemented
            assert(0);
            break;
        case EOSANEOS:
            dPdRho = ANEOSdPdRhoofRhoT(material->ANEOSmaterial, rho, T);
            break;
        default:
            assert(0);
    }

    return dPdRho;
}

/*
 * Calculates the derivative dPdT(rho,T) for a material
 */
double EOSdPdT(EOSMATERIAL *material, double rho, double T)
{
    double dPdT = 0;
    switch(material->matType)
    {
        case EOSIDEALGAS:
            // not implemented
            assert(0);
            break;
        case EOSTILLOTSON:
            // not implemented
            assert(0);
            break;
        case EOSANEOS:
            dPdT = ANEOSdPdTofRhoT(material->ANEOSmaterial, rho, T);
            break;
        default:
            assert(0);
    }

    return dPdT;
}

/*
 * Calculate the cold part of the internal energy, i.e., u(rho, T=0).
 */
double EOSUCold(EOSMATERIAL *material, double rho)
{
    double ucold = 0;
    switch(material->matType)
    {
        case EOSIDEALGAS:
            // ucold = 0
            break;
        case EOSTILLOTSON:
            ucold = tillColdULookup(material->tillmaterial, rho);
            break;
        case EOSANEOS:
            // For ANEOS temperatures below T_min cause problems
            ucold = ANEOSUofRhoT(material->ANEOSmaterial, rho, material->ANEOSmaterial->TAxis[0]);
            break;
        default:
            assert(0);
    }

    return ucold;
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
    double a, ua, Pa, b, ub, Pb, c = 0, uc = 0, Pc;
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

    b = 0.999*material2->rho0; // hard coded minimum
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
    while (2*(a-b)/(a+b) > 1e-10) {
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
#ifdef EOSLIB_VERBOSE
        fprintf(stderr, "Warning: Pressure is close to zero so no density correction is done.\n");
#endif
        return 1.0;
    }

    rho1 = EOSRhoofPT(material1, P, T);
    rho2 = EOSRhoofPT(material2, P, T);

    /* 
     * Limit correction in the expanded states because e.g. at phase transitions P(rho) is
     * monotonic and f_ij can become huge.
     */
    if ((rho1 <= 0.8*material1->rho0) || (rho2 <= 0.8*material2->rho0)) {
        return 1.0;
    }

    // If the density is unphysical return 1 so the density is not corrected.
    if ((rho1 <= 0.0) || (rho2 <= 0.0))
    {
        fprintf(stderr, "CalcWoolfsonCoeff: unphysical density (rho1= %g, rho2= %g).\n", rho1, rho2);
        //        return 1.0;
        return -1.0;
    }

#if 0
    // CR 10.12.2020: Debugging why some f_ij > 1e3
    if (rho1/rho2 > 1e3) {
        fprintf(stderr, "\n");
        fprintf(stderr, "EOSWoolfsonCoeff: ");
        fprintf(stderr, "iMat1= %i ", material1->iMat);
        fprintf(stderr, "iMat2= %i ", material2->iMat);
        fprintf(stderr, "P= %g ", P);
        fprintf(stderr, "T= %g ", T);
        fprintf(stderr, "rho1= %g ", rho1);
        fprintf(stderr, "rho2= %g ", rho2);
        fprintf(stderr, "\n");
    }
#endif
    return (rho1/rho2);
}
