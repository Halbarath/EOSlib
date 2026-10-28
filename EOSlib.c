/*  This file is part of EOSlib.
 *  Copyright (c) 2020-2026 Thomas Meier & Christian Reinhardt
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
 * This is a wrapper to be used in pkdgrav3, Gasoline and ballic.
 * It wraps around the Tillotson material library and the ANEOSmaterial library
 * but can be extended with additional material libraries.
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include "EOSlib.h"

/*
 * Fatal errors (see "Failure contract" in EOSlib.h): a message, assert(0)
 * (SIGABRT, so an MPI launcher kills the whole job) and exit(1) for NDEBUG
 * builds, for programming errors (unknown
 * material, function not implemented for the material) and for
 * initialisation failures. Lookup failures never come here; they return
 * EOS_FAIL_VALUE.
 */
#define EOS_FATAL(...) do { fprintf(stderr, __VA_ARGS__); assert(0); exit(1); } while (0)
#define EOS_UNKNOWN_MATERIAL(material) EOS_FATAL("%s was called for the unknown material %d.\n", __func__, (material)->iMat)
#define EOS_NOT_IMPLEMENTED(material) EOS_FATAL("%s is not implemented for material %d (%s).\n", __func__, (material)->iMat, (material)->MatString)

/*
 * Initialization of the material structures
 * iMat: Material number
 * dKpcUnit, dMsolUnit: Unit conversion factors
 * additional_data: EOS-specific options; a struct igeosParam for the ideal gas, NULL otherwise
 *
 * Required: all
 */
EOSMATERIAL *EOSinitMaterial(int iMat, double dKpcUnit, double dMsolUnit, const void * additional_data)
{
    EOSMATERIAL *material;
    material = (EOSMATERIAL *) calloc(1, sizeof(EOSMATERIAL));
    assert(material != NULL);
    material->iMat = iMat;
    material->bEntropyTableInit = EOS_FALSE;
    material->bEntropy = EOS_FALSE;

    /*
     * Unit system (the same rule in EOSlib and every sublibrary): both units > 0 and finite are code
     * units, both <= 0 is cgs; anything else (exactly one of them <= 0, NaN, inf) is a caller error.
     */
    if (!(((dKpcUnit > 0.0) && (dMsolUnit > 0.0) && isfinite(dKpcUnit) && isfinite(dMsolUnit)) || ((dKpcUnit <= 0.0) && (dMsolUnit <= 0.0)))) {
        fprintf(stderr, "EOSinitMaterial: invalid unit system dKpcUnit= %g dMsolUnit= %g (give both > 0 for code units or both <= 0 for cgs).\n", dKpcUnit, dMsolUnit);
        assert(0);
        exit(1);
    }

    /* Unit factors (code -> cgs), with the constants of the sublibraries. */
    if (dKpcUnit > 0.0 && dMsolUnit > 0.0) {
        const double MSOLG = 1.99e33;        /* solar mass in grams */
        const double GCGS = 6.67e-8;         /* G in cgs */
        const double KPCCM = 3.085678e21;    /* kiloparsec in centimeters */
        material->dErgPerGmUnit = GCGS*dMsolUnit*MSOLG/(dKpcUnit*KPCCM);
        material->dGmPerCcUnit = (dMsolUnit*MSOLG)/pow(dKpcUnit*KPCCM, 3.0);
    } else {
        material->dErgPerGmUnit = 1.0;
        material->dGmPerCcUnit = 1.0;
    }

    if (iMat == MAT_IDEALGAS)
    {
        if (additional_data == NULL) EOS_FATAL("EOSinitMaterial: the ideal gas (iMat= 0) needs a struct igeosParam as additional_data.\n");

        /* Pass additional parameter to the ideal gas EOS library. */
        struct igeosParam *Param;
        Param = (struct igeosParam *) additional_data;

        material->matType = EOSIDEALGAS;
        material->igeosmaterial = igeosInitMat(iMat, Param->dConstGamma, Param->dMeanMolMass, dKpcUnit, dMsolUnit);
        material->rho0 = material->igeosmaterial->rho0;
        material->rho0forInterfaceCorrection = 0.0; // Gas has no cutoff
        material->minSoundSpeed = 0;
        /* No entropy look up table required. */
        material->bEntropyTableInit = EOS_TRUE;
        material->bEntropy = EOS_TRUE;
        strcpy(material->MatString, "IDEAL GAS");
    } else if (iMat>=MAT_TILLOTSON_MIN && iMat<=MAT_TILLOTSON_MAX)
    {
#ifdef HAVE_TILLOTSON_H
        material->matType = EOSTILLOTSON;
        material->tillmaterial = tillInitMaterial(iMat, dKpcUnit, dMsolUnit);
        material->rho0 = material->tillmaterial->rho0;
        material->rho0forInterfaceCorrection = material->rho0;
        material->minSoundSpeed = sqrt(material->tillmaterial->A/material->tillmaterial->rho0);
        tilliMatString(material->tillmaterial, material->MatString);
        /*
         * Only the grid of the standard table (TILL_TABLE_* in tillotson.h: 1e-4..100 rho0,
         * 1 K..1e7 K at rho0) and its two edge isentropes, the cold curve and the top isentrope:
         * all that P, c, T(rho, u), u(rho, T), rho(P, T) and the table status need (77 kB, ~1 ms).
         * The full isentrope table (46 MB, ~0.6 s), which only EOSIsentropic needs, is built by
         * EOSinitIsentropicLookup. Up to 1.2.0 it was built here, for every material instance:
         * with one instance per thread 6 GB and over a minute of startup for 64 threads x 2
         * materials, also in codes that never call EOSIsentropic.
         */
        tillInitColdCurveDefault(material->tillmaterial);
        material->bEntropyTableInit = EOS_FALSE;
#else
        fprintf(stderr,"Tried to initialize a Tillotson material, but the Tillotson library is absent.\n");
        assert(0);
        exit(1);
#endif
    } else if (iMat>=MAT_ANEOS_MIN && iMat <=MAT_ANEOS_MAX)
    {
#ifdef HAVE_ANEOSMATERIAL_H
        material->matType = EOSANEOS;
        material->ANEOSmaterial = ANEOSinitMaterial(iMat, dKpcUnit, dMsolUnit);
        material->rho0 = ANEOSgetRho0(material->ANEOSmaterial);
        material->rho0forInterfaceCorrection = material->rho0;
        // The entropy look up table is initialized in ANEOSinitMaterial
        material->bEntropyTableInit = EOS_TRUE;
        material->bEntropy = EOS_TRUE;
        material->minSoundSpeed = ANEOSCofRhoT(material->ANEOSmaterial, material->rho0, ANEOSTmin(material->ANEOSmaterial));
        strcpy(material->MatString, material->ANEOSmaterial->matName);
#else
        fprintf(stderr,"Tried to initialize an ANEOSmaterial material, but the ANEOSmaterial library is absent.\n");
        assert(0);
        exit(1);
#endif
    } else if ((iMat >= MAT_REOS3_MIN) && (iMat <= MAT_REOS3_MAX)) {
#ifdef HAVE_REOS3_H
        material->matType = EOSREOS3;
        /* The last argument shifts u so that the whole table is positive (bSetRefU). */
        material->reos3material = reos3InitMaterial(iMat, dKpcUnit, dMsolUnit, 1);
        /*
         * rho0 of a REOS3 material is a convention (REOS3_RHO0_CGS, 1e-3 g/cc
         * unless outside the table), see the SCvH branch below.
         */
        material->rho0 = material->reos3material->rho0;
        material->rho0forInterfaceCorrection = 0.0; // Gas has no cutoff
        // So far we do not have entropy for REOS3.
        material->bEntropyTableInit = EOS_FALSE;
        // How do we define a minimum sound speed and reference density for H and He?
        material->minSoundSpeed = 0.0;
        snprintf(material->MatString, sizeof(material->MatString), "%s", material->reos3material->MatString);
#else
        fprintf(stderr,"Tried to initialize an reos3 material, but the reos3 library is absent.\n");
        assert(0);
        exit(1);
#endif
    } else if ((iMat >= MAT_SCVHEOS_MIN) && (iMat <= MAT_SCVHEOS_MAX)) {
#ifdef HAVE_SCVHEOS_H
        material->matType = EOSSCVHEOS;
        material->scvheosmaterial = scvheosInitMaterial(iMat, dKpcUnit, dMsolUnit);
        material->rho0 = material->scvheosmaterial->rho0;
        /*
         * rho0 of an SCvH material is a convention (SCVHEOS_RHO0_CGS, 1e-3 g/cc
         * unless outside the table): a gas has no reference density. It only
         * serves as the starting point of EOSSolveBC's lower bracket, which
         * shrinks by 0.9 up to EOS_SOLVEBC_MAX_EXPAND times, so interface
         * densities down to ~3e-8 g/cc are found from it. rho0forInterfaceCorrection
         * is 0 because the Woolfson cutoff 0.8*rho0 has no meaning for a gas.
         */
        material->rho0forInterfaceCorrection = 0.0; // Gas has no cutoff
        material->bEntropyTableInit = EOS_TRUE;
        material->bEntropy = EOS_TRUE;
        // How do we define a minimum sound speed and reference density for H and He?
        material->minSoundSpeed = 0.0;
        snprintf(material->MatString, sizeof(material->MatString), "%s", material->scvheosmaterial->MatString);
#else
        fprintf(stderr,"Tried to initialize an SCvH EOS material, but the SCvH EOS library is absent.\n");
        assert(0);
        exit(1);
#endif
    } else {
        fprintf(stderr, "EOSinitMaterial: iMat %i does not exist.\n",iMat);
        assert(0);
        exit(1);
    }

    /* Make sure that the flags are sent consistently. */
    if (material->bEntropy) assert(material->bEntropyTableInit);
    return material;
}

/*
 * Initialize the lookup tables needed for the isentropic evolution
 *
 * Required: pkdgrav3, gasoline (ISPH), fix initial conditions
 */
void EOSinitIsentropicLookup(EOSMATERIAL *material, const void * additional_data)
{
    if (material->bEntropyTableInit == EOS_TRUE) return;
    switch(material->matType)
    {
        case EOSIDEALGAS:
            // nothing to do, material->bEntropyTableInit is always EOS_TRUE
            break;
#ifdef HAVE_TILLOTSON_H
        case EOSTILLOTSON:
            /* The full isentrope table on the grid of the cold curve built in EOSinitMaterial. */
            tillInitLookupDefault(material->tillmaterial);
            material->bEntropyTableInit = EOS_TRUE;
            break;
#endif
#ifdef HAVE_ANEOSMATERIAL_H
        case EOSANEOS:
            // nothing to do, material->bEntropyTableInit is always EOS_TRUE
            break;
#endif
#ifdef HAVE_REOS3_H
        case EOSREOS3:
            // Not implemented yet: warn, do not fail.
            fprintf(stderr, "REOS3: Warning entropy lookup table not implemented yet.\n");
            break;
#endif
#ifdef HAVE_SCVHEOS_H
        case EOSSCVHEOS:
            // nothing to do, material->bEntropyTableInit is always EOS_TRUE
            break;
#endif
        default: EOS_UNKNOWN_MATERIAL(material);
    }
}

/*
 * Finalization of the material structures
 *
 * Required: all
 */
void EOSfinalizeMaterial(EOSMATERIAL *material)
{
    switch(material->matType)
    {
        case EOSIDEALGAS:
            igeosFinalizeMat(material->igeosmaterial);
            break;
#ifdef HAVE_TILLOTSON_H
        case EOSTILLOTSON:
            tillFinalizeMaterial(material->tillmaterial);
            break;
#endif
#ifdef HAVE_ANEOSMATERIAL_H
        case EOSANEOS:
            ANEOSfinalizeMaterial(material->ANEOSmaterial);
            break;
#endif
#ifdef HAVE_REOS3_H
        case EOSREOS3:
            reos3FinalizeMaterial(material->reos3material);
            break;
#endif
#ifdef HAVE_SCVHEOS_H
        case EOSSCVHEOS:
            scvheosFinalizeMaterial(material->scvheosmaterial);
            break;
#endif
        default: EOS_UNKNOWN_MATERIAL(material);
    }

    free(material);
}

/*
 * Print material data, e.g., EOS coefficients or the size of the EOS table.
 *
 * Required: ballic*
 */
void EOSPrintMat(EOSMATERIAL *material, FILE *fp)
{
    switch(material->matType)
    {
        case EOSIDEALGAS:
            igeosPrintMat(material->igeosmaterial, fp);
            break;
#ifdef HAVE_TILLOTSON_H
        case EOSTILLOTSON:
            tillPrintMat(material->tillmaterial, fp);
            break;
#endif
#ifdef HAVE_ANEOSMATERIAL_H
        case EOSANEOS:
            ANEOSPrintMat(material->ANEOSmaterial, fp);
            break;
#endif
#ifdef HAVE_REOS3_H
        case EOSREOS3:
            reos3PrintMat(material->reos3material, fp);
            break;
#endif
#ifdef HAVE_SCVHEOS_H
        case EOSSCVHEOS:
            scvheosPrintMat(material->scvheosmaterial, fp);
            break;
#endif
        default: EOS_UNKNOWN_MATERIAL(material);
    }
}

/*
 * Calculates the pressure P(rho,u) for a material
 *
 * Required: all
 */
double EOSPofRhoU(EOSMATERIAL *material, double rho, double u)
{
    double P = 0;

    switch(material->matType)
    {
        case EOSIDEALGAS:
            P = igeosPofRhoU(material->igeosmaterial, rho, u);
            break;
#ifdef HAVE_TILLOTSON_H
        case EOSTILLOTSON:
            P = tillPressure(material->tillmaterial, rho, u);
            break;
#endif
#ifdef HAVE_ANEOSMATERIAL_H
        case EOSANEOS:
            P = ANEOSPofRhoU(material->ANEOSmaterial, rho, u);
            break;
#endif
#ifdef HAVE_REOS3_H
        case EOSREOS3:
            P = reos3PofRhoU(material->reos3material, rho, u);
            break;
#endif
#ifdef HAVE_SCVHEOS_H
        case EOSSCVHEOS:
            P = scvheosPofRhoU(material->scvheosmaterial, rho, u);
            break;
#endif
        default: EOS_UNKNOWN_MATERIAL(material);
    }

    return P;
}

/*
 * Calculates the pressure P(rho,T) for a material
 *
 * Returns EOS_FAIL_VALUE if (rho, T) is outside the table of the material.
 * For Tillotson this is the case where u(rho, T) needs the cold curve and rho
 * lies outside the isentrope table (rho > rho_max): the failure value of
 * tillURhoTemp must not be passed on to tillPressure, which would clip the
 * resulting pressure to 0 and hide the failure.
 *
 * Required: pkdgrav3 (interface correction), tipsytools/tipsyprofile
 */
double EOSPofRhoT(EOSMATERIAL *material, double rho, double T)
{
    double P = 0;

    switch(material->matType)
    {
        case EOSIDEALGAS:
            P = igeosPofRhoT(material->igeosmaterial, rho, T);
            break;
#ifdef HAVE_TILLOTSON_H
        case EOSTILLOTSON:
        {
            /* The failure is already counted in nLookupFail by tillColdULookup. */
            double u = tillURhoTemp(material->tillmaterial, rho, T);
            P = tillFailed(u) ? EOS_FAIL_VALUE : tillPressure(material->tillmaterial, rho, u);
            break;
        }
#endif
#ifdef HAVE_ANEOSMATERIAL_H
        case EOSANEOS:
            P = ANEOSPofRhoT(material->ANEOSmaterial, rho, T);
            break;
#endif
#ifdef HAVE_REOS3_H
        case EOSREOS3:
            P = reos3PofRhoT(material->reos3material, rho, T);
            break;
#endif
#ifdef HAVE_SCVHEOS_H
        case EOSSCVHEOS:
            P = scvheosPofRhoT(material->scvheosmaterial, rho, T);
            break;
#endif
        default: EOS_UNKNOWN_MATERIAL(material);
    }

    return P;
}

/*
 * Calculates the sound speed c(rho,u) for a material
 *
 * Required: tipsytools
 */
double EOSCofRhoU(EOSMATERIAL *material, double rho, double u)
{
    double c = 0;

    switch(material->matType)
    {
        case EOSIDEALGAS:
            c = igeosCofRhoU(material->igeosmaterial, rho, u);
            break;
#ifdef HAVE_TILLOTSON_H
        case EOSTILLOTSON:
            tillPressureSound(material->tillmaterial, rho, u, &c);
            break;
#endif
#ifdef HAVE_ANEOSMATERIAL_H
        case EOSANEOS:
            c = ANEOSCofRhoU(material->ANEOSmaterial, rho, u);
            break;
#endif
#ifdef HAVE_REOS3_H
        case EOSREOS3:
            c = reos3CsofRhoU(material->reos3material, rho, u);
            break;
#endif
#ifdef HAVE_SCVHEOS_H
        case EOSSCVHEOS:
            c = scvheosCsofRhoU(material->scvheosmaterial, rho, u);
            break;
#endif
        default: EOS_UNKNOWN_MATERIAL(material);
    }

    return c;
}

/*
 * Calculates the pressure P(rho,u) and the sound speed c(rho,u) for a material
 * This should be used wherever both the pressure and the sound speed are needed for the same rho and u
 * as it is faster for some material models (esp. ANEOS)
 *
 * Required: gasoline
 */
double EOSPCofRhoU(EOSMATERIAL *material, double rho, double u, double *c)
{
    double P = 0;
    double T = 0;
    (void) T;   /* only the tabulated EOS fill it */

    switch(material->matType)
    {
        case EOSIDEALGAS:
            P = igeosPCofRhoU(material->igeosmaterial, rho, u, c);
            break;
#ifdef HAVE_TILLOTSON_H
        case EOSTILLOTSON:
            P = tillPressureSound(material->tillmaterial, rho, u, c);
            break;
#endif
#ifdef HAVE_ANEOSMATERIAL_H
        case EOSANEOS:
            P = ANEOSPCTofRhoU(material->ANEOSmaterial, rho, u, c, &T);
            break;
#endif
#ifdef HAVE_REOS3_H
        case EOSREOS3:
            P = reos3PCTofRhoU(material->reos3material, rho, u, c, &T);
            break;
#endif
#ifdef HAVE_SCVHEOS_H
        case EOSSCVHEOS:
            P = scvheosPCTofRhoU(material->scvheosmaterial, rho, u, c, &T);
            break;
#endif
        default: EOS_UNKNOWN_MATERIAL(material);
    }

    return P;
}

/*
 * Calculates the pressure P(rho,u), the sound speed c(rho,u) and the temperature T(rho,u) for a material
 * This should be used wherever all three values are are needed for the same rho and u
 * as it is faster for some material models (esp. ANEOS)
 *
 * Required: pkdgrav3
 */
double EOSPCTofRhoU(EOSMATERIAL *material, double rho, double u, double *c, double *T)
{
    double P = 0;

    switch(material->matType)
    {
        case EOSIDEALGAS:
            /* One call: a state below the cold curve u = 0 counted once (PC and T counted it twice). */
            P = igeosPCTofRhoU(material->igeosmaterial, rho, u, c, T);
            break;
#ifdef HAVE_TILLOTSON_H
        case EOSTILLOTSON:
            P = tillPressureSound(material->tillmaterial, rho, u, c);
            *T = tillTempRhoU(material->tillmaterial, rho, u);
            break;
#endif
#ifdef HAVE_ANEOSMATERIAL_H
        case EOSANEOS:
            P = ANEOSPCTofRhoU(material->ANEOSmaterial, rho, u, c, T);
            break;
#endif
#ifdef HAVE_REOS3_H
        case EOSREOS3:
            P = reos3PCTofRhoU(material->reos3material, rho, u, c, T);
            break;
#endif
#ifdef HAVE_SCVHEOS_H
        case EOSSCVHEOS:
            P = scvheosPCTofRhoU(material->scvheosmaterial, rho, u, c, T);
            break;
#endif
        default: EOS_UNKNOWN_MATERIAL(material);
    }

    return P;
}

/*
 * Calculates the internal energy u2(rho1, u1, rho2) after an isentropic evolution for a material
 *
 * Required: pkdgrav3 (ISPH), gasoline (ISPH), fix initial conditions
 */
double EOSIsentropic(EOSMATERIAL *material, double rho1, double u1, double rho2)
{
    double u2 = 0;
    switch(material->matType)
    {
        case EOSIDEALGAS:
            u2 = igeosIsentropicU(material->igeosmaterial, rho1, u1, rho2);
            break;
#ifdef HAVE_TILLOTSON_H
        case EOSTILLOTSON:
            /* The isentrope table is built by EOSinitIsentropicLookup, not by EOSinitMaterial: without it this is a programming error. */
            if (material->bEntropyTableInit != EOS_TRUE) {
                fprintf(stderr, "EOSIsentropic: the isentrope table of material %d (%s) is not built. "
                        "Call EOSinitIsentropicLookup after EOSinitMaterial for a Tillotson material.\n",
                        material->iMat, material->MatString);
                assert(0);
                exit(1);
            }
            u2 = tillLookupU(material->tillmaterial, rho1, u1, rho2, 0); // last argument is actually a particle number
            break;
#endif
#ifdef HAVE_ANEOSMATERIAL_H
        case EOSANEOS:
            u2 = ANEOSisentropicU(material->ANEOSmaterial, rho1, u1, rho2);
            break;
#endif
#ifdef HAVE_REOS3_H
        case EOSREOS3:
            EOS_NOT_IMPLEMENTED(material);   /* REOS3 has no entropy table */
            break;
#endif
#ifdef HAVE_SCVHEOS_H
        case EOSSCVHEOS:
            u2 = scvheosIsentropicU(material->scvheosmaterial, rho1, u1, rho2);
            break;
#endif
        default: EOS_UNKNOWN_MATERIAL(material);
    }

    return u2;
}

/*
 * Calculates the specific entropy S(rho, u) for a material.
 *
 * Returned in cgs (erg/g/K), NOT in code units; additive constant is
 * library-dependent. Returns EOS_FAIL_VALUE if not available.
 *
 * Required: none
 */
double EOSSofRhoU(EOSMATERIAL *material, double rho, double u)
{
    double S = 0;

    if (!material->bEntropy) return EOS_FAIL_VALUE;

    switch(material->matType)
    {
        case EOSIDEALGAS:
            S = igeosSofRhoU(material->igeosmaterial, rho, u);
            break;
#ifdef HAVE_TILLOTSON_H
        case EOSTILLOTSON:
            S = EOS_FAIL_VALUE;   /* no entropy for this EOS (bEntropy is EOS_FALSE, so this is not reached) */
            break;
#endif
#ifdef HAVE_ANEOSMATERIAL_H
        case EOSANEOS:
            S = ANEOSSofRhoU(material->ANEOSmaterial, rho, u);
            break;
#endif
#ifdef HAVE_REOS3_H
        case EOSREOS3:
            S = EOS_FAIL_VALUE;   /* no entropy for this EOS (bEntropy is EOS_FALSE, so this is not reached) */
            break;
#endif
#ifdef HAVE_SCVHEOS_H
        case EOSSCVHEOS:
            /* One cell lookup, counted once (T(rho, u) then s(rho, T) counted a density outside the table twice). */
            S = scvheosSofRhoU(material->scvheosmaterial, rho, u);
            break;
#endif
        default: EOS_UNKNOWN_MATERIAL(material);
    }

    return S;
}

/*
 * Calculates the specific entropy S(rho, T) for a material.
 *
 * Returned in cgs (erg/g/K), NOT in code units; additive constant is
 * library-dependent. Returns EOS_FAIL_VALUE if not available.
 *
 * Required: none
 */
double EOSSofRhoT(EOSMATERIAL *material, double rho, double T)
{
    double S = 0;

    if (!material->bEntropy) return EOS_FAIL_VALUE;

    switch(material->matType)
    {
        case EOSIDEALGAS:
            S = igeosSofRhoT(material->igeosmaterial, rho, T);
            break;
#ifdef HAVE_TILLOTSON_H
        case EOSTILLOTSON:
            S = EOS_FAIL_VALUE;   /* no entropy for this EOS (bEntropy is EOS_FALSE, so this is not reached) */
            break;
#endif
#ifdef HAVE_ANEOSMATERIAL_H
        case EOSANEOS:
            S = ANEOSSofRhoT(material->ANEOSmaterial, rho, T);
            break;
#endif
#ifdef HAVE_REOS3_H
        case EOSREOS3:
            S = EOS_FAIL_VALUE;   /* no entropy for this EOS (bEntropy is EOS_FALSE, so this is not reached) */
            break;
#endif
#ifdef HAVE_SCVHEOS_H
        case EOSSCVHEOS:
            S = scvheosSofRhoT(material->scvheosmaterial, rho, T);
            break;
#endif
        default: EOS_UNKNOWN_MATERIAL(material);
    }

    return S;
}

/*
 * Yield strength parameters of a material in code units, gathered once per evaluation. Returns the
 * yield strength model:
 *   -2: no strength model exists for this EOS (Tillotson); EOSYieldStrength is fatal
 *   -1: no yield stress (ANEOS material without a strength file)
 *    0: no limiter (gases, or a strength file with model 0)
 *    1: pressure-dependent yield strength of intact material (Collins et al. 2004)
 *    2: pressure-dependent yield strength of damaged material
 */
typedef struct {
    double Y0;    // Yield strength at 0 pressure
    double YM;    // Yield strength at infinite pressure
    double mui;   // coefficient of internal friction for intact material
    double mud;   // coefficient of internal friction for damaged material (-1: not defined)
    double xi;    // thermal softening parameter
} EOSYIELDPARAMS;

static int EOSYieldParams(EOSMATERIAL *material, EOSYIELDPARAMS *p)
{
    p->Y0 = 0.0;
    p->YM = 0.0;
    p->mui = 0.0;
    p->mud = -1.0;
    p->xi = 0.0;

    switch(material->matType)
    {
        case EOSIDEALGAS:
            return 0;   // No yield strength for gas type EOS
#ifdef HAVE_TILLOTSON_H
        case EOSTILLOTSON:
            return -2;  // Change this if/when yield is added to tillotson
#endif
#ifdef HAVE_ANEOSMATERIAL_H
        case EOSANEOS:
            return ANEOSYieldParameters(material->ANEOSmaterial, &p->Y0, &p->YM, &p->mui, &p->mud, &p->xi);
#endif
#ifdef HAVE_REOS3_H
        case EOSREOS3:
            return 0;   // No yield strength for gas type EOS
#endif
#ifdef HAVE_SCVHEOS_H
        case EOSSCVHEOS:
            return 0;   // No yield strength for gas type EOS
#endif
        default: EOS_UNKNOWN_MATERIAL(material);
    }
    return 0;   /* not reached */
}

/*
 * Return the yield strength model of the material (see EOSYieldParams for the codes).
 */
int EOSYieldModel(EOSMATERIAL *material){
    EOSYIELDPARAMS p;
    return EOSYieldParams(material, &p);
}

/*
 * Strength limiter from (rho, T, P) for the model and parameters from EOSYieldParams. The two
 * public entry points gather the parameters once (they did it twice in EOSYieldStrength: once in
 * EOSYieldModel and once here) and differ only in how T and P are obtained.
 */
static double EOSYieldStrengthCore(EOSMATERIAL *material, int yieldStrengthModel, const EOSYIELDPARAMS *p,
                                   double rho, double T, double P)
{
    double Y = 0.0;

    /*
     * Yield strength models (EOSYieldModel):
     * -2: No strength model for this EOS (Tillotson): fatal
     * -1: No yield stress
     *  0: No limiter, so the yield strength returned is irrelevant (not called).
     *  1: Pressure dependent yield strength for intact material (Collins et al. 2004)
     *  2: Pressure dependent yield strength for damaged material
     *
     * The pressure dependence is evaluated for P_eff = max(P, 0): in tension the
     * intact strength reduces to the cohesion Y0 and the damaged (cohesionless)
     * strength vanishes. This also keeps the intact-strength denominator
     * 1 + mui*P/(YM-Y0) strictly positive.
     */
    if (yieldStrengthModel == -2) EOS_NOT_IMPLEMENTED(material);   /* no strength model for Tillotson yet */

    /*
     * The parameters are validated when they are read (ANEOSReadStrengthParameters:
     * YM > Y0, mui > 0, xi > 0, mud >= 0 for model 2, model 0..2), so the checks
     * below cannot fire for a material from EOSinitMaterial.
     */
    if (yieldStrengthModel > 0) {
        /*
         * A failed pressure (EOS_FAIL_VALUE from a lookup outside the table, or a
         * failed P passed to EOSYieldStrengthofRhoTP) must not reach the tension
         * clamp below: fmax(-1e50, 0) = 0 turned it into the finite strength
         * Y0*f(T). EOSFailed also catches NaN. Not counted in nLookupFail: the
         * failed lookup that produced P was counted already.
         */
        if (EOSFailed(P)) {
            EOS_WARN("EOSYieldStrength: failed pressure P= %g for material %d (rho= %g T= %g).\n", P, material->iMat, rho, T);
            return EOS_FAIL_VALUE;
        }
        double Peff = fmax(P, 0.0);
        double Yi = p->Y0 + p->mui * Peff / (1.0 + p->mui * Peff / (p->YM - p->Y0));
        double Tmelt = 0.0;

        if (yieldStrengthModel == 1) {
            Y = Yi;
        } else if (yieldStrengthModel == 2) {
            if (p->mud < 0.0) EOS_FATAL("EOSYieldStrength: yieldStrengthModel 2 requested for material %d, but mud not defined.\n", material->iMat);
            Y = fmin(Yi, p->mud * Peff);
        } else {
            EOS_FATAL("EOSYieldStrength: yieldStrengthModel = %d not supported (material %d).\n", yieldStrengthModel, material->iMat);
        }

#ifdef HAVE_ANEOSMATERIAL_H
        /* Only ANEOS materials have a model > 0. */
        Tmelt = ANEOSTmeltofRho(material->ANEOSmaterial, rho);
#endif

        /*
         * Thermal softening (Ohnaka 1995), zero at and above the melt
         * temperature. T comes from a lookup (T_min fallback, never the
         * sentinel) and T_melt from the melt curve (sentinel for an invalid
         * rho): an invalid pair is a lookup failure, so return the sentinel
         * rather than exit.
         */
        if (!(T > 0.0) || !(Tmelt >= 0.0) || !isfinite(T)) {
            EOS_WARN("EOSYieldStrength: invalid T= %g or Tmelt= %g for material %d (rho= %g P= %g).\n", T, Tmelt, material->iMat, rho, P);
            return EOS_FAIL_VALUE;
        }
        Y *= fmax(tanh(p->xi * (Tmelt / T - 1.0)), 0.0);

        /* Guard against roundoff; the strength is never negative. */
        Y = fmax(Y, 0.0);
    }

    return Y;
}

/*
 * A state above the u column of the table (overheated material): T(rho, u) falls back to T_min by
 * contract and P, c to their values at (rho, T_min), so the strength was the full cold strength
 * Y0 f(T_min) of a material hotter than anything tabulated. Such a state is clamped to the top row
 * instead, T = T_max and P = P(rho, T_max), the convention of the isentropes for a state above the
 * table; the thermal softening then gives 0 wherever T_melt(rho) <= T_max.
 *
 * The fallback T is exactly TAxis[0], so only T == T_min needs the status test (ANEOSIsInTable,
 * one inversion of the u column; it agrees exactly with the lookup, also on the M-ANEOS columns
 * whose cold rows have du/dT < 0). A state in the table at T_min, or below the cold curve (cold
 * material, full strength), is left alone. Not counted again: the lookup that produced T counted
 * the fallback. Only ANEOS has strength models > 0.
 */
static void EOSYieldClampAboveTable(EOSMATERIAL *material, double rho, double u, double *pT, double *pP)
{
#ifdef HAVE_ANEOSMATERIAL_H
    if (material->matType == EOSANEOS) {
        ANEOSMATERIAL *m = material->ANEOSmaterial;
        if (*pT == m->TAxis[0] && !EOSFailed(*pP) && ANEOSIsInTable(m, rho, u) == ANEOS_OUTSIDE_TMAX) {
            *pT = m->TAxis[m->nT - 1];
            *pP = ANEOSPofRhoT(m, rho, *pT);
        }
    }
#else
    (void) material; (void) rho; (void) u; (void) pT; (void) pP;
#endif
}

/*
 * Calculate strength limiter from (rho, u). Derives T and P with one fused
 * lookup. A state above the table is evaluated on the top row (EOSYieldClampAboveTable).
 */
double EOSYieldStrength(EOSMATERIAL *material, double rho, double u)
{
    EOSYIELDPARAMS p;
    double T = 0.0, P = 0.0, c;
    int yieldStrengthModel = EOSYieldParams(material, &p);

    if (yieldStrengthModel > 0) {
        P = EOSPCTofRhoU(material, rho, u, &c, &T);
        EOSYieldClampAboveTable(material, rho, u, &T, &P);
    }
    return EOSYieldStrengthCore(material, yieldStrengthModel, &p, rho, T, P);
}

/*
 * Calculate strength limiter from (rho, T, P) supplied by the caller. Evaluated as given: a T from
 * the T_min fallback of a state above the table cannot be told from a cold state (use
 * EOSYieldStrength for a (rho, u) state).
 */
double EOSYieldStrengthofRhoTP(EOSMATERIAL *material, double rho, double T, double P)
{
    EOSYIELDPARAMS p;
    int yieldStrengthModel = EOSYieldParams(material, &p);

    return EOSYieldStrengthCore(material, yieldStrengthModel, &p, rho, T, P);
}

/*
 * Calculates the temperature T(rho,u) for a material
 *
 * Required: all
 */
double EOSTofRhoU(EOSMATERIAL *material, double rho, double u)
{
    double T = 0;
    switch(material->matType)
    {
        case EOSIDEALGAS:
            T = igeosTofRhoU(material->igeosmaterial, rho, u);
            break;
#ifdef HAVE_TILLOTSON_H
        case EOSTILLOTSON:
            T = tillTempRhoU(material->tillmaterial, rho, u);
            break;
#endif
#ifdef HAVE_ANEOSMATERIAL_H
        case EOSANEOS:
            T = ANEOSTofRhoU(material->ANEOSmaterial, rho, u);
            break;
#endif
#ifdef HAVE_REOS3_H
        case EOSREOS3:
            T = reos3TofRhoU(material->reos3material, rho, u);
            break;
#endif
#ifdef HAVE_SCVHEOS_H
        case EOSSCVHEOS:
            T = scvheosTofRhoU(material->scvheosmaterial, rho, u);
            break;
#endif
        default: EOS_UNKNOWN_MATERIAL(material);
    }

    return T;
}

/*
 * Calculates the internal energy u(rho,T) for a material
 *
 * Required: pkdgrav3, tipsytools
 */
double EOSUofRhoT(EOSMATERIAL *material, double rho, double T)
{
    double u = 0;
    switch(material->matType)
    {
        case EOSIDEALGAS:
            u = igeosUofRhoT(material->igeosmaterial, rho, T);
            break;
#ifdef HAVE_TILLOTSON_H
        case EOSTILLOTSON:
            u = tillURhoTemp(material->tillmaterial, rho, T);
            break;
#endif
#ifdef HAVE_ANEOSMATERIAL_H
        case EOSANEOS:
            u = ANEOSUofRhoT(material->ANEOSmaterial, rho, T);
            break;
#endif
#ifdef HAVE_REOS3_H
        case EOSREOS3:
            u = reos3UofRhoT(material->reos3material, rho, T);
            break;
#endif
#ifdef HAVE_SCVHEOS_H
        case EOSSCVHEOS:
            u = scvheosUofRhoT(material->scvheosmaterial, rho, T);
            break;
#endif
        default: EOS_UNKNOWN_MATERIAL(material);
    }

    return u;
}

/*
 * Calculates the density rho(p,T) for a material
 *
 * Required: internal, interface correction (gasoline,pkdgrav3)
 */
double EOSRhoofPT(EOSMATERIAL *material, double p, double T)
{
    double rho = 0;
    switch(material->matType)
    {
        case EOSIDEALGAS:
            rho = igeosRhoofPT(material->igeosmaterial, p, T);
            break;
#ifdef HAVE_TILLOTSON_H
        case EOSTILLOTSON:
            rho = tillRhoPTemp(material->tillmaterial, p, T);
            break;
#endif
#ifdef HAVE_ANEOSMATERIAL_H
        case EOSANEOS:
            rho = ANEOSRhoofPT(material->ANEOSmaterial, p, T);
            break;
#endif
#ifdef HAVE_REOS3_H
        case EOSREOS3:
            rho = reos3RhoofPT(material->reos3material, p, T);
            break;
#endif
#ifdef HAVE_SCVHEOS_H
        case EOSSCVHEOS:
            rho = scvheosRhoofPT(material->scvheosmaterial, p, T);
            break;
#endif
        default: EOS_UNKNOWN_MATERIAL(material);
    }

    return rho;
}

/*
 * Calculates the density rho(u,T) for a material
 *
 * Required: none
 */
double EOSRhoofUT(EOSMATERIAL *material, double u, double T)
{
    double rho = 0;
    switch(material->matType)
    {
        case EOSIDEALGAS:
            EOS_NOT_IMPLEMENTED(material);
            break;
#ifdef HAVE_TILLOTSON_H
        case EOSTILLOTSON:
            EOS_NOT_IMPLEMENTED(material);
            break;
#endif
#ifdef HAVE_ANEOSMATERIAL_H
        case EOSANEOS:
            rho = ANEOSRhoofUT(material->ANEOSmaterial, u, T);
            break;
#endif
#ifdef HAVE_REOS3_H
        case EOSREOS3:
            EOS_NOT_IMPLEMENTED(material);
            break;
#endif
#ifdef HAVE_SCVHEOS_H
        case EOSSCVHEOS:
            EOS_NOT_IMPLEMENTED(material);
            break;
#endif
        default: EOS_UNKNOWN_MATERIAL(material);
    }

    return rho;
}

/*
 * Tests if the given combination is below the cold curve, returns EOS_TRUE if so and EOS_FALSE if not.
 *
 * Returns EOS_FAIL if the cold curve cannot be evaluated at rho (outside the table; the failure is
 * counted in the sublibrary's nLookupFail) or u is not finite. Before, such a state was reported
 * as "not below", because u < EOS_FAIL_VALUE is false. Compare with == EOS_TRUE; a caller that
 * uses the truth value treats an undecidable state as below the cold curve, the safe side.
 *
 * For the tabulated EOS (ANEOS, REOS3, SCvH) "below the cold curve" is the table status of the
 * lookup, EOSIsInTable == EOS_OUTSIDE_VMIN: exactly the states that T(rho, u) cannot invert from
 * below and that fall back to T_min. It compared u with u(rho, T_min), the first table row, which
 * disagreed with EOSIsInTable (and the lookup) in two places:
 *   - M-ANEOS columns whose cold rows have du/dT < 0 (up to 18-58 K): states below u(rho, T_min)
 *     lie on the line at a slightly higher T, the lookup finds them without a fallback and
 *     EOSIsInTable reports EOS_SUCCESS, but this function said EOS_TRUE;
 *   - states within the edge tolerance below the first row, which the inversions clamp onto it.
 * Tillotson and the ideal gas have an analytic cold curve and keep u < u_cold(rho) (Tillotson
 * also below rho_min, where u_cold = 0 and EOSIsInTable reports EOS_OUTSIDE_RHOMIN).
 *
 * Required: ballic
 */
int EOSisbelowColdCurve(EOSMATERIAL *material, double rho, double u)
{
    /* u_cold first, also for the tabulated EOS: rho outside the table is counted there (EOSIsInTable is a status query and counts nothing). */
    double ucold = EOSUCold(material, rho);

    if (EOSFailed(ucold) || !isfinite(u)) {
        EOS_WARN("EOSisbelowColdCurve: cold curve not available (rho= %g u= %g u_cold= %g, material %d).\n", rho, u, ucold, material->iMat);
        return EOS_FAIL;
    }

    if (material->matType == EOSANEOS || material->matType == EOSREOS3 || material->matType == EOSSCVHEOS) {
        int st = EOSIsInTable(material, rho, u);
        if (st == EOS_OUTSIDE_VMIN) return EOS_TRUE;
        if (st == EOS_SUCCESS || st == EOS_OUTSIDE_VMAX) return EOS_FALSE;
        /* Not reached: rho is in the table (u_cold succeeded) and u is finite. */
        EOS_WARN("EOSisbelowColdCurve: table status %d at rho= %g u= %g (material %d).\n", st, rho, u, material->iMat);
        return EOS_FAIL;
    }
    return (u < ucold) ? EOS_TRUE : EOS_FALSE;
}

/*
 * Tests if the given (rho, u) lies inside the lookup table of the material.
 *
 * Returns EOS_SUCCESS if it does, otherwise one of EOS_OUTSIDE_RHOMIN,
 * EOS_OUTSIDE_RHOMAX, EOS_OUTSIDE_VMIN, EOS_OUTSIDE_VMAX indicating which
 * bound was violated, or EOS_FAIL if the sublibrary reports an unknown status.
 * No tolerance policy is applied here; use EOSIsUsable for that.
 *
 * Required: ballic, gasoline, tipsytools
 */
int EOSIsInTable(EOSMATERIAL *material, double rho, double u)
{
    int iret = EOS_FAIL;
    switch(material->matType)
    {
        case EOSIDEALGAS:
            /* No table; the cold curve is u = 0 (T_min = 0), as for Tillotson. */
            if (!(rho > 0.0)) return EOS_OUTSIDE_RHOMIN;
            if (!(u >= 0.0)) return EOS_OUTSIDE_VMIN;   /* also NaN */
            return EOS_SUCCESS;
#ifdef HAVE_TILLOTSON_H
        case EOSTILLOTSON:
            iret = tillIsInTable(material->tillmaterial, rho, u);
            if (iret == TILL_LOOKUP_SUCCESS) return EOS_SUCCESS;
            if (iret == TILL_LOOKUP_OUTSIDE_RHOMIN) return EOS_OUTSIDE_RHOMIN;
            if (iret == TILL_LOOKUP_OUTSIDE_RHOMAX) return EOS_OUTSIDE_RHOMAX;
            if (iret == TILL_LOOKUP_OUTSIDE_VMIN) return EOS_OUTSIDE_VMIN;
            if (iret == TILL_LOOKUP_OUTSIDE_VMAX) return EOS_OUTSIDE_VMAX;
            return EOS_FAIL;
#endif
#ifdef HAVE_ANEOSMATERIAL_H
        case EOSANEOS:
            iret = ANEOSIsInTable(material->ANEOSmaterial, rho, u);
            if (iret == ANEOS_SUCCESS) return EOS_SUCCESS;
            if (iret == ANEOS_OUTSIDE_RHOMIN) return EOS_OUTSIDE_RHOMIN;
            if (iret == ANEOS_OUTSIDE_RHOMAX) return EOS_OUTSIDE_RHOMAX;
            if (iret == ANEOS_OUTSIDE_TMIN) return EOS_OUTSIDE_VMIN;
            if (iret == ANEOS_OUTSIDE_TMAX) return EOS_OUTSIDE_VMAX;
            return EOS_FAIL;
#endif
#ifdef HAVE_REOS3_H
        case EOSREOS3:
            iret = reos3IsInExtrapLimit(material->reos3material, rho, u);
            if (iret == REOS3_SUCCESS) return EOS_SUCCESS;
            if (iret == REOS3_OUTSIDE_RHOMIN) return EOS_OUTSIDE_RHOMIN;
            if (iret == REOS3_OUTSIDE_RHOMAX) return EOS_OUTSIDE_RHOMAX;
            if (iret == REOS3_OUTSIDE_TMIN) return EOS_OUTSIDE_VMIN;
            if (iret == REOS3_OUTSIDE_TMAX) return EOS_OUTSIDE_VMAX;
            return EOS_FAIL;
#endif
#ifdef HAVE_SCVHEOS_H
        case EOSSCVHEOS:
            iret = scvheosIsInExtrapLimit(material->scvheosmaterial, rho, u);
            if (iret == SCVHEOS_SUCCESS) return EOS_SUCCESS;
            if (iret == SCVHEOS_OUTSIDE_RHOMIN) return EOS_OUTSIDE_RHOMIN;
            if (iret == SCVHEOS_OUTSIDE_RHOMAX) return EOS_OUTSIDE_RHOMAX;
            if (iret == SCVHEOS_OUTSIDE_TMIN) return EOS_OUTSIDE_VMIN;
            if (iret == SCVHEOS_OUTSIDE_TMAX) return EOS_OUTSIDE_VMAX;
            return EOS_FAIL;
#endif
        default: EOS_UNKNOWN_MATERIAL(material);
    }

    return iret;
}

/*
 * Policy layer on top of EOSIsInTable: can P, c and T be evaluated at (rho, u)?
 *
 * Analytic EOS (ideal gas, Tillotson) only need the lookup table for the
 * isentropic evolution, so u above the table is tolerated. Tabulated EOS
 * (ANEOS, REOS3, SCvH) cannot be evaluated outside the table at all.
 *
 * Tillotson below rho_min is usable as well: the cold curve is zero there (as
 * everywhere below rho0), so T(rho, u) = u/cv, and P and c are analytic; the
 * isentrope is integrated directly (tillLookupU). The state must still be valid:
 * rho > 0 and finite, u finite and not below the cold curve u = 0 (with the same
 * roundoff tolerance as tillIsInTable, so that the policy is continuous at rho_min).
 * Above rho_max the cold curve, u(rho, T) and T(rho, u) are not available, so
 * such states stay unusable.
 *
 * Required: pkdgrav3, gasoline, ballic
 */
int EOSIsUsable(EOSMATERIAL *material, double rho, double u)
{
    int iret = EOSIsInTable(material, rho, u);

    if (iret == EOS_SUCCESS) return EOS_TRUE;

    switch(material->matType)
    {
        case EOSIDEALGAS:
            /* Below the cold curve u = 0 (or rho <= 0): T, P and c are clamped, not usable. */
            return EOS_FALSE;
#ifdef HAVE_TILLOTSON_H
        case EOSTILLOTSON:
            /* Pressure is analytic; only the isentrope lookup needs u <= u_max. */
            if (iret == EOS_OUTSIDE_VMAX) return EOS_TRUE;
            /* Below rho_min (see above): the cold curve is u = 0. */
            if (iret == EOS_OUTSIDE_RHOMIN)
                return ((rho > 0.0) && isfinite(rho) && isfinite(u) &&
                        (u >= -tillUTol(material->tillmaterial, u))) ? EOS_TRUE : EOS_FALSE;
            return EOS_FALSE;
#endif
        default:
            /* Tabulated EOS: any out-of-table state is unusable. */
            return EOS_FALSE;
    }
}


/*
 * Calculates the derivative dPdrho(rho,u) for a material
 *
 * Required: ballic
 */
double EOSdPdRho(EOSMATERIAL *material, double rho, double u)
{
    double dPdRho = 0;
    switch(material->matType)
    {
        case EOSIDEALGAS:
            dPdRho = igeosdPdRho(material->igeosmaterial, rho, u);
            break;
#ifdef HAVE_TILLOTSON_H
        case EOSTILLOTSON:
            dPdRho = tilldPdrho(material->tillmaterial, rho, u);
            break;
#endif
#ifdef HAVE_ANEOSMATERIAL_H
        case EOSANEOS:
            dPdRho = ANEOSdPdRhoofRhoU(material->ANEOSmaterial, rho, u);
            break;
#endif
#ifdef HAVE_REOS3_H
        case EOSREOS3:
            dPdRho = reos3dPdRhoofRhoU(material->reos3material, rho, u);
            break;
#endif
#ifdef HAVE_SCVHEOS_H
        case EOSSCVHEOS:
            dPdRho = scvheosdPdRhoofRhoU(material->scvheosmaterial, rho, u);
            break;
#endif
        default: EOS_UNKNOWN_MATERIAL(material);
    }

    return dPdRho;
}

/*
 * Calculates the derivative dPdu(rho,u) for a material
 *
 * Required: ballic
 */
double EOSdPdU(EOSMATERIAL *material, double rho, double u)
{
    double dPdU = 0;
    switch(material->matType)
    {
        case EOSIDEALGAS:
            dPdU = igeosdPdU(material->igeosmaterial, rho, u);
            break;
#ifdef HAVE_TILLOTSON_H
        case EOSTILLOTSON:
            dPdU = tilldPdu(material->tillmaterial, rho, u);
            break;
#endif
#ifdef HAVE_ANEOSMATERIAL_H
        case EOSANEOS:
            dPdU = ANEOSdPdUofRhoU(material->ANEOSmaterial, rho, u);
            break;
#endif
#ifdef HAVE_REOS3_H
        case EOSREOS3:
            dPdU = reos3dPdUofRhoU(material->reos3material, rho, u);
            break;
#endif
#ifdef HAVE_SCVHEOS_H
        case EOSSCVHEOS:
            dPdU = scvheosdPdUofRhoU(material->scvheosmaterial, rho, u);
            break;
#endif
        default: EOS_UNKNOWN_MATERIAL(material);
    }

    return dPdU;
}

/*
 * Calculates the derivative dudrho(rho,u) = P/rho^2 for a material
 * (EOS_FAIL_VALUE where P(rho, u) fails: tabulated EOS outside the table).
 * The sentinel is passed on as it is, never divided by rho^2 (ANEOS and SCvH
 * did that up to 1.2.0; igeos and Tillotson have an analytic P).
 *
 * Required: ballic
 */
double EOSdUdRho(EOSMATERIAL *material, double rho, double u)
{
    double dUdRho = 0;
    switch(material->matType)
    {
        case EOSIDEALGAS:
            dUdRho = igeosPofRhoU(material->igeosmaterial, rho, u) / (rho * rho);
            break;
#ifdef HAVE_TILLOTSON_H
        case EOSTILLOTSON:
            dUdRho = tilldudrho(material->tillmaterial, rho, u);
            break;
#endif
#ifdef HAVE_ANEOSMATERIAL_H
        case EOSANEOS:
            dUdRho = ANEOSdUdRhoofRhoU(material->ANEOSmaterial, rho, u);
            break;
#endif
#ifdef HAVE_REOS3_H
        case EOSREOS3:
            dUdRho = reos3dUdRhoofRhoU(material->reos3material, rho, u);
            break;
#endif
#ifdef HAVE_SCVHEOS_H
        case EOSSCVHEOS:
            {
                /* No scvheosdUdRhoofRhoU: check the failed P here (see reos3dUdRhoofRhoU). */
                double P = scvheosPofRhoU(material->scvheosmaterial, rho, u);
                dUdRho = EOSFailed(P) ? EOS_FAIL_VALUE : P / (rho * rho);
            }
            break;
#endif
        default: EOS_UNKNOWN_MATERIAL(material);
    }

    return dUdRho;
}

/*
 * Calculates the derivative dPdrho(rho,T) for a material
 *
 * Required: none yet, maybe in future version of ballic
 */
double EOSdPdRhoofRhoT(EOSMATERIAL *material, double rho, double T)
{
    double dPdRho = 0;
    switch(material->matType)
    {
        case EOSIDEALGAS:
            EOS_NOT_IMPLEMENTED(material);
            break;
#ifdef HAVE_TILLOTSON_H
        case EOSTILLOTSON:
            EOS_NOT_IMPLEMENTED(material);
            break;
#endif
#ifdef HAVE_ANEOSMATERIAL_H
        case EOSANEOS:
            dPdRho = ANEOSdPdRhoofRhoT(material->ANEOSmaterial, rho, T);
            break;
#endif
#ifdef HAVE_REOS3_H
        case EOSREOS3:
            dPdRho = reos3dPdRhoofRhoT(material->reos3material, rho, T);
            break;
#endif
#ifdef HAVE_SCVHEOS_H
        case EOSSCVHEOS:
            dPdRho = scvheosdPdRhoofRhoT(material->scvheosmaterial, rho, T);
            break;
#endif
        default: EOS_UNKNOWN_MATERIAL(material);
    }

    return dPdRho;
}

/*
 * Calculates the derivative dPdT(rho,T) for a material
 *
 * Required: none yet, maybe in future version of ballic
 */
double EOSdPdT(EOSMATERIAL *material, double rho, double T)
{
    double dPdT = 0;
    switch(material->matType)
    {
        case EOSIDEALGAS:
            EOS_NOT_IMPLEMENTED(material);
            break;
#ifdef HAVE_TILLOTSON_H
        case EOSTILLOTSON:
            EOS_NOT_IMPLEMENTED(material);
            break;
#endif
#ifdef HAVE_ANEOSMATERIAL_H
        case EOSANEOS:
            dPdT = ANEOSdPdTofRhoT(material->ANEOSmaterial, rho, T);
            break;
#endif
#ifdef HAVE_REOS3_H
        case EOSREOS3:
            dPdT = reos3dPdTofRhoT(material->reos3material, rho, T);
            break;
#endif
#ifdef HAVE_SCVHEOS_H
        case EOSSCVHEOS:
            dPdT = scvheosdPdTofRhoT(material->scvheosmaterial, rho, T);
            break;
#endif
        default: EOS_UNKNOWN_MATERIAL(material);
    }

    return dPdT;
}

/*
 * Calculate the shear modulus Gamma(rho,T)
 * using a fixed guessed Poisson ratio nu
 */
double EOSGammaofRhoT(EOSMATERIAL *material, double rho, double T){
    double Gamma = 0.0;
    switch(material->matType)
    {
        case EOSIDEALGAS:
            // Gas has no shear modulus
            Gamma = 0.0;
            break;
#ifdef HAVE_TILLOTSON_H
        case EOSTILLOTSON:
            EOS_NOT_IMPLEMENTED(material);   /* no strength model for Tillotson yet */
            break;
#endif
#ifdef HAVE_ANEOSMATERIAL_H
        case EOSANEOS:
            Gamma = ANEOSGammaofRhoT(material->ANEOSmaterial, rho, T);
            break;
#endif
#ifdef HAVE_REOS3_H
        case EOSREOS3:
            // Gas has no shear modulus
            Gamma = 0.0;
            break;
#endif
#ifdef HAVE_SCVHEOS_H
        case EOSSCVHEOS:
            // Gas has no shear modulus
            Gamma = 0.0;
            break;
#endif
        default: EOS_UNKNOWN_MATERIAL(material);
    }
    return Gamma;
}

/*
 * Calculate the cold part of the internal energy, u(rho, T=0). For the
 * tabulated EOS this is u at the lowest tabulated temperature, i.e. the first
 * table row; where a table's cold rows are defective (u(T_min) above u at
 * higher T, the M-ANEOS tables below 18-58 K) it is not the energy minimum:
 * states slightly below it are inside the table (EOSIsInTable,
 * EOSisbelowColdCurve, which use the table status, not this value).
 *
 * Required: internal, ballic
 */
double EOSUCold(EOSMATERIAL *material, double rho)
{
    double ucold = 0;

    switch(material->matType)
    {
        case EOSIDEALGAS:
            // ucold = 0
            break;
#ifdef HAVE_TILLOTSON_H
        case EOSTILLOTSON:
            ucold = tillColdULookup(material->tillmaterial, rho);
            break;
#endif
#ifdef HAVE_ANEOSMATERIAL_H
        case EOSANEOS:
            // Lowest tabulated temperature: ANEOS cannot be evaluated below T_min.
            ucold = ANEOSUCold(material->ANEOSmaterial, rho);
            break;
#endif
#ifdef HAVE_REOS3_H
        case EOSREOS3:
            // Lowest tabulated temperature: the library cannot be evaluated below T_min.
            ucold = reos3UCold(material->reos3material, rho);
            break;
#endif
#ifdef HAVE_SCVHEOS_H
        case EOSSCVHEOS:
            // Lowest tabulated temperature: the library cannot be evaluated below T_min.
            ucold = scvheosUCold(material->scvheosmaterial, rho);
            break;
#endif
        default: EOS_UNKNOWN_MATERIAL(material);
    }

    return ucold;
}

/*
 * Evaluate u(rho,T) and P(rho,u) for a material and check that the state is
 * valid (inside the EOS table, finite). Returns EOS_SUCCESS or EOS_FAIL.
 */
static int EOSSolveBCEval(EOSMATERIAL *material, double rho, double T,
                          double *pu, double *pP)
{
    double u, P;

    if (!(rho > 0.0) || !isfinite(rho)) return EOS_FAIL;

    u = EOSUofRhoT(material, rho, T);
    if (EOSFailed(u) || !isfinite(u)) return EOS_FAIL;

    if (!EOSIsUsable(material, rho, u)) return EOS_FAIL;

    /* EOSFailed, not only isfinite: the sentinel -1e50 is finite (NaN is caught by both). */
    P = EOSPofRhoU(material, rho, u);
    if (EOSFailed(P) || !isfinite(P)) return EOS_FAIL;

    *pu = u;
    *pP = P;
    return EOS_SUCCESS;
}

/*
 * Move the bracket end *px of EOSSolveBC geometrically until P2 >= P (dir > 0, upper bracket,
 * factor EOS_SOLVEBC_EXPAND_FAC) or P2 <= P (dir < 0, lower bracket, factor 1/EOS_SOLVEBC_EXPAND_FAC),
 * keeping (*pu, *pP) the values at *px. A trial density where material 2 cannot be evaluated (outside
 * its table, or u(rho, T) not usable) does not end the search: the factor is replaced by its square
 * root, so the search closes in on the edge of the usable range and still brackets a root between
 * the last valid density and that edge. It fails when the factor is within EOS_SOLVEBC_TOL of 1 (the
 * edge is reached and P2 there is still on the wrong side of P) or after EOS_SOLVEBC_MAX_EXPAND
 * trials. The trials beyond the edge are lookups outside the table of material 2 and are counted in
 * its nLookupFail.
 *
 * Up to 1.2.0 the brackets moved by 1.1 and 0.9 per step (100 steps: a factor 1.4e4 up, 3.7e-5 down)
 * and the first trial outside the table ended the search.
 */
static int EOSSolveBCExpand(EOSMATERIAL *material, double T, double P, int dir,
                            double *px, double *pu, double *pP)
{
    double fac = (dir > 0) ? EOS_SOLVEBC_EXPAND_FAC : 1.0/EOS_SOLVEBC_EXPAND_FAC;
    int i;

    for (i = 0; (dir > 0) ? (*pP < P) : (*pP > P); i++) {
        double x, u, Px;

        if (i >= EOS_SOLVEBC_MAX_EXPAND) return EOS_FAIL;
        x = *px*fac;
        if (EOSSolveBCEval(material, x, T, &u, &Px) != EOS_SUCCESS) {
            fac = sqrt(fac);
            if (fabs(fac - 1.0) < EOS_SOLVEBC_TOL) return EOS_FAIL;
            continue;
        }
        *px = x;
        *pu = u;
        *pP = Px;
    }
    return EOS_SUCCESS;
}

/*
 * Given rho1, u1 (material 1) solve for rho2, u2 (material 2) at the interface
 * between two materials.
 *
 * The b.c. are:
 *
 * P1(rho1,u1) = P2(rho2,u2) and T1(rho1,u1) - Tjump = T2(rho2,u2)
 *
 * Since P = P(rho, T) and T2 is known we solve P2(rho2) - P1 = 0 by bisection,
 * assuming P2(rho, T2) is monotonically increasing in rho. The tolerances and
 * the discontinuity check are described in EOSlib.h.
 *
 * Returns EOS_SUCCESS, or EOS_FAIL (message under EOSLIB_VERBOSE) if not,
 * also at once if (rho1, u1) is not usable in material 1 (EOSIsUsable).
 * On failure *prho2 and *pu2 are untouched.
 */
int EOSSolveBC(EOSMATERIAL *material1, EOSMATERIAL *material2, double rho1, double u1,
               double *prho2, double *pu2, double Tjump)
{
    double P, T;
    double a, ua, Pa;   /* upper bracket: Pa >= P */
    double b, ub, Pb;   /* lower bracket: Pb <= P */
    double c, uc, Pc;
    double slope0 = 0.0; /* secant slope of P2(rho) once the bracket is EOS_SOLVEBC_SLOPE_WIDTH wide */
    int bHaveSlope = 0;
    int bHaveC = 0;      /* (c, uc) evaluated by the bisection */
    int i;

    /* State in material 1. */
    P = EOSPofRhoU(material1, rho1, u1);
    T = EOSTofRhoU(material1, rho1, u1) - Tjump;
    /*
     * A failed P (EOS_FAIL_VALUE, rho1 outside the table of material 1) is finite, so isfinite
     * alone let it through as the target pressure: every P2 > -1e50, the lower bracket shrank
     * EOS_SOLVEBC_MAX_EXPAND times (or until it left the table of material 2) and the solver
     * failed there with a message about material 2. EOSFailed also catches NaN.
     *
     * A state that T(rho, u) cannot invert (u above the column or below the cold curve of a
     * tabulated EOS, below the cold curve of Tillotson or the ideal gas) has a finite P and T,
     * but they are those of the fallback (rho1, T_min), not of (rho1, u1): the solver matched
     * the wrong target and reported success. EOSIsUsable is exactly the test whether P, c and T
     * are those of the state. It is checked after the lookups, which count the failure or the
     * fallback in nLookupFail (EOSIsUsable, a status query, counts nothing).
     *
     * T2 = T1 - Tjump = 0 is a valid target: the cold curve of Tillotson and the ideal gas
     * (T_min = 0), e.g. a cold compressed interface (it was rejected with T <= 0 up to 1.2.0).
     * A tabulated material 2 has T_min > 0, and EOSSolveBCEval fails cleanly there, because
     * u(rho, T2) of material 2 fails below its table. Only T2 < 0 is invalid input.
     */
    if (EOSFailed(P) || !isfinite(P) || !isfinite(T) || T < 0.0 ||
        EOSIsUsable(material1, rho1, u1) != EOS_TRUE) {
        EOS_WARN("EOSSolveBC: invalid state in material 1 "
                "(rho1= %g u1= %g P= %g T= %g, usable= %d).\n", rho1, u1, P, T, EOSIsUsable(material1, rho1, u1));
        return EOS_FAIL;
    }

    /*
     * Upper bracket: start at max(rho1, rho0 of material 2) and expand until P2 >= P. rho1 alone
     * (the denser component inside) failed at once where rho1 lies below the table of material 2:
     * a gas envelope (SCvH, REOS3 at 1e-6 g/cc) on a condensed core (ANEOS), whose interface
     * density is near the core's rho0, not near rho1. rho0 of material 2 is inside its table.
     */
    a = fmax(rho1, material2->rho0);
    if (EOSSolveBCEval(material2, a, T, &ua, &Pa) != EOS_SUCCESS) {
        EOS_WARN("EOSSolveBC: rho= %g T= %g not valid in material %d.\n",
                a, T, material2->iMat);
        return EOS_FAIL;
    }
    if (EOSSolveBCExpand(material2, T, P, +1, &a, &ua, &Pa) != EOS_SUCCESS) {
        EOS_WARN("EOSSolveBC: upper bracket not found in material %d "
                "(P= %g T= %g a= %g Pa= %g).\n", material2->iMat, P, T, a, Pa);
        return EOS_FAIL;
    }

    /*
     * Lower bracket: start just below the reference density, but never above
     * the upper bracket, and shrink until P2 <= P.
     */
    b = fmin(0.999*material2->rho0, 0.999*a);
    if (EOSSolveBCEval(material2, b, T, &ub, &Pb) != EOS_SUCCESS) {
        EOS_WARN("EOSSolveBC: rho= %g T= %g not valid in material %d.\n",
                b, T, material2->iMat);
        return EOS_FAIL;
    }
    if (EOSSolveBCExpand(material2, T, P, -1, &b, &ub, &Pb) != EOS_SUCCESS) {
        EOS_WARN("EOSSolveBC: lower bracket not found in material %d "
                "(P= %g T= %g b= %g Pb= %g).\n", material2->iMat, P, T, b, Pb);
        return EOS_FAIL;
    }

    /* The root must now be bracketed by (b, a). */
    if (!(b < a) || Pa < P || Pb > P) {
        EOS_WARN("EOSSolveBC: root not bracketed: P= %g T= %g "
                "a= %g ua= %g Pa= %g b= %g ub= %g Pb= %g\n",
                P, T, a, ua, Pa, b, ub, Pb);
        return EOS_FAIL;
    }

    /* Exact hits at the bracket ends (also covers the degenerate case a == b). */
    if (Pa == P) { *prho2 = a; *pu2 = ua; return EOS_SUCCESS; }
    if (Pb == P) { *prho2 = b; *pu2 = ub; return EOS_SUCCESS; }

    /*
     * Bisection. c is always set at least once because b < a here. The secant slope
     * (Pa - Pb)/(a - b) is recorded when the bracket first becomes narrower than
     * EOS_SOLVEBC_SLOPE_WIDTH; it is the reference for the discontinuity test below.
     */
    c = 0.5*(a + b);
    uc = ub;
    if (2.0*(a - b)/(a + b) <= EOS_SOLVEBC_SLOPE_WIDTH) {
        slope0 = (Pa - Pb)/(a - b);
        bHaveSlope = 1;
    }
    for (i = 0; 2.0*(a - b)/(a + b) > EOS_SOLVEBC_TOL; i++) {
        if (i >= EOS_SOLVEBC_MAX_BISECT) {
            EOS_WARN("EOSSolveBC: bisection did not converge "
                    "(a= %g b= %g).\n", a, b);
            return EOS_FAIL;
        }
        c = 0.5*(a + b);
        if (EOSSolveBCEval(material2, c, T, &uc, &Pc) != EOS_SUCCESS) {
            EOS_WARN("EOSSolveBC: invalid state during bisection "
                    "(rho= %g T= %g).\n", c, T);
            return EOS_FAIL;
        }
        bHaveC = 1;
        if (Pc < P) {
            b = c; ub = uc; Pb = Pc;
        } else {
            a = c; ua = uc; Pa = Pc;
        }
        if (!bHaveSlope && 2.0*(a - b)/(a + b) <= EOS_SOLVEBC_SLOPE_WIDTH) {
            slope0 = (Pa - Pb)/(a - b);
            bHaveSlope = 1;
        }
    }

    /*
     * Discontinuity test. On a continuous P2(rho, T) the pressure bracket closes
     * with the density bracket, Pa - Pb ~ (dP/drho) (a - b), so the secant slope
     * stays at the local dP/drho. Where P2 jumps by dP at the root (T(rho, u)
     * switches between roots of a non-monotonic u(T) column, e.g. at a melt
     * transition of a defective table), no density matches P and the slope grows
     * like dP/(a - b): by EOS_SOLVEBC_SLOPE_WIDTH/EOS_SOLVEBC_TOL = 1e4 between
     * the reference and the end. Fail if it grew by more than EOS_SOLVEBC_SLOPE_FAC.
     *
     * A test on the pressure gap alone would reject continuous stiff states: at
     * the end Pa - Pb ~ (dlnP/dlnrho) EOS_SOLVEBC_TOL P, which exceeds
     * EOS_SOLVEBC_PTOL P for K/P > 1e4 (rock below ~10 MPa). A gap below
     * EOS_SOLVEBC_PTOL P is accepted regardless of the slope, so jumps that small
     * pass as before, and roundoff in Pa - Pb at low K/P cannot cause a failure.
     * A root within EOS_SOLVEBC_SLOPE_WIDTH of a kink where dP/drho increases by
     * more than EOS_SOLVEBC_SLOPE_FAC fails this test; the one-sided slopes
     * below decide then.
     */
    if (Pa - Pb > EOS_SOLVEBC_PTOL * fabs(P) &&
        !(bHaveSlope && (Pa - Pb)/(a - b) <= EOS_SOLVEBC_SLOPE_FAC*slope0)) {
        /*
         * Second opinion for a root next to a kink. The reference slope is taken where the bracket
         * is EOS_SOLVEBC_SLOPE_WIDTH wide; if that bracket straddles a kink (P2 cut off at 0 on
         * one side, the stiff cold curve on the other, as Tillotson at rho0 for T -> 0) it is an
         * average over both sides and the final slope can exceed it by more than
         * EOS_SOLVEBC_SLOPE_FAC on a continuous P2. The one-sided secant slopes just outside the
         * final bracket, over EOS_SOLVEBC_PROBE times its width, tell the two apart: on a
         * continuous P2 the slope across the bracket is at most about the larger of them, at a
         * jump it is jump/(a - b), larger by orders of magnitude. Two more evaluations, only in
         * this case; if a probe cannot be evaluated (table edge) the first verdict stands.
         */
        const double h = EOS_SOLVEBC_PROBE*(a - b);
        double uL, PL, uR, PR;
        int bContinuous = 0;
        if (EOSSolveBCEval(material2, b - h, T, &uL, &PL) == EOS_SUCCESS &&
            EOSSolveBCEval(material2, a + h, T, &uR, &PR) == EOS_SUCCESS) {
            double sOut = fmax((Pb - PL)/h, (PR - Pa)/h);
            bContinuous = (sOut > 0.0 && (Pa - Pb)/(a - b) <= EOS_SOLVEBC_SLOPE_FAC*sOut);
        }
        if (!bContinuous) {
            EOS_WARN("EOSSolveBC: P2(rho, T) is discontinuous at the root "
                    "(P= %g T= %g rho= %g Pb= %g Pa= %g slope= %g slope0= %g).\n",
                    P, T, c, Pb, Pa, (Pa - Pb)/(a - b), slope0);
            return EOS_FAIL;
        }
    }

    /*
     * The returned (rho2, u2) must be a consistent pair. After the bisection (c, uc) are the values
     * of its last evaluation at c, so nothing is left to do (up to 1.2.0 that lookup was repeated
     * here, one evaluation of material 2 per call for nothing). Only if the loop never ran is uc
     * still ub, from another density: evaluate at c then. That cannot happen with the brackets
     * above, which start at least 0.1% apart (b <= 0.999 a), far wider than EOS_SOLVEBC_TOL; the
     * guard keeps the result consistent should the bracketing change.
     */
    if (!bHaveC && EOSSolveBCEval(material2, c, T, &uc, &Pc) != EOS_SUCCESS) return EOS_FAIL;

    *prho2 = c;
    *pu2 = uc;
    return EOS_SUCCESS;
}

/*
 * WOOLFSON_MIN_PRESSURE_CGS (a physical pressure) in the code units of the material.
 */
double EOSWoolfsonMinPressure(EOSMATERIAL *material)
{
    return WOOLFSON_MIN_PRESSURE_CGS/(material->dErgPerGmUnit*material->dGmPerCcUnit);
}

/*
 * Calculate the coefficient
 *
 * f_ij := rho_mat1(P, T)/rho_mat2(P, T)
 *
 * needed to correct the density at a material interface in gasoline.
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

    if (P < EOSWoolfsonMinPressure(material1)) {
        EOS_WARN("EOSWoolfsonCoeff: P= %g is close to zero, no density correction is done.\n", P);
        return 1.0;
    }

    rho1 = EOSRhoofPT(material1, P, T);
    rho2 = EOSRhoofPT(material2, P, T);

    /*
     * Unphysical density (rho <= 0, which includes the failure sentinel of
     * the inversions). This is deliberately NOT mapped to f_ij = 1: a failed
     * inversion means (P, T) is outside the range where the interface
     * correction is defined, and silently skipping the correction would hide
     * that. Return EOS_FAIL_VALUE so the caller can detect the failure and
     * abort or handle the particle explicitly. Callers must check for
     * f_ij < 0 (EOSFailed).
     */
    if ((rho1 <= 0.0) || (rho2 <= 0.0))
    {
        EOS_WARN("EOSWoolfsonCoeff: unphysical density (rho1= %g, rho2= %g) at P= %g T= %g.\n", rho1, rho2, P, T);
        return EOS_FAIL_VALUE;
    }

    /*
     * Limit the correction in expanded states: at phase transitions P(rho) is
     * flat (not strictly monotonic), so rho(P, T) is ill-conditioned and f_ij
     * can become huge.
     */
    if ((rho1 <= 0.8*material1->rho0) || (rho2 <= 0.8*material2->rho0)) {
        return 1.0;
    }

    return (rho1/rho2);
}