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

#ifndef EOSLIB_HINCLUDED
#define EOSLIB_HINCLUDED

#include <stdio.h>   /* FILE in EOSPrintMat */
 
#ifdef HAVE_TILLOTSON_H
    #include <tillotson.h>
#define TILL_VERSION_MAJOR_REQUIRED 3
#define TILL_VERSION_MINOR_REQUIRED 4
#define TILL_VERSION_PATCH_REQUIRED 4
#endif
#ifdef HAVE_ANEOSMATERIAL_H
    #include <ANEOSmaterial.h>
#define ANEOS_VERSION_MAJOR_REQUIRED 1
#define ANEOS_VERSION_MINOR_REQUIRED 0
#define ANEOS_VERSION_PATCH_REQUIRED 4
#endif
#ifdef HAVE_REOS3_H
    #include <reos3.h>
#define REOS3_VERSION_MAJOR_REQUIRED 1
#define REOS3_VERSION_MINOR_REQUIRED 0
#define REOS3_VERSION_PATCH_REQUIRED 1
#endif
#ifdef HAVE_SCVHEOS_H
    #include <scvheos.h>
#define SCVHEOS_VERSION_MAJOR_REQUIRED 1
#define SCVHEOS_VERSION_MINOR_REQUIRED 0
#define SCVHEOS_VERSION_PATCH_REQUIRED 1
#endif

#include "igeos.h"
#define IGEOS_VERSION_MAJOR_REQUIRED 1
#define IGEOS_VERSION_MINOR_REQUIRED 0
#define IGEOS_VERSION_PATCH_REQUIRED 2

#define EOS_VERSION_TEXT    "1.0.6"
#define EOS_VERSION_MAJOR   1
#define EOS_VERSION_MINOR   0
#define EOS_VERSION_PATCH   6

#define WOOLFSON_MIN_PRESSURE 1e-3

#define EOS_N_MATERIAL_MAX 130

#define EOSIDEALGAS  0
#define EOSTILLOTSON 1
#define EOSANEOS     2
#define EOSREOS3     3
#define EOSSCVHEOS   4

#define MAT_IDEALGAS 0
#define MAT_TILLOTSON_MIN 1
#define MAT_TILLOTSON_MAX 50
#define MAT_ANEOS_MIN 51
#define MAT_ANEOS_MAX 100
#define MAT_REOS3_MIN 101
#define MAT_REOS3_MAX 109
#define MAT_SCVHEOS_MIN 110
#define MAT_SCVHEOS_MAX 130

#define EOS_TRUE 1
#define EOS_FALSE 0

/*
 * Status codes. Note two conventions coexist in this library:
 *  - EOSIsInTable/EOSIsUsable use EOS_SUCCESS (1) / EOS_FAIL (0) / EOS_OUTSIDE_*.
 *  - EOSSolveBC follows the C convention: 0 on success, -1 on failure.
 */
#define EOS_SUCCESS 1
#define EOS_FAIL 0
#define EOS_OUTSIDE_RHOMIN -1
#define EOS_OUTSIDE_RHOMAX -2
#define EOS_OUTSIDE_VMIN -3
#define EOS_OUTSIDE_VMAX -4
 
// Pass additional parameters for some EOS
struct igeosParam
{ 
    /* Constant for ideal gas EOS. */
    double dConstGamma;
    double dMeanMolMass;	
};

typedef struct EOSmaterial
{
	int iMat; // Material number
	int matType; // Material type, e.g., 0: ideal gas, 1: Tillotson, 2: ANEOS
	double rho0; // reference density
	double rho0forInterfaceCorrection; // reference density for interface correction
	int bEntropyTableInit; // flag to signal if the entropy table is initialized
    int bEntropy; // flag to signal if the eos has entropy available (see EOSSofRhoU)
	double minSoundSpeed; // sound speed at reference values
    char MatString[256];
    IGEOSMAT *igeosmaterial; // Pointer to ideal gas material
#ifdef HAVE_TILLOTSON_H
	TILLMATERIAL *tillmaterial; // Pointer to tillotson material
#endif
#ifdef HAVE_ANEOSMATERIAL_H
	ANEOSMATERIAL *ANEOSmaterial; // Pointer to ANEOS material
#endif
#ifdef HAVE_REOS3_H
	REOS3MAT *reos3material; // Pointer to REOS3 material
#endif
#ifdef HAVE_SCVHEOS_H
	SCVHEOSMAT *scvheosmaterial; // Pointer to SCvH EOS material
#endif

} EOSMATERIAL;

#ifdef __cplusplus
extern "C" {
#endif
// Initialization and finalization

EOSMATERIAL *EOSinitMaterial(int iMat, double dKpcUnit, double dMsolUnit, const void * additional_data);
void EOSinitIsentropicLookup(EOSMATERIAL *material, const void * additional_data);
void EOSfinalizeMaterial(EOSMATERIAL *material);
void EOSPrintMat(EOSMATERIAL *material, FILE *fp);
int EOSCheckVersions();

// Access functions

// Main functions
double EOSPofRhoU(EOSMATERIAL *material, double rho, double u);
double EOSPofRhoT(EOSMATERIAL *material, double rho, double T);
double EOSCofRhoU(EOSMATERIAL *material, double rho, double u);
double EOSPCofRhoU(EOSMATERIAL *material, double rho, double u, double *c);
double EOSPCTofRhoU(EOSMATERIAL *material, double rho, double u, double *c, double *T);
double EOSIsentropic(EOSMATERIAL *material, double rho1, double u1, double rho2);

/*
 * Entropy functions.
 *
 * NOTE: Unlike all other quantities, the specific entropy is NOT converted to
 * code units. It is always returned in cgs (erg g^-1 K^-1), and the additive
 * constant is library-dependent, so only entropy differences within the same
 * material are meaningful. Both functions return -1e50 if the material has no
 * entropy available (bEntropy == EOS_FALSE) or the lookup fails.
 */
double EOSSofRhoU(EOSMATERIAL *material, double rho, double u);
double EOSSofRhoT(EOSMATERIAL *material, double rho, double T);

int EOSYieldModel(EOSMATERIAL *material);
double EOSYieldStrength(EOSMATERIAL *material, double rho, double u);
/*
 * Same as EOSYieldStrength, but takes T and P the caller already has (e.g.
 * from EOSPCTofRhoU) instead of re-deriving them from (rho, u).
 */
double EOSYieldStrengthofRhoTP(EOSMATERIAL *material, double rho, double T, double P);

// Inverse functions
double EOSTofRhoU(EOSMATERIAL *material, double rho, double u);
double EOSUofRhoT(EOSMATERIAL *material, double rho, double T);
double EOSRhoofPT(EOSMATERIAL *material, double p, double T);

double EOSRhoofUT(EOSMATERIAL *material, double u, double T);

int EOSisbelowColdCurve(EOSMATERIAL *material, double rho, double u);

/*
 * EOSIsInTable reports the raw table status of (rho, u) for the material:
 * EOS_SUCCESS, or one of the EOS_OUTSIDE_* codes, or EOS_FAIL. It never
 * applies a tolerance policy.
 *
 * EOSIsUsable applies the library policy: returns EOS_TRUE if P(rho, u),
 * c(rho, u) and T(rho, u) can be evaluated for this material even if the
 * state lies outside the lookup table (e.g. u > u_max for an analytic EOS such
 * as Tillotson), EOS_FALSE otherwise. Use this in the hydro code; use
 * EOSIsInTable when you need to know *why* a state is out of range or whether
 * isentropic lookups (EOSIsentropic) are available.
 */
int EOSIsInTable(EOSMATERIAL *material, double rho, double u);
int EOSIsUsable(EOSMATERIAL *material, double rho, double u);

// Derivatives
double EOSdPdRho(EOSMATERIAL *material, double rho, double u);
double EOSdPdU(EOSMATERIAL *material, double rho, double u);
double EOSdUdRho(EOSMATERIAL *material, double rho, double u);

double EOSdPdRhoofRhoT(EOSMATERIAL *material, double rho, double T);
double EOSdPdT(EOSMATERIAL *material, double rho, double T);

double EOSGammaofRhoT(EOSMATERIAL *material, double rho, double T);

double EOSUCold(EOSMATERIAL *material, double rho);

// Boundary condition solver (returns 0 on success, -1 on failure)
int EOSSolveBC(EOSMATERIAL *material1, EOSMATERIAL *material2, double rho1, double u1,
 double *prho2, double *pu2, double Tjump);
 
// Woolfson correction f_ij = rho_1(P,T)/rho_2(P,T). Returns -1.0 if the
// density inversion fails; callers must check for a negative result.
double EOSWoolfsonCoeff(EOSMATERIAL *material1, EOSMATERIAL *material2, double P, double T);
#ifdef __cplusplus
}
#endif
#endif
