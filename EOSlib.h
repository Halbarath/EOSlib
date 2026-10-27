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
 * Version history
 *   1.0.6  requires ANEOSmaterial 1.0.4; EOSSolveBC returns -1 where
 *          P2(rho, T) is discontinuous at the root instead of returning the
 *          location of the jump; solver tolerances exported from this header.
 *   1.0.7  requires scvheos 1.1.0: EOSPCofRhoU/EOSPCTofRhoU use the fused
 *          scvheosPCTofRhoU; EOSUCold no longer reaches into the SCvH
 *          tables; MatString set for SCvH materials; EOSSofRhoU returns
 *          -1e50 for SCvH when T(rho, u) fails; EOSIsInTable returns
 *          EOS_FAIL instead of a raw sublibrary code for REOS3 and SCvH.
 *   1.0.9  requires ANEOSmaterial 1.0.5, scvheos 1.2.0 and reos3 1.1.0, which
 *          share the interpBilinear core (closed table intervals with a
 *          roundoff tolerance at the edges, no GSL). EOSIsInTable treats
 *          the ANEOS table as closed with the same tolerance, so a state on
 *          the edge that the sublibrary evaluates is reported as inside.
 *   1.2.0  one failure contract across EOSlib and all sublibraries (see
 *          "Failure contract" below): requires tillotson 4.0.0, ANEOSmaterial
 *          1.2.0, reos3 1.3.0, scvheos 1.3.0, igeos 1.1.0. API changes:
 *          EOS_SUCCESS is 0 and EOS_FAIL is -99 (were 1 and 0), so test
 *          "== EOS_SUCCESS", never the truth value; EOSSolveBC returns
 *          EOS_SUCCESS/EOS_FAIL (0/-99, "!= 0" still means failure);
 *          EOSWoolfsonCoeff returns EOS_FAIL_VALUE (-1e50) instead of -1.0
 *          (still negative); EOSTofRhoU never fails (T_min of the table,
 *          0 for Tillotson); EOSYieldStrength returns EOS_FAIL_VALUE for an
 *          invalid T or T_melt instead of exiting; per-call diagnostics only
 *          with -DEOSLIB_VERBOSE; no GSL, no brent.o. The Tillotson isentrope
 *          table is the library's standard table (tillInitLookupDefault:
 *          1e-4..100 rho0, 1 K..1e7 K at rho0, 1200 x 1200), built once in
 *          EOSinitMaterial, the same physical table in every unit system
 *          (was 1e-4..200 and v_max = 1200 in code units, built twice).
 *          EOSSolveBC detects a discontinuous P2(rho, T) by the growth of the
 *          secant slope as the bracket closes (EOS_SOLVEBC_SLOPE_WIDTH,
 *          EOS_SOLVEBC_SLOPE_FAC) instead of the pressure gap alone, which
 *          rejected continuous stiff states (K/P > 1e4, rock below ~10 MPa).
 */

#ifndef EOSLIB_HINCLUDED
#define EOSLIB_HINCLUDED

#include <stdio.h>   /* FILE in EOSPrintMat */

#ifdef HAVE_TILLOTSON_H
    #include <tillotson.h>
#define TILL_VERSION_MAJOR_REQUIRED 4
#define TILL_VERSION_MINOR_REQUIRED 0
#define TILL_VERSION_PATCH_REQUIRED 0
#endif
#ifdef HAVE_ANEOSMATERIAL_H
    #include <ANEOSmaterial.h>
#define ANEOS_VERSION_MAJOR_REQUIRED 1
#define ANEOS_VERSION_MINOR_REQUIRED 2
#define ANEOS_VERSION_PATCH_REQUIRED 0
#endif
#ifdef HAVE_REOS3_H
    #include <reos3.h>
#define REOS3_VERSION_MAJOR_REQUIRED 1
#define REOS3_VERSION_MINOR_REQUIRED 3
#define REOS3_VERSION_PATCH_REQUIRED 0
#endif
#ifdef HAVE_SCVHEOS_H
    #include <scvheos.h>
#define SCVHEOS_VERSION_MAJOR_REQUIRED 1
#define SCVHEOS_VERSION_MINOR_REQUIRED 3
#define SCVHEOS_VERSION_PATCH_REQUIRED 0
#endif

#include "igeos.h"
#define IGEOS_VERSION_MAJOR_REQUIRED 1
#define IGEOS_VERSION_MINOR_REQUIRED 1
#define IGEOS_VERSION_PATCH_REQUIRED 0

#define EOS_VERSION_TEXT    "1.2.0"
#define EOS_VERSION_MAJOR   1
#define EOS_VERSION_MINOR   2
#define EOS_VERSION_PATCH   0

#define WOOLFSON_MIN_PRESSURE 1e-3

/* EOSSolveBC: bracket expansion limits and convergence tolerances. */
#define EOS_SOLVEBC_MAX_EXPAND  100   /* 1.1^100 ~ 1.4e4, 0.9^100 ~ 2.7e-5 */
#define EOS_SOLVEBC_MAX_BISECT  200
#define EOS_SOLVEBC_TOL         1e-10 /* relative, in density */
#define EOS_SOLVEBC_PTOL        1e-6  /* relative, in pressure: a gap below this always passes */
#define EOS_SOLVEBC_SLOPE_WIDTH 1e-6  /* relative bracket width of the reference secant slope */
#define EOS_SOLVEBC_SLOPE_FAC   10.0  /* the final secant slope may exceed the reference by this */

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
 * Status codes (the same values in every sublibrary): EOS_SUCCESS is 0,
 * EOS_FAIL a generic failure, EOS_OUTSIDE_* the reasons EOSIsInTable gives.
 * Always compare with == EOS_SUCCESS; the truth value of a status is not
 * meaningful (EOS_SUCCESS was 1 before 1.2.0).
 */
#define EOS_SUCCESS 0
#define EOS_FAIL  -99
#define EOS_OUTSIDE_RHOMIN -1
#define EOS_OUTSIDE_RHOMAX -2
#define EOS_OUTSIDE_VMIN -3
#define EOS_OUTSIDE_VMAX -4

/*
 * Failure contract (EOSlib and all sublibraries: igeos, tillotson,
 * ANEOSmaterial, reos3, scvheos).
 *
 * A failed evaluation never exits. Every function returning a double
 * returns exactly EOS_FAIL_VALUE (-1e50), in any unit system, when the
 * state is outside the table, an inversion has no solution, or the quantity
 * is not available for the material (e.g. the entropy); the fused lookups
 * set their output arguments to it as well. The one exception is T(rho, u):
 * EOSTofRhoU (and *T of EOSPCTofRhoU) never fails and returns the lowest
 * temperature of the table (T_min of the tabulated EOS, 0 for Tillotson)
 * when (rho, u) cannot be inverted; the hydro codes rely on this. Status
 * functions return EOS_FAIL or an EOS_OUTSIDE_* code.
 *
 * Every failure and every fallback is counted in the sublibrary's
 * nLookupFail (printed by EOSPrintMat) and reported on stderr only when the
 * libraries are compiled with -DEOSLIB_VERBOSE (EOS_WARN and the *_WARN of
 * each sublibrary). Test a returned value with EOSFailed(x), which is true
 * for the sentinel and for NaN.
 *
 * Fatal (message on stderr, then assert(0); exit(1)): initialisation errors
 * (unknown material, missing or defective table, version mismatch) and
 * calling a function that is not implemented for the material (EOSRhoofUT
 * for all but ANEOS, EOSIsentropic for REOS3, EOSdPdRhoofRhoT/EOSdPdT for the
 * analytic EOS, EOSGammaofRhoT/EOSYieldStrength for Tillotson); those are
 * programming errors, not lookup failures. The assert(0) comes first on
 * purpose: in an MPI run a plain exit() on one rank can leave the others
 * waiting in a collective and the job hangs, while SIGABRT makes the
 * launcher kill the whole job (the exit(1) is for -DNDEBUG builds). All
 * sublibraries do the same.
 */
#define EOS_FAIL_VALUE (-1e50)

static inline int EOSFailed(double x) {
    return !(x > -1e40);
}

/*
 * Diagnostic messages of the lookups and solvers, only with
 * -DEOSLIB_VERBOSE. The silent variant still type-checks the arguments but
 * compiles to nothing.
 */
#ifdef EOSLIB_VERBOSE
#define EOS_WARN(...) fprintf(stderr, __VA_ARGS__)
#else
#define EOS_WARN(...) do { if (0) fprintf(stderr, __VA_ARGS__); } while (0)
#endif

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
int EOSCheckVersions(void);   /* EOS_TRUE if every sublibrary is recent enough (fatal in EOSinitMaterial otherwise) */

// Access functions

// Main functions
double EOSPofRhoU(EOSMATERIAL *material, double rho, double u);
double EOSPofRhoT(EOSMATERIAL *material, double rho, double T);
double EOSCofRhoU(EOSMATERIAL *material, double rho, double u);
double EOSPCofRhoU(EOSMATERIAL *material, double rho, double u, double *c);
double EOSPCTofRhoU(EOSMATERIAL *material, double rho, double u, double *c, double *T);
/*
 * u2 on the isentrope through (rho1, u1) at rho2. If (rho1, u1) is below the
 * cold curve (the lowest tabulated temperature), all libraries follow the
 * coldest isentrope, i.e. return the cold curve at rho2 (counted in
 * nLookupFail); EOS_FAIL_VALUE only for an invalid input or rho2 outside a
 * tabulated EOS. Tillotson has an analytic P and integrates the ODE directly
 * where (rho1, u1) or rho2 is outside its table (counted as a fallback).
 * Not implemented (fatal) for REOS3.
 */
double EOSIsentropic(EOSMATERIAL *material, double rho1, double u1, double rho2);

/*
 * Entropy functions.
 *
 * NOTE: Unlike all other quantities, the specific entropy is NOT converted to
 * code units. It is always returned in cgs (erg g^-1 K^-1), and the additive
 * constant is library-dependent, so only entropy differences within the same
 * material are meaningful. Both functions return EOS_FAIL_VALUE if the
 * material has no entropy available (bEntropy == EOS_FALSE) or the lookup
 * fails.
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

/*
 * Boundary condition solver: given (rho1, u1) in material 1, find (rho2, u2)
 * in material 2 with equal pressure and T2 = T1 - Tjump. Converges to
 * EOS_SOLVEBC_TOL in density; the pressure residual is then at most
 * max(EOS_SOLVEBC_PTOL |P|, EOS_SOLVEBC_SLOPE_FAC (dP2/drho) EOS_SOLVEBC_TOL rho2),
 * i.e. ~ (dlnP/dlnrho) EOS_SOLVEBC_TOL relative for stiff, low-pressure states.
 * Returns EOS_SUCCESS (0), or EOS_FAIL on failure (no bracket, state outside the
 * table, or P2(rho, T2) discontinuous at the root, which happens where u(T)
 * of material 2 is not monotonic, e.g. at a melt transition of a defective
 * table), with a message under EOSLIB_VERBOSE. A discontinuity is detected by
 * the growth of the secant slope of P2(rho) as the bracket closes (see
 * EOSSolveBC in EOSlib.c), not by the pressure gap, which for a continuous but
 * stiff P2 can exceed EOS_SOLVEBC_PTOL |P|. On failure *prho2 and *pu2 are
 * untouched.
 */
int EOSSolveBC(EOSMATERIAL *material1, EOSMATERIAL *material2, double rho1, double u1,
 double *prho2, double *pu2, double Tjump);

/*
 * Woolfson correction f_ij = rho_1(P,T)/rho_2(P,T). Returns EOS_FAIL_VALUE
 * if a density inversion fails; callers must check for a negative result
 * (EOSFailed or f < 0).
 */
double EOSWoolfsonCoeff(EOSMATERIAL *material1, EOSMATERIAL *material2, double P, double T);
#ifdef __cplusplus
}
#endif
#endif