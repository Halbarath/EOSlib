/*
 * Header file for igeos.c.
 */
#ifndef IGEOS_HINCLUDED
#define IGEOS_HINCLUDED

#include <stdio.h>   /* FILE in igeosPrintMat */

/*
 * Version history
 *   1.1.0  EOSlib failure contract: IGEOS_FAIL = -99, IGEOS_FAIL_VALUE and
 *          igeosFailed(); igeosRhoofUT is fatal with a message (it is
 *          undefined for an ideal gas, not a lookup failure). igeosInitMat
 *          rejects an inconsistent unit pair (fatal): both units > 0 is code
 *          units, both <= 0 is cgs (exactly one <= 0 divided by zero).
 *          nLookupFail and IGEOS_WARN like the other libraries: every failure and
 *          fallback is counted once per call (see the failure contract below) and
 *          printed by igeosPrintMat.
 *          The cold curve is u = 0 (T_min = 0): a state below it (u < 0 or NaN)
 *          falls back to (rho, T = 0) like a state below the cold curve of a
 *          tabulated EOS, counted once: T = 0, P = 0, c = 0, dP/drho = 0, isentrope
 *          u2 = 0 (u < 0 gave T < 0, P < 0 and c = NaN). P(rho, T) and u(rho, T) return
 *          IGEOS_FAIL_VALUE for T < 0, rho(P, T) for T <= 0 or P < 0 (inf or a negative
 *          density before), the isentrope for rho1 or rho2 <= 0; all counted.
 *          New igeosPCTofRhoU (EOSPCTofRhoU counted a state twice through
 *          igeosPCofRhoU and igeosTofRhoU).
 *          igeosPrintMat prints rho0 and cv in cgs with the unit and the unit pair,
 *          like the tabulated libraries (it printed code-unit values unlabelled).
 *          New igeosdPdRhoofRhoT and igeosdPdTofRhoT (the analytic (rho, T) derivatives, for
 *          EOSdPdRhoofRhoT and EOSdPdT, which were fatal for the ideal gas).
 *          igeosInitMat rejects (fatal) dConstGamma <= 1 and dMeanMolMass <= 0 (or not finite): cv
 *          was inf, negative or NaN and every lookup returned garbage.
 *          An invalid density (rho <= 0, NaN, inf) fails like a density outside a table: P, c and
 *          the derivatives return IGEOS_FAIL_VALUE, T(rho, u) falls back to T = 0, all counted once
 *          (P was negative or 0 for rho <= 0, uncounted). New igeosUCold (0, or IGEOS_FAIL_VALUE for
 *          an invalid rho), for EOSUCold and EOSisbelowColdCurve.
 */
#define IGEOS_VERSION_TEXT    "1.1.0"
#define IGEOS_VERSION_MAJOR   1
#define IGEOS_VERSION_MINOR   1
#define IGEOS_VERSION_PATCH   0

/*
 * Status codes (EOSlib convention: 0 is success, -99 a generic failure,
 * -1..-4 are reserved for the OUTSIDE_* codes of the tabulated libraries).
 */
#define IGEOS_SUCCESS 0
#define IGEOS_FAIL  -99

/*
 * Failure contract (the EOSlib convention, shared by all EOS libraries).
 * The ideal gas is analytic; its cold curve is u = 0 and T_min = 0. A state
 * (rho, u) below the cold curve (u < 0, or NaN) falls back to (rho, T = 0),
 * like a state below the cold curve of a tabulated EOS: T = 0, P = 0, c = 0
 * (T(rho, u) never fails). An invalid density (rho <= 0, NaN, inf) is like a
 * density outside a table: every function that takes one returns
 * IGEOS_FAIL_VALUE (P, c, the derivatives, u_cold, the isentrope for rho1 or
 * rho2), except T(rho, u), which falls back to T = 0. IGEOS_FAIL_VALUE is also
 * returned by P(rho, T), u(rho, T) and the (rho, T) derivatives for T < 0 (below
 * T_min), by rho(P, T) for T <= 0 or P < 0, and by igeosSofRhoT/igeosSofRhoU for
 * T <= 0, where the entropy is undefined (all also for NaN). Every failure
 * and fallback is counted in Mat->nLookupFail, once per call, and reported on
 * stderr only with -DEOSLIB_VERBOSE (IGEOS_WARN). Test a returned value with
 * igeosFailed(x). Unimplemented functions (igeosRhoofUT) are fatal.
 */
#define IGEOS_FAIL_VALUE (-1e50)

static inline int igeosFailed(double x) {
    return !(x > -1e40);
}

/*
 * Diagnostic messages, only with -DEOSLIB_VERBOSE. The silent variant still
 * type-checks the arguments but compiles to nothing.
 */
#ifdef EOSLIB_VERBOSE
#define IGEOS_WARN(...) fprintf(stderr, __VA_ARGS__)
#else
#define IGEOS_WARN(...) do { if (0) fprintf(stderr, __VA_ARGS__); } while (0)
#endif

typedef struct igeosMat {
    /* Units. */
    double dKpcUnit;
    double dMsolUnit;
    double dGasConst;      /* k_B/m_H in code units (specific energy per K) */
    double dErgPerGmUnit;
    double dGmPerCcUnit;
    double dSecUnit;

    int iMat;

    /* Constant for ideal gas EOS. */
    double dConstGamma;
    double dMeanMolMass;

    /* The reference density required by ballic. */
    double rho0;

    /* Specific heat capacity. */
    double cv;

    /* Diagnostics: number of failed evaluations and fallbacks since init (see the failure contract above). */
    long nLookupFail;
} IGEOSMAT;

IGEOSMAT *igeosInitMat(int iMat, double dConstGamma, double dMeanMolMass, double dKpcUnit, double dMsolUnit);

int igeosFinalizeMat(IGEOSMAT *Mat);

void igeosPrintMat(IGEOSMAT *Mat, FILE *fp);

double igeosPofRhoU(IGEOSMAT *Mat, double rho, double u);
double igeosPofRhoT(IGEOSMAT *Mat, double rho, double T);
double igeosCofRhoU(IGEOSMAT *Mat, double rho, double u);
double igeosPCofRhoU(IGEOSMAT *Mat, double rho, double u, double *c);
double igeosPCTofRhoU(IGEOSMAT *Mat, double rho, double u, double *c, double *T);
double igeosIsentropicU(IGEOSMAT *Mat, double rho1, double u1, double rho2);
double igeosTofRhoU(IGEOSMAT *Mat, double rho, double u);
double igeosUofRhoT(IGEOSMAT *Mat, double rho, double T);
double igeosRhoofPT(IGEOSMAT *Mat, double P, double T);
double igeosRhoofUT(IGEOSMAT *Mat, double U, double T);
double igeosdPdRho(IGEOSMAT *Mat, double rho, double u);
double igeosdPdU(IGEOSMAT *Mat, double rho, double u);
/* dP/drho|T = (gamma-1) cv T and dP/dT|rho = (gamma-1) rho cv; IGEOS_FAIL_VALUE (counted) for T < 0 or NaN, like igeosPofRhoT. */
double igeosdPdRhoofRhoT(IGEOSMAT *Mat, double rho, double T);
double igeosdPdTofRhoT(IGEOSMAT *Mat, double rho, double T);
/* u(rho, T = 0) = 0, the cold curve; IGEOS_FAIL_VALUE (counted) for rho <= 0, NaN or inf. */
double igeosUCold(IGEOSMAT *Mat, double rho);
double igeosSofRhoU(IGEOSMAT *Mat, double rho, double u);
double igeosSofRhoT(IGEOSMAT *Mat, double rho, double T);

#endif