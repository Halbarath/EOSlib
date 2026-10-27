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
 * Test suite for EOSlib.
 *
 * Tests the wrapper layer: dispatch to the sublibraries, the ideal-gas EOS
 * against its closed form, the ANEOS path against the ANEOSmaterial library
 * directly, the table-status policy, the yield model, and the boundary
 * condition solver. Everything is an invariant; no reference values are
 * hard-coded. Usage:
 *
 *   testEOSmaterial [iMatA iMatB] [dKpcUnit dMsolUnit] [iMatSCVH] [iMatREOS3]
 *
 * Defaults: M-ANEOS iron (63) and forsterite (62), Earth-like code units,
 * SCvH 114 and REOS3 104. iMatA/iMatB must be ANEOS materials whose tables
 * are in the working directory; the Tillotson tests (dispatch, table status
 * and policy, interface solver between iron and granite) run only if
 * HAVE_TILLOTSON_H is defined and need no table files. With HAVE_SCVHEOS_H / HAVE_REOS3_H and a
 * non-zero iMatSCVH / iMatREOS3 the SCvH and REOS3 paths are tested too
 * (0 skips them), including the boundary condition solver for a gas envelope
 * on the denser ANEOS material; the gas tables must be in the working
 * directory. Exit status is the number of failed checks.
 *
 * The failure contract of the wrapper (EOS_FAIL_VALUE, the T_min fallback of
 * EOSTofRhoU, EOSSolveBC/EOSWoolfsonCoeff return codes) is tested per
 * material type; with -DEOSLIB_VERBOSE the libraries print a diagnostic for
 * every failed lookup, so messages are expected there.
 */
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "EOSlib.h"

/* ---------------------------------------------------------------- harness */

static int nChecks = 0;
static int nFails = 0;
static const char *currentTest = "";

#define CHECK(cond, ...) do { \
    nChecks++; \
    if (!(cond)) { \
        nFails++; \
        fprintf(stderr, "FAIL [%s] %s:%d: ", currentTest, __FILE__, __LINE__); \
        fprintf(stderr, __VA_ARGS__); \
        fprintf(stderr, "\n"); \
    } \
} while (0)

#define CHECK_REL(a, b, tol, ...) do { \
    double _a = (a), _b = (b); \
    double _den = fmax(fabs(_a), fabs(_b)); \
    double _err = (_den > 0.0) ? fabs(_a - _b) / _den : fabs(_a - _b); \
    nChecks++; \
    if (!(_err <= (tol))) { \
        nFails++; \
        fprintf(stderr, "FAIL [%s] %s:%d: ", currentTest, __FILE__, __LINE__); \
        fprintf(stderr, __VA_ARGS__); \
        fprintf(stderr, " (got %.15e vs %.15e, rel err %.3e > %.1e)\n", _a, _b, _err, (double)(tol)); \
    } \
} while (0)

/*
 * Pressure residual EOSSolveBC guarantees at a solved state (rho2, T2) of material e
 * (see EOSlib.h): EOS_SOLVEBC_PTOL |P1|, or EOS_SOLVEBC_SLOPE_FAC times the local
 * slope across the final density bracket (EOS_SOLVEBC_TOL rho2 wide). The slope is the
 * larger one-sided difference over the reference width, along the path the solver
 * uses (u(rho, T), then P(rho, u)); a factor 2 covers the midpoint and roundoff.
 */
static double solveBCPResidualBound(EOSMATERIAL *e, double rho2, double T2, double P1)
{
    double h = EOS_SOLVEBC_SLOPE_WIDTH * rho2;
    double P0 = EOSPofRhoU(e, rho2, EOSUofRhoT(e, rho2, T2));
    double Pp = EOSPofRhoU(e, rho2 + h, EOSUofRhoT(e, rho2 + h, T2));
    double Pm = EOSPofRhoU(e, rho2 - h, EOSUofRhoT(e, rho2 - h, T2));
    double slope = fmax(fabs(Pp - P0), fabs(P0 - Pm)) / h;
    if (!isfinite(slope)) slope = 0.0;
    return EOS_SOLVEBC_PTOL * fabs(P1) + 2.0 * EOS_SOLVEBC_SLOPE_FAC * slope * EOS_SOLVEBC_TOL * rho2;
}

#define TEST(name) do { currentTest = name; printf("  %-44s", name); fflush(stdout); } while (0)
#define END_TEST() do { printf("%s\n", nFails == failsBefore ? "ok" : "FAILED"); } while (0)

static int isSentinel(double v) { return EOSFailed(v); }

/* An energy strictly below the cold curve at rho, whatever the sign of
 * u_cold (the M-ANEOS tables have u_cold < 0 around the reference density). */
static double belowColdCurve(EOSMATERIAL *e, double rho)
{
    double ucold = EOSUCold(e, rho);
    return ucold - fabs(ucold) - 1.0;
}

/* ---------------------------------------------------------- ideal gas EOS */

static void testIdealGas(double dKpcUnit, double dMsolUnit)
{
    int failsBefore = nFails;
    TEST("ideal gas against closed form");
    struct igeosParam par = { 1.4, 1.0 };
    EOSMATERIAL *g = EOSinitMaterial(MAT_IDEALGAS, dKpcUnit, dMsolUnit, &par);
    double gm1 = par.dConstGamma - 1.0;

    CHECK(g->matType == EOSIDEALGAS, "matType");
    CHECK(g->bEntropy && g->bEntropyTableInit, "entropy flags");
    CHECK(EOSYieldModel(g) == 0, "yield model for gas");
    CHECK(EOSYieldStrength(g, 1.0, 1.0) == 0.0, "yield strength for gas");
    CHECK(EOSGammaofRhoT(g, 1.0, 300.0) == 0.0, "shear modulus for gas");

    double rhos[] = { 1e-6, 1e-3, 1.0, 50.0 };
    double us[]   = { 1e-4, 1e-1, 1.0, 1e3 };
    for (int a = 0; a < 4; a++) for (int b = 0; b < 4; b++) {
        double rho = rhos[a], u = us[b];
        double c, T;
        double P = EOSPCTofRhoU(g, rho, u, &c, &T);
        CHECK_REL(P, gm1 * rho * u, 1e-14, "P rho=%g u=%g", rho, u);
        CHECK_REL(c, sqrt(par.dConstGamma * gm1 * u), 1e-14, "c rho=%g u=%g", rho, u);
        CHECK_REL(EOSPofRhoU(g, rho, u), P, 1e-14, "PofRhoU vs PCT");
        CHECK_REL(EOSCofRhoU(g, rho, u), c, 1e-14, "CofRhoU vs PCT");
        CHECK_REL(EOSTofRhoU(g, rho, u), T, 1e-14, "TofRhoU vs PCT");
        CHECK_REL(EOSUofRhoT(g, rho, T), u, 1e-13, "u(rho,T(rho,u)) roundtrip");
        CHECK_REL(EOSPofRhoT(g, rho, T), P, 1e-13, "P(rho,T) vs P(rho,u)");
        CHECK_REL(EOSRhoofPT(g, P, T), rho, 1e-13, "rho(P,T) roundtrip");
        CHECK_REL(EOSdPdRho(g, rho, u), gm1 * u, 1e-14, "dP/drho");
        CHECK_REL(EOSdPdU(g, rho, u), gm1 * rho, 1e-14, "dP/du");
        CHECK_REL(EOSdUdRho(g, rho, u), P / (rho * rho), 1e-14, "du/drho");
        /* Isentrope: u2/u1 = (rho2/rho1)^(gamma-1), and entropy is conserved. */
        double rho2 = 2.0 * rho;
        double u2 = EOSIsentropic(g, rho, u, rho2);
        CHECK_REL(u2, u * pow(2.0, gm1), 1e-14, "isentropic u");
        CHECK_REL(EOSSofRhoU(g, rho2, u2), EOSSofRhoU(g, rho, u), 1e-12, "entropy conserved");
        CHECK(EOSIsInTable(g, rho, u) == EOS_SUCCESS, "ideal gas always in table");
        CHECK(EOSIsUsable(g, rho, u) == EOS_TRUE, "ideal gas always usable");
        CHECK(EOSUCold(g, rho) == 0.0, "ucold = 0 for gas");
    }
    EOSfinalizeMaterial(g);
    END_TEST();
}

/* ---------------------------------------------------------------- ANEOS */

#ifdef HAVE_ANEOSMATERIAL_H

/* 100 for a quick run (seconds); 1000 for a full sweep (minutes per material). */
#define NSAMPLE 1000

static double sampleRho(ANEOSMATERIAL *m, int k)
{
    double lo = log(m->rhoAxis[0]), hi = log(m->rhoAxis[m->nRho-1]);
    double f = 0.05 + 0.90 * (k + 0.5) / NSAMPLE;
    return exp(lo + f * (hi - lo)) * m->CGStoCodeUnitsforRho;
}

/* Log-spaced in T, from 300 K (above the defective cold rows of the M-ANEOS
 * tables) to the top of the table. */
static double sampleT(ANEOSMATERIAL *m, int k)
{
    double lo = log(fmax(m->TAxis[0], 300.0)), hi = log(m->TAxis[m->nT-1]);
    double f = 0.05 + 0.90 * (k + 0.5) / NSAMPLE;
    return exp(lo + f * (hi - lo));
}

static void testANEOSDispatch(EOSMATERIAL *e)
{
    int failsBefore = nFails;
    TEST("ANEOS dispatch matches ANEOSmaterial directly");
    ANEOSMATERIAL *m = e->ANEOSmaterial;
    CHECK(e->matType == EOSANEOS, "matType");
    CHECK(e->bEntropy && e->bEntropyTableInit, "entropy flags");
    CHECK_REL(e->rho0, ANEOSgetRho0(m), 1e-15, "rho0");
    CHECK(e->minSoundSpeed > 0.0, "minSoundSpeed = %g", e->minSoundSpeed);
    CHECK(strlen(e->MatString) > 0, "MatString empty");

    int nSkipped = 0;
    for (int a = 0; a < NSAMPLE; a++) for (int b = 0; b < NSAMPLE; b++) {
        double rho = sampleRho(m, a), T = sampleT(m, b);
        double u = EOSUofRhoT(e, rho, T);
        /* EOSUCold is the first table row. Where a defective cold row lies
         * above u(rho,T) (M-ANEOS forsterite above ~55 g/cc), EOSUCold,
         * EOSisbelowColdCurve and EOSIsInTable are wrong by construction. */
        if (u < EOSUCold(e, rho)) { nSkipped++; continue; }
        double c, Tf;
        double P = EOSPCTofRhoU(e, rho, u, &c, &Tf);
        CHECK_REL(u, ANEOSUofRhoT(m, rho, T), 0.0, "UofRhoT dispatch");
        /* The separate calls recompute the cell coordinate from T, so P and c
         * differ from the fused lookup by amplified roundoff. */
        CHECK_REL(P, ANEOSPofRhoU(m, rho, u), 1e-12, "PCTofRhoU P");
        CHECK_REL(c, ANEOSCofRhoU(m, rho, u), 1e-12, "PCTofRhoU c");
        /* T(rho,u) may return another root where u(T) is not monotonic (melt
         * transition, a table defect tested in ANEOSmaterial), so only require
         * that the returned T reproduces u. */
        CHECK_REL(EOSUofRhoT(e, rho, Tf), u, 1e-10, "u(rho, T(rho,u)) at rho=%g T=%g", rho, T);
        CHECK_REL(EOSPofRhoU(e, rho, u), P, 1e-12, "PofRhoU");
        CHECK_REL(EOSPofRhoT(e, rho, T), ANEOSPofRhoT(m, rho, T), 0.0, "PofRhoT");
        CHECK_REL(EOSCofRhoU(e, rho, u), c, 1e-12, "CofRhoU");
        CHECK_REL(EOSTofRhoU(e, rho, u), Tf, 0.0, "TofRhoU");
        double c2; double P2 = EOSPCofRhoU(e, rho, u, &c2);
        CHECK_REL(P2, P, 0.0, "PCofRhoU P"); CHECK_REL(c2, c, 0.0, "PCofRhoU c");
        CHECK_REL(EOSSofRhoU(e, rho, u), ANEOSSofRhoU(m, rho, u), 0.0, "SofRhoU");
        CHECK_REL(EOSSofRhoT(e, rho, T), ANEOSSofRhoT(m, rho, T), 0.0, "SofRhoT");
        CHECK_REL(EOSRhoofPT(e, P, T), ANEOSRhoofPT(m, P, T), 0.0, "RhoofPT");
        CHECK_REL(EOSdPdRho(e, rho, u), ANEOSdPdRhoofRhoU(m, rho, u), 0.0, "dPdRho");
        CHECK_REL(EOSdPdU(e, rho, u), ANEOSdPdUofRhoU(m, rho, u), 0.0, "dPdU");
        CHECK_REL(EOSdUdRho(e, rho, u), ANEOSdUdRhoofRhoU(m, rho, u), 0.0, "dUdRho");
        CHECK_REL(EOSdPdRhoofRhoT(e, rho, T), ANEOSdPdRhoofRhoT(m, rho, T), 0.0, "dPdRhoofRhoT");
        CHECK_REL(EOSdPdT(e, rho, T), ANEOSdPdTofRhoT(m, rho, T), 0.0, "dPdT");
        CHECK_REL(EOSUCold(e, rho), ANEOSUofRhoT(m, rho, m->TAxis[0]), 0.0, "UCold");
        CHECK(EOSisbelowColdCurve(e, rho, belowColdCurve(e, rho)) == 1, "below cold curve");
        CHECK(EOSisbelowColdCurve(e, rho, u) == 0, "above cold curve (u at T=%g)", T);
        CHECK(EOSIsInTable(e, rho, u) == EOS_SUCCESS, "in table at rho=%g T=%g", rho, T);
        CHECK(EOSIsUsable(e, rho, u) == EOS_TRUE, "usable at rho=%g T=%g", rho, T);
    }
    /* Entropy is in cgs: compare against the raw table at a node. */
    {
        int i = m->nRho / 2, j = m->nT / 2;
        double rho = m->rhoAxis[i] * m->CGStoCodeUnitsforRho;
        CHECK_REL(EOSSofRhoT(e, rho, m->TAxis[j]), m->sArray[i][j], 1e-13, "entropy returned in cgs");
    }
    CHECK(m->nLookupFail == 0, "nLookupFail = %ld after in-table dispatch", m->nLookupFail);
    if (nSkipped)
        printf("\n    note: %d of %d states skipped (u(rho,T) below the first table row, table defect)\n  %-44s",
               nSkipped, NSAMPLE * NSAMPLE, "");
    END_TEST();
}

static void testANEOSTableStatus(EOSMATERIAL *e)
{
    int failsBefore = nFails;
    TEST("EOSIsInTable / EOSIsUsable codes");
    ANEOSMATERIAL *m = e->ANEOSmaterial;
    double rhoOk = sampleRho(m, NSAMPLE / 2);
    double uOk = EOSUofRhoT(e, rhoOk, sampleT(m, NSAMPLE / 2));
    double rhoLo = 0.5 * m->rhoAxis[0] * m->CGStoCodeUnitsforRho;
    double rhoHi = 2.0 * m->rhoAxis[m->nRho-1] * m->CGStoCodeUnitsforRho;
    /* The table is closed at both ends with the relative tolerance
     * ANEOS_EDGE_TOL (1e-13): the nodes themselves are inside, even after the
     * code<->cgs round trip, and 1e-12 beyond them is outside. */
    double rhoMin = m->rhoAxis[0] * m->CGStoCodeUnitsforRho;
    double rhoMax = m->rhoAxis[m->nRho-1] * m->CGStoCodeUnitsforRho;
    double rhoEdge = m->rhoAxis[m->nRho-1] * (1.0 + 1e-12) * m->CGStoCodeUnitsforRho;
    double TMid = sampleT(m, NSAMPLE / 2);

    CHECK(EOSIsInTable(e, rhoLo, uOk) == EOS_OUTSIDE_RHOMIN, "below rhoMin");
    CHECK(EOSIsInTable(e, rhoHi, uOk) == EOS_OUTSIDE_RHOMAX, "above rhoMax");
    CHECK(EOSIsInTable(e, rhoEdge, uOk) == EOS_OUTSIDE_RHOMAX, "just above rhoMax");
    /* The edge columns are tested only where their first table row really is
     * the energy minimum: in a defective column (M-ANEOS forsterite above
     * ~55 g/cc) EOSUCold lies above u(rho, T) and EOSIsInTable reports
     * OUTSIDE_VMIN by construction, see EOSUCold. */
    {
        double uEdge = EOSUofRhoT(e, rhoMin, TMid);
        if (uEdge >= EOSUCold(e, rhoMin)) CHECK(EOSIsInTable(e, rhoMin, uEdge) == EOS_SUCCESS, "on rhoMin (closed interval)");
        else printf("\n    note: rhoMin column skipped (first table row above u(rho,T), table defect)\n  %-44s", "");
        uEdge = EOSUofRhoT(e, rhoMax, TMid);
        if (uEdge >= EOSUCold(e, rhoMax)) CHECK(EOSIsInTable(e, rhoMax, uEdge) == EOS_SUCCESS, "on rhoMax (closed interval)");
        else printf("\n    note: rhoMax column skipped (first table row above u(rho,T), table defect)\n  %-44s", "");
    }
    CHECK(EOSIsInTable(e, rhoOk, EOSUofRhoT(e, rhoOk, m->TAxis[0])) == EOS_SUCCESS, "on the cold curve");
    CHECK(EOSIsInTable(e, rhoOk, EOSUofRhoT(e, rhoOk, m->TAxis[m->nT-1])) == EOS_SUCCESS, "on the top row");
    CHECK(EOSIsInTable(e, -1.0, uOk) == EOS_OUTSIDE_RHOMIN, "negative rho");
    CHECK(EOSIsInTable(e, NAN, uOk) == EOS_OUTSIDE_RHOMIN, "NaN rho");
    CHECK(EOSIsInTable(e, rhoOk, NAN) == EOS_OUTSIDE_VMIN, "NaN u");
    CHECK(EOSIsInTable(e, rhoOk, belowColdCurve(e, rhoOk)) == EOS_OUTSIDE_VMIN, "below cold curve");
    CHECK(EOSIsInTable(e, rhoOk, 1e6 * uOk) == EOS_OUTSIDE_VMAX, "huge u");
    CHECK(EOSIsUsable(e, rhoLo, uOk) == EOS_FALSE, "tabulated EOS unusable below rhoMin");
    CHECK(EOSIsUsable(e, rhoOk, 1e6 * uOk) == EOS_FALSE, "tabulated EOS unusable above u_max");
    /* Consistency between the status and the actual lookup: anything EOSIsInTable
     * accepts must not produce a sentinel. */
    for (int a = 0; a < NSAMPLE; a++) for (int b = 0; b < NSAMPLE; b++) {
        double rho = sampleRho(m, a), u = EOSUofRhoT(e, rho, sampleT(m, b));
        if (EOSIsInTable(e, rho, u) == EOS_SUCCESS) {
            double c, T;
            double P = EOSPCTofRhoU(e, rho, u, &c, &T);
            CHECK(!isSentinel(P) && !isSentinel(c), "in-table state produced a sentinel at rho=%g", rho);
        }
    }
    END_TEST();
}

static void testANEOSYield(EOSMATERIAL *e)
{
    int failsBefore = nFails;
    TEST("yield strength model");
    ANEOSMATERIAL *m = e->ANEOSmaterial;
    int model = EOSYieldModel(e);
    CHECK(model == m->yieldStrengthModel, "EOSYieldModel = %d vs %d", model, m->yieldStrengthModel);

    if (model <= 0) {
        double rho = sampleRho(m, NSAMPLE / 2), u = EOSUofRhoT(e, rho, sampleT(m, 1));
        CHECK(EOSYieldStrength(e, rho, u) == 0.0, "no strength model but Y != 0");
        END_TEST();
        return;
    }

    double Y0, YM, mui, mud, xi;
    ANEOSYieldParameters(m, &Y0, &YM, &mui, &mud, &xi);

    for (int a = 1; a < NSAMPLE - 1; a++) {
        double rho = sampleRho(m, a);
        double Tcold = sampleT(m, 0);
        double Tm = ANEOSTmeltofRho(m, rho);
        double u = EOSUofRhoT(e, rho, Tcold);
        double c, T;
        double P = EOSPCTofRhoU(e, rho, u, &c, &T);

        double Y = EOSYieldStrength(e, rho, u);
        double Ytp = EOSYieldStrengthofRhoTP(e, rho, T, P);
        CHECK(Y >= 0.0 && isfinite(Y), "Y = %g at rho=%g", Y, rho);
        CHECK_REL(Y, Ytp, 1e-13, "Y(rho,u) vs Y(rho,T,P) at rho=%g", rho);

        /* Tension: strength at P < 0 equals strength at P = 0. */
        if (Tcold < Tm) {
            double Yzero = EOSYieldStrengthofRhoTP(e, rho, T, 0.0);
            double Yneg = EOSYieldStrengthofRhoTP(e, rho, T, -10.0 * fabs(P) - 1.0);
            CHECK_REL(Yneg, Yzero, 1e-14, "tension clamp at rho=%g", rho);
            /* Intact model: Y(P=0) = Y0 * thermal factor, and Y grows with P. */
            if (model == 1) {
                double thermal = fmax(tanh(xi * (Tm / T - 1.0)), 0.0);
                CHECK_REL(Yzero, Y0 * thermal, 1e-13, "Y(P=0) = Y0 * f(T) at rho=%g", rho);
                CHECK(EOSYieldStrengthofRhoTP(e, rho, T, 2.0 * fabs(P) + 1.0) >= Yzero, "Y increases with P");
            }
            /* Damaged model: Y <= intact Y and Y(P=0) = 0. */
            if (model == 2) CHECK(Yzero == 0.0, "damaged Y(P=0) = %g", Yzero);
        }
        /* At and above the melt temperature the strength vanishes. Needs a melt
         * temperature inside the table; the melt curve is 0 where the material
         * has no condensed phase (low-density vapour region). */
        if (Tm > m->TAxis[0] && Tm < m->TAxis[m->nT-1]) {
            double Tabove = fmin(1.5 * Tm, 0.99 * m->TAxis[m->nT-1]);
            double Pm = EOSPofRhoT(e, rho, Tabove);
            double Ym = EOSYieldStrengthofRhoTP(e, rho, Tabove, Pm);
            CHECK(Ym == 0.0, "Y above melt = %g at rho=%g", Ym, rho);
        }
        CHECK(EOSGammaofRhoT(e, rho, T) > 0.0, "shear modulus at rho=%g", rho);
    }
    END_TEST();
}

static void testSolveBC(EOSMATERIAL *inner, EOSMATERIAL *outer)
{
    int failsBefore = nFails;
    TEST("EOSSolveBC: P and T continuity at an interface");
    ANEOSMATERIAL *m = outer->ANEOSmaterial;
    int nSolved = 0, nSkipped = 0;
    for (int a = 2; a < NSAMPLE - 2; a++) for (int b = 1; b < NSAMPLE - 3; b++) {
        /* State in the outer (less dense) material. The solver is for condensed
         * interfaces (it brackets from rho1 upwards), so expanded or vapour
         * states are not meaningful inputs, and neither are pressures near
         * zero: below its saturation pressure a Maxwell-constructed table has
         * no state at all, so the same floor as for the Woolfson correction
         * applies. */
        double rho1 = sampleRho(m, a), T1 = sampleT(m, b);
        if (rho1 < 0.8 * outer->rho0) continue;
        double u1 = EOSUofRhoT(outer, rho1, T1);
        double P1 = EOSPofRhoU(outer, rho1, u1);
        if (P1 <= EOSWoolfsonMinPressure(outer)) continue;
        /* The solver derives its own T1 from (rho1, u1). Where T(rho,u) does
         * not reproduce T1 the column is defective (another root in broken
         * cold rows or a melt transition) and the solver is handed a
         * meaningless target; those states are skipped. */
        double T1s = EOSTofRhoU(outer, rho1, u1);
        if (fabs(T1s - T1) > 1e-8 * T1) { nSkipped++; continue; }
        double rho2, u2;
        for (int jump = 0; jump < 2; jump++) {
            double Tjump = jump ? 0.05 * T1 : 0.0;
            if (EOSSolveBC(outer, inner, rho1, u1, &rho2, &u2, Tjump) != 0) continue;
            nSolved++;
            /* The pressure residual the solver guarantees (see solveBCPResidualBound). */
            CHECK(fabs(EOSPofRhoU(inner, rho2, u2) - P1) <= solveBCPResidualBound(inner, rho2, T1 - Tjump, P1),
                  "P continuity rho1=%g T1=%g jump=%g: %.15e vs %.15e", rho1, T1, Tjump, EOSPofRhoU(inner, rho2, u2), P1);
            /* T continuity is checked through the forward lookup, since
             * T(rho2,u2) in the inner material may pick another root. */
            CHECK_REL(u2, EOSUofRhoT(inner, rho2, T1 - Tjump), 1e-10, "T continuity rho1=%g T1=%g jump=%g", rho1, T1, Tjump);
            CHECK(EOSIsUsable(inner, rho2, u2) == EOS_TRUE, "solved state not usable");
        }
    }
    CHECK(nSolved > 0, "EOSSolveBC never succeeded");
    /* Failure path: a negative target temperature. The library prints a
     * diagnostic to stderr for this; that message is expected. */
    {
        double rho1 = sampleRho(m, NSAMPLE - 1);
        double u1 = EOSUofRhoT(outer, rho1, sampleT(m, NSAMPLE - 1));
        double rho2 = -1.0, u2 = -1.0;
        int ret = EOSSolveBC(outer, inner, rho1, u1, &rho2, &u2, 10.0 * sampleT(m, NSAMPLE - 1));
        CHECK(ret == EOS_SUCCESS || ret == EOS_FAIL, "SolveBC: bad return %d", ret);
        if (ret == EOS_SUCCESS) CHECK(isfinite(rho2) && rho2 > 0.0, "SolveBC: success with rho2=%g", rho2);
        if (ret == EOS_FAIL) CHECK(rho2 == -1.0 && u2 == -1.0, "outputs touched on failure");
    }
    if (nSkipped) printf("\n    note: %d interface states skipped (T(rho,u) does not reproduce T, table defect)\n  %-44s", nSkipped, "");
    END_TEST();
}

static void testWoolfson(EOSMATERIAL *a, EOSMATERIAL *b)
{
    int failsBefore = nFails;
    TEST("EOSWoolfsonCoeff contract");
    ANEOSMATERIAL *m = a->ANEOSmaterial;
    double rho = sampleRho(m, NSAMPLE / 2), T = sampleT(m, NSAMPLE / 2);
    double P = EOSPofRhoT(a, rho, T);
    CHECK(EOSWoolfsonCoeff(a, a, P, T) == 1.0, "same material");
    CHECK(EOSWoolfsonCoeff(a, b, 0.5 * EOSWoolfsonMinPressure(a), T) == 1.0, "below minimum pressure");
    double f = EOSWoolfsonCoeff(a, b, P, T);
    CHECK(EOSFailed(f) || (f > 0.0 && isfinite(f)), "f_ij = %g", f);
    CHECK(EOSFailed(EOSWoolfsonCoeff(a, b, 1e6 * EOSPofRhoT(a, m->rhoAxis[m->nRho-1] * m->CGStoCodeUnitsforRho, T), T)), "f_ij above the table must be the sentinel");
    if (f > 0.0 && f != 1.0) {
        double r1 = EOSRhoofPT(a, P, T), r2 = EOSRhoofPT(b, P, T);
        CHECK_REL(f, r1 / r2, 1e-13, "f_ij = rho1/rho2");
    }
    END_TEST();
}

#endif /* HAVE_ANEOSMATERIAL_H */

/* ------------------------------------------------------------------ SCvH */

#ifdef HAVE_SCVHEOS_H

#define NSAMPLE_SCVH 300

static double scvhExp10(double x) { return exp(x * 2.302585092994045684); }

/* Log-spaced samples in the interior of the SCvH table, code units / Kelvin. */
static double scvhSampleRho(SCVHEOSMAT *s, int k)
{
    double f = 0.05 + 0.90 * (k + 0.5137) / NSAMPLE_SCVH;
    return scvhExp10(s->LogRhoMin + f * (s->LogRhoMax - s->LogRhoMin));
}

static double scvhSampleT(SCVHEOSMAT *s, int k)
{
    double f = 0.05 + 0.90 * (k + 0.5137) / NSAMPLE_SCVH;
    return scvhExp10(s->LogTMin + f * (s->LogTMax - s->LogTMin));
}

/* 1 if the table A increases along T in the rho column containing rho (both
 * node columns). Where it does not, T(rho, u) or T(rho, s) may return another
 * root (table defect), and the identities through them are not invariants. */
static int scvhMonotonic(SCVHEOSMAT *s, const double *A, double rho)
{
    double lr = log10(rho);
    int j = 0;
    while (j < s->nRho - 2 && s->dLogRhoAxis[j+1] <= lr) j++;
    for (int i = 0; i < s->nT - 1; i++) {
        if (!(A[(size_t)j * s->nT + i + 1] > A[(size_t)j * s->nT + i])) return 0;
        if (!(A[(size_t)(j+1) * s->nT + i + 1] > A[(size_t)(j+1) * s->nT + i])) return 0;
    }
    return 1;
}

static int scvhEnergyMonotonic(SCVHEOSMAT *s, double rho)  { return scvhMonotonic(s, s->dLogUArray, rho); }
static int scvhEntropyMonotonic(SCVHEOSMAT *s, double rho) { return scvhMonotonic(s, s->dLogSArray, rho); }

static void testSCVHDispatch(EOSMATERIAL *e)
{
    int failsBefore = nFails;
    TEST("SCvH dispatch matches scvheos directly");
    SCVHEOSMAT *s = e->scvheosmaterial;
    CHECK(e->matType == EOSSCVHEOS, "matType");
    CHECK(e->bEntropy && e->bEntropyTableInit, "entropy flags");
    CHECK_REL(e->rho0, s->rho0, 0.0, "rho0");
    CHECK(e->rho0forInterfaceCorrection == 0.0, "gas has no interface-correction cutoff");
    CHECK(strlen(e->MatString) > 0, "MatString empty");
    CHECK(EOSYieldModel(e) == 0, "yield model for gas");
    CHECK(EOSYieldStrength(e, scvhSampleRho(s, 1), 1.0) == 0.0, "yield strength for gas");
    CHECK(EOSGammaofRhoT(e, scvhSampleRho(s, 1), scvhSampleT(s, 1)) == 0.0, "shear modulus for gas");

    int nSkipped = 0, nSkippedS = 0;
    for (int a = 0; a < NSAMPLE_SCVH; a++) for (int b = 0; b < NSAMPLE_SCVH; b++) {
        double rho = scvhSampleRho(s, a), T = scvhSampleT(s, b);
        if (!scvhEnergyMonotonic(s, rho)) { nSkipped++; continue; }
        double u = EOSUofRhoT(e, rho, T);
        double c, Tf;
        double P = EOSPCTofRhoU(e, rho, u, &c, &Tf);
        CHECK_REL(u, scvheosUofRhoT(s, rho, T), 0.0, "UofRhoT dispatch");
        /* The inversion is closed-form on the bilinear interpolant, so the
         * fused lookup reproduces the (rho, T) lookups to roundoff. */
        CHECK_REL(P, scvheosPofRhoT(s, rho, T), 1e-12, "PCTofRhoU P at rho=%g T=%g", rho, T);
        CHECK_REL(c, scvheosCsofRhoT(s, rho, T), 1e-12, "PCTofRhoU c at rho=%g T=%g", rho, T);
        CHECK_REL(Tf, T, 1e-12, "PCTofRhoU T at rho=%g T=%g", rho, T);
        CHECK_REL(EOSPofRhoU(e, rho, u), P, 0.0, "PofRhoU");
        CHECK_REL(EOSCofRhoU(e, rho, u), c, 0.0, "CofRhoU");
        CHECK_REL(EOSTofRhoU(e, rho, u), scvheosTofRhoU(s, rho, u), 0.0, "TofRhoU");
        CHECK_REL(EOSPofRhoT(e, rho, T), scvheosPofRhoT(s, rho, T), 0.0, "PofRhoT");
        double c2; double P2 = EOSPCofRhoU(e, rho, u, &c2);
        CHECK_REL(P2, P, 0.0, "PCofRhoU P"); CHECK_REL(c2, c, 0.0, "PCofRhoU c");
        CHECK_REL(EOSSofRhoT(e, rho, T), scvheosSofRhoT(s, rho, T), 0.0, "SofRhoT");
        CHECK_REL(EOSSofRhoU(e, rho, u), scvheosSofRhoT(s, rho, T), 1e-12, "SofRhoU");
        CHECK_REL(EOSdPdRho(e, rho, u), scvheosdPdRhoofRhoU(s, rho, u), 0.0, "dPdRho");
        CHECK_REL(EOSdPdU(e, rho, u), scvheosdPdUofRhoU(s, rho, u), 0.0, "dPdU");
        CHECK_REL(EOSdUdRho(e, rho, u), P / (rho * rho), 1e-14, "dUdRho");
        CHECK_REL(EOSdPdRhoofRhoT(e, rho, T), scvheosdPdRhoofRhoT(s, rho, T), 0.0, "dPdRhoofRhoT");
        CHECK_REL(EOSdPdT(e, rho, T), scvheosdPdTofRhoT(s, rho, T), 0.0, "dPdT");
        /* EOSUCold evaluates u at pow(10, LogTMin); an ulp in T is 1e-16 in u. */
        CHECK_REL(EOSUCold(e, rho), scvheosUofRhoT(s, rho, scvhExp10(s->LogTMin)), 1e-13, "UCold");
        CHECK(EOSisbelowColdCurve(e, rho, 0.5 * EOSUCold(e, rho)) == 1, "below cold curve");
        CHECK(EOSisbelowColdCurve(e, rho, u) == 0, "above cold curve (u at T=%g)", T);
        CHECK(EOSIsInTable(e, rho, u) == EOS_SUCCESS, "in table at rho=%g T=%g", rho, T);
        CHECK(EOSIsUsable(e, rho, u) == EOS_TRUE, "usable at rho=%g T=%g", rho, T);
        /* Isentrope to the same density is the identity, where T(rho, s) is
         * unique (s(T) monotonic in this column). */
        if (scvhEntropyMonotonic(s, rho))
            CHECK_REL(EOSIsentropic(e, rho, u, rho), u, 1e-9, "isentrope to same rho");   /* s(T) round trip amplifies roundoff on flat stretches */
        else nSkippedS++;
    }
    /* Entropy is in cgs: compare against the raw table at a node. */
    {
        int j = s->nRho / 2, i = s->nT / 2;
        double rho = scvhExp10(s->dLogRhoAxis[j]), T = scvhExp10(s->dLogTAxis[i]);
        CHECK_REL(EOSSofRhoT(e, rho, T), scvhExp10(s->dLogSArray[(size_t)j * s->nT + i]), 1e-13, "entropy returned in cgs");
    }
    CHECK(s->nLookupFail == 0, "nLookupFail = %ld after in-table dispatch", s->nLookupFail);
    if (nSkipped || nSkippedS)
        printf("\n    note: %d of %d states skipped (u(T) not monotonic), %d isentrope identities skipped (s(T) not monotonic); table defects\n  %-44s",
               nSkipped, NSAMPLE_SCVH * NSAMPLE_SCVH, nSkippedS, "");
    END_TEST();
}

static void testSCVHTableStatus(EOSMATERIAL *e)
{
    int failsBefore = nFails;
    TEST("SCvH EOSIsInTable / EOSIsUsable codes");
    SCVHEOSMAT *s = e->scvheosmaterial;
    double rhoOk = scvhSampleRho(s, NSAMPLE_SCVH / 2);
    double uOk = EOSUofRhoT(e, rhoOk, scvhSampleT(s, NSAMPLE_SCVH / 2));
    double rhoLo = 0.5 * scvhExp10(s->LogRhoMin), rhoHi = 2.0 * scvhExp10(s->LogRhoMax);
    double uHi = 2.0 * EOSUofRhoT(e, rhoOk, scvhExp10(s->LogTMax));

    CHECK(EOSIsInTable(e, rhoLo, uOk) == EOS_OUTSIDE_RHOMIN, "below rhoMin");
    CHECK(EOSIsInTable(e, rhoHi, uOk) == EOS_OUTSIDE_RHOMAX, "above rhoMax");
    CHECK(EOSIsInTable(e, rhoOk, 0.5 * EOSUCold(e, rhoOk)) == EOS_OUTSIDE_VMIN, "below cold curve");
    CHECK(EOSIsInTable(e, rhoOk, uHi) == EOS_OUTSIDE_VMAX, "above top row");
    CHECK(EOSIsInTable(e, rhoOk, -1.0) == EOS_OUTSIDE_VMIN, "negative u");
    CHECK(EOSIsInTable(e, -1.0, uOk) == EOS_OUTSIDE_RHOMIN, "negative rho");
    CHECK(EOSIsUsable(e, rhoLo, uOk) == EOS_FALSE, "tabulated EOS unusable below rhoMin");
    CHECK(EOSIsUsable(e, rhoOk, uHi) == EOS_FALSE, "tabulated EOS unusable above u_max");
    CHECK(EOSIsUsable(e, rhoOk, -1.0) == EOS_FALSE, "tabulated EOS unusable for u < 0");
    /* Anything EOSIsInTable accepts must evaluate without a sentinel. */
    for (int a = 0; a < NSAMPLE_SCVH; a++) for (int b = 0; b < NSAMPLE_SCVH; b++) {
        double rho = scvhSampleRho(s, a);
        if (!scvhEnergyMonotonic(s, rho)) continue;
        double u = EOSUofRhoT(e, rho, scvhSampleT(s, b));
        if (EOSIsInTable(e, rho, u) == EOS_SUCCESS) {
            double c, T;
            double P = EOSPCTofRhoU(e, rho, u, &c, &T);
            CHECK(!isSentinel(P) && !isSentinel(c) && !isSentinel(T), "in-table state produced a sentinel at rho=%g", rho);
        }
    }
    END_TEST();
}

#ifdef HAVE_ANEOSMATERIAL_H
/* Gas envelope (SCvH, outer) on a condensed core (ANEOS, inner). */
static void testSCVHSolveBC(EOSMATERIAL *inner, EOSMATERIAL *outer)
{
    int failsBefore = nFails;
    TEST("EOSSolveBC: SCvH envelope on ANEOS core");
    SCVHEOSMAT *s = outer->scvheosmaterial;
    int nSolved = 0;
    for (int a = 2; a < NSAMPLE_SCVH - 2; a += 4) for (int b = 1; b < NSAMPLE_SCVH - 3; b += 4) {
        /* State at the base of the envelope. The solver brackets from rho1
         * upwards in the inner material, so the envelope state must be denser
         * than the core's low-density table edge, and its pressure must exceed
         * the interface-correction floor like in the ANEOS test. */
        double rho1 = scvhSampleRho(s, a), T1 = scvhSampleT(s, b);
        if (!scvhEnergyMonotonic(s, rho1)) continue;
        double u1 = EOSUofRhoT(outer, rho1, T1);
        double P1 = EOSPofRhoU(outer, rho1, u1);
        if (P1 <= EOSWoolfsonMinPressure(outer)) continue;
        double rho2, u2;
        if (EOSSolveBC(outer, inner, rho1, u1, &rho2, &u2, 0.0) != 0) continue;
        nSolved++;
        CHECK(fabs(EOSPofRhoU(inner, rho2, u2) - P1) <= solveBCPResidualBound(inner, rho2, T1, P1),
              "P continuity rho1=%g T1=%g: %.15e vs %.15e", rho1, T1, EOSPofRhoU(inner, rho2, u2), P1);
        CHECK_REL(u2, EOSUofRhoT(inner, rho2, T1), 1e-10, "T continuity rho1=%g T1=%g", rho1, T1);
        CHECK(EOSIsUsable(inner, rho2, u2) == EOS_TRUE, "solved state not usable");
    }
    if (nSolved == 0) printf("\n    note: no envelope state had a solution in the core table (table ranges do not overlap)\n  %-44s", "");
    END_TEST();
}
#endif

#endif /* HAVE_SCVHEOS_H */

/* ----------------------------------------------------------------- REOS3 */

#ifdef HAVE_REOS3_H

#define NSAMPLE_REOS3 300

static double reos3Exp10(double x) { return exp(x * 2.302585092994045684); }

/* Log-spaced samples in the interior of the REOS3 TABLE (not the
 * extrapolation region), code units / Kelvin. */
static double reos3SampleRho(REOS3MAT *r, int k)
{
    double f = 0.05 + 0.90 * (k + 0.5137) / NSAMPLE_REOS3;
    return reos3Exp10(r->dLogRhoAxis[0] + f * (r->LogRhoMax - r->dLogRhoAxis[0]));
}

static double reos3SampleT(REOS3MAT *r, int k)
{
    double f = 0.05 + 0.90 * (k + 0.5137) / NSAMPLE_REOS3;
    return reos3Exp10(r->LogTMin + f * (r->LogTMax - r->LogTMin));
}

/* Samples in the ideal-gas extrapolation region below the first rho node. */
static int reos3HasExtrapolation(REOS3MAT *r) { return r->LogRhoMin < r->dLogRhoAxis[0]; }

static double reos3SampleRhoExtrap(REOS3MAT *r, int k)
{
    double f = 0.05 + 0.90 * (k + 0.5137) / NSAMPLE_REOS3;
    return reos3Exp10(r->LogRhoMin + f * (r->dLogRhoAxis[0] - r->LogRhoMin));
}

/* 1 if u(T) increases along both node columns of the rho cell containing rho
 * (the first column for the extrapolation region). Elsewhere T(rho, u) fails
 * by contract (table defect) and the identities through it are not invariants. */
static int reos3EnergyMonotonic(REOS3MAT *r, double rho)
{
    double lr = log10(rho);
    int j = 0;
    while (j < r->nRho - 2 && r->dLogRhoAxis[j+1] <= lr) j++;
    for (int i = 0; i < r->nT - 1; i++) {
        if (!(r->dLogUArray[(size_t)j * r->nT + i + 1] > r->dLogUArray[(size_t)j * r->nT + i])) return 0;
        if (!(r->dLogUArray[(size_t)(j+1) * r->nT + i + 1] > r->dLogUArray[(size_t)(j+1) * r->nT + i])) return 0;
    }
    return 1;
}

/* Cell containing lt on the T axis, and the amplification 1/dlogu_cell of
 * roundoff in the T(rho, u) inversion there: u carries the uRef offset, so
 * where the cell's thermal energy is a tiny fraction of u the round trips
 * are only accurate to eps*(1+|logu|)/dlogu. */
static double reos3InversionTol(REOS3MAT *r, double rho, double T)
{
    double lr = log10(rho), lt = log10(T);
    int j = 0, i = 0;
    while (j < r->nRho - 2 && r->dLogRhoAxis[j+1] <= lr) j++;
    while (i < r->nT - 2 && r->dLogTAxis[i+1] <= lt) i++;
    double w = (lr - r->dLogRhoAxis[j]) / (r->dLogRhoAxis[j+1] - r->dLogRhoAxis[j]);
    if (w < 0.0) w = 0.0;   /* extrapolation region: first column */
    const double *c0 = r->dLogUArray + (size_t)j * r->nT, *c1 = c0 + r->nT;
    double lu0 = (1.0 - w) * c0[i] + w * c1[i], lu1 = (1.0 - w) * c0[i+1] + w * c1[i+1];
    double dlu = fabs(lu1 - lu0);
    if (!(dlu > 0.0)) return INFINITY;
    return 1e-12 + 64.0 * DBL_EPSILON * (1.0 + fmax(fabs(lu0), fabs(lu1))) / dlu;
}

static void testREOS3Dispatch(EOSMATERIAL *e)
{
    int failsBefore = nFails;
    TEST("REOS3 dispatch matches reos3 directly");
    REOS3MAT *r = e->reos3material;
    CHECK(e->matType == EOSREOS3, "matType");
    CHECK(!e->bEntropy && !e->bEntropyTableInit, "REOS3 has no entropy");
    CHECK_REL(e->rho0, r->rho0, 0.0, "rho0");
    CHECK(e->rho0forInterfaceCorrection == 0.0, "gas has no interface-correction cutoff");
    CHECK(strlen(e->MatString) > 0, "MatString empty");
    CHECK(EOSYieldModel(e) == 0, "yield model for gas");
    CHECK(EOSYieldStrength(e, reos3SampleRho(r, 1), 1.0) == 0.0, "yield strength for gas");
    CHECK(EOSGammaofRhoT(e, reos3SampleRho(r, 1), reos3SampleT(r, 1)) == 0.0, "shear modulus for gas");
    CHECK(EOSSofRhoU(e, reos3SampleRho(r, 1), 1.0) == EOS_FAIL_VALUE, "entropy sentinel (rho,u)");
    CHECK(EOSSofRhoT(e, reos3SampleRho(r, 1), reos3SampleT(r, 1)) == EOS_FAIL_VALUE, "entropy sentinel (rho,T)");
    /* EOSinitIsentropicLookup only warns for REOS3 (not implemented). */
    EOSinitIsentropicLookup(e, NULL);
    CHECK(!e->bEntropyTableInit, "isentropic lookup must stay uninitialised");

    int nSkipped = 0;
    for (int a = 0; a < NSAMPLE_REOS3; a++) for (int b = 0; b < NSAMPLE_REOS3; b++) {
        double rho = reos3SampleRho(r, a), T = reos3SampleT(r, b);
        if (!reos3EnergyMonotonic(r, rho)) { nSkipped++; continue; }
        double tol = reos3InversionTol(r, rho, T);
        double u = EOSUofRhoT(e, rho, T);
        double c, Tf;
        double P = EOSPCTofRhoU(e, rho, u, &c, &Tf);
        CHECK_REL(u, reos3UofRhoT(r, rho, T), 0.0, "UofRhoT dispatch");
        /* Closed-form inversion: the fused lookup reproduces the (rho, T)
         * lookups to the (amplified) roundoff of the inversion. */
        CHECK_REL(P, reos3PofRhoT(r, rho, T), 100.0 * tol, "PCTofRhoU P at rho=%g T=%g", rho, T);
        CHECK_REL(c, reos3CsofRhoT(r, rho, T), 100.0 * tol, "PCTofRhoU c at rho=%g T=%g", rho, T);
        CHECK_REL(Tf, T, tol, "PCTofRhoU T at rho=%g T=%g", rho, T);
        CHECK_REL(EOSPofRhoU(e, rho, u), P, 0.0, "PofRhoU");
        CHECK_REL(EOSCofRhoU(e, rho, u), c, 0.0, "CofRhoU");
        CHECK_REL(EOSTofRhoU(e, rho, u), reos3TofRhoU(r, rho, u), 0.0, "TofRhoU");
        CHECK_REL(EOSPofRhoT(e, rho, T), reos3PofRhoT(r, rho, T), 0.0, "PofRhoT");
        double c2; double P2 = EOSPCofRhoU(e, rho, u, &c2);
        CHECK_REL(P2, P, 0.0, "PCofRhoU P"); CHECK_REL(c2, c, 0.0, "PCofRhoU c");
        CHECK_REL(EOSdPdRho(e, rho, u), reos3dPdRhoofRhoU(r, rho, u), 0.0, "dPdRho");
        CHECK_REL(EOSdPdU(e, rho, u), reos3dPdUofRhoU(r, rho, u), 0.0, "dPdU");
        CHECK_REL(EOSdUdRho(e, rho, u), P / (rho * rho), 1e-14, "dUdRho");
        CHECK_REL(EOSdPdRhoofRhoT(e, rho, T), reos3dPdRhoofRhoT(r, rho, T), 0.0, "dPdRhoofRhoT");
        CHECK_REL(EOSdPdT(e, rho, T), reos3dPdTofRhoT(r, rho, T), 0.0, "dPdT");
        /* EOSUCold evaluates u at pow(10, LogTMin); an ulp in T is 1e-16 in u. */
        CHECK_REL(EOSUCold(e, rho), reos3UofRhoT(r, rho, reos3Exp10(r->LogTMin)), 1e-13, "UCold");
        CHECK(EOSisbelowColdCurve(e, rho, belowColdCurve(e, rho)) == 1, "below cold curve");
        CHECK(EOSisbelowColdCurve(e, rho, u) == 0, "above cold curve (u at T=%g)", T);
        CHECK(EOSIsInTable(e, rho, u) == EOS_SUCCESS, "in table at rho=%g T=%g", rho, T);
        CHECK(EOSIsUsable(e, rho, u) == EOS_TRUE, "usable at rho=%g T=%g", rho, T);
    }
    /* rho(P, T) dispatch and round trip, where the row is monotonic (the
     * library fails otherwise, by contract). */
    for (int a = 0; a < NSAMPLE_REOS3; a += 7) for (int b = 0; b < NSAMPLE_REOS3; b += 7) {
        double rho = reos3SampleRho(r, a), T = reos3SampleT(r, b);
        double P = EOSPofRhoT(e, rho, T);
        if (reos3IsInExtrapLimit(r, rho, EOSUofRhoT(e, rho, T)) != REOS3_SUCCESS) continue;
        double lt = log10(T);
        int i = 0; while (i < r->nT - 2 && r->dLogTAxis[i+1] <= lt) i++;
        int mono = 1;
        for (int j = 0; j < r->nRho - 1 && mono; j++)
            if (!(r->dLogPArray[(size_t)(j+1) * r->nT + i] > r->dLogPArray[(size_t)j * r->nT + i]) ||
                !(r->dLogPArray[(size_t)(j+1) * r->nT + i + 1] > r->dLogPArray[(size_t)j * r->nT + i + 1])) mono = 0;
        if (!mono) continue;
        CHECK_REL(EOSRhoofPT(e, P, T), reos3RhoofPT(r, P, T), 0.0, "RhoofPT dispatch");
        CHECK_REL(EOSRhoofPT(e, P, T), rho, 1e-12, "rho(P,T) roundtrip at rho=%g T=%g", rho, T);
    }
    CHECK(r->nLookupFail == 0, "nLookupFail = %ld after in-table dispatch", r->nLookupFail);
    if (nSkipped)
        printf("\n    note: %d of %d states skipped (u(T) not monotonic, table defect)\n  %-44s",
               nSkipped, NSAMPLE_REOS3 * NSAMPLE_REOS3, "");
    END_TEST();
}

static void testREOS3Extrapolation(EOSMATERIAL *e)
{
    int failsBefore = nFails;
    TEST("REOS3 ideal-gas extrapolation through EOSlib");
    REOS3MAT *r = e->reos3material;
    if (!reos3HasExtrapolation(r)) {
        printf("\n    note: extrapolation disabled (REOS3_RHOMIN_CGS at or above the first node)\n  %-44s", "");
        END_TEST();
        return;
    }
    /* Below the first rho node EOSlib must report the state as inside the
     * lookup limits (EOSIsInTable delegates to reos3IsInExtrapLimit) and
     * evaluate it as the isothermal ideal gas the library defines. */
    double rho0 = reos3Exp10(r->dLogRhoAxis[0]);
    int nSkipped = 0;
    for (int a = 0; a < NSAMPLE_REOS3; a += 3) for (int b = 0; b < NSAMPLE_REOS3; b += 3) {
        double rho = reos3SampleRhoExtrap(r, a), T = reos3SampleT(r, b);
        double u = EOSUofRhoT(e, rho, T);
        CHECK_REL(u, EOSUofRhoT(e, rho0, T), 1e-12, "u of the first column at rho=%g T=%g", rho, T);
        CHECK_REL(EOSPofRhoT(e, rho, T), EOSPofRhoT(e, rho0, T) * rho / rho0, 1e-12, "P scales with rho at rho=%g T=%g", rho, T);
        CHECK(EOSIsInTable(e, rho, u) == EOS_SUCCESS, "extrapolation region is inside the lookup limits at rho=%g", rho);
        CHECK(EOSIsUsable(e, rho, u) == EOS_TRUE, "extrapolation region usable at rho=%g", rho);
        if (!reos3EnergyMonotonic(r, rho)) { nSkipped++; continue; }
        double c, Tf;
        double P = EOSPCTofRhoU(e, rho, u, &c, &Tf);
        CHECK(!isSentinel(P) && !isSentinel(c) && !isSentinel(Tf), "fused lookup in the extrapolation region at rho=%g", rho);
        CHECK_REL(Tf, T, reos3InversionTol(r, rho, T), "T(rho,u) in the extrapolation region at rho=%g T=%g", rho, T);
        CHECK_REL(EOSdPdRho(e, rho, u), P / rho, 1e-12, "dP/drho|u = P/rho in the extrapolation region at rho=%g", rho);
    }
    if (nSkipped) printf("\n    note: %d states skipped (u(T) not monotonic in the first column)\n  %-44s", nSkipped, "");
    END_TEST();
}

static void testREOS3TableStatus(EOSMATERIAL *e)
{
    int failsBefore = nFails;
    TEST("REOS3 EOSIsInTable / EOSIsUsable codes");
    REOS3MAT *r = e->reos3material;
    double rhoOk = reos3SampleRho(r, NSAMPLE_REOS3 / 2);
    double uOk = EOSUofRhoT(e, rhoOk, reos3SampleT(r, NSAMPLE_REOS3 / 2));
    double rhoLo = 0.5 * reos3Exp10(r->LogRhoMin), rhoHi = 2.0 * reos3Exp10(r->LogRhoMax);
    double uMin = EOSUofRhoT(e, rhoOk, reos3Exp10(r->LogTMin)), uMax = EOSUofRhoT(e, rhoOk, reos3Exp10(r->LogTMax));
    double uLo = uMin - fabs(uMin) - 1.0, uHi = uMax + fabs(uMax) + 1.0;

    CHECK(EOSIsInTable(e, rhoOk, uOk) == EOS_SUCCESS, "in table");
    CHECK(EOSIsInTable(e, reos3Exp10(r->LogRhoMin), uOk) == EOS_SUCCESS, "on the lookup limit");
    CHECK(EOSIsInTable(e, reos3Exp10(r->LogRhoMax), EOSUofRhoT(e, reos3Exp10(r->LogRhoMax), reos3SampleT(r, NSAMPLE_REOS3 / 2))) == EOS_SUCCESS, "on rhoMax");
    CHECK(EOSIsInTable(e, rhoOk, uMin) == EOS_SUCCESS, "on the cold curve");
    CHECK(EOSIsInTable(e, rhoOk, uMax) == EOS_SUCCESS, "on the top row");
    CHECK(EOSIsInTable(e, rhoLo, uOk) == EOS_OUTSIDE_RHOMIN, "below the lookup limit");
    CHECK(EOSIsInTable(e, rhoHi, uOk) == EOS_OUTSIDE_RHOMAX, "above rhoMax");
    CHECK(EOSIsInTable(e, rhoOk, uLo) == EOS_OUTSIDE_VMIN, "below cold curve");
    CHECK(EOSIsInTable(e, rhoOk, uHi) == EOS_OUTSIDE_VMAX, "above top row");
    CHECK(EOSIsInTable(e, -1.0, uOk) == EOS_OUTSIDE_RHOMIN, "negative rho");
    CHECK(EOSIsInTable(e, rhoOk, NAN) == EOS_OUTSIDE_VMIN, "NaN u");
    CHECK(EOSIsUsable(e, rhoLo, uOk) == EOS_FALSE, "tabulated EOS unusable below the lookup limit");
    CHECK(EOSIsUsable(e, rhoOk, uHi) == EOS_FALSE, "tabulated EOS unusable above u_max");
    CHECK(EOSIsUsable(e, rhoOk, uLo) == EOS_FALSE, "tabulated EOS unusable below the cold curve");
    /* Anything EOSIsInTable accepts must evaluate without a sentinel. */
    for (int a = 0; a < NSAMPLE_REOS3; a++) for (int b = 0; b < NSAMPLE_REOS3; b++) {
        double rho = reos3SampleRho(r, a);
        if (!reos3EnergyMonotonic(r, rho)) continue;
        double u = EOSUofRhoT(e, rho, reos3SampleT(r, b));
        if (EOSIsInTable(e, rho, u) == EOS_SUCCESS) {
            double c, T;
            double P = EOSPCTofRhoU(e, rho, u, &c, &T);
            CHECK(!isSentinel(P) && !isSentinel(c) && !isSentinel(T), "in-table state produced a sentinel at rho=%g", rho);
        } else {
            CHECK(0, "in-table state rejected at rho=%g u=%g", rho, u);
        }
    }
    END_TEST();
}

#ifdef HAVE_ANEOSMATERIAL_H
/* Gas envelope (REOS3, outer) on a condensed core (ANEOS, inner). */
static void testREOS3SolveBC(EOSMATERIAL *inner, EOSMATERIAL *outer)
{
    int failsBefore = nFails;
    TEST("EOSSolveBC: REOS3 envelope on ANEOS core");
    REOS3MAT *r = outer->reos3material;
    int nSolved = 0;
    for (int a = 2; a < NSAMPLE_REOS3 - 2; a += 4) for (int b = 1; b < NSAMPLE_REOS3 - 3; b += 4) {
        /* Same conditions as the SCvH envelope test: the solver brackets from
         * rho1 upwards in the inner material, and P must exceed the
         * interface-correction floor. */
        double rho1 = reos3SampleRho(r, a), T1 = reos3SampleT(r, b);
        if (!reos3EnergyMonotonic(r, rho1)) continue;
        double u1 = EOSUofRhoT(outer, rho1, T1);
        double P1 = EOSPofRhoU(outer, rho1, u1);
        if (P1 <= EOSWoolfsonMinPressure(outer)) continue;
        double rho2, u2;
        if (EOSSolveBC(outer, inner, rho1, u1, &rho2, &u2, 0.0) != 0) continue;
        nSolved++;
        CHECK(fabs(EOSPofRhoU(inner, rho2, u2) - P1) <= solveBCPResidualBound(inner, rho2, T1, P1),
              "P continuity rho1=%g T1=%g: %.15e vs %.15e", rho1, T1, EOSPofRhoU(inner, rho2, u2), P1);
        CHECK_REL(u2, EOSUofRhoT(inner, rho2, T1), 1e-10, "T continuity rho1=%g T1=%g", rho1, T1);
        CHECK(EOSIsUsable(inner, rho2, u2) == EOS_TRUE, "solved state not usable");
    }
    if (nSolved == 0) printf("\n    note: no envelope state had a solution in the core table (table ranges do not overlap)\n  %-44s", "");
    END_TEST();
}
#endif

#endif /* HAVE_REOS3_H */

/* ------------------------------------------------------ failure contract */

/*
 * The wrapper's failure contract for a tabulated material (ANEOS, SCvH,
 * REOS3): outside the table every lookup returns EOS_FAIL_VALUE, except
 * EOSTofRhoU, which falls back to the lowest tabulated temperature; the
 * fused lookup sets its outputs the same way; the status function names the
 * reason; everything is counted in the sublibrary's nLookupFail. rhoHi is a
 * density above the table, TMin the lowest tabulated temperature, (rhoOk,
 * uOk) a state inside the table and uHi an energy above its column.
 */
static void testFailureContract(EOSMATERIAL *e, double rhoOk, double uOk, double uHi, double rhoHi, double TMin, long *pnLookupFail)
{
    int failsBefore = nFails;
    TEST("failure contract through the wrapper");
    long nf0 = *pnLookupFail;
    double c = 1.0, T = 1.0, P;

    CHECK(EOSIsInTable(e, rhoHi, uOk) == EOS_OUTSIDE_RHOMAX, "status above rhoMax");
    CHECK(EOSIsInTable(e, rhoOk, uHi) == EOS_OUTSIDE_VMAX, "status above the column");
    CHECK(EOSIsUsable(e, rhoHi, uOk) == EOS_FALSE, "usable above rhoMax");
    CHECK(EOSFailed(EOSPofRhoU(e, rhoHi, uOk)), "P above rhoMax");
    CHECK(EOSFailed(EOSCofRhoU(e, rhoHi, uOk)), "c above rhoMax");
    CHECK(EOSFailed(EOSUofRhoT(e, rhoHi, TMin)), "u(rho, T) above rhoMax");
    CHECK(EOSFailed(EOSPofRhoT(e, rhoHi, TMin)), "P(rho, T) above rhoMax");
    CHECK(EOSFailed(EOSdPdRho(e, rhoHi, uOk)), "dP/drho above rhoMax");
    CHECK(EOSFailed(EOSdPdU(e, rhoHi, uOk)), "dP/du above rhoMax");
    CHECK_REL(EOSTofRhoU(e, rhoHi, uOk), TMin, 1e-12, "T(rho, u) above rhoMax falls back to T_min");
    CHECK_REL(EOSTofRhoU(e, rhoOk, uHi), TMin, 1e-12, "T(rho, u) above the column falls back to T_min");
    CHECK_REL(EOSPofRhoU(e, rhoOk, uHi), EOSPofRhoT(e, rhoOk, TMin), 1e-12, "P(rho, u) above the column = P(rho, T_min)");
    P = EOSPCTofRhoU(e, rhoHi, uOk, &c, &T);
    CHECK(EOSFailed(P) && EOSFailed(c), "fused lookup above rhoMax: P=%g c=%g", P, c);
    CHECK_REL(T, TMin, 1e-12, "fused lookup above rhoMax: T");
    P = EOSPCTofRhoU(e, rhoOk, uHi, &c, &T);
    CHECK(!EOSFailed(P) && !EOSFailed(c), "fused lookup above the column: P=%g c=%g", P, c);
    CHECK_REL(T, TMin, 1e-12, "fused lookup above the column: T");
    if (e->bEntropy) CHECK(EOSFailed(EOSSofRhoU(e, rhoHi, uOk)), "s above rhoMax");
    CHECK(*pnLookupFail > nf0, "failures not counted (%ld -> %ld)", nf0, *pnLookupFail);
    *pnLookupFail = nf0;
    END_TEST();
}

/* ------------------------------------------------------------- Tillotson */

#ifdef HAVE_TILLOTSON_H
#define NSAMPLE_TILL 300

/* Log-spaced in rho over the lookup table (code units). */
static double tillSampleRho(TILLMATERIAL *m, int k)
{
    double lo = log(m->rhomin), hi = log(m->rhomax);
    double f = 0.02 + 0.96 * (k + 0.5) / NSAMPLE_TILL;
    return exp(lo + f * (hi - lo));
}

/* Log-spaced in T from 10 K to 1e5 K: u = u_cold + cv*T stays far below the
 * top of the isentrope table (vmax) for every material. */
static double tillSampleT(int k)
{
    double lo = log(10.0), hi = log(1e5);
    double f = 0.02 + 0.96 * (k + 0.5) / NSAMPLE_TILL;
    return exp(lo + f * (hi - lo));
}

static void testTillotsonDispatch(EOSMATERIAL *e)
{
    int failsBefore = nFails;
    TEST("Tillotson dispatch matches tillotson directly");
    TILLMATERIAL *m = e->tillmaterial;
    CHECK(e->matType == EOSTILLOTSON, "matType");
    CHECK(!e->bEntropy && e->bEntropyTableInit, "entropy flags");
    CHECK_REL(e->rho0, m->rho0, 0.0, "rho0");
    CHECK_REL(e->rho0forInterfaceCorrection, m->rho0, 0.0, "interface-correction cutoff");
    CHECK_REL(e->minSoundSpeed, sqrt(m->A / m->rho0), 1e-15, "minSoundSpeed");
    CHECK(strlen(e->MatString) > 0, "MatString empty");
    CHECK(EOSYieldModel(e) == -2, "yield model %d (no strength for Tillotson)", EOSYieldModel(e));
    CHECK(EOSSofRhoU(e, m->rho0, 1.0) == EOS_FAIL_VALUE, "entropy sentinel (rho,u)");
    CHECK(EOSSofRhoT(e, m->rho0, 300.0) == EOS_FAIL_VALUE, "entropy sentinel (rho,T)");
    long nf0 = m->nLookupFail;

    int nRhoPT = 0, nInTable = 0, nAbove = 0;
    for (int a = 0; a < NSAMPLE_TILL; a++) for (int b = 0; b < NSAMPLE_TILL; b++) {
        double rho = tillSampleRho(m, a), T = tillSampleT(b);
        double u = EOSUofRhoT(e, rho, T);
        double c, Tf;
        double P = EOSPCTofRhoU(e, rho, u, &c, &Tf);
        CHECK_REL(u, tillURhoTemp(m, rho, T), 0.0, "UofRhoT dispatch");
        CHECK_REL(P, tillPressure(m, rho, u), 0.0, "PCTofRhoU P");
        CHECK_REL(EOSPofRhoU(e, rho, u), P, 0.0, "PofRhoU");
        CHECK_REL(EOSPofRhoT(e, rho, T), P, 0.0, "PofRhoT");
        CHECK_REL(EOSCofRhoU(e, rho, u), c, 0.0, "CofRhoU");
        CHECK(c >= e->minSoundSpeed * (1.0 - 1e-12) && isfinite(c), "c = %g below the floor at rho=%g T=%g", c, rho, T);
        CHECK_REL(Tf, tillTempRhoU(m, rho, u), 0.0, "TofRhoU");
        CHECK_REL(EOSTofRhoU(e, rho, u), Tf, 0.0, "TofRhoU dispatch");
        /* u(rho, T) = u_cold + cv*T is exactly invertible up to roundoff in u_cold + cv*T. */
        CHECK_REL(Tf, T, 1e-8, "T(rho, u(rho,T)) round trip at rho=%g T=%g", rho, T);
        double c2; double P2 = EOSPCofRhoU(e, rho, u, &c2);
        CHECK_REL(P2, P, 0.0, "PCofRhoU P"); CHECK_REL(c2, c, 0.0, "PCofRhoU c");
        CHECK_REL(EOSdPdRho(e, rho, u), tilldPdrho(m, rho, u), 0.0, "dPdRho");
        CHECK_REL(EOSdPdU(e, rho, u), tilldPdu(m, rho, u), 0.0, "dPdU");
        CHECK_REL(EOSdUdRho(e, rho, u), tilldudrho(m, rho, u), 0.0, "dUdRho");
        CHECK_REL(EOSUCold(e, rho), tillColdULookup(m, rho), 0.0, "UCold");
        CHECK(EOSUCold(e, rho) >= 0.0 || rho >= m->rho0, "u_cold < 0 in the expanded region at rho=%g", rho);
        CHECK(EOSisbelowColdCurve(e, rho, belowColdCurve(e, rho)) == 1, "below cold curve");
        CHECK(EOSisbelowColdCurve(e, rho, u) == 0, "above cold curve (u at T=%g)", T);
        /* The isentrope table ends at s_max; in the expanded region that is a
         * low u, so hot expanded states lie above it. They stay usable. */
        int st = EOSIsInTable(e, rho, u);
        CHECK(st == EOS_SUCCESS || st == EOS_OUTSIDE_VMAX, "status %d at rho=%g T=%g", st, rho, T);
        CHECK(EOSIsUsable(e, rho, u) == EOS_TRUE, "usable at rho=%g T=%g", rho, T);
        if (st == EOS_SUCCESS) {
            /* Isentropic lookup: identity at the same density, dispatch otherwise. */
            CHECK_REL(EOSIsentropic(e, rho, u, rho), u, 1e-8, "isentrope to same rho at rho=%g T=%g", rho, T);
            if (a + 1 < NSAMPLE_TILL) {
                double rho2 = tillSampleRho(m, a + 1);
                CHECK_REL(EOSIsentropic(e, rho, u, rho2), tillLookupU(m, rho, u, rho2, 0), 0.0, "isentrope dispatch");
            }
            nInTable++;
        } else nAbove++;
        /* rho(P, T): unique in the condensed region (rho >= rho0) at P > 0. */
        if (rho >= m->rho0 && P > 0.0) {
            CHECK_REL(EOSRhoofPT(e, P, T), tillRhoPTemp(m, P, T), 0.0, "RhoofPT dispatch");
            CHECK_REL(EOSRhoofPT(e, P, T), rho, 1e-8, "rho(P(rho,T),T) at rho=%g T=%g", rho, T);
            nRhoPT++;
        }
    }
    CHECK(nRhoPT > 0, "no rho(P,T) round trip tested");
    CHECK(nInTable > 0, "no state inside the isentrope table");
    CHECK(m->nLookupFail == nf0, "nLookupFail %ld -> %ld after in-table dispatch", nf0, m->nLookupFail);
    if (nAbove) printf("\n    note: %d of %d states above the isentrope table (hot expanded states), usable but no isentrope\n  %-44s", nAbove, nInTable + nAbove, "");
    END_TEST();
}

static void testTillotsonTableStatus(EOSMATERIAL *e)
{
    int failsBefore = nFails;
    TEST("Tillotson EOSIsInTable / EOSIsUsable codes");
    TILLMATERIAL *m = e->tillmaterial;
    double rhoOk = 2.0 * m->rho0;
    double uOk = EOSUofRhoT(e, rhoOk, 1000.0);
    double ucold = EOSUCold(e, rhoOk);
    double uHuge = 1e4 * EOSUofRhoT(e, rhoOk, 1e6);   /* far above the isentrope table */
    double rhoLo = 0.5 * m->rhomin, rhoHi = 2.0 * m->rhomax;
    long nf0 = m->nLookupFail;

    CHECK(rhoOk < m->rhomax, "table does not cover 2 rho0");
    CHECK(EOSIsInTable(e, rhoOk, uOk) == EOS_SUCCESS, "in table");
    CHECK(EOSIsInTable(e, m->rhomin, EOSUCold(e, m->rhomin)) == EOS_SUCCESS, "on rhomin (closed interval)");
    CHECK(EOSIsInTable(e, m->rhomax, EOSUCold(e, m->rhomax)) == EOS_SUCCESS, "on rhomax (closed interval)");
    CHECK(EOSIsInTable(e, rhoOk, ucold) == EOS_SUCCESS, "on the cold curve");
    CHECK(EOSIsInTable(e, rhoLo, uOk) == EOS_OUTSIDE_RHOMIN, "below rhomin");
    CHECK(EOSIsInTable(e, rhoHi, uOk) == EOS_OUTSIDE_RHOMAX, "above rhomax");
    CHECK(EOSIsInTable(e, -1.0, uOk) == EOS_OUTSIDE_RHOMIN, "negative rho");
    CHECK(EOSIsInTable(e, NAN, uOk) == EOS_OUTSIDE_RHOMIN, "NaN rho");
    CHECK(EOSIsInTable(e, rhoOk, 0.5 * ucold) == EOS_OUTSIDE_VMIN, "below the cold curve");
    CHECK(EOSIsInTable(e, rhoOk, uHuge) == EOS_OUTSIDE_VMAX, "above the table");
    /* A non-finite u is no state at all: EOS_FAIL, not OUTSIDE_VMAX (which EOSIsUsable tolerates). */
    CHECK(EOSIsInTable(e, rhoOk, NAN) == EOS_FAIL, "NaN u");
    CHECK(EOSIsInTable(e, rhoOk, INFINITY) == EOS_FAIL, "+inf u");

    /* Policy: the pressure is analytic, only u above the table is tolerated. */
    CHECK(EOSIsUsable(e, rhoOk, uOk) == EOS_TRUE, "usable in table");
    CHECK(EOSIsUsable(e, rhoOk, uHuge) == EOS_TRUE, "usable above the table (analytic P)");
    CHECK(EOSIsUsable(e, rhoOk, 0.5 * ucold) == EOS_FALSE, "not usable below the cold curve");
    CHECK(EOSIsUsable(e, rhoLo, uOk) == EOS_FALSE, "not usable below rhomin");
    CHECK(EOSIsUsable(e, rhoHi, uOk) == EOS_FALSE, "not usable above rhomax");
    CHECK(EOSIsUsable(e, rhoOk, NAN) == EOS_FALSE, "not usable for NaN u");
    CHECK(EOSIsUsable(e, rhoOk, INFINITY) == EOS_FALSE, "not usable for +inf u");
    CHECK(isfinite(EOSPofRhoU(e, rhoOk, uHuge)), "P above the table is finite");
    CHECK(isfinite(EOSCofRhoU(e, rhoOk, uHuge)), "c above the table is finite");
    CHECK(EOSisbelowColdCurve(e, rhoOk, 0.5 * ucold) == 1, "below cold curve");
    CHECK(EOSisbelowColdCurve(e, rhoOk, uOk) == 0, "above cold curve");

    /* Failure contract: outside the table in rho the cold curve is unknown;
     * the isentrope lookup integrates directly instead (a counted fallback). */
    long nf = m->nLookupFail;
    CHECK(EOSFailed(EOSUofRhoT(e, rhoHi, 300.0)), "u(rho, T) above rhomax");
    /* P(rho, T) needs u(rho, T): the failure must reach the caller, not be clipped to P = 0. */
    CHECK(EOSFailed(EOSPofRhoT(e, rhoHi, 300.0)), "P(rho, T) above rhomax");
    CHECK(EOSFailed(EOSUCold(e, rhoHi)), "u_cold above rhomax");
    CHECK(EOSTofRhoU(e, rhoHi, uOk) == 0.0, "T(rho, u) above rhomax falls back to T = 0");
    CHECK(EOSUCold(e, rhoLo) == 0.0, "u_cold below rhomin is zero (rho < rho0)");
    CHECK(EOSFailed(EOSRhoofPT(e, -1.0, 300.0)), "rho(P < 0, T)");
    /* Below the cold curve the isentrope is the cold curve (the ANEOS/SCvH convention), counted. */
    CHECK(EOSIsentropic(e, rhoOk, 0.5 * ucold, 1.5 * rhoOk) == EOSUCold(e, 1.5 * rhoOk), "isentrope from below the cold curve clamps to u_cold(rho2)");
    CHECK(EOSFailed(EOSIsentropic(e, rhoOk, NAN, rhoOk)), "isentrope from a NaN u");
    CHECK(m->nLookupFail == nf + 7, "failures not counted (%ld -> %ld)", nf, m->nLookupFail);
    CHECK(isfinite(EOSIsentropic(e, rhoHi, uOk, rhoOk)), "isentrope from above rhomax integrates directly");
    CHECK(isfinite(EOSIsentropic(e, rhoOk, uOk, rhoHi)), "isentrope to above rhomax integrates directly");
    CHECK(m->nLookupFail == nf + 9, "direct integration counted as a fallback (%ld -> %ld)", nf, m->nLookupFail);
    CHECK(EOSFailed(EOSIsentropic(e, rhoOk, uOk, NAN)), "isentrope to a NaN rho2");
    CHECK(EOSFailed(EOSIsentropic(e, rhoOk, uOk, INFINITY)), "isentrope to an infinite rho2");
    CHECK(EOSFailed(EOSIsentropic(e, NAN, uOk, rhoOk)), "isentrope from a NaN rho1");
    CHECK(m->nLookupFail == nf + 12, "invalid densities counted (%ld -> %ld)", nf, m->nLookupFail);
    /* T(rho, u) below the cold curve falls back to T = 0 (contract), counted, never negative. */
    CHECK(EOSTofRhoU(e, rhoOk, 0.5 * ucold) == 0.0, "T(rho, u) below the cold curve falls back to T = 0");
    {
        double c = -1.0, T = -1.0;
        double P = EOSPCTofRhoU(e, rhoOk, 0.5 * ucold, &c, &T);
        CHECK(T == 0.0 && isfinite(P) && isfinite(c), "fused lookup below the cold curve: T=%g (expected 0)", T);
    }
    CHECK(EOSTofRhoU(e, rhoOk, ucold) == 0.0, "T(rho, u) on the cold curve is 0");
    CHECK(m->nLookupFail == nf + 14, "T fallback below the cold curve counted (%ld -> %ld)", nf, m->nLookupFail);
    /* u < 0 (below the cold curve) at the former singularity w0 = 0 of the Tillotson Gamma: finite. */
    {
        double uneg = -m->u0 * (rhoOk / m->rho0) * (rhoOk / m->rho0);
        double c = -1.0, T = -1.0;
        double P = EOSPCTofRhoU(e, rhoOk, uneg, &c, &T);
        CHECK(isfinite(P) && P >= 0.0 && isfinite(c) && c > 0.0, "P, c at u = -u0 eta^2: P=%g c=%g", P, c);
        CHECK(T == 0.0, "T at u = -u0 eta^2 is 0, not negative: T=%g", T);
        CHECK(isfinite(EOSPofRhoU(e, 2.0 * rhoOk, 2.0 * uneg)), "P at u < -u0 eta^2");
    }
    CHECK(EOSFailed(EOSUofRhoT(e, -1.0, 300.0)), "u(rho < 0, T)");
    {
        double c = 1.0, T = 1.0;
        double P = EOSPCTofRhoU(e, rhoHi, uOk, &c, &T);
        CHECK(isfinite(P) && isfinite(c) && T == 0.0, "fused lookup above rhomax: P=%g c=%g T=%g (analytic P, c; T = 0)", P, c, T);
    }
    m->nLookupFail = nf0;
    END_TEST();
}

/* Interface between two Tillotson materials: the denser one is the inner one. */
static void testTillotsonSolveBC(EOSMATERIAL *inner, EOSMATERIAL *outer)
{
    int failsBefore = nFails;
    TEST("EOSSolveBC / EOSWoolfsonCoeff: two Tillotson materials");
    TILLMATERIAL *m = outer->tillmaterial;
    int nSolved = 0, nWoolfson = 0;
    for (int a = 0; a < NSAMPLE_TILL; a++) for (int b = 0; b < NSAMPLE_TILL; b++) {
        double rho1 = tillSampleRho(m, a), T1 = tillSampleT(b);
        if (rho1 < outer->rho0 || rho1 > 3.0 * outer->rho0) continue;   /* condensed interface states */
        double u1 = EOSUofRhoT(outer, rho1, T1);
        double P1 = EOSPofRhoU(outer, rho1, u1);
        if (P1 <= EOSWoolfsonMinPressure(outer)) continue;
        double rho2, u2;
        for (int jump = 0; jump < 2; jump++) {
            double Tjump = jump ? 0.05 * T1 : 0.0;
            int ret = EOSSolveBC(outer, inner, rho1, u1, &rho2, &u2, Tjump);
            CHECK(ret == EOS_SUCCESS || ret == EOS_FAIL, "SolveBC: bad return %d", ret);
            if (ret != EOS_SUCCESS) continue;
            nSolved++;
            CHECK(fabs(EOSPofRhoU(inner, rho2, u2) - P1) <= solveBCPResidualBound(inner, rho2, T1 - Tjump, P1),
                  "P continuity rho1=%g T1=%g jump=%g: %.15e vs %.15e", rho1, T1, Tjump, EOSPofRhoU(inner, rho2, u2), P1);
            CHECK_REL(EOSTofRhoU(inner, rho2, u2), T1 - Tjump, 1e-8, "T continuity rho1=%g T1=%g jump=%g", rho1, T1, Tjump);
            CHECK(EOSIsUsable(inner, rho2, u2) == EOS_TRUE, "solved state not usable");
        }
        /* Woolfson correction: the ratio of the two inversions, or 1 in the expanded region. */
        double f = EOSWoolfsonCoeff(outer, inner, P1, T1);
        CHECK(EOSFailed(f) || (f > 0.0 && isfinite(f)), "f_ij = %g at P=%g T=%g", f, P1, T1);
        if (!EOSFailed(f) && f != 1.0) {
            double r1 = EOSRhoofPT(outer, P1, T1), r2 = EOSRhoofPT(inner, P1, T1);
            CHECK_REL(f, r1 / r2, 1e-14, "f_ij = rho1/rho2 at P=%g T=%g", P1, T1);
            nWoolfson++;
        }
    }
    CHECK(nSolved > 0, "EOSSolveBC never succeeded");
    CHECK(nWoolfson > 0, "EOSWoolfsonCoeff never corrected");
    CHECK(EOSWoolfsonCoeff(outer, outer, 1.0, 300.0) == 1.0, "same material");
    CHECK(EOSWoolfsonCoeff(outer, inner, 0.5 * EOSWoolfsonMinPressure(outer), 300.0) == 1.0, "below minimum pressure");
    /* Failure path: a negative target temperature. */
    {
        double rho1 = 1.5 * outer->rho0, u1 = EOSUofRhoT(outer, rho1, 1000.0);
        double rho2 = -1.0, u2 = -1.0;
        int ret = EOSSolveBC(outer, inner, rho1, u1, &rho2, &u2, 1e4);
        CHECK(ret == EOS_FAIL, "SolveBC with T < 0 must fail (ret=%d)", ret);
        CHECK(rho2 == -1.0 && u2 == -1.0, "outputs touched on failure");
    }
    END_TEST();
}

/*
 * Stiff, low-pressure interface: granite on iron at small compression and very low T, where
 * K/P is 1e4..1e7. The pressure gap of the final density bracket, ~ (K/P) EOS_SOLVEBC_TOL P,
 * then exceeds EOS_SOLVEBC_PTOL P although P2(rho, T) is continuous; the solver must succeed
 * (it rejected these states as discontinuous up to EOSlib 1.2.0 before the slope test).
 */
static void testTillotsonSolveBCStiff(EOSMATERIAL *inner, EOSMATERIAL *outer)
{
    int failsBefore = nFails;
    TEST("EOSSolveBC: stiff low-pressure interface");
    static const double dmu[] = {0.0, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3};
    static const double Ts[] = {1e-3, 1e-2, 0.1, 1.0, 10.0};
    int nStiff = 0;
    for (size_t i = 0; i < sizeof(dmu)/sizeof(dmu[0]); i++) for (size_t j = 0; j < sizeof(Ts)/sizeof(Ts[0]); j++) {
        double rho1 = outer->rho0 * (1.0 + dmu[i]), T1 = Ts[j];
        double u1 = EOSUofRhoT(outer, rho1, T1);
        double P1 = EOSPofRhoU(outer, rho1, u1);
        CHECK(P1 > 0.0, "P1 = %g at rho1=%g T1=%g", P1, rho1, T1);
        double rho2 = -1.0, u2 = -1.0;
        int ret = EOSSolveBC(outer, inner, rho1, u1, &rho2, &u2, 0.0);
        CHECK(ret == EOS_SUCCESS, "SolveBC failed for a continuous P2: rho1=%g T1=%g P1=%g (ret=%d)", rho1, T1, P1, ret);
        if (ret != EOS_SUCCESS) continue;
        double P2 = EOSPofRhoU(inner, rho2, u2);
        /* Stiffness of the inner material at the solution, K/P = dlnP/dlnrho. */
        double h = 1e-6 * rho2;
        double K = rho2 * (EOSPofRhoT(inner, rho2 + h, T1) - EOSPofRhoT(inner, rho2 - h, T1)) / (2.0 * h);
        if (K > 1e4 * P1) nStiff++;
        CHECK(fabs(P2 - P1) <= solveBCPResidualBound(inner, rho2, T1, P1),
              "P continuity rho1=%g T1=%g K/P=%g: %.15e vs %.15e", rho1, T1, K / P1, P2, P1);
        CHECK_REL(EOSTofRhoU(inner, rho2, u2), T1, 1e-8, "T continuity rho1=%g T1=%g", rho1, T1);
        CHECK(EOSIsUsable(inner, rho2, u2) == EOS_TRUE, "solved state not usable");
    }
    /* The test is only meaningful if it reaches the regime the old gap test rejected. */
    CHECK(nStiff >= 10, "only %d interface states with K/P > 1e4", nStiff);
    END_TEST();
}

/*
 * A genuine discontinuity must still be detected. With TILL_PRESS_MELOSH the cold term of
 * expanded cold states is cut below 0.8 rho0, so for a hot enough (but still cold-region)
 * state P(rho, T) of granite jumps from P+ > 0 to 0 at 0.8 rho0. A target pressure P+/2,
 * taken from an iron state at the same T, has no root: EOSSolveBC must fail. Without the
 * flag the same P(rho, T) is continuous and the solver must succeed.
 */
static void testTillotsonSolveBCJump(EOSMATERIAL *iron, EOSMATERIAL *granite)
{
    int failsBefore = nFails;
    TEST("EOSSolveBC: discontinuous P2(rho, T)");
    TILLMATERIAL *g = granite->tillmaterial;
    /* u = cv T below us everywhere below rho0 (the cold curve is zero there). */
    double T = 0.5 * g->us / g->cv;
    double rhoJ = 0.8 * granite->rho0;
    double Pplus = EOSPofRhoT(granite, rhoJ * (1.0 + 1e-9), T);
    double Pminus = EOSPofRhoT(granite, rhoJ * (1.0 - 1e-9), T);
    CHECK(Pplus > 0.0, "P(0.8 rho0+) = %g must be positive for this test", Pplus);
#ifdef TILL_PRESS_MELOSH
    CHECK(Pminus == 0.0, "P(0.8 rho0-) = %g, expected the Melosh cut", Pminus);
#else
    CHECK_REL(Pminus, Pplus, 1e-6, "P(rho, T) continuous at 0.8 rho0");
#endif
    /* P(rho, T) increases on (0.8, 1) rho0, so the only crossing of P+/2 is at the jump. */
    {
        double Pprev = Pplus;
        for (int k = 1; k <= 100; k++) {
            double rho = rhoJ * (1.0 + 0.2499 * k / 100.0);
            double Pk = EOSPofRhoT(granite, rho, T);
            CHECK(Pk > Pprev, "P(rho, T) not increasing at rho=%g: %g after %g", rho, Pk, Pprev);
            Pprev = Pk;
        }
    }
    double P = 0.5 * Pplus;
    double rho1 = EOSRhoofPT(iron, P, T);
    CHECK(!EOSFailed(rho1) && rho1 > granite->rho0, "iron rho(P=%g, T=%g) = %g", P, T, rho1);
    if (!EOSFailed(rho1)) {
        double u1 = EOSUofRhoT(iron, rho1, T);
        CHECK_REL(EOSPofRhoU(iron, rho1, u1), P, 1e-8, "iron state P");
        double rho2 = -1.0, u2 = -1.0;
        int ret = EOSSolveBC(iron, granite, rho1, u1, &rho2, &u2, 0.0);
#ifdef TILL_PRESS_MELOSH
        CHECK(ret == EOS_FAIL, "SolveBC across a pressure jump must fail (ret=%d, rho2=%g)", ret, rho2);
        CHECK(rho2 == -1.0 && u2 == -1.0, "outputs touched on failure");
#else
        CHECK(ret == EOS_SUCCESS, "SolveBC on a continuous P2 failed (ret=%d)", ret);
        if (ret == EOS_SUCCESS) {
            CHECK(rho2 < granite->rho0, "rho2 = %g should be expanded", rho2);
            CHECK(fabs(EOSPofRhoU(granite, rho2, u2) - P) <= solveBCPResidualBound(granite, rho2, T, P),
                  "P continuity: %.15e vs %.15e", EOSPofRhoU(granite, rho2, u2), P);
        }
#endif
    }
    END_TEST();
}

/*
 * The Woolfson pressure threshold is a physical pressure: the same physical (P, T) gives the same
 * f_ij in every unit system, no correction just below WOOLFSON_MIN_PRESSURE_CGS and the correction
 * just above it. Up to EOSlib 1.2.0 it was 1e-3 in code units (1.44 GPa in Earth units, 1e-4 Pa
 * in cgs). The unit factors of EOSMATERIAL must agree with those of the sublibraries.
 */
static void testWoolfsonUnits(void)
{
    int failsBefore = nFails;
    TEST("EOSWoolfsonCoeff threshold in every unit system");
    const double units[3][2] = {
        {0.0, 0.0},                           /* cgs */
        {2.06701e-13, 4.80438e-08},           /* Earth radius and mass */
        {1.3e8 / 3.085678e21, 1.0},           /* one solar mass in 1300 km */
    };
    const double T = 300.0;
    double fAbove[3];
    for (int iu = 0; iu < 3; iu++) {
        EOSMATERIAL *granite = EOSinitMaterial(1, units[iu][0], units[iu][1], NULL);
        EOSMATERIAL *iron = EOSinitMaterial(2, units[iu][0], units[iu][1], NULL);
        struct igeosParam gp = {5.0 / 3.0, 1.0};
        EOSMATERIAL *gas = EOSinitMaterial(0, units[iu][0], units[iu][1], &gp);
        const double Punit = granite->dErgPerGmUnit * granite->dGmPerCcUnit;
        CHECK_REL(granite->dErgPerGmUnit, granite->tillmaterial->dErgPerGmUnit, 1e-15, "dErgPerGmUnit vs tillotson (units %d)", iu);
        CHECK_REL(granite->dGmPerCcUnit, granite->tillmaterial->dGmPerCcUnit, 1e-15, "dGmPerCcUnit vs tillotson (units %d)", iu);
        CHECK_REL(gas->dErgPerGmUnit, gas->igeosmaterial->dErgPerGmUnit, 1e-15, "dErgPerGmUnit vs igeos (units %d)", iu);
        CHECK_REL(gas->dGmPerCcUnit, gas->igeosmaterial->dGmPerCcUnit, 1e-15, "dGmPerCcUnit vs igeos (units %d)", iu);
        CHECK_REL(EOSWoolfsonMinPressure(granite) * Punit, WOOLFSON_MIN_PRESSURE_CGS, 1e-14, "threshold in cgs (units %d)", iu);
        CHECK(EOSWoolfsonCoeff(granite, iron, WOOLFSON_MIN_PRESSURE_CGS * (1.0 - 1e-6) / Punit, T) == 1.0, "no correction just below the threshold (units %d)", iu);
        fAbove[iu] = EOSWoolfsonCoeff(granite, iron, WOOLFSON_MIN_PRESSURE_CGS * (1.0 + 1e-6) / Punit, T);
        CHECK(!EOSFailed(fAbove[iu]) && fAbove[iu] > 0.0 && fAbove[iu] < 1.0, "correction just above the threshold (units %d): f= %g", iu, fAbove[iu]);
        /* The physical f at a high pressure, also independent of the units. */
        CHECK_REL(EOSWoolfsonCoeff(granite, iron, 1e11 / Punit, T),
                  EOSRhoofPT(granite, 1e11 / Punit, T) / EOSRhoofPT(iron, 1e11 / Punit, T), 1e-14, "f = rho1/rho2 at 10 GPa (units %d)", iu);
        EOSfinalizeMaterial(granite);
        EOSfinalizeMaterial(iron);
        EOSfinalizeMaterial(gas);
    }
    for (int iu = 1; iu < 3; iu++) CHECK_REL(fAbove[iu], fAbove[0], 1e-8, "f above the threshold in unit system %d", iu);
    END_TEST();
}

static void testTillotson(double dKpcUnit, double dMsolUnit)
{
    /* Iron (2) and granite (1): the two Tillotson materials with a unique
     * rho(P, T) in the condensed region; ice and water are not (intermediate region). */
    EOSMATERIAL *iron = EOSinitMaterial(2, dKpcUnit, dMsolUnit, NULL);
    EOSMATERIAL *granite = EOSinitMaterial(1, dKpcUnit, dMsolUnit, NULL);
    printf("  T1: %s  T2: %s\n", granite->MatString, iron->MatString);
    testTillotsonDispatch(granite);
    testTillotsonDispatch(iron);
    testTillotsonTableStatus(granite);
    testTillotsonTableStatus(iron);
    CHECK(iron->rho0 > granite->rho0, "iron must be denser than granite");
    testTillotsonSolveBC(iron, granite);
    testTillotsonSolveBCStiff(iron, granite);
    testTillotsonSolveBCJump(iron, granite);
    EOSfinalizeMaterial(iron);
    EOSfinalizeMaterial(granite);
    testWoolfsonUnits();
}
#endif

/* ------------------------------------------------------------------- main */

int main(int argc, char *argv[])
{
    int iMatA = 63, iMatB = 62;
    int iMatSCVH = 114;                /* 0: SCvH tests skipped */
    int iMatREOS3 = 104;               /* 0: REOS3 tests skipped */
    double dKpcUnit = 2.06701e-13;   /* Earth radius */
    double dMsolUnit = 4.80438e-08;  /* Earth mass */

    if (argc > 2) { iMatA = atoi(argv[1]); iMatB = atoi(argv[2]); }
    if (argc > 4) { dKpcUnit = atof(argv[3]); dMsolUnit = atof(argv[4]); }
    if (argc > 5) iMatSCVH = atoi(argv[5]);
    if (argc > 6) iMatREOS3 = atoi(argv[6]);

    printf("EOSlib %s test suite: materials %d and %d, dKpcUnit = %g, dMsolUnit = %g\n",
           EOS_VERSION_TEXT, iMatA, iMatB, dKpcUnit, dMsolUnit);

    {
        int failsBefore = nFails;
        TEST("library versions");
        CHECK(EOSCheckVersions() == EOS_TRUE, "EOSCheckVersions failed");
        END_TEST();
    }

    testIdealGas(dKpcUnit, dMsolUnit);

#ifdef HAVE_ANEOSMATERIAL_H
    EOSMATERIAL *A = EOSinitMaterial(iMatA, dKpcUnit, dMsolUnit, NULL);
    EOSMATERIAL *B = EOSinitMaterial(iMatB, dKpcUnit, dMsolUnit, NULL);
    EOSinitIsentropicLookup(A, NULL);
    EOSinitIsentropicLookup(B, NULL);
    printf("  A: %s  B: %s\n", A->MatString, B->MatString);

    testANEOSDispatch(A);
    testANEOSDispatch(B);
    testANEOSTableStatus(A);
    testANEOSTableStatus(B);
    testANEOSYield(A);
    testANEOSYield(B);
    {
        ANEOSMATERIAL *m = A->ANEOSmaterial;
        double rhoOk = sampleRho(m, NSAMPLE / 2), TMin = m->TAxis[0];
        double uOk = EOSUofRhoT(A, rhoOk, sampleT(m, NSAMPLE / 2));
        double uHi = 2.0 * fabs(EOSUofRhoT(A, rhoOk, m->TAxis[m->nT-1])) + 1.0;
        testFailureContract(A, rhoOk, uOk, uHi, 2.0 * m->rhoAxis[m->nRho-1] * m->CGStoCodeUnitsforRho, TMin, &m->nLookupFail);
    }
    /* Convention: the denser material is the inner one. */
    if (A->rho0 >= B->rho0) { testSolveBC(A, B); testWoolfson(B, A); }
    else                    { testSolveBC(B, A); testWoolfson(A, B); }

#ifdef HAVE_SCVHEOS_H
    if (iMatSCVH > 0) {
        EOSMATERIAL *S = EOSinitMaterial(iMatSCVH, dKpcUnit, dMsolUnit, NULL);
        printf("  S: %s\n", S->MatString);
        testSCVHDispatch(S);
        testSCVHTableStatus(S);
        {
            SCVHEOSMAT *s = S->scvheosmaterial;
            double rhoOk = scvhSampleRho(s, 1), TMin = pow(10.0, s->LogTMin);
            double uOk = EOSUofRhoT(S, rhoOk, scvhSampleT(s, 1));
            double uHi = 2.0 * EOSUofRhoT(S, rhoOk, pow(10.0, s->LogTMax));
            testFailureContract(S, rhoOk, uOk, uHi, 2.0 * pow(10.0, s->LogRhoMax), TMin, &s->nLookupFail);
        }
        if (A->rho0 >= B->rho0) testSCVHSolveBC(A, S); else testSCVHSolveBC(B, S);
        EOSfinalizeMaterial(S);
    }
#endif

#ifdef HAVE_REOS3_H
    if (iMatREOS3 > 0) {
        EOSMATERIAL *R = EOSinitMaterial(iMatREOS3, dKpcUnit, dMsolUnit, NULL);
        printf("  R: %s\n", R->MatString);
        testREOS3Dispatch(R);
        testREOS3Extrapolation(R);
        testREOS3TableStatus(R);
        {
            REOS3MAT *r = R->reos3material;
            double rhoOk = reos3SampleRho(r, 1), TMin = pow(10.0, r->LogTMin);
            double uOk = EOSUofRhoT(R, rhoOk, reos3SampleT(r, 1));
            double uHi = 2.0 * EOSUofRhoT(R, rhoOk, pow(10.0, r->LogTMax));
            testFailureContract(R, rhoOk, uOk, uHi, 2.0 * pow(10.0, r->LogRhoMax), TMin, &r->nLookupFail);
        }
        if (A->rho0 >= B->rho0) testREOS3SolveBC(A, R); else testREOS3SolveBC(B, R);
        EOSfinalizeMaterial(R);
    }
#endif

    EOSfinalizeMaterial(A);
    EOSfinalizeMaterial(B);
#else
    printf("  (ANEOS tests skipped: HAVE_ANEOSMATERIAL_H not defined)\n");
#endif

#ifdef HAVE_TILLOTSON_H
    testTillotson(dKpcUnit, dMsolUnit);
#endif

    printf("%d checks, %d failed\n", nChecks, nFails);
    return nFails;
}