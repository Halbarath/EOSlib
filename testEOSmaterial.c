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
 * are in the working directory; the Tillotson tests run only if
 * HAVE_TILLOTSON_H is defined. With HAVE_SCVHEOS_H / HAVE_REOS3_H and a
 * non-zero iMatSCVH / iMatREOS3 the SCvH and REOS3 paths are tested too
 * (0 skips them), including the boundary condition solver for a gas envelope
 * on the denser ANEOS material; the gas tables must be in the working
 * directory. Exit status is the number of failed checks.
 *
 * The REOS3 tests only visit states inside the lookup limits (the library
 * exits on out-of-range input unless built with -DREOS3_NOEXIT, which this
 * test does not need).
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

#define TEST(name) do { currentTest = name; printf("  %-44s", name); fflush(stdout); } while (0)
#define END_TEST() do { printf("%s\n", nFails == failsBefore ? "ok" : "FAILED"); } while (0)

static int isSentinel(double v) { return v < -1e40; }

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
        if (P1 <= WOOLFSON_MIN_PRESSURE) continue;
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
            /* The solver guarantees the pressure bracket closed to EOS_SOLVEBC_PTOL. */
            CHECK_REL(EOSPofRhoU(inner, rho2, u2), P1, EOS_SOLVEBC_PTOL, "P continuity rho1=%g T1=%g jump=%g", rho1, T1, Tjump);
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
        CHECK(ret == -1 || (isfinite(rho2) && rho2 > 0.0), "SolveBC: bad return %d with rho2=%g", ret, rho2);
        if (ret == -1) CHECK(rho2 == -1.0 && u2 == -1.0, "outputs touched on failure");
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
    CHECK(EOSWoolfsonCoeff(a, b, 0.5 * WOOLFSON_MIN_PRESSURE, T) == 1.0, "below minimum pressure");
    double f = EOSWoolfsonCoeff(a, b, P, T);
    CHECK(f == -1.0 || (f > 0.0 && isfinite(f)), "f_ij = %g", f);
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
        if (P1 <= WOOLFSON_MIN_PRESSURE) continue;
        double rho2, u2;
        if (EOSSolveBC(outer, inner, rho1, u1, &rho2, &u2, 0.0) != 0) continue;
        nSolved++;
        CHECK_REL(EOSPofRhoU(inner, rho2, u2), P1, EOS_SOLVEBC_PTOL, "P continuity rho1=%g T1=%g", rho1, T1);
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
    CHECK(EOSSofRhoU(e, reos3SampleRho(r, 1), 1.0) == -1e50, "entropy sentinel (rho,u)");
    CHECK(EOSSofRhoT(e, reos3SampleRho(r, 1), reos3SampleT(r, 1)) == -1e50, "entropy sentinel (rho,T)");
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
        if (P1 <= WOOLFSON_MIN_PRESSURE) continue;
        double rho2, u2;
        if (EOSSolveBC(outer, inner, rho1, u1, &rho2, &u2, 0.0) != 0) continue;
        nSolved++;
        CHECK_REL(EOSPofRhoU(inner, rho2, u2), P1, EOS_SOLVEBC_PTOL, "P continuity rho1=%g T1=%g", rho1, T1);
        CHECK_REL(u2, EOSUofRhoT(inner, rho2, T1), 1e-10, "T continuity rho1=%g T1=%g", rho1, T1);
        CHECK(EOSIsUsable(inner, rho2, u2) == EOS_TRUE, "solved state not usable");
    }
    if (nSolved == 0) printf("\n    note: no envelope state had a solution in the core table (table ranges do not overlap)\n  %-44s", "");
    END_TEST();
}
#endif

#endif /* HAVE_REOS3_H */

/* ------------------------------------------------------------- Tillotson */

#ifdef HAVE_TILLOTSON_H
static void testTillotsonPolicy(double dKpcUnit, double dMsolUnit)
{
    int failsBefore = nFails;
    TEST("Tillotson: usable above the lookup table");
    EOSMATERIAL *t = EOSinitMaterial(1, dKpcUnit, dMsolUnit, NULL);   /* granite */
    CHECK(t->matType == EOSTILLOTSON, "matType");
    CHECK(!t->bEntropy, "Tillotson has no entropy");
    CHECK(EOSSofRhoU(t, t->rho0, 1.0) == -1e50, "entropy sentinel");
    double rho = t->rho0;
    double u = 1e4 * EOSUofRhoT(t, rho, 300.0);   /* far above any lookup table */
    int st = EOSIsInTable(t, rho, u);
    if (st == EOS_OUTSIDE_VMAX) CHECK(EOSIsUsable(t, rho, u) == EOS_TRUE, "OUTSIDE_VMAX must be usable for Tillotson");
    CHECK(isfinite(EOSPofRhoU(t, rho, u)), "P above table is finite");
    EOSfinalizeMaterial(t);
    END_TEST();
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
        CHECK(EOSCheckVersions() == 1, "EOSCheckVersions failed");
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
    /* Convention: the denser material is the inner one. */
    if (A->rho0 >= B->rho0) { testSolveBC(A, B); testWoolfson(B, A); }
    else                    { testSolveBC(B, A); testWoolfson(A, B); }

#ifdef HAVE_SCVHEOS_H
    if (iMatSCVH > 0) {
        EOSMATERIAL *S = EOSinitMaterial(iMatSCVH, dKpcUnit, dMsolUnit, NULL);
        printf("  S: %s\n", S->MatString);
        testSCVHDispatch(S);
        testSCVHTableStatus(S);
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
    testTillotsonPolicy(dKpcUnit, dMsolUnit);
#endif

    printf("%d checks, %d failed\n", nChecks, nFails);
    return nFails;
}