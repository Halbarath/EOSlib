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
 *   testEOSmaterial [iMatA iMatB] [dKpcUnit dMsolUnit]
 *
 * Defaults: ANEOS dunite (54) and ANEOS iron (55), Earth-like code units.
 * iMatA/iMatB must be ANEOS materials whose tables are in the working
 * directory; the Tillotson tests run only if HAVE_TILLOTSON_H is defined.
 * Exit status is the number of failed checks.
 */
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
    /* Just above the exclusive upper bound, with enough margin to survive the
     * code<->cgs round trip (the exact node is one ulp ambiguous). */
    double rhoEdge = m->rhoAxis[m->nRho-1] * (1.0 + 1e-12) * m->CGStoCodeUnitsforRho;

    CHECK(EOSIsInTable(e, rhoLo, uOk) == EOS_OUTSIDE_RHOMIN, "below rhoMin");
    CHECK(EOSIsInTable(e, rhoHi, uOk) == EOS_OUTSIDE_RHOMAX, "above rhoMax");
    CHECK(EOSIsInTable(e, rhoEdge, uOk) == EOS_OUTSIDE_RHOMAX, "just above rhoMax");
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
    double dKpcUnit = 2.06701e-13;   /* Earth radius */
    double dMsolUnit = 4.80438e-08;  /* Earth mass */

    if (argc > 2) { iMatA = atoi(argv[1]); iMatB = atoi(argv[2]); }
    if (argc > 4) { dKpcUnit = atof(argv[3]); dMsolUnit = atof(argv[4]); }

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
