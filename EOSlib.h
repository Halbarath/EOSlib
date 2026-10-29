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
 *   1.2.0  Compared with 1.0.9. Requires tillotson 4.0.0, ANEOSmaterial 1.2.0, reos3 1.3.0,
 *          scvheos 1.3.0 and igeos 1.1.0, each with the functions of this round named below
 *          (raise the *_VERSION_*_REQUIRED values below when the sublibrary versions are bumped).
 *
 *          API changes (callers must check):
 *          - EOS_SUCCESS is 0 and EOS_FAIL is -99 (were 1 and 0): test "== EOS_SUCCESS", never
 *            the truth value. EOSSolveBC returns EOS_SUCCESS/EOS_FAIL ("!= 0" still means
 *            failure). EOSWoolfsonCoeff returns EOS_FAIL_VALUE (-1e50) instead of -1.0 (still
 *            negative).
 *          - EOSCheckVersions is removed: the sublibrary versions are checked at compile time
 *            (#error in this header; callers drop the call). Its static "checked" flag in
 *            EOSinitMaterial was also unsynchronised shared state between threads.
 *          - Tillotson: EOSinitMaterial builds only the cold curve and the top isentrope of the
 *            standard table (77 kB, ~1 ms, tillInitColdCurveDefault); the full isentrope table
 *            (46 MB, ~0.6 s) is built by EOSinitIsentropicLookup, and EOSIsentropic is fatal
 *            without it (callers that use EOSIsentropic must add the call; bEntropyTableInit is
 *            EOS_FALSE until then). The table is the library's standard table (TILL_TABLE_*:
 *            rho from 1e-4 rho0 to 100.58 rho0, the first node at or above 100 rho0 so that rho0
 *            is a node; 1 K..1e7 K at rho0; 1200 x 1200), the same physical table for every
 *            material and unit system. Before 1.2.0 the table was fixed in code units, rho in
 *            1e-4..200 and v = u(rho0) up to 1200, and was built twice in EOSinitMaterial for every
 *            instance. Its physical range therefore depended on the material and the unit system:
 *            with dKpcUnit = 2.06701e-13, dMsolUnit = 4.80438e-08 (code density 0.3685 g/cc) it
 *            covered 1.4e-5..27 rho0 for granite and 4.7e-6..9.4 rho0 for iron, up to 1.5e6 K at
 *            rho0 for granite; in cgs 3.7e-5..74 rho0 for granite but only up to 1.5e-4 K. Results
 *            of states inside both tables differ by the interpolation error of the two grids. The
 *            material with the cold curve only gives bit-identical results to the full table.
 *          - EOS_N_MATERIAL_MAX is the size of an array indexed by iMat, MAT_SCVHEOS_MAX + 1 = 131
 *            (it was 130, the largest id: such an array overflowed by one), like
 *            TILL_N_MATERIAL_MAX. Loops over ids must use iMat < EOS_N_MATERIAL_MAX.
 *          - The Woolfson pressure threshold is a physical pressure, WOOLFSON_MIN_PRESSURE_CGS =
 *            1.44 GPa, converted per material (EOSWoolfsonMinPressure); WOOLFSON_MIN_PRESSURE
 *            (1e-3 in code units: 1.44 GPa in Earth units but 1e-4 Pa in cgs) is gone.
 *          - One rule for the unit pair in EOSlib and every sublibrary: both dKpcUnit and
 *            dMsolUnit > 0 (and finite) is code units, both <= 0 is cgs, anything else is fatal
 *            (igeos divided by zero, Tillotson only asserted, ANEOS, SCvH and REOS3 silently used
 *            cgs). EOSMATERIAL records the unit factors dErgPerGmUnit and dGmPerCcUnit.
 *          - No GSL, no brent.o.
 *
 *          Failure contract (see "Failure contract" below), one for EOSlib and all sublibraries:
 *          - A failed evaluation never exits; it returns EOS_FAIL_VALUE. EOSTofRhoU never fails
 *            (T_min of the table, 0 for Tillotson and the ideal gas). EOSYieldStrength returns
 *            EOS_FAIL_VALUE for an invalid T or T_melt instead of exiting.
 *          - nLookupFail counts every failure and fallback once per call: ANEOS P, c, s, phase
 *            and the (rho, u) derivatives counted a density outside the table twice, the ANEOS
 *            and SCvH isentropes up to three times, EOSSofRhoU for SCvH twice (it uses the new
 *            scvheosSofRhoU). EOSPrintMat prints nLookupFail for every library (only reos3 did).
 *          - The Tillotson exception is named: P, c and the (rho, u) derivatives are analytic and
 *            evaluated at the given (rho, u) also below the cold curve and outside the table,
 *            uncounted, for a valid density.
 *          - Invalid density (rho <= 0, NaN, inf) for the analytic EOS: Tillotson and the ideal
 *            gas fail like a density outside a table (P, c and the derivatives EOS_FAIL_VALUE,
 *            T = 0, counted once; EOSUCold of the gas EOS_FAIL_VALUE via the new igeosUCold,
 *            EOSisbelowColdCurve EOS_FAIL; EOSIsInTable for the gas EOS_OUTSIDE_RHOMIN, or
 *            EOS_OUTSIDE_RHOMAX for inf). The gas returned P < 0 for rho < 0, EOSdUdRho = NaN at
 *            rho = 0 and EOSUCold = 0 for any rho; Tillotson returned NaN, uncounted.
 *            EOSPCTofRhoU for Tillotson counts an invalid rho once (it counted it twice).
 *          - EOSdUdRho returns EOS_FAIL_VALUE where P(rho, u) fails; for ANEOS and SCvH it
 *            returned -1e50/rho^2 (not the sentinel; for rho > 1e5 not even EOSFailed), for the
 *            gas NaN at rho = 0.
 *          - Output: per-call diagnostics only with -DEOSLIB_VERBOSE, and EOSinitMaterial and
 *            EOSinitIsentropicLookup are silent without it (reos3 printed the table name, uRef and
 *            the unit conversion, reos3 and scvheos their monotonicity reports, sound speed
 *            repairs and the note that rho0 lies outside the table, EOSinitIsentropicLookup
 *            "not implemented" for REOS3 on every call; with one instance per thread all of it
 *            appeared 64 times). Errors before a fatal exit are still printed.
 *          - EOSPrintMat prints every dimensional value in cgs with its unit (g/cc, K, erg/g,
 *            ...) for every library: igeos and Tillotson printed code units without a label,
 *            ANEOS, SCvH and REOS3 bare numbers, SCvH its table boundaries as log10 values.
 *
 *          Ideal gas (igeos):
 *          - igeosInitMat rejects (fatal) dConstGamma <= 1 and dMeanMolMass <= 0 or not finite.
 *          - u < 0 lies below its cold curve u = 0 and falls back to T = 0, P = 0, c = 0
 *            (counted), like a state below the cold curve of a table; it gave T < 0, P < 0 and
 *            c = NaN. EOSIsInTable reports EOS_OUTSIDE_VMIN there and EOSIsUsable is false.
 *          - EOSRhoofPT returns EOS_FAIL_VALUE for T <= 0 or P < 0 (inf before). EOSPCTofRhoU
 *            uses the new igeosPCTofRhoU (one count per call).
 *          - EOSdPdRhoofRhoT and EOSdPdT: (gamma-1) cv T and (gamma-1) rho cv (igeosdPdRhoofRhoT,
 *            igeosdPdTofRhoT); they were fatal (not implemented).
 *
 *          Table status and cold curve:
 *          - EOSlib no longer reads table internals of the sublibraries: EOSIsInTable uses
 *            ANEOSIsInTable, EOSUCold uses ANEOSUCold/reos3UCold/scvheosUCold/igeosUCold, the
 *            init uses ANEOSTmin. ANEOSmaterial, scvheos and reos3 vendor interpBilinear 1.0.3;
 *            their headers check at compile time that the copy in effect in this translation
 *            unit has their type layout (INTERP_ABI_VERSION), so EOSlib and the libraries cannot
 *            disagree on the embedded INTERPGRID.
 *          - EOSIsInTable/EOSIsUsable for ANEOS agree exactly with the lookups (ANEOSIsInTable
 *            runs the same inversion): about a quarter of the states near the cold curve of the
 *            M-ANEOS tables, which the lookup finds, were reported EOS_OUTSIDE_VMIN and unusable.
 *          - EOSIsUsable accepts Tillotson states below rho_min (rho > 0, u >= 0 within
 *            roundoff): P, c and T are defined there and the isentrope is integrated directly.
 *          - EOSisbelowColdCurve returns EOS_FAIL where the cold curve cannot be evaluated or u
 *            is not finite (it answered EOS_FALSE). For ANEOS, REOS3 and SCvH it is the table
 *            status of the lookup (EOSIsInTable == EOS_OUTSIDE_VMIN): comparing u with
 *            u(rho, T_min) said "below" for M-ANEOS states below the first row that the lookup
 *            finds without a fallback (columns with du/dT < 0 at the cold rows), and for states
 *            within the edge tolerance of the first row. For Tillotson it is
 *            tillIsBelowColdCurve, u < u_cold - tillUTol(u), the tolerance of tillIsInTable,
 *            tillTempRhoU and EOSIsUsable (u < u_cold said "below" within roundoff of the cold
 *            curve). For the ideal gas u < 0.
 *
 *          EOSSolveBC:
 *          - A discontinuous P2(rho, T) is detected by the growth of the secant slope as the
 *            bracket closes (EOS_SOLVEBC_SLOPE_WIDTH, EOS_SOLVEBC_SLOPE_FAC) instead of the
 *            pressure gap alone, which rejected continuous stiff states (K/P > 1e4, rock below
 *            ~10 MPa). Where the slope test fails, the one-sided secant slopes just outside the
 *            final bracket give a second opinion (EOS_SOLVEBC_PROBE): a root within
 *            EOS_SOLVEBC_SLOPE_WIDTH of a kink of a continuous P2 (Tillotson at rho0 near T = 0:
 *            P cut off at 0 below, the stiff cold curve above) was rejected whenever the
 *            reference bracket straddled the kink.
 *          - It fails at once for a failed state in material 1 (P = EOS_FAIL_VALUE; isfinite(P)
 *            let the sentinel through as the target) and for a state that is not usable in
 *            material 1 (EOSIsUsable: above the column or below the cold curve P and T are those
 *            of the fallback, and the solver matched that wrong target and reported success).
 *          - It accepts a target temperature T2 = T1 - Tjump = 0 (the cold curve of Tillotson and
 *            the ideal gas); only T2 < 0 is rejected.
 *          - Bracket search: where material 2 cannot be evaluated at rho1 < rho0 of material 2,
 *            the upper bracket starts at that rho0 (it failed at once for a gas envelope, SCvH or
 *            REOS3, below the density range of a condensed core, ANEOS). Where the steps of 1.1
 *            and 0.9 do not find a bracket, up to EOS_SOLVEBC_MAX_COARSE steps of a factor 2
 *            follow (EOS_SOLVEBC_EXPAND_FAC), and a trial beyond the table of material 2 refines
 *            the step instead of failing. Wherever the old search succeeded, the brackets and so
 *            the result are bit-identical (factor-2 steps from the start ended on another root
 *            of a non-monotonic P2, M-ANEOS iron on forsterite, and ballic's integration no
 *            longer converged).
 *          - The last bisection lookup is not repeated at the end (one evaluation of material 2
 *            per call less, bit-identical).
 *
 *          EOSIsentropic:
 *          - For ANEOS and SCvH the clamp to the table edge the isentrope leaves through is
 *            documented (ANEOS clamped to the cold curve on both sides). A state above the u
 *            column starts from the top row s(rho1, T_max), not from the T_min fallback, which put
 *            u2 on the cold curve at rho2. For rho1 outside the table ANEOS returns
 *            EOS_FAIL_VALUE, like SCvH (ANEOSisentropicU returned the cold curve at rho2).
 *          - Tillotson still integrates the ODE outside its table: its pressure is analytic, so
 *            the isentrope exists there.
 *
 *          Strength:
 *          - EOSYieldModel documents -2 (Tillotson: no strength model). EOSYieldStrength gathers
 *            the yield parameters once per call (EOSYieldModel and the core each called
 *            ANEOSYieldParameters).
 *          - EOSYieldStrength and EOSYieldStrengthofRhoTP return EOS_FAIL_VALUE for a failed
 *            pressure (models 1 and 2): the tension clamp max(P, 0) turned P = EOS_FAIL_VALUE into
 *            P_eff = 0 and a finite Y0 f(T) for a state outside the table.
 *          - EOSYieldStrength evaluates a state above the u column (overheated material) on the
 *            top row T_max, P(rho, T_max), i.e. without strength wherever T_melt <= T_max; with
 *            the T_min fallback it had the full cold strength Y0 f(T_min). The test costs one
 *            inversion of the u column only when T == T_min. EOSYieldStrengthofRhoTP evaluates the
 *            T, P it is given (from T alone the fallback cannot be recognised).
 *
 *          Documentation only: rho0forInterfaceCorrection is not read by EOSlib (EOSWoolfsonCoeff
 *          cuts at 0.8*rho0 for every material, the gases included).
 *   1.3.0  Compared with 1.2.0. Requires ANEOSmaterial 1.3.0 (ANEOSTmax), scvheos 1.3.1 and reos3
 *          1.3.1 (the status of non-finite input), all three with interpBilinear 1.0.4; tillotson
 *          4.0.0 and igeos 1.1.0 as before.
 *          - The Woolfson pressure threshold WOOLFSON_MIN_PRESSURE_CGS is 3.684e6 dyn/cm^2
 *            (0.37 MPa): 1e-3 in the pkdgrav3/ballic units dKpcUnit = 2.06701e-13, dMsolUnit =
 *            4.80438e-08 (9.56e25 g, not the Earth mass), i.e. the threshold that runs in these
 *            units had before 1.2.0 (to roundoff). 1.2.0 used 1.44 GPa, the old code-unit threshold
 *            of true Earth units (R_E = 6.371e8 cm, M_E = 5.972e27 g), which raised it by a factor
 *            3900 for the pkdgrav3 units. Runs in cgs (1e-4 Pa before 1.2.0) and in true Earth units
 *            use 0.37 MPa as well.
 *          - EOSIsInTable reports non-finite input by the side it lies on, the same for every
 *            library: rho <= 0 and NaN EOS_OUTSIDE_RHOMIN, rho = +inf EOS_OUTSIDE_RHOMAX, u = NaN
 *            and -inf EOS_OUTSIDE_VMIN, u = +inf EOS_OUTSIDE_VMAX. EOSIsUsable is EOS_FALSE for all
 *            of them. Before: the ideal gas reported u = +inf as EOS_SUCCESS (and usable), SCvH and
 *            REOS3 reported +inf as RHOMIN and VMIN, Tillotson a non-finite u as EOS_FAIL, and ANEOS
 *            u = +-inf as EOS_SUCCESS (found at T_min by the node tolerance of the tolerant inversion,
 *            fixed in interpBilinear 1.0.4); T(rho, +-inf) fell back to T_min uncounted for ANEOS.
 *            EOS_FAIL is left for an unknown sublibrary status.
 *          - EOSYieldStrength clamps a state above the ANEOS table with ANEOSTmin and ANEOSTmax: it
 *            read the T axis of the table, contrary to the 1.2.0 entry (EOSlib no longer reads
 *            sublibrary internals). Same results.
 *          - EOSSolveBC: each evaluation of material 2 skips EOSIsUsable for ANEOS and Tillotson,
 *            where a u from a successful u(rho, T) is always usable (testANEOSmaterial and
 *            testtillotson check it): for ANEOS one of the two u-column inversions per step.
 *            Bit-identical.
 *          - EOSisbelowColdCurve for Tillotson does one cold-curve lookup per call (it did two,
 *            EOSUCold and the one inside tillIsBelowColdCurve). Same results and counts.
 *   1.4.0  Compared with 1.3.0. Requires reos3 1.4.0 (reos3InitIsentropicLookup,
 *          reos3IsentropicU); the other sublibraries as before.
 *          - EOSIsentropic for REOS3 (it was fatal, "not implemented"): REOS3 has no entropy, so
 *            reos3 integrates du/dlogrho = P(rho, u)/rho of its interpolant into an isentrope table,
 *            as Tillotson does (see "Isentrope table" in reos3.h). EOSinitIsentropicLookup builds it
 *            (~2 s, 16 MB for REOS3_H_HE) and sets bEntropyTableInit; EOSIsentropic without it is
 *            fatal, as for Tillotson (callers that use EOSIsentropic with REOS3 must add the call).
 *            The clamps are those of ANEOS and SCvH (below the cold curve: the isentrope through
 *            (rho1, T_min); above the top row: through (rho1, T_max); leaving the table: its edge
 *            at rho2; counted once). bEntropy stays EOS_FALSE (EOSSofRhoU/EOSSofRhoT unchanged).
 *          - EOSinitIsentropicLookup no longer prints a note for REOS3.
 *          - EOS_VERSION_PATCH was 1 for EOS_VERSION_TEXT "1.3.0"; both are 1.4.0 now.
 */

#ifndef EOSLIB_HINCLUDED
#define EOSLIB_HINCLUDED

#include <stdio.h>   /* FILE in EOSPrintMat */

/*
 * Minimum versions of the sublibraries. They are compile-time constants on both sides, so they
 * are checked here with #error (as the sublibraries check interpBilinear) instead of at run time
 * (EOSCheckVersions did that before EOSlib 1.2.0 and has been removed).
 */
#ifdef HAVE_TILLOTSON_H
    #include <tillotson.h>
#define TILL_VERSION_MAJOR_REQUIRED 4
#define TILL_VERSION_MINOR_REQUIRED 0
#define TILL_VERSION_PATCH_REQUIRED 1
#if (TILL_VERSION_MAJOR * 10000 + TILL_VERSION_MINOR * 100 + TILL_VERSION_PATCH) < \
    (TILL_VERSION_MAJOR_REQUIRED * 10000 + TILL_VERSION_MINOR_REQUIRED * 100 + TILL_VERSION_PATCH_REQUIRED)
#error "EOSlib: the Tillotson library (tillotson.h) is older than the version EOSlib requires (TILL_VERSION_*_REQUIRED)"
#endif
#endif
#ifdef HAVE_ANEOSMATERIAL_H
    #include <ANEOSmaterial.h>
#define ANEOS_VERSION_MAJOR_REQUIRED 1
#define ANEOS_VERSION_MINOR_REQUIRED 3
#define ANEOS_VERSION_PATCH_REQUIRED 0
#if (ANEOS_VERSION_MAJOR * 10000 + ANEOS_VERSION_MINOR * 100 + ANEOS_VERSION_PATCH) < \
    (ANEOS_VERSION_MAJOR_REQUIRED * 10000 + ANEOS_VERSION_MINOR_REQUIRED * 100 + ANEOS_VERSION_PATCH_REQUIRED)
#error "EOSlib: the ANEOS library (ANEOSmaterial.h) is older than the version EOSlib requires (ANEOS_VERSION_*_REQUIRED)"
#endif
#endif
#ifdef HAVE_REOS3_H
    #include <reos3.h>
#define REOS3_VERSION_MAJOR_REQUIRED 1
#define REOS3_VERSION_MINOR_REQUIRED 4
#define REOS3_VERSION_PATCH_REQUIRED 0
#if (REOS3_VERSION_MAJOR * 10000 + REOS3_VERSION_MINOR * 100 + REOS3_VERSION_PATCH) < \
    (REOS3_VERSION_MAJOR_REQUIRED * 10000 + REOS3_VERSION_MINOR_REQUIRED * 100 + REOS3_VERSION_PATCH_REQUIRED)
#error "EOSlib: the REOS3 library (reos3.h) is older than the version EOSlib requires (REOS3_VERSION_*_REQUIRED)"
#endif
#endif
#ifdef HAVE_SCVHEOS_H
    #include <scvheos.h>
#define SCVHEOS_VERSION_MAJOR_REQUIRED 1
#define SCVHEOS_VERSION_MINOR_REQUIRED 3
#define SCVHEOS_VERSION_PATCH_REQUIRED 1
#if (SCVHEOS_VERSION_MAJOR * 10000 + SCVHEOS_VERSION_MINOR * 100 + SCVHEOS_VERSION_PATCH) < \
    (SCVHEOS_VERSION_MAJOR_REQUIRED * 10000 + SCVHEOS_VERSION_MINOR_REQUIRED * 100 + SCVHEOS_VERSION_PATCH_REQUIRED)
#error "EOSlib: the SCvH library (scvheos.h) is older than the version EOSlib requires (SCVHEOS_VERSION_*_REQUIRED)"
#endif
#endif

#include "igeos.h"
#define IGEOS_VERSION_MAJOR_REQUIRED 1
#define IGEOS_VERSION_MINOR_REQUIRED 2
#define IGEOS_VERSION_PATCH_REQUIRED 0
#if (IGEOS_VERSION_MAJOR * 10000 + IGEOS_VERSION_MINOR * 100 + IGEOS_VERSION_PATCH) < \
    (IGEOS_VERSION_MAJOR_REQUIRED * 10000 + IGEOS_VERSION_MINOR_REQUIRED * 100 + IGEOS_VERSION_PATCH_REQUIRED)
#error "EOSlib: the ideal gas library (igeos.h) is older than the version EOSlib requires (IGEOS_VERSION_*_REQUIRED)"
#endif

#define EOS_VERSION_TEXT    "1.4.0"
#define EOS_VERSION_MAJOR   1
#define EOS_VERSION_MINOR   4
#define EOS_VERSION_PATCH   0

/*
 * EOSWoolfsonCoeff does no correction below this pressure (cgs, dyn/cm^2 = 0.1 Pa). It is a
 * physical pressure, converted with the unit system of the material (EOSWoolfsonMinPressure).
 * Before EOSlib 1.2.0 the threshold was WOOLFSON_MIN_PRESSURE = 1e-3 in code units. The value
 * here is that threshold in the pkdgrav3/ballic units dKpcUnit = 2.06701e-13 (Earth radius),
 * dMsolUnit = 4.80438e-08 (9.56e25 g, not the Earth mass): energy unit 9.998e9 erg/g, density
 * unit 0.3685 g/cc, pressure unit 3.684e9 dyn/cm^2, so 1e-3 = 3.684e6 dyn/cm^2 = 0.37 MPa and
 * runs in these units behave as before 1.2.0. In other unit systems the old code-unit threshold
 * meant a different pressure (1e-4 Pa in cgs, 1.44 GPa in true Earth units R_E, M_E); all of them
 * use 0.37 MPa now. EOSlib 1.2.0 used 1.44 GPa (see the version history).
 */
#define WOOLFSON_MIN_PRESSURE_CGS 3.684121632893488e6

/* EOSSolveBC: bracket expansion limits and convergence tolerances. */
#define EOS_SOLVEBC_MAX_EXPAND  100   /* bracket steps of 1.1 up, 0.9 down (as before 1.2.0): 1.1^100 ~ 1.4e4, 0.9^100 ~ 2.7e-5 */
#define EOS_SOLVEBC_EXPAND_FAC  2.0   /* then steps of x2 up, /2 down ... */
#define EOS_SOLVEBC_MAX_COARSE  100   /* ... at most this many more (2^100 ~ 1e30); trials refined at a table edge count too */
#define EOS_SOLVEBC_MAX_BISECT  200
#define EOS_SOLVEBC_TOL         1e-10 /* relative, in density */
#define EOS_SOLVEBC_PTOL        1e-6  /* relative, in pressure: a gap below this always passes */
#define EOS_SOLVEBC_SLOPE_WIDTH 1e-6  /* relative bracket width of the reference secant slope */
#define EOS_SOLVEBC_SLOPE_FAC   10.0  /* the final secant slope may exceed the reference by this */
#define EOS_SOLVEBC_PROBE       4.0   /* width of the one-sided slopes outside the final bracket, in bracket widths */

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

/*
 * Size of an array indexed by the material id (iMat = 0 .. MAT_SCVHEOS_MAX), i.e. the largest id
 * plus one, the convention of TILL_N_MATERIAL_MAX (10 for ids up to 7). It was 130, the largest id
 * itself, so that an array of EOS_N_MATERIAL_MAX entries indexed by iMat overflowed by one for the
 * SCvH id 130. Loop over ids with iMat < EOS_N_MATERIAL_MAX.
 */
#define EOS_N_MATERIAL_MAX (MAT_SCVHEOS_MAX + 1)
#if EOS_N_MATERIAL_MAX <= MAT_TILLOTSON_MAX || EOS_N_MATERIAL_MAX <= MAT_ANEOS_MAX || \
    EOS_N_MATERIAL_MAX <= MAT_REOS3_MAX || EOS_N_MATERIAL_MAX <= MAT_SCVHEOS_MAX
#error "EOSlib: EOS_N_MATERIAL_MAX must exceed every material id (it sizes arrays indexed by iMat)"
#endif
#if defined(ANEOS_N_MATERIAL_MAX) && (ANEOS_N_MATERIAL_MAX <= MAT_ANEOS_MAX)
#error "EOSlib: ANEOS_N_MATERIAL_MAX must exceed MAT_ANEOS_MAX (it sizes arrays indexed by the ANEOS iMat)"
#endif
#if defined(TILL_N_MATERIAL_MAX) && (TILL_N_MATERIAL_MAX <= FORSTERITE)
#error "EOSlib: TILL_N_MATERIAL_MAX must exceed the largest Tillotson material id"
#endif

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
 * temperature of the table (T_min of the tabulated EOS, 0 for Tillotson and
 * the ideal gas) when (rho, u) cannot be inverted; the hydro codes rely on
 * this. A state below the cold curve (u < u(rho, T_min); u < 0 for the ideal
 * gas) is evaluated at (rho, T_min) by every (rho, u) function. An invalid
 * density (rho <= 0, NaN, inf) is outside every table, also of the analytic
 * EOS (ideal gas, Tillotson): EOS_FAIL_VALUE, T(rho, u) = T_min, EOSIsInTable
 * EOS_OUTSIDE_RHOMIN (EOS_OUTSIDE_RHOMAX for inf), EOSisbelowColdCurve
 * EOS_FAIL. A non-finite u at a valid density is outside as well: EOSIsInTable
 * EOS_OUTSIDE_VMIN for NaN and -inf, EOS_OUTSIDE_VMAX for +inf, and T(rho, u)
 * falls back to T_min. Status functions return EOS_FAIL or an EOS_OUTSIDE_* code.
 *
 * Exception, Tillotson: its pressure is an analytic function of (rho, u),
 * so EOSPofRhoU, EOSCofRhoU, EOSPCofRhoU, the P and c of EOSPCTofRhoU,
 * EOSdPdRho, EOSdPdU and EOSdUdRho are evaluated at the given (rho, u)
 * everywhere: below the cold curve (u < 0 enters the Grueneisen term as
 * Gamma = a + b), below rho_min and above rho_max. For a valid density
 * (rho > 0 and finite) they never fail, never fall back and are not counted;
 * an invalid density (rho <= 0, NaN, inf) returns EOS_FAIL_VALUE, counted,
 * as for every library. Only the functions that need the cold curve
 * fall back or fail there (T(rho, u) gives 0, counted; u(rho, T), P(rho, T),
 * rho(P, T) and u_cold give EOS_FAIL_VALUE above rho_max, counted), so below
 * the cold curve P(rho, u) is not P(rho, T(rho, u)). Whether P, c and T of a
 * state belong together is what EOSIsUsable reports: EOS_FALSE below the
 * cold curve and above rho_max, for Tillotson as for every other library.
 *
 * Every failure and every fallback is counted in the sublibrary's
 * nLookupFail (printed by EOSPrintMat for every library), once per call: a
 * call that fails or falls back counts exactly 1, however many of its
 * internal steps do (a density outside the table is one bad state, not a
 * T fallback plus a failed lookup). It is reported on stderr only when the
 * libraries are compiled with -DEOSLIB_VERBOSE (EOS_WARN and the *_WARN of
 * each sublibrary). Test a returned value with EOSFailed(x), which is true
 * for the sentinel and for NaN.
 *
 * Fatal (message on stderr, then assert(0); exit(1)): initialisation errors
 * (unknown material, missing or defective table, version mismatch) and
 * calling a function that is not implemented for the material (EOSRhoofUT
 * for all but ANEOS, EOSdPdRhoofRhoT/EOSdPdT for Tillotson,
 * EOSGammaofRhoT/EOSYieldStrength for Tillotson) or EOSIsentropic for
 * Tillotson and REOS3 without EOSinitIsentropicLookup; those are
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
    /*
     * Intended reference density for the interface correction: rho0 for Tillotson and ANEOS, 0 for
     * the gases (no natural cutoff). EOSlib never reads it: EOSWoolfsonCoeff cuts at 0.8*rho0 for
     * every material, the gases included (8e-4 g/cc for SCvH and REOS3 with the default rho0 of
     * 1e-3 g/cc). It is kept for callers that read it.
     */
    double rho0forInterfaceCorrection;
    int bEntropyTableInit; // flag to signal if the entropy table is initialized
    int bEntropy; // flag to signal if the eos has entropy available (see EOSSofRhoU)
    double minSoundSpeed; // sound speed at reference values
    /* Unit system: code energy per mass -> erg/g and code density -> g/cc (1 in cgs). */
    double dErgPerGmUnit;
    double dGmPerCcUnit;
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
/*
 * Builds what EOSIsentropic needs beyond EOSinitMaterial: for Tillotson the full isentrope table
 * (46 MB, ~0.6 s per material instance; EOSinitMaterial builds only its cold curve and top
 * isentrope), for REOS3 its isentrope table (reos3InitIsentropicLookup, ~16 MB and ~2 s for
 * REOS3_H_HE: REOS3 has no entropy, the isentropes are integrated from P(rho, u) as for
 * Tillotson). Call it once per material instance that uses EOSIsentropic; a second call does
 * nothing. The ideal gas, ANEOS and SCvH need nothing (bEntropyTableInit is EOS_TRUE from
 * EOSinitMaterial). Sets bEntropyTableInit for Tillotson and REOS3; bEntropy stays EOS_FALSE for
 * both (no entropy values, EOSSofRhoU/EOSSofRhoT give EOS_FAIL_VALUE).
 */
void EOSinitIsentropicLookup(EOSMATERIAL *material, const void * additional_data);
void EOSfinalizeMaterial(EOSMATERIAL *material);
void EOSPrintMat(EOSMATERIAL *material, FILE *fp);

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
 * nLookupFail). If it is above the u column of a tabulated EOS (ANEOS, SCvH),
 * the isentrope starts from the top row, s(rho1, T_max) (counted); the energy
 * above the top row is lost. Where the isentrope leaves a tabulated EOS at
 * rho2, u2 is clamped to the edge it leaves through: the cold curve u(rho2,
 * T_min) below, the top row u(rho2, T_max) above (counted as well);
 * EOS_FAIL_VALUE only for an invalid input or rho1 or rho2 outside a
 * tabulated EOS (no entropy at rho1, no s(T) line at rho2). Tillotson has an analytic P and integrates the ODE directly
 * where (rho1, u1) or rho2 is outside its table (counted as a fallback); it
 * needs EOSinitIsentropicLookup first (fatal without). REOS3 has no entropy:
 * its isentropes are those of du/dlogrho = P(rho, u)/rho, integrated into a
 * table by EOSinitIsentropicLookup (fatal without), with the same clamps as
 * ANEOS and SCvH; rho1 or rho2 outside the lookup limits (the table and its
 * ideal-gas extrapolation) gives EOS_FAIL_VALUE.
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

/*
 * Yield strength model: -2 no strength model for this EOS (Tillotson, EOSYieldStrength is
 * fatal), -1 no yield stress, 0 no limiter, 1 intact, 2 damaged material (see EOSlib.c).
 */
int EOSYieldModel(EOSMATERIAL *material);
/*
 * Yield strength for models 1 and 2 (0 for models -1 and 0). Returns EOS_FAIL_VALUE if the
 * pressure lookup fails (rho outside the table) or T, T_melt is invalid; test with EOSFailed.
 * A state above the u column (overheated, T(rho, u) falls back to T_min) is evaluated on the top
 * row, T = T_max and P = P(rho, T_max), so it has no strength wherever T_melt <= T_max; a state
 * below the cold curve has the strength of the cold curve.
 */
double EOSYieldStrength(EOSMATERIAL *material, double rho, double u);
/*
 * Strength at (rho, T, P) as given. A T from the T_min fallback of a state above the table cannot
 * be told from a cold state and gets the full cold strength: for a (rho, u) state use
 * EOSYieldStrength. A failed P (EOSFailed(P), NaN included) gives EOS_FAIL_VALUE; it is not
 * counted again in nLookupFail.
 */
double EOSYieldStrengthofRhoTP(EOSMATERIAL *material, double rho, double T, double P);

// Inverse functions
double EOSTofRhoU(EOSMATERIAL *material, double rho, double u);
double EOSUofRhoT(EOSMATERIAL *material, double rho, double T);
double EOSRhoofPT(EOSMATERIAL *material, double p, double T);

double EOSRhoofUT(EOSMATERIAL *material, double u, double T);

/*
 * EOS_TRUE if (rho, u) is below the cold curve, EOS_FALSE if not, EOS_FAIL if u_cold(rho) cannot
 * be evaluated (rho outside the table) or u is not finite. Compare with == EOS_TRUE. For the
 * tabulated EOS "below" is the table status (EOSIsInTable == EOS_OUTSIDE_VMIN, the states that
 * fall back to T_min), which agrees with the lookup also where the cold rows of a table are
 * defective (M-ANEOS); for Tillotson it is u < u_cold(rho) - tillUTol(u), the roundoff tolerance of
 * its table status and T(rho, u), and for the ideal gas u < 0. An invalid density (rho <= 0,
 * NaN, inf) gives EOS_FAIL for every material.
 */
int EOSisbelowColdCurve(EOSMATERIAL *material, double rho, double u);

/*
 * EOSIsInTable reports the raw table status of (rho, u) for the material:
 * EOS_SUCCESS, or one of the EOS_OUTSIDE_* codes, or EOS_FAIL (an unknown
 * sublibrary status). It never applies a tolerance policy. Non-finite input
 * gets the code of the side it lies on (rho <= 0, NaN: RHOMIN; rho = +inf:
 * RHOMAX; u = NaN, -inf: VMIN; u = +inf: VMAX), for every library.
 *
 * EOSIsUsable applies the library policy: returns EOS_TRUE if P(rho, u),
 * c(rho, u) and T(rho, u) can be evaluated for this material even if the
 * state lies outside the lookup table (for Tillotson: a finite u > u_max, and
 * rho < rho_min with a finite u >= 0, where the cold curve is zero; for the
 * ideal gas every finite rho > 0 and u >= 0), EOS_FALSE otherwise. Use this in the hydro code; use
 * EOSIsInTable when you need to know *why* a state is out of range or whether
 * isentropic lookups (EOSIsentropic) are available.
 */
int EOSIsInTable(EOSMATERIAL *material, double rho, double u);
int EOSIsUsable(EOSMATERIAL *material, double rho, double u);

// Derivatives
double EOSdPdRho(EOSMATERIAL *material, double rho, double u);
double EOSdPdU(EOSMATERIAL *material, double rho, double u);
/* du/drho|_s = P/rho^2; EOS_FAIL_VALUE where P(rho, u) fails. */
double EOSdUdRho(EOSMATERIAL *material, double rho, double u);

/*
 * dP/drho at constant T and dP/dT at constant rho. Analytic for the ideal gas ((gamma-1) cv T and
 * (gamma-1) rho cv, EOS_FAIL_VALUE for T < 0), from the tables for ANEOS, REOS3 and SCvH; not
 * implemented (fatal) for Tillotson.
 */
double EOSdPdRhoofRhoT(EOSMATERIAL *material, double rho, double T);
double EOSdPdT(EOSMATERIAL *material, double rho, double T);

double EOSGammaofRhoT(EOSMATERIAL *material, double rho, double T);

double EOSUCold(EOSMATERIAL *material, double rho);

/*
 * Boundary condition solver: given (rho1, u1) in material 1, find (rho2, u2)
 * in material 2 with equal pressure and T2 = T1 - Tjump >= 0 (T2 = 0 is the
 * cold curve of Tillotson and the ideal gas; below T_min of a tabulated
 * material 2 the solver fails). Converges to
 * EOS_SOLVEBC_TOL in density; the pressure residual is then at most
 * max(EOS_SOLVEBC_PTOL |P|, EOS_SOLVEBC_SLOPE_FAC (dP2/drho) EOS_SOLVEBC_TOL rho2),
 * i.e. ~ (dlnP/dlnrho) EOS_SOLVEBC_TOL relative for stiff, low-pressure states.
 * Returns EOS_SUCCESS (0), or EOS_FAIL on failure (no bracket, state outside the
 * table, or P2(rho, T2) discontinuous at the root, which happens where u(T)
 * of material 2 is not monotonic, e.g. at a melt transition of a defective
 * table), with a message under EOSLIB_VERBOSE. A discontinuity is detected by
 * the growth of the secant slope of P2(rho) as the bracket closes (see
 * EOSSolveBC in EOSlib.c), not by the pressure gap, which for a continuous but
 * stiff P2 can exceed EOS_SOLVEBC_PTOL |P|. A state (rho1, u1) that is not
 * usable in material 1 (EOSIsUsable: rho1 outside the table, or u1 outside the
 * u(T) column, where T(rho, u) falls back to T_min and P would be that of the
 * fallback) fails at once, without evaluating material 2. The upper bracket
 * starts at rho1 (at rho0 of material 2 if rho1 is below it and cannot be
 * evaluated there), the lower one just below rho0 (and below the upper
 * bracket). Both move by 1.1 and 0.9 per step as before 1.2.0, so wherever that
 * search succeeds the result is the same; where it fails they go on by
 * EOS_SOLVEBC_EXPAND_FAC per step and close in on the edge of the table of
 * material 2 instead of failing at the first trial beyond it (those trials count
 * in its nLookupFail), so a gas envelope far below the density range of a
 * condensed core solves.
 * On failure *prho2
 * and *pu2 are untouched.
 */
int EOSSolveBC(EOSMATERIAL *material1, EOSMATERIAL *material2, double rho1, double u1,
 double *prho2, double *pu2, double Tjump);

/*
 * Woolfson correction f_ij = rho_1(P,T)/rho_2(P,T). Returns 1 for the same
 * material and for P < EOSWoolfsonMinPressure(material1), EOS_FAIL_VALUE if a
 * density inversion fails; callers must check for a negative result
 * (EOSFailed or f < 0). Both materials must use the same unit system.
 * Returns 1 also where rho_1 <= 0.8*rho0 of material 1 or rho_2 <= 0.8*rho0
 * of material 2 (expanded states, ill-conditioned inversion), for every
 * material including the gases; rho0forInterfaceCorrection is not used.
 */
double EOSWoolfsonCoeff(EOSMATERIAL *material1, EOSMATERIAL *material2, double P, double T);
/* WOOLFSON_MIN_PRESSURE_CGS in the code units of the material. */
double EOSWoolfsonMinPressure(EOSMATERIAL *material);
#ifdef __cplusplus
}
#endif
#endif