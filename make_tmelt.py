#!/usr/bin/env python3
"""
Build a T_melt(rho) table for an EOS material from an analytic melting curve
T_melt(P), using the EOSlib python wrapper (eoslib.py). Works for any material
EOSlib can init, not only ANEOS.

  1. evaluate T_melt(P) on a log pressure grid
  2. invert to rho(P, T_melt) with EOSRhoofPT (brentq on P(rho, T) as fallback)
  3. interpolate T_melt onto the table's rho axis (axes.in) or a log rho grid
  4. set T_melt = 0 below cutoff*rho0 and write <out> in one of two binary formats
     (native byte order, as the C libraries read it):

       aneos      int32 nrho, float64 T_melt[nrho]
                  read by ANEOSmaterial as <table>_tmelt. The densities are not in the
                  file: they are the rho axis of the ANEOS table, so --axes (the axes.in
                  of that table) is required and nrho must be its nRho.
       tillotson  int32 n, float64 rho[n] (g/cc), float64 T_melt[n] (K)
                  read by tillotson (>= 4.1.0) as tillotson_<name>.in_tmelt, e.g.
                  tillotson_granite.in_tmelt. The file carries its own density axis
                  (--axes or the log grid), positive and strictly increasing.

     --format auto (default) picks the format from --imat: 1..50 Tillotson,
     51..100 ANEOS (the EOSlib material ranges).

Run from the directory with the EOS table files, as for the C tools.

Note: EOSlib reads the strength files when the material is initialised, and a
strength file with yieldStrengthModel 1 or 2 without its _tmelt file is fatal.
Build the melt curve before the strength file is in place (or in another
directory: a Tillotson material needs no table files).

Examples:
  python make_tmelt.py --curve gonzalez2016 --imat 62 --axes axes.in --out MANEOStable_forsterite.in_tmelt
  python make_tmelt.py --curve senft2009 --imat 6 --rho-min 1e-4 --rho-max 1e3 --nrho 1000 --out tillotson_dunite.in_tmelt
"""
import argparse
import numpy as np
from scipy import optimize
import eoslib

# --- melting curves, P in GPa, T in K
def melting_curve_senft2009(P):
    """Senft & Stewart (2009), https://www.sciencedirect.com/science/article/pii/S0012821X09005172"""
    return np.where(P < 0.7, 2000.0, 2000.0 + (P - 0.7) * 85.7)

def melting_curve_gonzalez2016(P):
    """Gonzalez-Cataldo et al. (2016), https://www.nature.com/articles/srep26537"""
    return 1968.5 + 478.7 * P ** 0.426

def melting_curve_duncan2018(P):
    """Duncan et al. (2018), https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2018GL078182"""
    return np.where(P <= 10.0, -4.877 * P ** 2 + 120.2 * P + 1361.0,
           np.where(P <= 23.0, -1.323 * (P - 10.0) ** 2 + 38.18 * (P - 10.0) + 2075.0,
                    77.75 * (P - 23.0) + 2348.0))

curves = {"senft2009": melting_curve_senft2009,
          "gonzalez2016": melting_curve_gonzalez2016,
          "duncan2018": melting_curve_duncan2018}

# EOSlib material number ranges (EOSlib.h)
MAT_TILLOTSON = (1, 50)
MAT_ANEOS = (51, 100)

# --- arguments
parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("--curve", required=True, choices=sorted(curves), help="melting curve T(P)")
parser.add_argument("--imat", type=int, required=True, help="EOSlib material number")
parser.add_argument("--out", required=True, help="output file")
parser.add_argument("--format", choices=["auto", "aneos", "tillotson"], default="auto",
                    help="file format (default auto: from --imat, 1..50 tillotson, 51..100 aneos)")
parser.add_argument("--axes", help="axes.in of an ANEOS table: use its rho axis (required for --format aneos)")
parser.add_argument("--rho-min", type=float, default=1e-4, help="g/cc, rho grid without --axes")
parser.add_argument("--rho-max", type=float, default=1e3, help="g/cc, rho grid without --axes")
parser.add_argument("--nrho", type=int, default=1000, help="rho grid size without --axes")
parser.add_argument("--p-min", type=float, default=1e-10, help="GPa, sampling range of T(P)")
parser.add_argument("--p-max", type=float, default=1e6, help="GPa, sampling range of T(P)")
parser.add_argument("--np", type=int, default=1000000, help="number of pressure samples")
parser.add_argument("--cutoff", type=float, default=0.9, help="T_melt = 0 below cutoff*rho0 (default 0.9)")
parser.add_argument("--gamma", type=float, help="ideal gas only")
parser.add_argument("--mu", type=float, help="ideal gas only")
parser.add_argument("--save-curve", help="also write the inverted curve as text: T [K] P [GPa] rho [g/cc]")
parser.add_argument("--ascii", help="also write rho [g/cc], T_melt [K] of the final table as text")
args = parser.parse_args()

# --- output format (checked before the material is initialised)
fmt = args.format
if fmt == "auto":
    if MAT_TILLOTSON[0] <= args.imat <= MAT_TILLOTSON[1]:
        fmt = "tillotson"
    elif MAT_ANEOS[0] <= args.imat <= MAT_ANEOS[1]:
        fmt = "aneos"
    else:
        parser.error("--format auto: material %d is neither Tillotson (%d..%d) nor ANEOS (%d..%d); "
                     "give --format" % ((args.imat,) + MAT_TILLOTSON + MAT_ANEOS))
if fmt == "aneos" and not args.axes:
    # The ANEOS file has no density axis: T_melt[i] belongs to rhoAxis[i] of the table. A log grid
    # is rejected at load (nRho mismatch) or, if nrho happens to match, read on the wrong densities.
    parser.error("--format aneos needs --axes (the axes.in of the ANEOS table)")
print("format:", fmt)

GPA = 1e10   # dyn/cm^2, EOSlib is used in cgs here

# --- material (cgs units)
mat = eoslib.Material(args.imat, 0.0, 0.0, gamma=args.gamma, mu=args.mu)
rho0 = mat.rho0
print("material:", mat, " rho0 = %g g/cc" % rho0)

# --- T(P) and rho(P, T)
P = np.geomspace(args.p_min, args.p_max, args.np)
T = curves[args.curve](P)

rho = mat.RhoofPT(P * GPA, T)
bad = ~(np.isfinite(rho) & (rho > 0.0))
if bad.any():
    # fall back to a root finder on P(rho, T), as the old calcPressureRhoT version did
    print("RhoofPT failed for %d points, using brentq on P(rho, T)" % bad.sum())
    for i in np.flatnonzero(bad):
        f = lambda r: P[i] * GPA - mat.PofRhoT(r, T[i])
        try:
            rho[i] = optimize.root_scalar(f, method="brentq", bracket=(1e-25, 1e3)).root
        except ValueError:
            rho[i] = np.nan
ok = np.isfinite(rho) & (rho > 0.0)
if not ok.all():
    print("dropping %d of %d points" % ((~ok).sum(), ok.size))
P, T, rho = P[ok], T[ok], rho[ok]
if rho.size < 2:
    raise SystemExit("error: fewer than two points of the melting curve could be inverted")
print("inverted curve: rho = %g .. %g g/cc, T = %g .. %g K" % (rho.min(), rho.max(), T.min(), T.max()))

if args.save_curve:
    header = "{:>14s}{:>15s}{:>15s}".format("T [K]", "P [GPa]", "rho [g/cm^3]")
    np.savetxt(args.save_curve, np.column_stack([T, P, rho]), fmt="%15.7E", header=header, comments="#")
    print("wrote", args.save_curve)

# --- T_melt on the density axis
if args.axes:
    axes = np.loadtxt(args.axes)
    nrho = int(axes[0])
    rho_axis = axes[2:nrho + 2]
else:
    nrho = args.nrho
    rho_axis = np.geomspace(args.rho_min, args.rho_max, nrho)

order = np.argsort(rho)
T_melt = np.interp(rho_axis, rho[order], T[order], left=0.0, right=T[order][-1])
T_melt[rho_axis < args.cutoff * rho0] = 0.0

# --- checks of the readers (ANEOSReadMeltCurve, tillReadMeltCurve), so that a bad file fails here and not at load
rho_axis = np.asarray(rho_axis, dtype=np.float64)
T_melt = np.asarray(T_melt, dtype=np.float64)
if rho_axis.size != nrho:
    raise SystemExit("error: the rho axis has %d values, expected nrho = %d" % (rho_axis.size, nrho))
if not (np.all(np.isfinite(T_melt)) and np.all(T_melt >= 0.0)):
    raise SystemExit("error: T_melt is negative or not finite")
if fmt == "tillotson":
    if nrho < 2:
        raise SystemExit("error: the Tillotson melt curve needs at least two nodes (nrho = %d)" % nrho)
    if not (np.all(np.isfinite(rho_axis)) and np.all(rho_axis > 0.0) and np.all(np.diff(rho_axis) > 0.0)):
        raise SystemExit("error: the rho axis must be positive, finite and strictly increasing")

# --- output
with open(args.out, "wb") as f:
    f.write(np.int32(nrho).tobytes())
    if fmt == "tillotson":
        f.write(rho_axis.tobytes())      # g/cc
    f.write(T_melt.tobytes())            # K
print("wrote %s (%s format, nrho = %d)" % (args.out, fmt, nrho))

if args.ascii:
    np.savetxt(args.ascii, np.column_stack([rho_axis, T_melt]), header="rho [g/cc]  T_melt [K]")
    print("wrote", args.ascii)