#!/usr/bin/env python3
"""
Build a T_melt(rho) table for an EOS material from an analytic melting curve
T_melt(P), using the EOSlib python wrapper (eoslib.py). Works for any material
EOSlib can init, not only ANEOS.

  1. evaluate T_melt(P) on a log pressure grid
  2. invert to rho(P, T_melt) with EOSRhoofPT (brentq on P(rho, T) as fallback)
  3. interpolate T_melt onto the table's rho axis (axes.in) or a log rho grid
  4. set T_melt = 0 below cutoff*rho0 and write <out>: int32 nrho, float64 T_melt[nrho]

Run from the directory with the EOS table files, as for the C tools.

Examples:
  python make_tmelt.py --curve gonzalez2016 --imat 52 --axes axes.in --out MANEOStable_dunite.in_tmelt
  python make_tmelt.py --curve senft2009 --imat 2 --rho-min 1e-4 --rho-max 1e3 --nrho 1000 --out till_dunite_tmelt
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

# --- arguments
parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("--curve", required=True, choices=sorted(curves), help="melting curve T(P)")
parser.add_argument("--imat", type=int, required=True, help="EOSlib material number")
parser.add_argument("--out", required=True, help="output file (int32 nrho + float64 T_melt)")
parser.add_argument("--axes", help="axes.in of an ANEOS table: use its rho axis")
parser.add_argument("--rho-min", type=float, default=1e-4, help="g/cc, rho grid without --axes")
parser.add_argument("--rho-max", type=float, default=1e3, help="g/cc, rho grid without --axes")
parser.add_argument("--nrho", type=int, default=1000, help="rho grid size without --axes")
parser.add_argument("--p-min", type=float, default=1e-10, help="GPa, sampling range of T(P)")
parser.add_argument("--p-max", type=float, default=1e4, help="GPa, sampling range of T(P)")
parser.add_argument("--np", type=int, default=10000, help="number of pressure samples")
parser.add_argument("--cutoff", type=float, default=0.9, help="T_melt = 0 below cutoff*rho0 (default 0.9)")
parser.add_argument("--gamma", type=float, help="ideal gas only")
parser.add_argument("--mu", type=float, help="ideal gas only")
parser.add_argument("--save-curve", help="also write the inverted curve as text: T [K] P [GPa] rho [g/cc]")
parser.add_argument("--ascii", help="also write rho [g/cc], T_melt [K] of the final table as text")
args = parser.parse_args()

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

# --- output
with open(args.out, "wb") as f:
    f.write(np.int32(nrho).tobytes())
    f.write(T_melt.astype(np.float64).tobytes())
print("wrote %s (nrho = %d)" % (args.out, nrho))

if args.ascii:
    np.savetxt(args.ascii, np.column_stack([rho_axis, T_melt]), header="rho [g/cc]  T_melt [K]")
    print("wrote", args.ascii)