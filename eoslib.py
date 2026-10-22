"""
Python wrapper around EOSlib (ctypes).

Build the shared library once with

    make libEOSlib

in the EOSlib directory; this module loads libEOSlib.so from the directory
it lives in (override with the EOSLIB_PATH environment variable). Run Python
from the directory that contains the EOS table files, exactly as for the
compiled tools.

The material tables are loaded once in EOSinitMaterial and stay in memory for
the lifetime of the Material object.

Example:

    import eoslib
    iron = eoslib.Material(52, dKpcUnit, dMsolUnit)
    P, c, T = iron.PCTofRhoU(rho, u)
    P = iron.PofRhoU(rho_array, u_array)      # NumPy arrays work too

Error behaviour: the C library terminates the process (assert/exit) on
unknown materials or unimplemented functions. Those cases cannot be turned
into Python exceptions; avoid calling, e.g., GammaofRhoT on Tillotson.
Sentinel return values (-1e50 for entropy, -1.0 for the Woolfson coefficient)
are mapped to NaN.
"""
import ctypes as ct
import os
import sys

import numpy as np

# --- library handle -----------------------------------------------------------

_path = os.environ.get("EOSLIB_PATH",
                       os.path.join(os.path.dirname(os.path.abspath(__file__)), "libEOSlib.so"))
_lib = ct.CDLL(_path)

# --- constants from EOSlib.h --------------------------------------------------

EOSIDEALGAS = 0
EOSTILLOTSON = 1
EOSANEOS = 2
EOSREOS3 = 3
EOSSCVHEOS = 4
MAT_TYPE_NAMES = {0: "ideal gas", 1: "Tillotson", 2: "ANEOS", 3: "REOS3", 4: "SCvH"}

MAT_IDEALGAS = 0
MAT_TILLOTSON_MIN, MAT_TILLOTSON_MAX = 1, 50
MAT_ANEOS_MIN, MAT_ANEOS_MAX = 51, 100
MAT_REOS3_MIN, MAT_REOS3_MAX = 101, 109
MAT_SCVHEOS_MIN, MAT_SCVHEOS_MAX = 110, 130

EOS_TRUE, EOS_FALSE = 1, 0
EOS_SUCCESS, EOS_FAIL = 1, 0
EOS_OUTSIDE_RHOMIN = -1
EOS_OUTSIDE_RHOMAX = -2
EOS_OUTSIDE_VMIN = -3
EOS_OUTSIDE_VMAX = -4
STATUS_NAMES = {1: "in table", 0: "fail", -1: "rho < rho_min", -2: "rho > rho_max",
                -3: "u < u_min", -4: "u > u_max"}

_S_FAIL = -1e50  # entropy sentinel

# --- C structures -------------------------------------------------------------


class IgeosParam(ct.Structure):
    _fields_ = [("dConstGamma", ct.c_double), ("dMeanMolMass", ct.c_double)]


class _EOSMaterialHead(ct.Structure):
    """Leading fields of EOSMATERIAL. Everything after MatString is #ifdef-
    dependent pointers, so only this prefix is declared (read-only use)."""
    _fields_ = [("iMat", ct.c_int),
                ("matType", ct.c_int),
                ("rho0", ct.c_double),
                ("rho0forInterfaceCorrection", ct.c_double),
                ("bEntropyTableInit", ct.c_int),
                ("bEntropy", ct.c_int),
                ("minSoundSpeed", ct.c_double),
                ("MatString", ct.c_char * 256)]


# --- function signatures ------------------------------------------------------

_d, _i, _p, _dp = ct.c_double, ct.c_int, ct.c_void_p, ct.POINTER(ct.c_double)

_sigs = {
    "EOSinitMaterial": (_p, [_i, _d, _d, _p]),
    "EOSinitIsentropicLookup": (None, [_p, _p]),
    "EOSfinalizeMaterial": (None, [_p]),
    "EOSPrintMat": (None, [_p, _p]),
    "EOSCheckVersions": (_i, []),
    "EOSPofRhoU": (_d, [_p, _d, _d]),
    "EOSPofRhoT": (_d, [_p, _d, _d]),
    "EOSCofRhoU": (_d, [_p, _d, _d]),
    "EOSPCofRhoU": (_d, [_p, _d, _d, _dp]),
    "EOSPCTofRhoU": (_d, [_p, _d, _d, _dp, _dp]),
    "EOSIsentropic": (_d, [_p, _d, _d, _d]),
    "EOSSofRhoU": (_d, [_p, _d, _d]),
    "EOSSofRhoT": (_d, [_p, _d, _d]),
    "EOSYieldModel": (_i, [_p]),
    "EOSYieldStrength": (_d, [_p, _d, _d]),
    "EOSYieldStrengthofRhoTP": (_d, [_p, _d, _d, _d]),
    "EOSTofRhoU": (_d, [_p, _d, _d]),
    "EOSUofRhoT": (_d, [_p, _d, _d]),
    "EOSRhoofPT": (_d, [_p, _d, _d]),
    "EOSRhoofUT": (_d, [_p, _d, _d]),
    "EOSisbelowColdCurve": (_i, [_p, _d, _d]),
    "EOSIsInTable": (_i, [_p, _d, _d]),
    "EOSIsUsable": (_i, [_p, _d, _d]),
    "EOSdPdRho": (_d, [_p, _d, _d]),
    "EOSdPdU": (_d, [_p, _d, _d]),
    "EOSdUdRho": (_d, [_p, _d, _d]),
    "EOSdPdRhoofRhoT": (_d, [_p, _d, _d]),
    "EOSdPdT": (_d, [_p, _d, _d]),
    "EOSGammaofRhoT": (_d, [_p, _d, _d]),
    "EOSUCold": (_d, [_p, _d]),
    "EOSSolveBC": (_i, [_p, _p, _d, _d, _dp, _dp, _d]),
    "EOSWoolfsonCoeff": (_d, [_p, _p, _d, _d]),
}
for _name, (_res, _args) in _sigs.items():
    _f = getattr(_lib, _name)
    _f.restype = _res
    _f.argtypes = _args


def check_versions():
    """True if all sublibraries have the required versions."""
    return bool(_lib.EOSCheckVersions())


def _c_stdout():
    """FILE* of the C stdout (for EOSPrintMat)."""
    libc = ct.CDLL(None)
    for sym in ("stdout", "__stdoutp"):  # glibc, macOS
        try:
            return ct.c_void_p.in_dll(libc, sym)
        except ValueError:
            pass
    raise RuntimeError("cannot locate C stdout")


# --- material -----------------------------------------------------------------


class Material:
    """One EOS material. Scalar arguments return floats; NumPy arrays (or
    anything broadcastable) are evaluated element-wise."""

    def __init__(self, iMat, dKpcUnit, dMsolUnit, gamma=None, mu=None):
        """iMat: material number. For the ideal gas (iMat=0) gamma and mu
        (mean molecular mass) are required."""
        self._mat = None
        param = None
        if iMat == MAT_IDEALGAS:
            if gamma is None or mu is None:
                raise ValueError("ideal gas needs gamma and mu")
            self._param = IgeosParam(gamma, mu)  # must outlive the call
            param = ct.cast(ct.pointer(self._param), ct.c_void_p)
        self._mat = _lib.EOSinitMaterial(iMat, dKpcUnit, dMsolUnit, param)
        if not self._mat:
            raise RuntimeError("EOSinitMaterial returned NULL")
        self._head = ct.cast(self._mat, ct.POINTER(_EOSMaterialHead)).contents
        self.dKpcUnit = dKpcUnit
        self.dMsolUnit = dMsolUnit

    def __del__(self):
        if getattr(self, "_mat", None):
            _lib.EOSfinalizeMaterial(self._mat)
            self._mat = None

    def __repr__(self):
        return "Material(iMat=%d, %s, '%s')" % (self.iMat, MAT_TYPE_NAMES.get(self.matType, "?"),
                                                 self.MatString)

    # -- read-only struct fields
    @property
    def iMat(self):
        return self._head.iMat

    @property
    def matType(self):
        return self._head.matType

    @property
    def rho0(self):
        return self._head.rho0

    @property
    def rho0forInterfaceCorrection(self):
        return self._head.rho0forInterfaceCorrection

    @property
    def bEntropy(self):
        return bool(self._head.bEntropy)

    @property
    def bEntropyTableInit(self):
        return bool(self._head.bEntropyTableInit)

    @property
    def minSoundSpeed(self):
        return self._head.minSoundSpeed

    @property
    def MatString(self):
        return self._head.MatString.decode(errors="replace")

    # -- setup / info
    def initIsentropicLookup(self):
        _lib.EOSinitIsentropicLookup(self._mat, None)

    def printMat(self):
        sys.stdout.flush()
        _lib.EOSPrintMat(self._mat, _c_stdout())

    def yieldModel(self):
        return _lib.EOSYieldModel(self._mat)

    # -- element-wise evaluation helper
    def _ev(self, cfunc, *args, nout=1, otypes=None):
        if all(np.ndim(a) == 0 for a in args):
            return cfunc(self._mat, *[float(a) for a in args])
        f = np.vectorize(lambda *a: cfunc(self._mat, *a),
                         otypes=otypes or [float] * nout)
        return f(*args)

    # -- main functions
    def PofRhoU(self, rho, u):
        return self._ev(_lib.EOSPofRhoU, rho, u)

    def PofRhoT(self, rho, T):
        return self._ev(_lib.EOSPofRhoT, rho, T)

    def CofRhoU(self, rho, u):
        return self._ev(_lib.EOSCofRhoU, rho, u)

    def PCofRhoU(self, rho, u):
        """Returns (P, c)."""
        def f(mat, rho, u):
            c = ct.c_double()
            P = _lib.EOSPCofRhoU(mat, rho, u, ct.byref(c))
            return P, c.value
        return self._ev(f, rho, u, nout=2)

    def PCTofRhoU(self, rho, u):
        """Returns (P, c, T)."""
        def f(mat, rho, u):
            c, T = ct.c_double(), ct.c_double()
            P = _lib.EOSPCTofRhoU(mat, rho, u, ct.byref(c), ct.byref(T))
            return P, c.value, T.value
        return self._ev(f, rho, u, nout=3)

    def Isentropic(self, rho1, u1, rho2):
        """u2 after isentropic evolution from (rho1, u1) to rho2."""
        return self._ev(_lib.EOSIsentropic, rho1, u1, rho2)

    def SofRhoU(self, rho, u):
        """Specific entropy in cgs (erg/g/K); NaN if unavailable."""
        S = self._ev(_lib.EOSSofRhoU, rho, u)
        return np.where(np.asarray(S) <= _S_FAIL, np.nan, S)[()]

    def SofRhoT(self, rho, T):
        """Specific entropy in cgs (erg/g/K); NaN if unavailable."""
        S = self._ev(_lib.EOSSofRhoT, rho, T)
        return np.where(np.asarray(S) <= _S_FAIL, np.nan, S)[()]

    def YieldStrength(self, rho, u):
        return self._ev(_lib.EOSYieldStrength, rho, u)

    def YieldStrengthofRhoTP(self, rho, T, P):
        return self._ev(_lib.EOSYieldStrengthofRhoTP, rho, T, P)

    # -- inverse functions
    def TofRhoU(self, rho, u):
        return self._ev(_lib.EOSTofRhoU, rho, u)

    def UofRhoT(self, rho, T):
        return self._ev(_lib.EOSUofRhoT, rho, T)

    def RhoofPT(self, P, T):
        return self._ev(_lib.EOSRhoofPT, P, T)

    def RhoofUT(self, u, T):
        """ANEOS only (other materials abort the process)."""
        return self._ev(_lib.EOSRhoofUT, u, T)

    # -- table checks
    def isbelowColdCurve(self, rho, u):
        return self._ev(_lib.EOSisbelowColdCurve, rho, u, otypes=[bool])

    def isInTable(self, rho, u):
        """Raw status: EOS_SUCCESS or one of the EOS_OUTSIDE_* codes (see STATUS_NAMES)."""
        return self._ev(_lib.EOSIsInTable, rho, u, otypes=[int])

    def isUsable(self, rho, u):
        return self._ev(_lib.EOSIsUsable, rho, u, otypes=[bool])

    # -- derivatives
    def dPdRho(self, rho, u):
        return self._ev(_lib.EOSdPdRho, rho, u)

    def dPdU(self, rho, u):
        return self._ev(_lib.EOSdPdU, rho, u)

    def dUdRho(self, rho, u):
        return self._ev(_lib.EOSdUdRho, rho, u)

    def dPdRhoofRhoT(self, rho, T):
        """Tabulated EOS only."""
        return self._ev(_lib.EOSdPdRhoofRhoT, rho, T)

    def dPdT(self, rho, T):
        """Tabulated EOS only."""
        return self._ev(_lib.EOSdPdT, rho, T)

    def GammaofRhoT(self, rho, T):
        """Shear modulus; ANEOS only (Tillotson aborts the process)."""
        return self._ev(_lib.EOSGammaofRhoT, rho, T)

    def UCold(self, rho):
        return self._ev(_lib.EOSUCold, rho)


# --- two-material functions ---------------------------------------------------


def solveBC(material1, material2, rho1, u1, Tjump=0.0):
    """Boundary condition at a material interface: given (rho1, u1) in
    material1 return (rho2, u2) in material2 with equal pressure and
    T2 = T1 - Tjump. Raises RuntimeError on failure (the C library prints the
    reason to stderr)."""
    rho2, u2 = ct.c_double(), ct.c_double()
    ret = _lib.EOSSolveBC(material1._mat, material2._mat, float(rho1), float(u1),
                          ct.byref(rho2), ct.byref(u2), float(Tjump))
    if ret != 0:
        raise RuntimeError("EOSSolveBC failed for rho1=%g u1=%g (iMat %d -> %d)"
                           % (rho1, u1, material1.iMat, material2.iMat))
    return rho2.value, u2.value


def woolfsonCoeff(material1, material2, P, T):
    """f_ij = rho_1(P,T)/rho_2(P,T); NaN where the density inversion fails."""
    def f(P, T):
        return _lib.EOSWoolfsonCoeff(material1._mat, material2._mat, P, T)
    if np.ndim(P) == 0 and np.ndim(T) == 0:
        r = f(float(P), float(T))
        return np.nan if r < 0.0 else r
    r = np.vectorize(f, otypes=[float])(P, T)
    return np.where(r < 0.0, np.nan, r)