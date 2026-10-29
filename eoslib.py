"""
Python wrapper around EOSlib (ctypes), matching EOSlib.h 1.4.0.

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

Units: both dKpcUnit and dMsolUnit > 0 (and finite) are code units, both <= 0
is cgs; anything else raises ValueError (the C library would abort).

Failure contract (EOSlib.h): a failed lookup returns EOS_FAIL_VALUE (-1e50)
from every double-valued function. With fail_to_nan = True (the default) the
wrapper maps it (and NaN) to NaN, for scalars and arrays alike; set
eoslib.fail_to_nan = False to get the raw sentinel (test it with failed()).
T(rho, u) never fails (it falls back to T_min).

Fatal errors: the C library terminates the process (assert/exit) on unknown
materials, missing tables and functions that are not implemented for a
material. The wrapper checks the cases it can decide in Python (material id
range, unit system, ideal gas parameters, not-implemented functions such as
GammaofRhoT on Tillotson) and raises ValueError/NotImplementedError instead.
A missing table file or a sublibrary that was not compiled in still aborts.

Isentropes: Tillotson and REOS3 need EOSinitIsentropicLookup before
EOSIsentropic (fatal without). The wrapper builds the table on the first call
of Material.Isentropic, or at construction with isentropic=True.

Diagnostic messages of the C library (lookup failures, EOSSolveBC) are only
printed if it was compiled with -DEOSLIB_VERBOSE.
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

# Version of EOSlib.h this wrapper was written for (the C library exports no
# version symbol; the sublibrary versions are checked at compile time).
EOS_VERSION_TEXT = "1.4.0"
EOS_VERSION_MAJOR, EOS_VERSION_MINOR, EOS_VERSION_PATCH = 1, 4, 0

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
EOS_N_MATERIAL_MAX = MAT_SCVHEOS_MAX + 1   # size of an array indexed by iMat

EOS_TRUE, EOS_FALSE = 1, 0
EOS_SUCCESS, EOS_FAIL = 0, -99
EOS_OUTSIDE_RHOMIN = -1
EOS_OUTSIDE_RHOMAX = -2
EOS_OUTSIDE_VMIN = -3
EOS_OUTSIDE_VMAX = -4
STATUS_NAMES = {EOS_SUCCESS: "in table", EOS_FAIL: "fail",
                EOS_OUTSIDE_RHOMIN: "rho < rho_min", EOS_OUTSIDE_RHOMAX: "rho > rho_max",
                EOS_OUTSIDE_VMIN: "u < u_min", EOS_OUTSIDE_VMAX: "u > u_max"}

EOS_FAIL_VALUE = -1e50

# Yield strength models (EOSYieldModel)
YIELD_MODEL_NAMES = {-2: "no strength model for this EOS", -1: "no yield stress",
                     0: "no limiter", 1: "intact (Collins et al. 2004)", 2: "damaged"}

# EOSWoolfsonCoeff does no correction below this pressure (dyn/cm^2)
WOOLFSON_MIN_PRESSURE_CGS = 3.684121632893488e6

# EOSSolveBC tolerances
EOS_SOLVEBC_MAX_EXPAND = 100
EOS_SOLVEBC_EXPAND_FAC = 2.0
EOS_SOLVEBC_MAX_COARSE = 100
EOS_SOLVEBC_MAX_BISECT = 200
EOS_SOLVEBC_TOL = 1e-10
EOS_SOLVEBC_PTOL = 1e-6
EOS_SOLVEBC_SLOPE_WIDTH = 1e-6
EOS_SOLVEBC_SLOPE_FAC = 10.0
EOS_SOLVEBC_PROBE = 4.0

# Map EOS_FAIL_VALUE (and NaN) to NaN in returned values.
fail_to_nan = True


def failed(x):
    """EOSFailed: True for the sentinel EOS_FAIL_VALUE and for NaN (element-wise)."""
    r = ~(np.asarray(x, dtype=float) > -1e40)
    return bool(r) if r.ndim == 0 else r


def _clean(x):
    if np.ndim(x) == 0:
        x = float(x)
        return x if x > -1e40 else np.nan
    a = np.asarray(x, dtype=float)
    return np.where(a > -1e40, a, np.nan)


# --- C structures -------------------------------------------------------------


class IgeosParam(ct.Structure):
    _fields_ = [("dConstGamma", ct.c_double), ("dMeanMolMass", ct.c_double)]


class _EOSMaterialHead(ct.Structure):
    """Leading fields of EOSMATERIAL (EOSlib.h 1.4.0). Everything after
    MatString is #ifdef-dependent pointers, so only this prefix is declared
    (read-only use)."""
    _fields_ = [("iMat", ct.c_int),
                ("matType", ct.c_int),
                ("rho0", ct.c_double),
                ("rho0forInterfaceCorrection", ct.c_double),
                ("bEntropyTableInit", ct.c_int),
                ("bEntropy", ct.c_int),
                ("minSoundSpeed", ct.c_double),
                ("dErgPerGmUnit", ct.c_double),
                ("dGmPerCcUnit", ct.c_double),
                ("MatString", ct.c_char * 256)]


# --- function signatures ------------------------------------------------------

_d, _i, _p, _dp = ct.c_double, ct.c_int, ct.c_void_p, ct.POINTER(ct.c_double)

_sigs = {
    "EOSinitMaterial": (_p, [_i, _d, _d, _p]),
    "EOSinitIsentropicLookup": (None, [_p, _p]),
    "EOSfinalizeMaterial": (None, [_p]),
    "EOSPrintMat": (None, [_p, _p]),
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
    "EOSWoolfsonMinPressure": (_d, [_p]),
}
for _name, (_res, _args) in _sigs.items():
    _f = getattr(_lib, _name)
    _f.restype = _res
    _f.argtypes = _args


def _c_stdout():
    """FILE* of the C stdout (for EOSPrintMat)."""
    libc = ct.CDLL(None)
    for sym in ("stdout", "__stdoutp"):  # glibc, macOS
        try:
            return ct.c_void_p.in_dll(libc, sym)
        except ValueError:
            pass
    raise RuntimeError("cannot locate C stdout")


def _valid_iMat(iMat):
    return (iMat == MAT_IDEALGAS
            or MAT_TILLOTSON_MIN <= iMat <= MAT_TILLOTSON_MAX
            or MAT_ANEOS_MIN <= iMat <= MAT_ANEOS_MAX
            or MAT_REOS3_MIN <= iMat <= MAT_REOS3_MAX
            or MAT_SCVHEOS_MIN <= iMat <= MAT_SCVHEOS_MAX)


def _valid_units(dKpcUnit, dMsolUnit):
    code = dKpcUnit > 0.0 and dMsolUnit > 0.0 and np.isfinite(dKpcUnit) and np.isfinite(dMsolUnit)
    cgs = dKpcUnit <= 0.0 and dMsolUnit <= 0.0
    return code or cgs


def _vec(func, args, nout=1, otypes=None, clean=True):
    """Evaluate func(*args) element-wise. Scalars give Python scalars."""
    if all(np.ndim(a) == 0 for a in args):
        r = func(*[float(a) for a in args])
        if otypes is not None:
            r = otypes[0](r)
    else:
        arrs = [np.asarray(a, dtype=float) for a in args]
        r = np.vectorize(func, otypes=otypes or [float] * nout)(*arrs)
    if not clean or not fail_to_nan:
        return r
    if nout == 1:
        return _clean(r)
    return tuple(_clean(x) for x in r)


# --- material -----------------------------------------------------------------


class Material:
    """One EOS material. Scalar arguments return floats; NumPy arrays (or
    anything broadcastable) are evaluated element-wise."""

    def __init__(self, iMat, dKpcUnit, dMsolUnit, gamma=None, mu=None, isentropic=False):
        """iMat: material number. For the ideal gas (iMat=0) gamma (> 1) and
        mu (mean molecular mass, > 0) are required. isentropic=True builds the
        isentrope table at once (Tillotson, REOS3); otherwise it is built on
        the first call of Isentropic."""
        self._mat = None
        iMat = int(iMat)
        dKpcUnit = float(dKpcUnit)
        dMsolUnit = float(dMsolUnit)
        if not _valid_iMat(iMat):
            raise ValueError("iMat %d does not exist" % iMat)
        if not _valid_units(dKpcUnit, dMsolUnit):
            raise ValueError("invalid unit system dKpcUnit=%g dMsolUnit=%g (give both > 0 for "
                             "code units or both <= 0 for cgs)" % (dKpcUnit, dMsolUnit))
        param = None
        if iMat == MAT_IDEALGAS:
            if gamma is None or mu is None:
                raise ValueError("ideal gas needs gamma and mu")
            if not (gamma > 1.0 and np.isfinite(gamma)):
                raise ValueError("ideal gas needs gamma > 1 (got %g)" % gamma)
            if not (mu > 0.0 and np.isfinite(mu)):
                raise ValueError("ideal gas needs mu > 0 (got %g)" % mu)
            self._param = IgeosParam(gamma, mu)  # must outlive the call
            param = ct.cast(ct.pointer(self._param), ct.c_void_p)
        self._mat = _lib.EOSinitMaterial(iMat, dKpcUnit, dMsolUnit, param)
        if not self._mat:
            raise RuntimeError("EOSinitMaterial returned NULL")
        self._head = ct.cast(self._mat, ct.POINTER(_EOSMaterialHead)).contents
        self.dKpcUnit = dKpcUnit
        self.dMsolUnit = dMsolUnit
        self._yieldModel = _lib.EOSYieldModel(self._mat)
        if isentropic:
            self.initIsentropicLookup()

    def __del__(self):
        if getattr(self, "_mat", None) and _lib is not None:
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
        """Not read by EOSlib (EOSWoolfsonCoeff cuts at 0.8*rho0 for every material)."""
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
    def dErgPerGmUnit(self):
        """Code energy per mass -> erg/g (1 in cgs)."""
        return self._head.dErgPerGmUnit

    @property
    def dGmPerCcUnit(self):
        """Code density -> g/cc (1 in cgs)."""
        return self._head.dGmPerCcUnit

    @property
    def MatString(self):
        return self._head.MatString.decode(errors="replace")

    # -- internal helpers
    def _ev(self, cfunc, *args, nout=1, otypes=None, clean=True):
        mat = self._mat
        return _vec(lambda *a: cfunc(mat, *a), args, nout=nout, otypes=otypes, clean=clean)

    def _require(self, allowed, name):
        if self.matType not in allowed:
            raise NotImplementedError("%s is not implemented for material %d (%s)"
                                      % (name, self.iMat, MAT_TYPE_NAMES.get(self.matType, "?")))

    # -- setup / info
    def initIsentropicLookup(self):
        """Build the isentrope table (Tillotson ~46 MB/0.6 s, REOS3 ~16 MB/2 s);
        a no-op for the other EOS and on a second call."""
        _lib.EOSinitIsentropicLookup(self._mat, None)

    def printMat(self):
        """Print the material data (all dimensional values in cgs) and nLookupFail."""
        sys.stdout.flush()
        _lib.EOSPrintMat(self._mat, _c_stdout())

    def yieldModel(self):
        """-2 no strength model (Tillotson), -1 no yield stress, 0 no limiter,
        1 intact, 2 damaged (see YIELD_MODEL_NAMES)."""
        return self._yieldModel

    def woolfsonMinPressure(self):
        """WOOLFSON_MIN_PRESSURE_CGS in the code units of the material."""
        return _lib.EOSWoolfsonMinPressure(self._mat)

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
        """Returns (P, c, T). T never fails (falls back to T_min)."""
        def f(mat, rho, u):
            c, T = ct.c_double(), ct.c_double()
            P = _lib.EOSPCTofRhoU(mat, rho, u, ct.byref(c), ct.byref(T))
            return P, c.value, T.value
        return self._ev(f, rho, u, nout=3)

    def Isentropic(self, rho1, u1, rho2):
        """u2 after isentropic evolution from (rho1, u1) to rho2. Builds the
        isentrope table first if needed (Tillotson, REOS3)."""
        if not self._head.bEntropyTableInit:
            self.initIsentropicLookup()
        return self._ev(_lib.EOSIsentropic, rho1, u1, rho2)

    def SofRhoU(self, rho, u):
        """Specific entropy in cgs (erg/g/K), additive constant library-dependent;
        NaN if unavailable (Tillotson, REOS3) or the lookup fails."""
        S = self._ev(_lib.EOSSofRhoU, rho, u, clean=False)
        return _clean(S)

    def SofRhoT(self, rho, T):
        """Specific entropy in cgs (erg/g/K); NaN if unavailable."""
        S = self._ev(_lib.EOSSofRhoT, rho, T, clean=False)
        return _clean(S)

    def YieldStrength(self, rho, u):
        """Strength limiter from (rho, u); a state above the table is evaluated
        on the top row T_max."""
        if self._yieldModel == -2:
            raise NotImplementedError("no strength model for material %d (%s)"
                                      % (self.iMat, self.MatString))
        return self._ev(_lib.EOSYieldStrength, rho, u)

    def YieldStrengthofRhoTP(self, rho, T, P):
        """Strength limiter at (rho, T, P) as given."""
        if self._yieldModel == -2:
            raise NotImplementedError("no strength model for material %d (%s)"
                                      % (self.iMat, self.MatString))
        return self._ev(_lib.EOSYieldStrengthofRhoTP, rho, T, P)

    # -- inverse functions
    def TofRhoU(self, rho, u):
        """Never fails: T_min of the table (0 for Tillotson and the ideal gas)
        where (rho, u) cannot be inverted."""
        return self._ev(_lib.EOSTofRhoU, rho, u)

    def UofRhoT(self, rho, T):
        return self._ev(_lib.EOSUofRhoT, rho, T)

    def RhoofPT(self, P, T):
        return self._ev(_lib.EOSRhoofPT, P, T)

    def RhoofUT(self, u, T):
        """ANEOS only."""
        self._require((EOSANEOS,), "RhoofUT")
        return self._ev(_lib.EOSRhoofUT, u, T)

    # -- table checks
    def isbelowColdCurve(self, rho, u):
        """Raw status: EOS_TRUE, EOS_FALSE, or EOS_FAIL where the cold curve
        cannot be evaluated (rho outside the table) or u is not finite.
        Compare with == EOS_TRUE."""
        return self._ev(_lib.EOSisbelowColdCurve, rho, u, otypes=[int], clean=False)

    def isInTable(self, rho, u):
        """Raw status: EOS_SUCCESS (0) or one of the EOS_OUTSIDE_* codes, or
        EOS_FAIL (see STATUS_NAMES)."""
        return self._ev(_lib.EOSIsInTable, rho, u, otypes=[int], clean=False)

    def isUsable(self, rho, u):
        return self._ev(_lib.EOSIsUsable, rho, u, otypes=[bool], clean=False)

    # -- derivatives
    def dPdRho(self, rho, u):
        return self._ev(_lib.EOSdPdRho, rho, u)

    def dPdU(self, rho, u):
        return self._ev(_lib.EOSdPdU, rho, u)

    def dUdRho(self, rho, u):
        """P/rho^2; NaN where P(rho, u) fails."""
        return self._ev(_lib.EOSdUdRho, rho, u)

    def dPdRhoofRhoT(self, rho, T):
        """Ideal gas and tabulated EOS (not Tillotson)."""
        self._require((EOSIDEALGAS, EOSANEOS, EOSREOS3, EOSSCVHEOS), "dPdRhoofRhoT")
        return self._ev(_lib.EOSdPdRhoofRhoT, rho, T)

    def dPdT(self, rho, T):
        """Ideal gas and tabulated EOS (not Tillotson)."""
        self._require((EOSIDEALGAS, EOSANEOS, EOSREOS3, EOSSCVHEOS), "dPdT")
        return self._ev(_lib.EOSdPdT, rho, T)

    def GammaofRhoT(self, rho, T):
        """Shear modulus; ANEOS (0 for the gases, not Tillotson)."""
        self._require((EOSIDEALGAS, EOSANEOS, EOSREOS3, EOSSCVHEOS), "GammaofRhoT")
        return self._ev(_lib.EOSGammaofRhoT, rho, T)

    def UCold(self, rho):
        """u at the lowest tabulated temperature; NaN outside the table."""
        return self._ev(_lib.EOSUCold, rho)


# --- two-material functions ---------------------------------------------------


def _check_same_units(material1, material2):
    if (material1.dErgPerGmUnit != material2.dErgPerGmUnit
            or material1.dGmPerCcUnit != material2.dGmPerCcUnit):
        raise ValueError("materials %d and %d use different unit systems"
                         % (material1.iMat, material2.iMat))


def solveBC(material1, material2, rho1, u1, Tjump=0.0):
    """Boundary condition at a material interface: given (rho1, u1) in
    material1 return (rho2, u2) in material2 with equal pressure and
    T2 = T1 - Tjump >= 0. Raises RuntimeError on failure (the C library prints
    the reason to stderr only if compiled with -DEOSLIB_VERBOSE)."""
    _check_same_units(material1, material2)
    rho2, u2 = ct.c_double(), ct.c_double()
    ret = _lib.EOSSolveBC(material1._mat, material2._mat, float(rho1), float(u1),
                          ct.byref(rho2), ct.byref(u2), float(Tjump))
    if ret != EOS_SUCCESS:
        raise RuntimeError("EOSSolveBC failed for rho1=%g u1=%g Tjump=%g (iMat %d -> %d)"
                           % (rho1, u1, Tjump, material1.iMat, material2.iMat))
    return rho2.value, u2.value


def woolfsonCoeff(material1, material2, P, T):
    """f_ij = rho_1(P,T)/rho_2(P,T). 1 for the same material, below
    material1.woolfsonMinPressure() and in expanded states (rho <= 0.8 rho0);
    NaN where a density inversion fails (EOS_FAIL_VALUE with fail_to_nan = False)."""
    _check_same_units(material1, material2)
    m1, m2 = material1._mat, material2._mat
    return _vec(lambda P, T: _lib.EOSWoolfsonCoeff(m1, m2, P, T), (P, T))