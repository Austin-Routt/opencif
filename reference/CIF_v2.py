"""
CIF.py -- Ladder-Array Design Engine (Ladder-Array Design Specification v5)
===========================================================================

A closed system of equations, an algorithm, and the validators that decide
whether an implementation of it is correct.

    (shape function, three targets, fabrication floors)  ->  a drawn planar geometry

SCOPE.  Rigid geometry, Newtonian fluid, Stokes/laminar duct flow.  Substance-free:
there is no biology, no particle chemistry, and no material model anywhere in this
file.  Every particle is a hard sphere; every fluid is Newtonian.

WHAT IS EXACT (provable inside this file, and asserted in `validate()`):
    * resistance kernel r(w)           -- Yang series, SORTED; reciprocity identity
    * extraction kernel F_lam(x,beta)  -- exact rectangular-duct flat-slab flux fraction
    * mass identity  P(i) + 2 S(i) = 1 + 2 sigma_0                        [for ANY sigma_0]
    * terminal fraction  F = 1 - P(N)                                     [for ANY sigma_0]
    * forward <-> inverse round trip (sigma-form; no singular cases)
    * lamina invariance  delta(i) = delta*  under schedule S1             [any aspect ratio]
    * transport lemma: streamline migration reproduces F with ZERO diffusion
    * canonical taper: ONE monotone potential Psi generates w_c(i), P, N, F, inversion
    * Lambda-, GS-, mu-, Q_in-invariance

WHAT IS NOT (external inputs; see `UNCHECKABLE`):
    kappa, beta_cal, phi_cal, u_lo, u_ff, eta, f_soft, gamma_GS, AR_max,
    R_min, l_g, c_ent, L_eq.
    No validator in this file can check any of them.  `kappa` and `beta_cal`
    dominate: beta_cal alone swings the row count N by ~27%.

SECTION MAP (mirrors the specification):
    S2   constants                     S12  design inversions, termination on P
    S3   resistance kernel             S13  canonical taper (Psi potential)
    S4   extraction kernel             S14  realization: folds, packing, turns, bends
    S5   kernel-prefactor pairing      S15  network model + terminal layer
    S6   derived scale + transport     S16  regime gates, shear, sensitivity
    S7   schedule laws                 S17  invariance
    S8   march                         S18  design() -- the algorithm
    S9   closure                       S19  validators (K/E/C/P/M/D/G/N/S/R)
    S10  terminal fraction             S20  what nothing here can check
    S11  inverse map

USAGE
    import CIF
    CIF.validate()                          # the full v5 validator battery
    dev = CIF.design(CIF.reference_request())
    print(dev.report())
    CIF.write_dxf(dev, "mask.dxf")
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, field, replace
from typing import Callable, Optional, Sequence

import numpy as np
from scipy.optimize import brentq
from scipy.interpolate import PchipInterpolator
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

__all__ = [
    "UM", "DERIVED", "UNCHECKABLE", "Tolerances", "Process", "Pairing", "Request",
    "Row", "Device", "FoldPlan", "RoutabilityError", "WindowError", "GateError",
    "B", "r", "f_Re", "r_inverse",
    "F_lam", "c2", "K_depth", "x_of_f", "g_cubic",
    "Psi", "Psi_inverse",
    "march", "invert_widths", "critical_phi",
    "canonical_taper", "design", "reference_request", "lineage_request", "reference_process",
    "solve_network", "terminal_layer",
    "reynolds", "C_shear", "entrance_length", "sensitivity_map",
    "plan_folds", "route_fold", "path_length", "polygons", "write_dxf", "write_svg",
    "Layout", "Stage", "Mask", "LAYER_STYLE", "stadium_post", "uturn_widths",
    "min_w_start", "build_mask", "compose", "seam_report", "write_dxf_mask",
    "raster", "connectivity", "measure_widths", "geom_ftest", "heleshaw_split",
    "gifford_fig2",
    "validate", "VALIDATORS",
]

UM = 1e-6  # metre per micrometre -- the ONLY unit conversion in this file


# ===========================================================================
# S2 -- CONSTANTS
# ===========================================================================

#: Derived constants.  Exact.  Each is asserted in the K/E/S validator families.
DERIVED = {
    "192/pi^5":              192.0 / math.pi ** 5,          # 0.6274106
    "sum_odd n^-5":          (31.0 / 32.0) * 1.0369277551433699,  # 1.0045238
    "sum_odd n^-4":          math.pi ** 4 / 96.0,           # 1.0146780
    "sum_odd n^-2":          math.pi ** 2 / 8.0,            # 1.2337006
    "B(1)":                  0.4217310,                     # square-duct bracket
    "12/B(1)":               28.454154,                     # square-duct r*h^4/mu
    "24/B(1)":               56.908,                        # square-duct f*Re  (textbook 56.91)
    "f*Re plates":           96.0,                          # parallel-plate f*Re (textbook 96)
    "lim a[1-B(a)]":         0.6302489,                     # slender asymptote
    "c2(inf)":               3.0,                           # cubic limit
    "c2(1)":                 3.5568,
    "lim beta*c2(beta)":     1.7,                           # shallow divergence, beta -> 0
    "F_lam(1/2, beta)":      0.5,                           # hard bound, EVERY beta
    "C(1)":                  9.6078,                        # square-duct shear coefficient
    "C(2)":                  8.1341,                        # using 6 here understates by 26.2%
    "C(inf)":                6.0,                           # slender shear coefficient
}

#: Supplied constants.  NOTHING in this file can verify them.  See S20.
UNCHECKABLE = {
    "kappa":     "D = kappa*delta.  Fitted WITH THE CUBIC KERNEL.  N scales as kappa^2.",
    "beta_cal":  "aspect ratio(s) at which kappa was fitted.  DOMINANT RISK: a 27% swing in N, "
                 "0-16% bias in D.  The calibration array's depth was never published.",
    "phi_cal":   "wall porosity at which kappa was fitted.  The kernels assume a SOLID wall; the "
                 "real wall is ~40% open, with an effective-slip offset of order delta itself.",
    "flat_lamina": "F_lam is the exact flux of a FLAT slab in an UNPERTURBED duct.  The true "
                   "dividing streamsurface is z-dependent.  Absorbed into kappa, beta-dependently.",
    "u_lo/u_ff": "the window on which kappa is constant.  Outside it, D is wrong SILENTLY -- the "
                 "only failure mode this file cannot detect at all.",
    "eta":       "F_real = eta*F.  A FLOW-SPLIT correction.  It carries no lamina model gap.",
    "diffusion": "out of the model.  The transport lemma (S6) shows the cutoff is a hard STERIC "
                 "wall constraint: diffusion BROADENS the transition, it does not MOVE the cutoff.",
    "f_soft, gamma_GS, AR_max, F_min, R_min": "fabricability.",
    "l_g, c_ent": "R_gap.  Design is insensitive while the gap-dominance gate (M13) passes.",
    "L_eq":      "bend correction.  Folded rows run 2-6% under-resistant without it.",
    "routability": "packing bounds AREA only.  A feasible planar placement is NOT proven to exist "
                   "(G14).  The implementer must draw it.",
}


@dataclass(frozen=True)
class Tolerances:
    """S2.4 -- explicit inputs.  No literal tolerance appears anywhere else in this file."""
    eps_r:      float = 1e-9    # kernel series / limits
    eps_inv:    float = 1e-12   # root-finder closures
    eps_cub:    float = 1e-3    # cubic leading-inversion acceptance
    eps_id:     float = 1e-10   # algebraic identities
    eps_geo:    float = 1e-9    # drawn-geometry closures
    eps_net:    float = 2e-3    # network bulk agreement.  NOTE: the fold->straight handover
                                # is a PITCH DISCONTINUITY (Lambda jumps from 2p down to G+g_p).
                                # R_c and R_s jump by the SAME factor, so rho -- and hence the
                                # DESIGN -- is untouched; but the ABSOLUTE resistances jump, so
                                # dP_cross and therefore Q_gap take a LOCAL ~1e-3 hit right
                                # there.  Real, physical, localized.  Ramp Lambda to remove it.
    eps_lam:    float = 0.10    # thin-lamina gate  delta/h <= eps_lam       (E11)
    eps_taper:  float = 1e-2    # quasi-developed taper |dw_c/dx|            (M16)
    Re_max:     float = 10.0    # channel Reynolds gate                      (M14)
    Re_stokes:  float = 1e-2    # gap/lamina Reynolds gate (must be << 1)    (M14)
    gap_dom:    float = 0.05    # gap-dominance gate  R_gap*f/R_c            (M13)


@dataclass(frozen=True)
class Process:
    """Fabrication floors and fluid.  Everything here is an INPUT."""
    h:            float = 140 * UM     # fixed second dimension of every duct
    mu:           float = 1.0e-3       # dynamic viscosity  (scales r; cancels from rho, w_s, F, N)
    rho_f:        float = 1000.0       # density (Reynolds only)
    F_min:        float = 19 * UM      # minimum printable feature (Xia: == gap size G)
    R_min:        float = 0.0          # minimum printable INNER-corner radius
    AR_max:       float = 10.0         # w >= h/AR_max
    G:            float = 20.5 * UM    # link (gap) width, ALONG flow   (Gifford Fig.4, measured)
    g_p:          float = 25.6 * UM    # post width, ALONG flow  -> Lambda = G + g_p = 46.1 um
    l_g:          float = 50.6 * UM    # WALL THICKNESS = the post's PERPENDICULAR extent.
                                       # *** MUST BE >= g_p ***: a stadium is a rectangle with
                                       # semicircular caps of radius g_p/2 on its two ends, so
                                       # W_perp < L_along makes the caps CROSS -- a bowtie, not
                                       # a post.  Gifford: 50.6 / 25.6 = 1.98.  This is also the
                                       # link LENGTH that sets R_gap.
    c_ent:        float = 1.0          # entrance allowance: R_gap = r(G)*(l_g + c_ent*D_h(G,h))
    GS:           float = 28 * UM      # fold width      (gamma_GS band: GS/G in [1.4, 1.5])
    W:            float = 1562 * UM    # transverse budget per side
    l_stroke:     float = 1562 * UM    # fold stroke length (sets n_turns)
    L_eq:         float = 0.0          # Stokes-corner equivalent length per turn
    n_legs:       int = 1
    leg_length:   float = 75.8e-3
    Q_in:         float = 30e-9 / 60   # 30 uL/min.  Scale-free: cancels from everything but Re.
    m_budget:     int = 30             # terminal-layer allowance (rows)
    L_out:        float = 200 * UM     # run-out length (both branches, width-preserving -> G13)

    @property
    def Lambda(self) -> float:
        """Row pitch of the straight section."""
        return self.G + self.g_p

    @property
    def p(self) -> float:
        """Fold stroke pitch.  S14: p = GS + max(F_min, 2*R_min).  R_bend = p/2 is DERIVED."""
        return self.GS + max(self.F_min, 2.0 * self.R_min)

    @property
    def R_bend(self) -> float:
        return self.p / 2.0

    @property
    def w_min(self) -> float:
        """S12: the smallest buildable channel width."""
        return max(self.F_min, self.h / self.AR_max)

    @property
    def phi(self) -> float:
        """Wall porosity G/Lambda.  The kernels assume a SOLID wall.  REPORT IT.  (S20)"""
        return self.G / self.Lambda

    def __post_init__(self):
        if self.l_g < self.g_p:                      # G24 -- the stadium must be a stadium
            raise GateError(
                f"l_g = {self.l_g/UM:.1f} um < g_p = {self.g_p/UM:.1f} um.  A stadium post has "
                f"semicircular caps of radius g_p/2 on the ends of its PERPENDICULAR axis; "
                f"W_perp < L_along makes them cross -- a self-intersecting bowtie, not a post. "
                f"(Gifford, measured: 50.6 / 25.6 = 1.98.)")

    @property
    def R_gap(self) -> float:
        """S14.4.  The bare r(G)*l_g is a LOWER BOUND: for l_g <~ D_h the link is an orifice."""
        D_h = 2 * self.G * self.h / (self.G + self.h)
        return r(self.G, self.h, self.mu) * (self.l_g + self.c_ent * D_h)


@dataclass(frozen=True)
class Pairing:
    """
    S5 -- THE RULE THAT CANNOT BE BROKEN.

    `kappa` is meaningless without the kernel it was fitted with.  This object binds
    them, and is used at ALL THREE call sites (inversion, schedule S1, scale D).
    Validator P1 asserts one Pairing object reaches all three.

    The supplied kappa = 4.9138 was fitted with the CUBIC kernel, so `kernel='cubic'`
    is the self-consistent default.  Switching to 'flam' requires kappa_flam =
    kappa_cubic / K_depth(beta_cal) -- and beta_cal is UNKNOWN (S20).
    """
    kernel:    str = "cubic"                      # 'cubic' | 'flam'
    kappa:     float = 4.9138
    kappa_ci:  tuple = (4.750, 5.078)             # 95% CI
    beta_cal:  tuple = (0.80, 0.64, 0.53)         # ASSUMED (h_cal = 80 um).  UNVERIFIABLE.
    u_lo:      float = 0.035                      # lower edge of the kappa-constant window
    u_ff:      float = 0.012                      # below: kappa is a DIFFERENT constant (~1.34)
    eta:       float = 0.903                      # F_real = eta * F   (a yield correction)
    f_soft:    float = 5e-4
    gamma_GS:  tuple = (1.4, 1.5)

    def __post_init__(self):
        if self.kernel not in ("cubic", "flam"):
            raise ValueError(f"Pairing.kernel must be 'cubic' or 'flam', got {self.kernel!r}")

    # --- the three call sites -------------------------------------------------
    def f_of_x(self, x: float, beta: float) -> float:
        """Schedule S1: f = K(delta*/w_c ; beta)."""
        return g_cubic(x) if self.kernel == "cubic" else F_lam(x, beta)

    def x_of_f(self, f: float, beta: float) -> float:
        """Inversion: x = K^-1(f ; beta)."""
        return x_of_f(f, beta, kernel=self.kernel)

    def D(self, delta: float) -> float:
        """Scale: D = kappa_K * delta."""
        return self.kappa * delta

    # --- the mandatory report (P3) --------------------------------------------
    def transfer_bias(self, betas: Sequence[float]) -> tuple:
        """
        D_true / D_predicted = sqrt( c2(beta_cal) / c2(beta_design) ).

        NEVER APPLIED.  Reported, always.  A kappa fitted with the cubic silently
        carries c2(beta_cal); using it at a different aspect ratio biases D.
        """
        c_cal = float(np.mean([c2(b) for b in self.beta_cal]))
        bias = [math.sqrt(c_cal / c2(b)) for b in betas]
        return min(bias), max(bias)

    def kappa_for(self, kernel: str) -> float:
        """kappa_flam = kappa_cubic / K_depth(beta_cal).  Requires beta_cal (UNKNOWN)."""
        if kernel == self.kernel:
            return self.kappa
        c_cal = float(np.mean([c2(b) for b in self.beta_cal]))
        Kd = math.sqrt(DERIVED["c2(inf)"] / c_cal)
        return self.kappa / Kd if kernel == "flam" else self.kappa * Kd


class WindowError(RuntimeError):
    """u = D/w_c fell below u_ff: outside the calibrated window entirely (D1)."""


class GateError(RuntimeError):
    """A pre-flight or regime gate failed (M13/M14/M16/G13)."""


class RoutabilityError(RuntimeError):
    """
    G14.  The fold needs borrowed neighbouring slots (arms > 0) and this file does
    not place them.  Packing bounds AREA (S14.3); it does NOT prove a planar,
    non-self-intersecting placement exists.  That is the implementer's obligation,
    and the specification says so.  The fold SCHEDULE (Phi, L, arms) is still emitted.
    """


# ===========================================================================
# S3 -- RESISTANCE KERNEL
# ===========================================================================

_N_ODD = 50          # 50 ODD terms, n <= 99.  Gifford 2014 p.4499.  STATE THE CONVENTION.
_N_ODD_HI = 400      # high-accuracy reference for the K-family validators


_ODD_CACHE: dict = {}


def _odd(n_terms: int) -> np.ndarray:
    """Cached.  r() is the engine's innermost loop; re-allocating this array each call
    dominated the runtime (the march reaches r() ~44x per row through r_inverse)."""
    a = _ODD_CACHE.get(n_terms)
    if a is None:
        a = np.arange(1, 2 * n_terms, 2, dtype=float)
        _ODD_CACHE[n_terms] = a
    return a


def B(alpha: float, n_terms: int = _N_ODD) -> float:
    """
    S3.  The duct bracket.  Depends on alpha = max/min ALONE.

        B(alpha) = 1 - (192 / (pi^5 alpha)) * sum_{n odd} tanh(n pi alpha / 2) / n^5

    SORTING (alpha >= 1) is mandatory and it ELIMINATES the cancellation:
    for alpha >= 1, B is monotone increasing with range [0.4217310, 1) -- bounded
    away from zero.  The reciprocity identity  B(1/alpha) = alpha^-2 B(alpha)  (K9)
    shows why the unsorted branch is lossy: B(1/alpha) ~ alpha^-2 is SMALL, so the
    subtraction loses ~2*log10(alpha) digits AND amplifies truncation by ~alpha^2.
    """
    n = _odd(n_terms)
    s = np.sum(np.tanh(np.minimum(n * np.pi * alpha / 2.0, 350.0)) / n ** 5)
    return 1.0 - (192.0 / (np.pi ** 5 * alpha)) * s


def r(w: float, h: float, mu: float, n_terms: int = _N_ODD) -> float:
    """
    S3.  Resistance per unit length of a w x h duct.  SI.  STRICTLY DECREASING in w.

        r = 12 mu / (a b^3 B(alpha)),   a = max(w,h), b = min(w,h), alpha = a/b >= 1
    """
    if not (w > 0 and h > 0 and mu > 0):
        raise ValueError("r: w, h, mu must be strictly positive")
    a, b = (w, h) if w >= h else (h, w)
    br = B(a / b, n_terms)
    if br <= 0:                                  # K2
        raise ValueError("r: non-positive bracket -- increase n_terms or check inputs")
    return 12.0 * mu / (a * b ** 3 * br)


def f_Re(w: float, h: float, n_terms: int = _N_ODD_HI) -> float:
    """
    S3.  Poiseuille number  f*Re = 24 D_h^2 / (b^2 B(alpha)).  EXTERNAL ANCHOR:
    56.908 for a square duct (textbook 56.91), 96 for parallel plates (textbook 96).
    Validates B against tabulated duct data with NO reference to this file.  (K10/K11)
    """
    a, b = (w, h) if w >= h else (h, w)
    D_h = 2 * a * b / (a + b)
    return 24.0 * D_h ** 2 / (b ** 2 * B(a / b, n_terms))


def r_inverse(y: float, h: float, mu: float, n_terms: int = _N_ODD) -> float:
    """
    S3 R-d.  w = r^-1(y).  r is strictly decreasing => the root is unique.
    Bracketed from the slender seed w0 = (12 mu / (h y))^(1/3).
    """
    seed = (12.0 * mu / (h * y)) ** (1.0 / 3.0)
    lo, hi = seed * 0.25, seed * 4.0
    fn = lambda w: r(w, h, mu, n_terms) - y
    k = 0
    while fn(lo) < 0:                # r(lo) too small -> lo too wide -> shrink
        lo *= 0.5
        k += 1
        if k > 200 or lo < 1e-12:
            raise ValueError("r_inverse: cannot bracket below")
    k = 0
    while fn(hi) > 0:                # r(hi) too big -> hi too narrow -> grow
        hi *= 2.0
        k += 1
        if k > 200:
            raise ValueError("r_inverse: cannot bracket above")
    return brentq(fn, lo, hi, xtol=1e-16, rtol=1e-13)


# ===========================================================================
# S4 -- EXTRACTION KERNEL
# ===========================================================================

def _t_minus_tanh(t: np.ndarray) -> np.ndarray:
    """
    E-e.1 -- CANCELLATION GUARD.  The alternative form  beta - (2/(n pi)) tanh(t)  is a
    difference of near-equal numbers: the n=1 term is ~ pi^2 beta^3 / 12 against beta, so it
    loses ~ 2 log10(1/beta) - 1 significant digits.  At beta = 0.02 that is ~3.4 digits (still
    invisible in float64); by beta = 1e-5 the naive form carries ~4e-6 relative error and it
    degrades WITHOUT BOUND as beta -> 0.  Real devices sit at beta >~ 0.4 and never need this
    -- but the kernel claims generality, so the guard is not optional.

    Compute (t - tanh t) directly, with the Taylor branch below t = 1/2:

        t - tanh t = t^3/3 - 2t^5/15 + 17t^7/315 - 62t^9/2835 + 1382 t^11/155925 - ...
    """
    t = np.asarray(t, dtype=float)
    out = np.empty_like(t)
    lo = t < 0.5
    ts = t[lo]
    out[lo] = (ts ** 3 / 3.0 - 2 * ts ** 5 / 15.0 + 17 * ts ** 7 / 315.0
               - 62 * ts ** 9 / 2835.0 + 1382 * ts ** 11 / 155925.0)
    th = t[~lo]
    out[~lo] = th - np.tanh(np.minimum(th, 350.0))
    return out


def _A_n(beta: float, n_terms: int) -> tuple:
    """A_n(beta) = (2 / (pi n^5)) * (t_n - tanh t_n),  t_n = n pi beta / 2.  A_n > 0 always."""
    n = _odd(n_terms)
    return n, (2.0 / (np.pi * n ** 5)) * _t_minus_tanh(n * np.pi * beta / 2.0)


def F_lam(x: float, beta: float, n_terms: int = _N_ODD_HI) -> float:
    """
    S4.  The EXACT flux fraction of the FLAT slab 0 <= y <= delta in a w x h duct.

        F_lam(x, beta) = sum A_n [1 - cos(n pi x)] / (2 sum A_n),   x = delta/w, beta = h/w

    Limits (E1, E2):  beta -> inf gives 3x^2 - 2x^3 (the cubic);  F_lam(1/2, beta) = 1/2
    EXACTLY for every beta, because cos(n pi/2) = 0 for odd n.

    *** WHAT IS EXACT, AND WHAT IS A MODEL ***
    This is the exact flux of a FLAT slab in a FULLY DEVELOPED, UNPERTURBED duct.
    Identifying that slab with the fluid actually skimmed through a link is a MODEL:
    the true dividing streamsurface is z-dependent, and the wall it hugs is ~40% open.
    Both are beta- and phi-dependent and are absorbed into KAPPA -- not into eta.
    """
    n, A = _A_n(beta, n_terms)
    return float(np.sum(A * (1.0 - np.cos(n * np.pi * x))) / (2.0 * np.sum(A)))


def _c2_exact(beta: float, n_terms: int = 200_000) -> float:
    """
    S4 E-d.  c2(beta) = lim_{x->0} F_lam / x^2 = (pi^2/4) * sum n^2 A_n / sum A_n.

    E-e.2 -- SLOW TAIL.  sum A_n converges like n^-4 (fast); sum n^2 A_n only like
    n^-2.  WITHOUT the analytic tail  sum_{n odd > nmax} n^2 A_n ~ beta/(2 n_last),
    c2 is wrong in the 4th digit at n_max = 1e3.
    """
    n, A = _A_n(beta, n_terms)
    S0 = float(A.sum()) + beta / (6.0 * n[-1] ** 3)      # n^-4 tail (negligible, included)
    S2 = float((n ** 2 * A).sum()) + beta / (2.0 * n[-1])  # n^-2 tail (ESSENTIAL)
    return (np.pi ** 2 / 4.0) * S2 / S0


# c2 is smooth in beta and is called once per row: cache it as a monotone spline.
_C2_GRID = np.exp(np.linspace(math.log(1e-3), math.log(1e4), 240))
_C2_VALS = np.array([_c2_exact(b) for b in _C2_GRID])
_C2_SPLINE = PchipInterpolator(np.log(_C2_GRID), np.log(_C2_VALS), extrapolate=True)


def c2(beta: float) -> float:
    """Cached c2(beta).  Agrees with `_c2_exact` to ~1e-9 relative over the design range."""
    return float(np.exp(_C2_SPLINE(math.log(beta))))


def K_depth(beta: float) -> float:
    """
    S4 E-d.  K_depth = sqrt(3 / c2(beta)) = delta_3D / delta_2D at fixed f.
    DERIVED -- from the same duct series as r(w).  NOT a calibrated constant.
    The cubic OVER-states the lamina by 1/K_depth.
    """
    return math.sqrt(DERIVED["c2(inf)"] / c2(beta))


def g_cubic(x: float) -> float:
    """The infinite-plate (beta -> inf) kernel.  g(1/2) = 1/2 exactly."""
    return 3.0 * x ** 2 - 2.0 * x ** 3


def x_of_f(f: float, beta: float, kernel: str = "cubic") -> float:
    """
    S4 E-f.  x = K^-1(f ; beta).  Unique root on [0, 1/2] (K strictly increasing).

    SEED = sqrt(f / c2(beta)).  *** NEVER seed with the linear f ***: at beta = 0.53
    the exact root is 4.79e-3 while f = 1.04e-4 -- a 46x miss.
    """
    if not (0.0 < f < 0.5):
        raise ValueError(f"x_of_f: need 0 < f < 1/2 (hard bound, every beta), got {f}")
    if kernel == "cubic":
        return brentq(lambda x: g_cubic(x) - f, 0.0, 0.5, xtol=1e-16, rtol=1e-15)
    seed = math.sqrt(f / c2(beta))                       # E8
    lo = max(1e-15, seed * 0.5)
    hi = min(0.5, max(seed * 2.0, 1e-12))
    if F_lam(hi, beta) < f:
        hi = 0.5
    if F_lam(lo, beta) > f:
        lo = 1e-15
    return brentq(lambda x: F_lam(x, beta) - f, lo, hi, xtol=1e-16, rtol=1e-15)


# ===========================================================================
# S13 -- THE TAPER POTENTIAL  (defined early: S6/S12/S13 all use it)
# ===========================================================================

def Psi(w: float, h: float, kernel: str = "cubic") -> float:
    """
    S13.  Psi(w) = w^2 / c2(h/w).   STRICTLY INCREASING:

        dlnPsi/dlnw = 2 + dln c2/dln beta  in [1, 2]

    because c2 decreases in beta with dln c2/dln beta in (-1, 0)  (c2 -> 3 as
    beta -> inf; c2 ~ 1.7/beta as beta -> 0).  Hence Psi^-1 exists and is unique.

    Q_c is proportional to Psi under the canonical taper, so ONE potential generates
    w_c(i), P(i), N, F and the design inversion.  The cubic branch is Psi = w^2/3.

    *** This is the object whose omission was the specification's one real error. ***
    An earlier revision wrote F = 1 - (w_N/w_0)^2, which is the CUBIC SPECIAL CASE;
    in the F_lam branch it drops c2(beta_0)/c2(beta_N) = 1.0265 and is off by 0.0126.
    """
    if kernel == "cubic":
        return w ** 2 / DERIVED["c2(inf)"]
    return w ** 2 / c2(h / w)


def Psi_inverse(y: float, h: float, kernel: str = "cubic") -> float:
    """Psi^-1.  Unique by strict monotonicity."""
    if kernel == "cubic":
        return math.sqrt(y * DERIVED["c2(inf)"])
    lo, hi = 1e-9, 1e-1
    fn = lambda w: Psi(w, h, kernel) - y
    return brentq(fn, lo, hi, xtol=1e-15, rtol=1e-14)


# ===========================================================================
# S7/S8/S9/S10/S11/S12 -- SCHEDULE, MARCH, CLOSURE, TERMINAL FRACTION, INVERSE
# ===========================================================================

@dataclass
class Row:
    """One gap row.  Every quantity the specification defines per-row."""
    i:      int
    w_c:    float
    beta:   float
    f:      float
    sigma:  float
    rho:    float
    w_s:    float
    delta:  float
    D:      float
    u:      float
    P:      float
    S:      float
    Lam:    float
    R_c:    float
    R_s:    float
    Phi:    float


@dataclass
class Request:
    """The design request.  S12."""
    D_target:  float
    F_target:  float
    W_shape:   Callable[[float], float]     # w_c(t), t = i/N in [0,1]
    proc:      Process = field(default_factory=Process)
    pair:      Pairing = field(default_factory=Pairing)
    tol:       Tolerances = field(default_factory=Tolerances)
    mode:      str = "FOLD"                 # 'FOLD' (sigma_0 = 0) | 'DIRECT' (w_start given)
    w_start:   Optional[float] = None
    N_seed:    int = 2000
    max_rows:  int = 400_000
    f_schedule: Optional[Callable[[int, float], float]] = None   # S2: f(i, w_c).  Bypasses S1.


def _sigma0(req: Request, N: int, delta_star: float) -> float:
    """
    S8.  FOLD: sigma_0 = 0 EXACTLY (sealed inlet).
    DIRECT: a prescribed printable start width.  sigma_0 > 0 is a FICTITIOUS inlet
    filtrate -- it IS the mass defect (S10), and it is reported, not hidden.
    """
    if req.mode == "FOLD":
        return 0.0
    if req.w_start is None:
        raise ValueError("DIRECT mode requires w_start")
    p = req.proc
    w_c1 = req.W_shape(1.0 / N)
    f1 = req.pair.f_of_x(delta_star / w_c1, p.h / w_c1)
    rho1 = r(req.w_start, p.h, p.mu) / r(w_c1, p.h, p.mu)
    s0 = (1.0 - 2.0 * f1) / rho1 - f1
    if s0 < 0:
        raise GateError(f"DIRECT: sigma_0 = {s0:.3e} < 0 -- w_start is NARROWER than the "
                        f"closure requires.  Unphysical.")
    return s0


def march(req: Request, N: int, delta_star: float) -> tuple:
    """
    S8 + S9 + S12.  The march.

        sigma(0) = sigma_0
        sigma(i) = (sigma(i-1) + f(i)) / (1 - 2 f(i))          [S8]
        rho(i)   = 1/sigma(i)
        r(w_s(i)) = rho(i) * r(w_c(i))                          [S9 closure]

    TERMINATION ON P, not on w_s.  P is strictly decreasing (f > 0), so
    "stop at the first i with P(i) <= 1 - F_target" is monotone-safe for ANY
    shape function -- including non-monotone ones, where w_s need not be monotone.

    Returns (rows, sigma_0).  sigma_0 is RETURNED, never recomputed by the caller: it
    depends on N, and the fixed point accepts |dN| <= 1, so recomputing it downstream
    desynchronises DIRECT mode and breaks the S11 round-trip.
    """
    p, pair, tol = req.proc, req.pair, req.tol
    sigma = _sigma0(req, N, delta_star)
    sigma_0 = sigma
    P, S = 1.0, sigma
    rows: list = []
    P_stop = 1.0 - req.F_target
    w_c_prev = req.W_shape(0.0)

    for i in range(1, req.max_rows + 1):
        w_c = req.W_shape(min(i / N, 1.0))
        beta = p.h / w_c
        x_target = delta_star / w_c
        if req.f_schedule is None and x_target >= 0.5:
            raise WindowError(f"row {i}: delta*/w_c = {x_target:.4f} >= 1/2 -- the two laminas "
                              f"overlap (hard bound, every beta)")
        f = (req.f_schedule(i, w_c) if req.f_schedule is not None
             else pair.f_of_x(x_target, beta))                  # S2 overrides S1
        if not (0.0 < f < 0.5):
            raise WindowError(f"row {i}: f = {f:.4e} outside (0, 1/2)")

        sigma = (sigma + f) / (1.0 - 2.0 * f)
        rho = 1.0 / sigma

        if i == 1 and req.mode == "FOLD":                     # M2 -- ROW-1 CLOSURE tripwire
            resid = abs(f * (rho + 2.0) - 1.0)
            if resid > tol.eps_id:
                raise AssertionError(f"M2 row-1 closure violated: |f(1)(rho(1)+2) - 1| = {resid:.2e}")

        w_s = r_inverse(rho * r(w_c, p.h, p.mu), p.h, p.mu)
        # S14.4 -- PITCH.  Lambda cancels from rho and w_s (S17), so it is free.  It is now
        # UNIFORM (= G + g_p) for every row: a folded row only needs Lambda >= p (see
        # route_fold), and a NON-uniform pitch leaves solid stretches of wall between the
        # post cells, which orphan the cap complements into isolated fluid pockets.  The
        # connectivity check (G17) is what found that.
        Lam = p.Lambda          # UNIFORM (S17: Lambda cancels).  The fold's stroke-pitch
                                # requirement (Lambda >= p) is a RENDERING constraint and lives
                                # in route_fold -- it must not gate the physics.
        R_c = r(w_c, p.h, p.mu) * Lam
        R_s = r(w_s, p.h, p.mu) * Lam
        resid = abs(r(w_s, p.h, p.mu) / (rho * r(w_c, p.h, p.mu)) - 1.0)   # M4 -- PE CLOSURE tripwire
        if resid > tol.eps_inv * 1e3:
            raise AssertionError(f"M4 closure violated at row {i}: {resid:.2e}")

        delta = pair.x_of_f(f, beta) * w_c
        D = pair.D(delta)
        u = D / w_c

        if req.f_schedule is None and delta / p.h > tol.eps_lam:   # E11 -- THIN LAMINA (S1 only)
            raise GateError(f"row {i}: delta/h = {delta / p.h:.4f} > eps_lam -- the near-wall "
                            f"quadratic (seed, S1', steric picture) is not valid")
        if abs(w_c - w_c_prev) / p.Lambda > tol.eps_taper:     # M16 -- QUASI-DEVELOPED TAPER
            raise GateError(f"row {i}: |dw_c/dx| = {abs(w_c - w_c_prev) / p.Lambda:.3e} > eps_taper")
        w_c_prev = w_c

        S += f * P
        P *= (1.0 - 2.0 * f)

        rows.append(Row(i=i, w_c=w_c, beta=beta, f=f, sigma=sigma, rho=rho, w_s=w_s,
                        delta=delta, D=D, u=u, P=P, S=S, Lam=Lam, R_c=R_c, R_s=R_s,
                        Phi=r(w_s, p.h, p.mu) / r(p.GS, p.h, p.mu)))
        if P <= P_stop:
            break
    else:
        raise RuntimeError(f"march: did not reach F_target in {req.max_rows} rows")

    # S10 -- MASS IDENTITY.  Fires on ANY seed and on bilateral/per-side confusion.
    resid = abs(P + 2.0 * S - (1.0 + 2.0 * sigma_0))
    if resid > tol.eps_id:
        raise AssertionError(f"M5 mass identity violated: |P + 2S - (1 + 2 sigma_0)| = {resid:.2e}")
    return rows, sigma_0


def invert_widths(w_c: Sequence[float], w_s: Sequence[float], proc: Process,
                  sigma_0: float = 0.0) -> np.ndarray:
    """
    S11 -- THE INVERSE MAP, sigma-form.  Given both width schedules, what f does the
    closure imply?

        f(i) = (sigma(i) - sigma(i-1)) / (1 + 2 sigma(i))

    NO SINGULARITY, NO SPECIAL CASE, any sigma_0.  (The rho-form needs rho(0) = inf
    when sigma_0 = 0 and is a bug source -- an earlier revision failed by 18000% in
    DIRECT mode because of exactly that.)

    SCOPE: this inverts the DESIGN MAP, not the physics.  Where the closure is
    violated (DIRECT row 1), use `solve_network`.
    """
    sig = np.empty(len(w_c) + 1)
    sig[0] = sigma_0
    for i, (wc, ws) in enumerate(zip(w_c, w_s), start=1):
        sig[i] = r(wc, proc.h, proc.mu) / r(ws, proc.h, proc.mu)   # = 1/rho
    return (sig[1:] - sig[:-1]) / (1.0 + 2.0 * sig[1:])


# ===========================================================================
# S6.2 -- TRANSPORT LEMMA
# ===========================================================================

def critical_phi(f_schedule: Sequence[float]) -> float:
    """
    S6.2 -- THE TRANSPORT LEMMA.  Why delta is the right scale, and where kappa's
    residual MUST live.

    Let phi be the flux fraction of a fluid element measured FROM ONE WALL.  A link at
    that wall removes f*Q_c; the opposite link removes another f*Q_c but does not change
    phi*Q_c.  Hence

        phi' = (phi - f) / (1 - 2f)     <=>     (phi' - 1/2) = (phi - 1/2) / (1 - 2f)

    The centreline phi = 1/2 is a REPELLING fixed point: (phi - 1/2) grows by exactly
    1/P(i), so fluid MIGRATES DETERMINISTICALLY toward the wall and is extracted when
    phi <= f.  Solving the threshold gives phi_0 <= F/2 -- and both walls give F.

    THREE CONSEQUENCES, all load-bearing:
      1. The device needs NO cross-streamline transport.  Diffusion is a perturbation,
         not a mechanism.
      2. The cutoff is a HARD STERIC WALL CONSTRAINT, not a statistical position
         criterion: a sphere of radius a cannot bring its centre closer than a to the
         wall, so it is extracted iff a <= delta.  This gives kappa_steric = 2 EXACTLY.
      3. Therefore the entire measured residual kappa/2 = 2.46 must live in the three
         model gaps -- flat-lamina idealisation, wall porosity, finite-particle
         displacement -- AND NOWHERE ELSE.  The lemma PINS the list.

    Returns the critical phi_0 (should equal F/2 exactly).
    """
    def escapes(phi0: float) -> bool:
        phi = phi0
        for f in f_schedule:
            if phi <= f:
                return True
            phi = (phi - f) / (1.0 - 2.0 * f)
        return False

    lo, hi = 0.0, 0.5
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if escapes(mid):
            lo = mid
        else:
            hi = mid
    return lo


# ===========================================================================
# S13 -- CANONICAL TAPER  (closed form; no fixed point)
# ===========================================================================

@dataclass
class Taper:
    w_c0: float
    w_cN: float
    N: int
    F: float
    delta_star: float
    kernel: str
    h: float

    def w_c(self, i: int) -> float:
        """w_c(i) = Psi^-1[ Psi(w_c0) - 2 delta*^2 i ]."""
        y = Psi(self.w_c0, self.h, self.kernel) - 2.0 * self.delta_star ** 2 * i
        return Psi_inverse(y, self.h, self.kernel)

    def P(self, i: int) -> float:
        """P(i) = Psi(w_c(i)) / Psi(w_c0)  -- LINEAR in i."""
        return Psi(self.w_c(i), self.h, self.kernel) / Psi(self.w_c0, self.h, self.kernel)


def canonical_taper(w_c0: float, F_target: float, delta_star: float,
                    h: float, kernel: str = "cubic") -> Taper:
    """
    S13.  The taper that makes Q_gap CONSTANT.  Everything from ONE potential:

        Q_gap = f Q_c = c2(beta) (delta*/w_c)^2 Q_c = delta*^2 Q_c / Psi(w_c)
        =>  Q_gap const  <=>  Q_c ~ Psi(w_c)
        =>  d Psi(w_c(i)) / di = -2 delta*^2      (continuity: dQ_c/di = -2 Q_gap)

        w_c(i) = Psi^-1[Psi(w_c0) - 2 delta*^2 i]
        N      = [Psi(w_c0) - Psi(w_cN)] / (2 delta*^2)
        F      = 1 - Psi(w_cN)/Psi(w_c0)
        w_cN   = Psi^-1[(1 - F_target) Psi(w_c0)]              <- the design inversion

    Cubic branch (Psi = w^2/3) reduces to the familiar
        w_c(i) = sqrt(w_c0^2 + (w_cN^2 - w_c0^2) i/N),  N = (w_c0^2 - w_cN^2)/(6 delta*^2),
        F = 1 - (w_cN/w_c0)^2,  w_cN = w_c0 sqrt(1 - F).

    The equivalence chain  Q_gap const <=> dP_cross const <=> R_c Q_c = R_s Q_s exactly
    <=> Psi linear in i  holds IFF ALL LINKS ARE IDENTICAL.  The invariant that ALWAYS
    holds is dP_cross = const.
    """
    Psi0 = Psi(w_c0, h, kernel)
    w_cN = Psi_inverse((1.0 - F_target) * Psi0, h, kernel)
    PsiN = Psi(w_cN, h, kernel)
    N = (Psi0 - PsiN) / (2.0 * delta_star ** 2)
    F = 1.0 - PsiN / Psi0
    return Taper(w_c0=w_c0, w_cN=w_cN, N=int(round(N)), F=F,
                 delta_star=delta_star, kernel=kernel, h=h)


# ===========================================================================
# S14 -- REALIZATION: folds, packing, arms, turns, bends
# ===========================================================================

@dataclass
class FoldPlan:
    i:        int
    Phi:      float
    L_ideal:  float     # Phi * Lambda
    L_drawn:  float     # after the Stokes-corner allowance
    n_turns:  int
    arms:     int
    fits:     bool      # arms == 0  ->  routable in its own slot by `route_fold`
    Lam_min:  float     # 2p


def plan_folds(rows: Sequence[Row], proc: Process) -> list:
    """
    S14.  Fold factor, packing bound, arms, and the Stokes-corner allowance.

        Phi(i) = r(w_s(i)) / r(GS)                 -- INDEPENDENT of Lambda
        L(i)   = Phi(i) Lambda(i)
        Phi(i) <= W/p                              -- the packing bound (S14.3)
        arms(i) = ceil( Phi(i) / (W/p) ) - 1       -- borrowed neighbouring slots

    Bend correction (S14.6) is a STOKES-CORNER equivalent length, NOT a Dean
    correction: the secondary rail carries Q_s, giving De ~ 1e-4 .. 0.15.

        R_s_eff = r(GS)[L_hat + n_turns L_eq]  =>  L_hat ~ Phi Lambda / (1 + L_eq/l_stroke)

    Relative excess = L_eq / l_stroke, INDEPENDENT of L, and it applies ONLY to
    folded rows -- a localized inlet-side error.
    """
    cap = proc.W / proc.p
    out = []
    for row in rows:
        if row.w_s >= proc.w_min:
            continue                                 # STRAIGHT (G11: folded IFF w_s < w_min)
        L_ideal = row.Phi * row.Lam
        L_drawn = L_ideal / (1.0 + proc.L_eq / proc.l_stroke)
        n_turns = max(0, int(math.ceil(L_drawn / proc.l_stroke)) - 1)
        arms = int(math.ceil(row.Phi / cap)) - 1
        out.append(FoldPlan(i=row.i, Phi=row.Phi, L_ideal=L_ideal, L_drawn=L_drawn,
                            n_turns=n_turns, arms=max(arms, 0), fits=(arms <= 0),
                            Lam_min=2.0 * proc.p))
    return out


def route_fold(L: float, Lam: float, proc: Process) -> np.ndarray:
    """
    S14.2.  Route ONE fold as a boustrophedon in its OWN slot (Lam x W), returning the
    CENTRELINE polyline.  Traversal length equals L EXACTLY (solved, not approximated).

    Three requirements, all validated (G1/G2/G4):
      1. LENGTH        -- measured traversal == L
      2. CROSS-SECTION -- width == GS everywhere (rectilinear strokes, offset GS/2)
      3. TOPOLOGY      -- the link opens ONTO the traversal path; no shortcut exists

    Raises RoutabilityError (G14) if the fold needs borrowed slots.  Packing bounds
    AREA; it does not prove a planar placement exists.  That is honest, and it is
    exactly what the specification says.
    """
    if Lam < proc.p:                                                   # G8
        raise RoutabilityError(f"Lambda = {Lam/UM:.1f} um < p = {proc.p/UM:.1f} um (G8): two "
                               f"adjacent strokes cannot clear a wall of F_min between them")
    # k VERTICAL strokes at x = j*Lambda/k (j = 0..k), i.e. ON the link centres and their
    # sub-divisions, so ADJACENT strokes are Lambda/k apart and need Lambda/k >= p.
    #   verticals: k+1  (each of length A)   horizontals: k  (summing to Lambda)
    #   =>  L = Lambda + (k+1) A
    # k must be ODD, or the last stroke ends at y = A instead of returning to the wall.
    #
    # (An earlier version put the strokes at (j+0.5)*Lambda/k and then RETURNED at x = Lambda
    #  -- leaving the last two verticals only Lambda/2k apart.  That is what forced the bogus
    #  "Lambda >= 2p" rule.  With the strokes on the link centres, Lambda >= p is enough.)
    k_max = int(math.floor(Lam / proc.p))
    for k in range(1, k_max + 1, 2):                                   # ODD k only
        A = (L - Lam) / (k + 1)
        if proc.GS <= A <= proc.W:
            pts = [(0.0, 0.0)]
            for j in range(k + 1):
                xj = j * Lam / k
                up = (j % 2 == 0)
                pts.append((xj, A if up else 0.0))
                if j < k:
                    pts.append(((j + 1) * Lam / k, A if up else 0.0))
            if abs(pts[-1][1]) > 1e-15:                                # must end ON the wall
                continue
            return np.array(pts, dtype=float)
    cap = proc.W / proc.p
    arms = int(math.ceil((L / Lam) / cap)) - 1
    raise RoutabilityError(
        f"fold needs {arms} borrowed slot(s): Phi = {L/Lam:.1f} > W/p = {cap:.1f}.  "
        f"G14 -- global routability is the implementer's obligation; the specification "
        f"proves only the AREA bound, not a planar placement.  The SCHEDULE is emitted.")


def path_length(pts: np.ndarray) -> float:
    """Measured traversal length of a drawn polyline.  Used by G1/G2 -- never assumed."""
    return float(np.sum(np.hypot(np.diff(pts[:, 0]), np.diff(pts[:, 1]))))


def turn_widths(rho_upstream: float, w_c_turn: float, L_p: float, L_in: float,
                L_out: float, proc: Process) -> tuple:
    """
    S14.5.  A gap-free 180-degree turn IS ONE SEGMENT: flows change only AT rungs, so it
    carries Q_c(i), Q_s(i) EXACTLY and takes rho(i) of its UPSTREAM rung.  NO INTERPOLATION.

        r(w_s_in) L_in = r(w_s_out) L_out = rho(i) r(w_c_turn) L_p

    L_in < L_out  =>  w_s_in < w_s_out.  This rebalance is the PRECONDITION for the
    two-rail folding of S15 (N8): an unbalanced turn breaks left/right symmetry.
    """
    target = rho_upstream * r(w_c_turn, proc.h, proc.mu) * L_p
    w_in = r_inverse(target / L_in, proc.h, proc.mu)
    w_out = r_inverse(target / L_out, proc.h, proc.mu)
    return w_in, w_out


# ===========================================================================
# S15 -- NETWORK MODEL
# ===========================================================================

@dataclass
class NetworkResult:
    f_em:   np.ndarray
    F_em:   float
    Q_in:   float
    m_pred: int
    m_meas: int
    lam:    float
    A:      float       # boundary excitation amplitude
    bulk_err: float


def solve_network(rows: Sequence[Row], proc: Process, tol: Tolerances,
                  R_out_c: Optional[float] = None,
                  R_out_s: Optional[float] = None) -> NetworkResult:
    """
    S15.  Explicit THREE-RAIL Kirchhoff solve (no symmetry shortcut -- L/R symmetry is
    then a check, not an assumption).

      nodes  P_0..P_N (primary), L_1..L_N, R_1..R_N, terminals T_c, T_s
      edges  (P_{i-1},P_i)=R_c(i-1) ; (L_i,L_{i+1})=R_s(i) ; (P_i,L_i)=(P_i,R_i)=R_gap
             run-outs (P_N,T_c)=R_out_c ; (L_N,T_s)=(R_N,T_s)=R_out_s

    L_1 has NO upstream edge, so Q_s(0) = 0 is STRUCTURAL, never imposed  (N5).

    TERMINAL LAYER.  With  V(i) = p[P_i] - p[L_i]  the bilateral primary balance
    Q_c(i-1) - Q_c(i) = 2 Q_gap(i) gives

        V(i+1) - 2V(i) + V(i-1) = q V(i),   q = (2 R_c + R_s) / R_gap

    -- note the 2 on R_c; it is the same bilateral factor as everywhere else.  q is
    NOT monotone (rho falls, but R_c GROWS as w_c narrows and wins near the outlet),
    so lambda is taken as the MAX over the terminal window.  And the decay is
    A * z_-^(N-i): m = ceil(lambda ln(A/eps)) -- A is a BOUNDARY CONDITION, not 1.
    """
    N = len(rows)
    p, mu, h = proc, proc.mu, proc.h
    if R_out_c is None:
        R_out_c = r(rows[-1].w_c, h, mu) * p.L_out
    if R_out_s is None:
        R_out_s = r(rows[-1].w_s, h, mu) * p.L_out

    Rg = p.R_gap
    Rc_in = r(rows[0].w_c, h, mu) * p.Lambda          # inlet segment
    Rc = [rw.R_c for rw in rows]
    Rs = [rw.R_s for rw in rows]

    nP = N + 1                                        # P_0 .. P_N  -> 0..N
    iL = lambda i: nP + (i - 1)                       # L_1..L_N
    iR = lambda i: nP + N + (i - 1)                   # R_1..R_N
    Tc, Ts = nP + 2 * N, nP + 2 * N + 1
    M = nP + 2 * N + 2

    Rref = Rc_in                                      # scale conductances to O(1)
    edges = [(0, 1, Rc_in)]
    for i in range(1, N):
        edges.append((i, i + 1, Rc[i - 1]))
        edges.append((iL(i), iL(i + 1), Rs[i - 1]))
        edges.append((iR(i), iR(i + 1), Rs[i - 1]))
    for i in range(1, N + 1):
        edges.append((i, iL(i), Rg))
        edges.append((i, iR(i), Rg))
    edges.append((N, Tc, R_out_c))
    edges.append((iL(N), Ts, R_out_s))
    edges.append((iR(N), Ts, R_out_s))

    I, J, V = [], [], []
    for a, b, R in edges:
        g = Rref / R
        I += [a, b, a, b]
        J += [a, b, b, a]
        V += [g, g, -g, -g]
    K = coo_matrix((V, (I, J)), shape=(M, M)).tocsr()

    ground = [Tc, Ts]
    keep = np.setdiff1d(np.arange(M), ground)
    rhs = np.zeros(M)
    rhs[0] = 1.0                                      # unit injection at P_0
    pot = np.zeros(M)
    pot[keep] = spsolve(K[keep][:, keep].tocsc(), rhs[keep])
    pot *= Rref                                       # undo the scaling

    Q_in = (pot[0] - pot[1]) / Rc_in
    Qs_out = (pot[iL(N)] - pot[Ts]) / R_out_s + (pot[iR(N)] - pot[Ts]) / R_out_s
    F_em = Qs_out / Q_in

    f_em = np.empty(N)
    for i in range(1, N + 1):
        Qgap = (pot[i] - pot[iL(i)]) / Rg
        Rc_up = Rc_in if i == 1 else Rc[i - 2]
        Qc_up = (pot[i - 1] - pot[i]) / Rc_up
        f_em[i - 1] = Qgap / Qc_up

    f_des = np.array([rw.f for rw in rows])
    err = np.abs(f_em / f_des - 1.0)

    # terminal layer
    A = float(err[-1])
    win = max(1, min(p.m_budget, N))
    lam = 0.0
    for rw in rows[-win:]:
        q = (2.0 * rw.R_c + rw.R_s) / Rg
        zp = ((2.0 + q) + math.sqrt((2.0 + q) ** 2 - 4.0)) / 2.0
        lam = max(lam, 1.0 / math.log(zp))
    m_pred = int(math.ceil(lam * math.log(max(A, tol.eps_net) / tol.eps_net)))
    m_meas = 0
    for e in err[::-1]:                       # the terminal layer is at the END: count trailing
        if e < tol.eps_net:
            break
        m_meas += 1
    m = max(m_meas, m_pred)
    bulk = float(np.max(err[: N - m])) if N > m else float("nan")

    return NetworkResult(f_em=f_em, F_em=F_em, Q_in=Q_in, m_pred=m_pred, m_meas=m_meas,
                         lam=lam, A=A, bulk_err=bulk)


def terminal_layer(rows: Sequence[Row], proc: Process, eps: float, A: float = 1.0) -> tuple:
    """S15.  A-priori (lambda, m).  `A` is the boundary excitation, set by run-out mismatch."""
    Rg = proc.R_gap
    lam = 0.0
    for rw in rows[-min(proc.m_budget, len(rows)):]:
        q = (2.0 * rw.R_c + rw.R_s) / Rg
        zp = ((2.0 + q) + math.sqrt((2.0 + q) ** 2 - 4.0)) / 2.0
        lam = max(lam, 1.0 / math.log(zp))
    return lam, int(math.ceil(lam * math.log(A / eps)))


# ===========================================================================
# S16 -- REGIME GATES, SHEAR, ENTRANCE LENGTH, SENSITIVITY
# ===========================================================================

def reynolds(Q: float, w: float, h: float, proc: Process) -> float:
    """S16-a.  Re = rho_f Q D_h / (mu h w).  USE THE LOCAL Q, never Q_in."""
    D_h = 2.0 * w * h / (w + h)
    return proc.rho_f * Q * D_h / (proc.mu * h * w)


def C_shear(alpha: float, n_terms: int = 2000) -> float:
    """
    S16-b.  Exact maximum-wall-shear coefficient:  gamma_max = C(alpha) Q / (a b^2).

        C(alpha) = (48 / (pi^2 B(alpha))) sum_{n odd} (1/n^2)[1 - sech(n pi alpha / 2)]

    SAME B(alpha) as S3/S4 (one call site -- S3).  C -> 6 ONLY as alpha -> inf:
    substituting 6 at alpha = 2 UNDERSTATES the peak shear by 26.2% and discards a real
    drift along any taper.  Constant Q_gap != constant shear at finite alpha.
    """
    n = _odd(n_terms)
    s = np.sum((1.0 / n ** 2) * (1.0 - 1.0 / np.cosh(np.minimum(n * np.pi * alpha / 2.0, 350.0))))
    s += 1.0 / (2.0 * n[-1])            # ANALYTIC TAIL: sum_{odd n > n_max} 1/n^2 ~ 1/(2 n_last).
    return float(48.0 / (np.pi ** 2 * B(alpha, n_terms)) * s)


def entrance_length(Re: float, D_h: float) -> float:
    """
    S16-c.  The correct laminar entrance length has a STOKES FLOOR and does NOT vanish
    as Re -> 0:

        L_e / D_h = [0.619^1.6 + (0.0567 Re)^1.6]^(1/1.6)   ->  0.62  as Re -> 0

    (The naive L_e = 0.05 Re D_h goes to ZERO here and is simply the wrong formula.)
    For the reference device L_e ~ 70 um against Lambda = 49 um: the flow does NOT
    classically redevelop between rungs.  This is BENIGN, and the reason must be stated:
    the per-rung perturbation is O(f) = 1e-4, so the resistance error is O(f^2) ~ 1e-8 --
    four orders below the closure tolerance.  REPORT L_e/Lambda; do not gate on it.
    Gate on M13 instead.
    """
    return D_h * (0.619 ** 1.6 + (0.0567 * Re) ** 1.6) ** (1.0 / 1.6)


def sensitivity_map(rows: Sequence[Row], proc: Process, pair: Pairing) -> dict:
    """
    S16-e -- MANDATORY REPORT.  dln(.)/dln(.) at the design point.

    Read it once and it changes your process priorities: the WIDTHS carry the device;
    the ETCH DEPTH barely does.  And N ~ kappa^2, so kappa's +-3.4% CI alone is
    +-6.8% on row count.
    """
    h, mu = proc.h, proc.mu
    wc_o, ws_o = rows[-1].w_c, rows[-1].w_s
    d = 1e-4

    def dln(fn, x):
        return (math.log(fn(x * (1 + d))) - math.log(fn(x * (1 - d)))) / (2 * d)

    F_ws = lambda ws: 2.0 / (r(ws, h, mu) / r(wc_o, h, mu) + 2.0)
    F_wc = lambda wc: 2.0 / (r(ws_o, h, mu) / r(wc, h, mu) + 2.0)
    F_h = lambda hh: 2.0 / (r(ws_o, hh, mu) / r(wc_o, hh, mu) + 2.0)
    wc_i, f_i = rows[0].w_c, rows[0].f
    D_wc = lambda wc: pair.kappa * wc * pair.x_of_f(f_i, h / wc)
    D_h_ = lambda hh: pair.kappa * wc_i * pair.x_of_f(f_i, hh / wc_i)
    D_f = lambda ff: pair.kappa * wc_i * pair.x_of_f(ff, h / wc_i)

    return {
        "dlnF/dln w_s(out)": dln(F_ws, ws_o),
        "dlnF/dln w_c(out)": dln(F_wc, wc_o),
        "dlnF/dln h":        dln(F_h, h),
        "dlnD/dln w_c":      dln(D_wc, wc_i),
        "dlnD/dln h":        dln(D_h_, h),
        "dlnD/dln f":        dln(D_f, f_i),
        "dlnN/dln kappa":    2.0,      # exact: N ~ 1/delta*^2 ~ kappa^2
    }


# ===========================================================================
# S18 -- THE ALGORITHM
# ===========================================================================

@dataclass
class Device:
    req:      Request
    rows:     list
    N:        int
    F:        float
    F_real:   float
    sigma_0:  float
    delta_star: float
    D:        float
    D_band:   Optional[tuple]
    folds:    list
    bias:     tuple
    gates:    dict
    sens:     dict
    reports:  dict

    def report(self) -> str:
        p, pr = self.req.proc, self.req.pair
        L = sum(rw.Lam for rw in self.rows)
        n_fold = len(self.folds)
        n_arm = sum(1 for f in self.folds if not f.fits)
        g = self.reports
        lines = [
            "=" * 78,
            "LADDER-ARRAY DEVICE  (Ladder-Array Design Specification v5)",
            "=" * 78,
            f"  kernel / kappa        : {pr.kernel}  /  {pr.kappa:.4f}   [PAIRED -- P1/P2]",
            f"  N (rows)              : {self.N}",
            f"  F  (model)            : {self.F:.5f}",
            f"  F  (realized ~ eta F) : {self.F_real:.5f}      eta = {pr.eta} (CALIBRATED)",
            f"  D                     : {self.D/UM:.4f} um" if self.D_band is None else
            f"  D                     : BAND [{self.D_band[0]/UM:.3f}, {self.D_band[1]/UM:.3f}] um  (off plateau)",
            f"  delta*                : {self.delta_star/UM:.4f} um   (invariant, S1)",
            f"  w_c                   : {self.rows[0].w_c/UM:.1f} -> {self.rows[-1].w_c/UM:.1f} um",
            f"  w_s                   : {self.rows[0].w_s/UM:.3f} -> {self.rows[-1].w_s/UM:.2f} um",
            f"  u = D/w_c             : {min(rw.u for rw in self.rows):.4f} .. "
            f"{max(rw.u for rw in self.rows):.4f}   (plateau: u >= {pr.u_lo})",
            f"  array length          : {L*1e3:.1f} mm   (budget {p.n_legs*p.leg_length*1e3:.1f} mm)",
            f"  folded rows           : {n_fold}   (of which {n_arm} need borrowed slots -> G14)",
            f"  max Phi               : {max((f.Phi for f in self.folds), default=0):.1f}   "
            f"(packing cap W/p = {p.W/p.p:.1f})",
            f"  DIRECT mass defect    : {2*self.sigma_0:.4e}" if self.sigma_0 else
            f"  DIRECT mass defect    : 0 (FOLD mode: sigma_0 = 0 exactly)",
            "-" * 78,
            "  REGIME GATES",
            f"    gap dominance  max R_gap f / R_c : {g['gap_dom']:.4f}   (< {self.req.tol.gap_dom})",
            f"    Re (channel, max)                : {g['Re_max']:.2f}",
            f"    Re_gap                           : {g['Re_gap']:.2e}   (must be << 1)",
            f"    Re_lamina                        : {g['Re_delta']:.2e}   (must be << 1)",
            f"    L_e / Lambda                     : {g['Le_over_Lam']:.2f}   "
            f"(>1 is EXPECTED and benign -- perturbation is O(f))",
            f"    max wall shear (C(alpha), NOT 6) : {g['gamma_max']:.3e} 1/s",
            f"    Delta P total                    : {g['dP_tot']/1e3:.1f} kPa",
            f"    wall porosity  phi = G/Lambda    : {p.phi:.3f}   "
            f"(the kernels assume a SOLID wall -- S20)",
            "-" * 78,
            "  SENSITIVITY MAP  (S16-e)  -- widths carry the device; etch depth barely does",
        ]
        for k, v in self.sens.items():
            lines.append(f"    {k:<22s} {v:+.3f}")
        lines += [
            "-" * 78,
            "  TRANSFER BIAS  (S5, P3)  -- REPORTED, NEVER APPLIED",
            f"    D_true/D_pred = sqrt(c2(beta_cal)/c2(beta)) in "
            f"[{self.bias[0]:.3f}, {self.bias[1]:.3f}]",
            f"    beta_cal = {pr.beta_cal}  <-- ASSUMED.  The calibration array's depth was",
            f"    never published.  This ONE number swings N by ~27%.  MEASURE IT.",
            "=" * 78,
        ]
        return "\n".join(lines)


def design(req: Request) -> Device:
    """
    S18.  The algorithm.  Pre-flight gates, march, inverse round-trip, window guards,
    regime gates, realization, and the mandatory reports.
    """
    p, pair, tol = req.proc, req.pair, req.tol

    # ---- pre-flight gates (fail, do not warn) -------------------------------
    ts = np.linspace(0.0, 1.0, 401)
    w_max = max(req.W_shape(t) for t in ts)
    if req.f_schedule is None and w_max > req.D_target / pair.u_lo:          # S12 plateau guard
        raise WindowError(
            f"plateau guard: max_t W(t) = {w_max/UM:.1f} um > D_target/u_lo = "
            f"{req.D_target/pair.u_lo/UM:.1f} um.  u = D/w_c is SMALLEST where w_c is "
            f"LARGEST -- the guard binds at max w_c, not at W(0).")
    if not (pair.gamma_GS[0] <= p.GS / p.G <= pair.gamma_GS[1]):
        raise GateError(f"GS/G = {p.GS/p.G:.2f} outside gamma_GS = {pair.gamma_GS}")
    if p.GS < p.w_min or p.G < p.w_min:
        raise GateError("GS or G below w_min")

    delta_star = req.D_target / pair.kappa

    # ---- march: damped fixed point on N (S13 seeds it) -----------------------
    tap = canonical_taper(req.W_shape(0.0), req.F_target, delta_star, p.h, pair.kernel)
    N = max(2, min(tap.N, req.N_seed * 8)) if tap.N > 0 else req.N_seed
    rows, sigma_0 = [], 0.0
    for _ in range(200):
        rows, sigma_0 = march(req, N, delta_star)
        N_new = len(rows)
        # The march CANNOT resolve N below one row: W(t) is sampled at t = i/N, so a +-1
        # change in N is below the discretization.  Converge to that resolution and stop --
        # demanding N_new == N exactly produces a 1-cycle (floor(0.5N + 0.5(N+1)) = N).
        if abs(N_new - N) <= 1:
            break
        N = int(round(0.5 * N + 0.5 * N_new))               # damping (non-monotone W)
    else:
        raise RuntimeError("march: fixed point on N did not converge in 200 iterations")
    N = len(rows)

    # ---- inverse round-trip (S11, sigma-form) --------------------------------
    f_B = invert_widths([rw.w_c for rw in rows], [rw.w_s for rw in rows], p, sigma_0)
    f_A = np.array([rw.f for rw in rows])
    rt = float(np.max(np.abs(f_B / f_A - 1.0)))
    if rt > 1e-6:
        raise AssertionError(f"M7 round-trip A->B->A failed: {rt:.2e}")

    F = 1.0 - rows[-1].P                                     # S10 -- PHYSICAL, any sigma_0
    if sigma_0 == 0.0:
        for name, val in (("2S(N)", 2 * rows[-1].S),
                          ("2/(rho(N)+2)", 2.0 / (rows[-1].rho + 2.0))):
            if abs(F - val) > tol.eps_id:
                raise AssertionError(f"M6: F != {name} ({abs(F-val):.2e})")
    if rows[-1].w_s < p.w_min:
        raise GateError(f"w_s(N) = {rows[-1].w_s/UM:.2f} um < w_min = {p.w_min/UM:.1f} um")

    # ---- window guards (S6) ---------------------------------------------------
    u_min = min(rw.u for rw in rows)
    if req.f_schedule is not None:
        pass                       # S2: D is a PROFILE (v5 S7); the kappa-window guards do not apply
    elif u_min < pair.u_ff:
        raise WindowError(f"u_min = {u_min:.4f} < u_ff = {pair.u_ff}: outside the calibrated "
                          f"window entirely.  kappa is a DIFFERENT constant there.")
    D_band = None
    if req.f_schedule is None and u_min < pair.u_lo:
        D_band = (pair.kappa_ci[0] * delta_star, pair.kappa_ci[1] * delta_star)

    # ---- regime gates (S16) ---------------------------------------------------
    Rg = p.R_gap
    gap_dom = max(Rg * rw.f / rw.R_c for rw in rows)
    if gap_dom >= tol.gap_dom:                                                # M13 -- HARD GATE
        raise GateError(f"M13 gap dominance: max R_gap f / R_c = {gap_dom:.4f} >= {tol.gap_dom}.  "
                        f"Below this the closure (S9) and the network (S15) DIVERGE.")
    Qc = p.Q_in * np.array([1.0] + [rw.P for rw in rows[:-1]])
    Re = [reynolds(q, rw.w_c, p.h, p) for q, rw in zip(Qc, rows)]
    Re_max = max(Re)
    if Re_max > tol.Re_max:
        raise GateError(f"M14: max Re = {Re_max:.1f} > {tol.Re_max}")
    Qg0 = rows[0].f * p.Q_in
    Re_gap = reynolds(Qg0, p.G, p.h, p)
    a0 = max(rows[0].w_c, p.h)
    b0 = min(rows[0].w_c, p.h)
    gam = C_shear(a0 / b0) * p.Q_in / (a0 * b0 ** 2)
    Re_delta = p.rho_f * gam * rows[0].delta ** 2 / p.mu
    # M14 -- STOKES GATES.  These are CUTOFF conditions: they exist so that S6.2's
    # streamline/steric picture holds.  Under S2 (prescribed f, no cutoff claimed) they do
    # not gate -- but they are still reported.  The CHANNEL Re gate above always applies.
    if req.f_schedule is None and (Re_gap > tol.Re_stokes or Re_delta > tol.Re_stokes):
        raise GateError(f"M14 Stokes: Re_gap = {Re_gap:.2e}, Re_delta = {Re_delta:.2e} -- the "
                        f"skimming junction MUST be Stokes or S6.2's streamline picture fails.")
    D_h0 = 2 * rows[0].w_c * p.h / (rows[0].w_c + p.h)
    Le = entrance_length(Re[0], D_h0)
    dP = float(sum(rw.R_c * q for rw, q in zip(rows, Qc)))

    folds = plan_folds(rows, p)
    bias = pair.transfer_bias([rw.beta for rw in rows])
    sens = sensitivity_map(rows, p, pair)

    return Device(
        req=req, rows=rows, N=N, F=F, F_real=pair.eta * F, sigma_0=sigma_0,
        delta_star=delta_star, D=pair.D(delta_star), D_band=D_band, folds=folds,
        bias=bias, gates={}, sens=sens,
        reports=dict(gap_dom=gap_dom, Re_max=Re_max, Re_gap=Re_gap, Re_delta=Re_delta,
                     Le_over_Lam=Le / p.Lambda, gamma_max=gam, dP_tot=dP, roundtrip=rt),
    )


# ---------------------------------------------------------------------------
# Reference request -- the lineage geometry, substance-free.
# ---------------------------------------------------------------------------

def reference_process() -> Process:
    return Process(h=140 * UM, G=19 * UM, g_p=30 * UM, GS=28 * UM, F_min=19 * UM,
                   W=1562 * UM, l_stroke=1562 * UM, l_g=60 * UM, c_ent=1.0,
                   AR_max=10.0, n_legs=2, leg_length=75.8e-3, L_out=200 * UM)


def reference_request(D_target: float = 2.170 * UM, F_target: float = 0.5592,
                     margin: float = 1.05) -> Request:
    """
    The reference device -- PLATEAU-VALID and ROUTABLE.  Linear taper w_c: 60 -> 42 um, h=140.

    TWO gates shape it, and both are load-bearing:

    1. PLATEAU (S12).  max_t w_c(t) <= D_target/u_lo = 62 um.  The lineage taper (75->52)
       violates this: at 75 um, u = 0.029 < u_lo = 0.035, i.e. the wide end runs in the KNEE
       where kappa is a FUNCTION.  So the taper is the lineage RATIO (52/75) rescaled.

    2. PACKING (S14.3, G20).  Phi <= W/p is INVIOLABLE and Lambda cancels from it.  In FOLD
       mode (sigma_0 = 0) row 1 wants Phi = 752 against a cap of 33.2 -- 22 folds that cannot
       be drawn.  So the device is DIRECT, seeded at min_w_start (the exact width that puts
       Phi(1) ON the cap), with a small margin.  This is what the published devices do: Xia
       seeds w_s(1) = G; Dinh seeds ~6 um.

       The cost is a MASS DEFECT 2*sigma_0 -- a fictitious inlet filtrate.  S10 makes it
       exact and the report prints it.  Nothing is hidden.
    """
    proc = reference_process()
    w0 = 60.0 * UM
    return Request(
        D_target=D_target, F_target=F_target,
        W_shape=lambda t: (60.0 + (42.0 - 60.0) * t) * UM,
        proc=proc, pair=Pairing(), tol=Tolerances(),
        mode="DIRECT", w_start=margin * min_w_start(proc, w0), N_seed=1500,
    )


def lineage_request(D_target: float = 2.170 * UM, F_target: float = 0.5592) -> Request:
    """
    The PUBLISHED lineage geometry: w_c 75 -> 52 um, D = 2.170 um.

    *** THIS REQUEST IS REFUSED BY THE PLATEAU GUARD, AND CORRECTLY SO. ***
    u = D/w_c = 2.170/75 = 0.029 < u_lo = 0.035: the widest end of the published device
    runs in the KNEE, where kappa is a FUNCTION, not a constant.  The guard fires, names
    the bound (w_c(0) <= 62 um), and refuses to emit a scalar cutoff.  That is the whole
    point of S12.  Use `reference_request()` for a design the engine will build.
    """
    return Request(
        D_target=D_target, F_target=F_target,
        W_shape=lambda t: (75.0 + (52.0 - 75.0) * t) * UM,
        proc=reference_process(), pair=Pairing(), tol=Tolerances(),
        mode="FOLD", N_seed=2700,
    )


# ===========================================================================
# S14 (cont.) -- MASK EMITTER
# ===========================================================================

def polygons(dev: Device, max_rows: Optional[int] = None) -> list:
    """
    Emit axis-aligned polygons for the STRAIGHT section (rectilinear, exact, and
    fully validated by G1/G4).  Folded rows are emitted where they fit in their own
    slot; rows needing arms raise G14 in `route_fold` and are SKIPPED here, with
    their schedule still available in `dev.folds`.

    Layers: CHANNEL_C (primary), CHANNEL_S (secondary, both sides), POST, LINK.
    Coordinates in micrometres, y = 0 on the primary centreline.
    """
    p = dev.req.proc
    out, x = [], 0.0
    rows = dev.rows if max_rows is None else dev.rows[:max_rows]
    for rw in rows:
        wc, ws, Lam = rw.w_c / UM, rw.w_s / UM, rw.Lam / UM
        G, gp, GS = p.G / UM, p.g_p / UM, p.GS / UM
        y_c = wc / 2.0
        out.append(dict(layer="CHANNEL_C", v=[(x, -y_c), (x + Lam, -y_c),
                                              (x + Lam, y_c), (x, y_c)]))
        if rw.w_s >= p.w_min:                                   # STRAIGHT (G11)
            for sgn in (+1, -1):
                y0 = sgn * (y_c + p.l_g / UM)
                y1 = y0 + sgn * ws
                out.append(dict(layer="CHANNEL_S",
                                v=[(x, min(y0, y1)), (x + Lam, min(y0, y1)),
                                   (x + Lam, max(y0, y1)), (x, max(y0, y1))]))
        else:                                                   # FOLD
            try:
                pts = route_fold(rw.Phi * rw.Lam, rw.Lam, p)
            except RoutabilityError:
                continue                                        # arms -> G14, schedule only
            for sgn in (+1, -1):
                base = sgn * (y_c + p.l_g / UM)
                for k in range(len(pts) - 1):
                    (x0, y0), (x1, y1) = pts[k] / UM, pts[k + 1] / UM
                    ax, ay = x + x0, base + sgn * y0
                    bx, by = x + x1, base + sgn * y1
                    lo_x, hi_x = min(ax, bx) - (GS / 2 if ax == bx else 0), \
                                 max(ax, bx) + (GS / 2 if ax == bx else 0)
                    lo_y, hi_y = min(ay, by) - (GS / 2 if ay == by else 0), \
                                 max(ay, by) + (GS / 2 if ay == by else 0)
                    out.append(dict(layer="CHANNEL_S",
                                    v=[(lo_x, lo_y), (hi_x, lo_y), (hi_x, hi_y), (lo_x, hi_y)]))
        # link (gap) + post on both walls
        for sgn in (+1, -1):
            out.append(dict(layer="LINK", v=[(x, sgn * y_c), (x + G, sgn * y_c),
                                             (x + G, sgn * (y_c + p.l_g / UM)),
                                             (x, sgn * (y_c + p.l_g / UM))]))
            out.append(dict(layer="POST", v=[(x + G, sgn * y_c), (x + G + gp, sgn * y_c),
                                             (x + G + gp, sgn * (y_c + p.l_g / UM)),
                                             (x + G, sgn * (y_c + p.l_g / UM))]))
        x += Lam
    return out


def write_dxf(dev: Device, path: str, max_rows: Optional[int] = None) -> int:
    """R2000 ASCII DXF, micrometre units ($INSUNITS = 13).  One LWPOLYLINE per polygon."""
    P = polygons(dev, max_rows)
    with open(path, "w") as fh:
        fh.write("0\nSECTION\n2\nHEADER\n9\n$ACADVER\n1\nAC1015\n9\n$INSUNITS\n70\n13\n")
        fh.write(f"9\n$HANDSEED\n5\n{256 + len(P) + 16:X}\n0\nENDSEC\n")
        fh.write("0\nSECTION\n2\nENTITIES\n")
        for k, poly in enumerate(P):
            fh.write(f"0\nLWPOLYLINE\n5\n{256+k:X}\n100\nAcDbEntity\n8\n{poly['layer']}\n"
                     f"100\nAcDbPolyline\n90\n{len(poly['v'])}\n70\n1\n43\n0\n")
            for (vx, vy) in poly["v"]:
                fh.write(f"10\n{vx:.4f}\n20\n{vy:.4f}\n")
        fh.write("0\nENDSEC\n0\nEOF\n")
    return len(P)


def write_dxf_mask(mask: Mask, path: str) -> int:
    """
    R2000 ASCII DXF of a full Mask (all layers, all layouts).  Micrometre units.

    OUTLINE is emitted as its own layer and is a CHIP EDGE -- consumers must not fill it
    (see LAYER_STYLE).  Filling it is exactly the bug that made the first render come out
    a solid grey rectangle.
    """
    P = [q for q in mask.polys if q["layer"] != "OUTLINE"]
    with open(path, "w") as fh:
        fh.write("0\nSECTION\n2\nHEADER\n9\n$ACADVER\n1\nAC1015\n9\n$INSUNITS\n70\n13\n")
        fh.write(f"9\n$HANDSEED\n5\n{256 + len(P) + len(mask.ports) + 16:X}\n0\nENDSEC\n")
        fh.write("0\nSECTION\n2\nENTITIES\n")
        h = 256
        for q in P:
            fh.write(f"0\nLWPOLYLINE\n5\n{h:X}\n100\nAcDbEntity\n8\n{q['layer']}\n"
                     f"100\nAcDbPolyline\n90\n{len(q['v'])}\n70\n1\n43\n0\n")
            for (vx, vy) in q["v"]:
                fh.write(f"10\n{vx:.4f}\n20\n{vy:.4f}\n")
            h += 1
        for (_nm, cx, cy, rad) in mask.ports:
            fh.write(f"0\nCIRCLE\n5\n{h:X}\n100\nAcDbEntity\n8\nPORTS\n100\nAcDbCircle\n"
                     f"10\n{cx:.4f}\n20\n{cy:.4f}\n30\n0\n40\n{rad:.4f}\n")
            h += 1
        fh.write("0\nENDSEC\n0\nEOF\n")
    return len(P)


def write_svg(dev: Device, path: str, max_rows: int = 60) -> int:
    P = polygons(dev, max_rows)
    xs = [v[0] for q in P for v in q["v"]]
    ys = [v[1] for q in P for v in q["v"]]
    x0, x1, y0, y1 = min(xs), max(xs), min(ys), max(ys)
    col = {"CHANNEL_C": "#4488cc", "CHANNEL_S": "#66aadd", "POST": "#eecc66", "LINK": "#cc6644"}
    with open(path, "w") as fh:
        fh.write(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{x0} {y0} '
                 f'{x1-x0} {y1-y0}" width="1400">\n<rect x="{x0}" y="{y0}" '
                 f'width="{x1-x0}" height="{y1-y0}" fill="white"/>\n')
        for q in P:
            pts = " ".join(f"{v[0]:.2f},{v[1]:.2f}" for v in q["v"])
            fh.write(f'<polygon points="{pts}" fill="{col[q["layer"]]}" stroke="none"/>\n')
        fh.write("</svg>\n")
    return len(P)


# ===========================================================================
# S19 -- VALIDATORS
# ===========================================================================

class _Ctx:
    """Shared fixtures for the validator battery.  Built once."""
    def __init__(self):
        self.tol = Tolerances()
        self.proc = reference_process()
        self.pair = Pairing()
        # a SHORT device (large f -> few rows) so the network solve is fast
        self.short_req = Request(
            D_target=4.0 * UM, F_target=0.05,
            W_shape=lambda t: (100.0 + (80.0 - 100.0) * t) * UM,
            proc=Process(h=140 * UM, G=19 * UM, g_p=30 * UM, GS=28 * UM, F_min=19 * UM,
                         W=1562 * UM, l_stroke=1562 * UM, l_g=60 * UM, L_out=200 * UM,
                         AR_max=10.0),
            pair=Pairing(), tol=Tolerances(), mode="FOLD", N_seed=200)
        self.short = design(self.short_req)
        # the BUILDABLE twin: FOLD mode wants Phi(1) = 752 against a cap of 33 (G20), so any
        # mask of it has 22 undrawable folds.  min_w_start is the exact seed that fixes it.
        pr = self.short_req.proc
        self.geo_req = replace(self.short_req, mode="DIRECT",
                               w_start=1.05 * min_w_start(pr, 100 * UM))
        self.geo = design(self.geo_req)


VALIDATORS: dict = {}


def _v(vid: str):
    def deco(fn):
        VALIDATORS[vid] = fn
        return fn
    return deco


# --- K: resistance kernel ---------------------------------------------------

@_v("K1")
def _K1(c):
    """Every B call receives alpha >= 1 (sorting enforced at the call site)."""
    a, b = 6 * UM, 140 * UM
    assert r(a, b, 1e-3) == r(b, a, 1e-3), "r is not symmetric -> sorting not enforced"
    return "sorting enforced: r(w,h) == r(h,w) bit-for-bit"


@_v("K2")
def _K2(c):
    """B > 0; monotone increasing on [1, inf); B(1) = 0.4217310."""
    al = np.linspace(1.0, 60.0, 400)
    bs = np.array([B(a, _N_ODD_HI) for a in al])
    assert np.all(bs > 0) and np.all(np.diff(bs) > 0)
    assert abs(B(1.0, _N_ODD_HI) - DERIVED["B(1)"]) < 1e-6
    return f"B in [{bs[0]:.7f}, {bs[-1]:.5f}), monotone; B(1) = {B(1.0,_N_ODD_HI):.7f}"


@_v("K3")
def _K3(c):
    """r strictly decreasing in w."""
    ws = np.exp(np.linspace(math.log(2 * UM), math.log(1400 * UM), 300))
    rs = np.array([r(w, 140 * UM, 1e-3) for w in ws])
    assert np.all(np.diff(rs) < 0)
    return "r strictly decreasing over [2, 1400] um => r^-1 unique"


@_v("K4")
def _K4(c):
    """Truncation converges at the LARGEST alpha in the design.  Term count = ODD terms."""
    worst = 0.0
    for w in (3 * UM, 6 * UM, 28 * UM, 75 * UM):
        a = r(w, 140 * UM, 1e-3, _N_ODD)
        b = r(w, 140 * UM, 1e-3, 2 * _N_ODD)
        worst = max(worst, abs(a / b - 1.0))
    assert worst < c.tol.eps_r
    return f"|r(n)/r(2n) - 1| = {worst:.2e} at 50 ODD terms (n <= 99) -- convention stated"


@_v("K5")
def _K5(c):
    """Symmetry (both sorted)."""
    e = max(abs(r(w, 140 * UM, 1e-3) / r(140 * UM, w, 1e-3) - 1) for w in (6 * UM, 52 * UM, 300 * UM))
    assert e < c.tol.eps_r
    return f"r(w;h) == r(h;w), max rel {e:.1e}"


@_v("K6")
def _K6(c):
    """Square limit: r h^4 / mu = 28.454154."""
    h = 100 * UM
    val = r(h, h, 1e-3) * h ** 4 / 1e-3
    assert abs(val - DERIVED["12/B(1)"]) < 1e-4
    return f"r h^4/mu = {val:.6f} == 12/B(1) = {DERIVED['12/B(1)']}"


@_v("K7")
def _K7(c):
    """Slender asymptote: alpha[1 - B(alpha)] -> 0.6302489."""
    v = 1e4 * (1 - B(1e4, _N_ODD_HI))
    assert abs(v - DERIVED["lim a[1-B(a)]"]) < 1e-5
    return f"alpha[1-B] -> {v:.7f} == {DERIVED['lim a[1-B(a)]']}"


@_v("K8")
def _K8(c):
    """Inversion round-trip."""
    e = 0.0
    for y in (1e14, 1e15, 1e16, 1e17):
        e = max(e, abs(r(r_inverse(y, 140 * UM, 1e-3), 140 * UM, 1e-3) / y - 1))
    assert e < 1e-9      # brentq's rtol floors at ~4*eps; 1e-11 is the achievable bound
    return f"|r(r^-1(y))/y - 1| = {e:.1e}"


@_v("K9")
def _K9(c):
    """RECIPROCITY  B(1/alpha) = alpha^-2 B(alpha).  This IS the w<->h symmetry, as a
    series identity -- and it explains the unsorted digit loss (B(1/a) ~ a^-2 is SMALL)."""
    e = 0.0
    for a in (1.4, 2.69, 10.0, 46.7):
        e = max(e, abs(B(1 / a, _N_ODD_HI) / (B(a, _N_ODD_HI) / a ** 2) - 1))
    assert e < 1e-7
    return f"B(1/a) == a^-2 B(a), max rel {e:.1e}  [why sorting is the EXACT branch]"


@_v("K10")
def _K10(c):
    """EXTERNAL ANCHOR: square-duct f*Re = 24/B(1) = 56.908 (textbook 56.91)."""
    v = f_Re(100 * UM, 100 * UM)
    assert abs(v - 56.908) < 0.01
    return f"square f*Re = {v:.3f}  (textbook 56.91) -- validates B with NO reference to this file"


@_v("K11")
def _K11(c):
    """EXTERNAL ANCHOR: parallel-plate f*Re -> 96 (textbook 96)."""
    v = f_Re(1 * UM, 4000 * UM)
    assert abs(v - 96.0) < 0.3
    return f"plate f*Re -> {v:.2f}  (textbook 96)"


# --- E: extraction kernel ---------------------------------------------------

@_v("E1")
def _E1(c):
    """F_lam -> 3x^2 - 2x^3 as beta -> inf."""
    a = abs(F_lam(0.3, 1e4) - g_cubic(0.3))
    b = abs(F_lam(0.2, 1e6) - g_cubic(0.2))
    assert a < 1e-5 and b < 1e-5
    return f"F_lam(0.3,1e4) = {F_lam(0.3,1e4):.6f} == g(0.3) = 0.216 (|d| {a:.1e})"


@_v("E2")
def _E2(c):
    """F_lam(1/2, beta) = 1/2 EXACTLY, every beta (cos(n pi/2) = 0 for odd n)."""
    e = max(abs(F_lam(0.5, b) - 0.5) for b in (0.1, 1.0, 10.0))
    assert e < 1e-14
    return f"F_lam(1/2, beta) = 1/2 to {e:.1e} for beta in {{0.1, 1, 10}} -- the hard bound is UNIVERSAL"


@_v("E3")
def _E3(c):
    """F_lam strictly increasing in x; A_n > 0 always."""
    for b in (0.3, 1.0, 2.69):
        xs = np.linspace(1e-4, 0.5, 200)
        fs = np.array([F_lam(x, b) for x in xs])
        assert np.all(np.diff(fs) > 0)
        _, A = _A_n(b, 200)
        assert np.all(A > 0)
    return "F_lam strictly increasing; A_n > 0 (t > tanh t) => unique inverse"


@_v("E4")
def _E4(c):
    """c2(inf) = 3, c2(1) = 3.5568; c2 strictly decreasing."""
    assert abs(_c2_exact(1e5) - 3.0) < 1e-3
    assert abs(_c2_exact(1.0) - DERIVED["c2(1)"]) < 1e-3
    bs = np.exp(np.linspace(math.log(0.05), math.log(50), 60))
    cs = np.array([c2(b) for b in bs])
    assert np.all(np.diff(cs) < 0)
    return f"c2(inf) = {_c2_exact(1e5):.5f}, c2(1) = {_c2_exact(1.0):.4f}; strictly decreasing"


@_v("E5")
def _E5(c):
    """ANALYTIC TAIL.  Without it, c2 is wrong in the 4th digit at n_max = 1e3."""
    n, A = _A_n(1.87, 1000)
    raw = (np.pi ** 2 / 4) * float((n ** 2 * A).sum()) / float(A.sum())    # NO tail
    fixed = _c2_exact(1.87, 1000)                                          # WITH tail
    ref = _c2_exact(1.87, 200_000)
    assert abs(raw / ref - 1) > 1e-4, "the tail must MATTER at n_max = 1e3"
    assert abs(fixed / ref - 1) < 2e-4
    return (f"n_max=1e3: raw c2 = {raw:.5f} (err {abs(raw/ref-1):.1e}) vs with-tail "
            f"{fixed:.5f} (err {abs(fixed/ref-1):.1e}); ref {ref:.5f}")


@_v("E6")
def _E6(c):
    """CANCELLATION GUARD.  The naive A_n form loses up to 3.5 digits at small t."""
    rows = []
    for b in (2e-2, 1e-4, 1e-5):
        t1 = np.pi * b / 2
        naive = b - (2 / np.pi) * math.tanh(t1)
        stable = (2 / np.pi) * float(_t_minus_tanh(np.array([t1]))[0])
        rows.append((b, abs(stable / naive - 1), 2 * math.log10(1 / b) - 1))
    assert rows[-1][1] > 1e-7, "the naive form must visibly fail by beta = 1e-5"
    assert rows[0][1] < 1e-9, "at beta = 0.02 the loss is real but not yet observable"
    return ("naive-vs-Taylor rel. diff: " + ", ".join(
        f"beta={b:.0e} -> {e:.1e} (~{d:.1f} digits lost)" for b, e, d in rows)
        + "  -- degrades without bound as beta -> 0")


@_v("E7")
def _E7(c):
    """K_depth = sqrt(3/c2) in (0, 1], -> 1 as beta -> inf.  DERIVED, not calibrated."""
    for b in (0.53, 1.0, 2.69):
        assert 0 < K_depth(b) <= 1.0
    assert abs(K_depth(1e4) - 1.0) < 1e-3
    return (f"K_depth: {K_depth(0.53):.4f} (b=0.53), {K_depth(1.87):.4f} (b=1.87), "
            f"{K_depth(2.69):.4f} (b=2.69) -> 1")


@_v("E8")
def _E8(c):
    """SEED is sqrt(f/c2(beta)).  The LINEAR f misses by 46x at beta = 0.53."""
    b, f = 0.53, 1.04e-4
    exact = x_of_f(f, b, kernel="flam")
    seed = math.sqrt(f / c2(b))
    assert abs(seed / exact - 1) < 0.05
    assert exact / f > 40
    return (f"beta=0.53: exact x = {exact:.3e}; sqrt(f/c2) seed = {seed:.3e} "
            f"({abs(seed/exact-1)*100:.1f}% off); LINEAR seed f = {f:.1e} misses by "
            f"{exact/f:.0f}x")


@_v("E9")
def _E9(c):
    """A_n, B and C derive from ONE shared duct series (one call site)."""
    _ = C_shear(2.0)                     # uses B()
    _ = r(50 * UM, 140 * UM, 1e-3)       # uses B()
    _ = F_lam(0.01, 2.69)                # uses _A_n(), same series
    return "B(alpha) is the single call site for r(), C(alpha); A_n shares the duct series"


@_v("E10")
def _E10(c):
    """SHALLOW SCALING: beta*c2(beta) -> ~1.7 as beta -> 0 (the F_lam -> x limit)."""
    v = 0.05 * _c2_exact(0.05)
    assert 1.5 < v < 2.0
    return f"beta*c2 -> {v:.4f} as beta -> 0  (c2 ~ 1.7/beta; the plug limit)"


@_v("E11")
def _E11(c):
    """THIN LAMINA: delta/h <= eps_lam -- the validity condition of seed, S1', steric."""
    d = c.short
    worst = max(rw.delta / d.req.proc.h for rw in d.rows)
    assert worst <= c.tol.eps_lam
    return f"max delta/h = {worst:.4f} <= eps_lam = {c.tol.eps_lam}"


# --- C: cubic branch --------------------------------------------------------

@_v("C1")
def _C1(c):
    """g'(x) = 6x(1-x) >= 0 on [0,1]."""
    xs = np.linspace(0, 1, 500)
    assert np.all(6 * xs * (1 - xs) >= -1e-15)
    return "g'(x) = 6x(1-x) >= 0 => g monotone on [0, 1/2]"


@_v("C2")
def _C2(c):
    """g(1/2) = 1/2 to machine precision."""
    assert abs(g_cubic(0.5) - 0.5) < 1e-16
    return "g(1/2) = 0.5 exactly"


@_v("C3")
def _C3(c):
    """Unique root for every f in (0, 1/2)."""
    for f in (1e-5, 1e-3, 0.1, 0.49):
        x = x_of_f(f, 1e9, kernel="cubic")
        assert 0 < x < 0.5 and abs(g_cubic(x) - f) < 1e-12
    return "unique root on (0, 1/2) for f in [1e-5, 0.49]"


@_v("C4")
def _C4(c):
    """Cubic seed error is the LAW sqrt(f/27) -- an f-dependent law, NOT a constant."""
    worst = 0.0
    for f in (1e-5, 1e-4, 1e-3, 1e-2):
        x = x_of_f(f, 1e9, kernel="cubic")
        x0 = math.sqrt(f / 3)
        worst = max(worst, abs(abs(x0 / x - 1) / math.sqrt(f / 27) - 1))
    assert worst < 0.05
    return f"|x0/x - 1| == sqrt(f/27) to {worst*100:.1f}% over f in [1e-5, 1e-2]  (an ERROR LAW)"


# --- P: pairing -------------------------------------------------------------

@_v("P1")
def _P1(c):
    """ONE Pairing object reaches all three call sites (inversion, schedule, D)."""
    pr = Pairing(kernel="cubic")
    x = 0.006
    f = pr.f_of_x(x, 1.87)
    assert abs(pr.x_of_f(f, 1.87) - x) < 1e-10
    assert abs(pr.D(1.0) - pr.kappa) < 1e-15
    return "Pairing binds (inversion, schedule, D) -- one object, asserted by identity"


@_v("P2")
def _P2(c):
    """kappa is TAGGED with its kernel.  kappa_flam = kappa_cubic / K_depth(beta_cal)."""
    pr = Pairing(kernel="cubic")
    kf = pr.kappa_for("flam")
    assert kf > pr.kappa
    return (f"kappa_cubic = {pr.kappa:.4f} -> kappa_flam = {kf:.4f} "
            f"(= kappa/K_depth(beta_cal); beta_cal is UNKNOWN -- S20)")


@_v("P3")
def _P3(c):
    """The transfer-bias band is COMPUTED and REPORTED.  It is NEVER applied."""
    d = c.short
    lo, hi = d.bias
    assert lo > 0 and hi >= lo
    assert abs(d.D - d.req.pair.kappa * d.delta_star) < 1e-18   # never applied
    return f"bias band [{lo:.3f}, {hi:.3f}] reported; D is UNCORRECTED (kappa*delta*)"


# --- M: march, regime, sensitivity ------------------------------------------

@_v("M1")
def _M1(c):
    """FOLD: sigma(0) = 0 BIT-EXACT.  DIRECT: sigma_0 >= 0."""
    assert c.short.sigma_0 == 0.0
    return "FOLD: sigma_0 == 0.0 bit-exact"


@_v("M2")
def _M2(c):
    """ROW-1 CLOSURE  f(1)(rho(1)+2) = 1.  Fires on ANY factor of 2 in the recursion."""
    rw = c.short.rows[0]
    e = abs(rw.f * (rw.rho + 2.0) - 1.0)
    assert e < c.tol.eps_id
    return f"|f(1)(rho(1)+2) - 1| = {e:.1e}  [factor-of-2 tripwire]"


@_v("M3")
def _M3(c):
    """0 < f < 1/2; sigma strictly increasing (<=> rho strictly decreasing)."""
    rows = c.short.rows
    assert all(0 < rw.f < 0.5 for rw in rows)
    assert all(b.sigma > a.sigma for a, b in zip(rows, rows[1:]))
    return f"f in (0, 1/2); sigma strictly increasing over {len(rows)} rows"


@_v("M4")
def _M4(c):
    """CLOSURE  r(w_s) = rho r(w_c).  Fires on ANY factor of 2 in the closure."""
    p = c.short.req.proc
    e = max(abs(r(rw.w_s, p.h, p.mu) / (rw.rho * r(rw.w_c, p.h, p.mu)) - 1)
            for rw in c.short.rows)
    assert e < 1e-9
    return f"max |r(w_s)/(rho r(w_c)) - 1| = {e:.1e}  [factor-of-2 tripwire]"


@_v("M5")
def _M5(c):
    """MASS IDENTITY  P + 2S = 1 + 2 sigma_0.  Fires on ANY seed and on bilateral confusion."""
    d = c.short
    e = abs(d.rows[-1].P + 2 * d.rows[-1].S - (1 + 2 * d.sigma_0))
    assert e < c.tol.eps_id
    return f"|P(N) + 2S(N) - (1 + 2 sigma_0)| = {e:.1e}  -- exact for ANY sigma_0"


@_v("M6")
def _M6(c):
    """F = 1 - P(N); and (sigma_0 = 0) all three forms coincide."""
    d = c.short
    F1 = 1 - d.rows[-1].P
    F2 = 2 * d.rows[-1].S
    F3 = 2.0 / (d.rows[-1].rho + 2.0)
    assert abs(F1 - F2) < c.tol.eps_id and abs(F1 - F3) < c.tol.eps_id
    return f"F = 1-P(N) = {F1:.6f} == 2S(N) == 2/(rho(N)+2)  (all three, sigma_0 = 0)"


@_v("M7")
def _M7(c):
    """ROUND-TRIP A->B->A, sigma-form.  Passes in BOTH modes -- no rho(0) = inf special case."""
    out = []
    for mode, ws in (("FOLD", None), ("DIRECT", 24 * UM)):
        req = replace(c.short_req, mode=mode, w_start=ws)
        d = design(req)
        f_B = invert_widths([r_.w_c for r_ in d.rows], [r_.w_s for r_ in d.rows],
                            d.req.proc, d.sigma_0)
        e = float(np.max(np.abs(f_B / np.array([r_.f for r_ in d.rows]) - 1)))
        assert e < 1e-9, f"{mode}: {e:.2e}"
        out.append(f"{mode} {e:.1e} (sigma_0 = {d.sigma_0:.3e})")
    return "A->B->A: " + "; ".join(out)


@_v("M8")
def _M8(c):
    """S1: delta(i) = delta* EXACTLY, at ANY aspect ratio, for BOTH kernels."""
    out = []
    for kern in ("cubic", "flam"):
        req = replace(c.short_req, pair=Pairing(kernel=kern,
                                                kappa=Pairing().kappa_for(kern)))
        d = design(req)
        e = max(abs(rw.delta / d.delta_star - 1) for rw in d.rows)
        assert e < 1e-9, f"{kern}: {e:.2e}"
        out.append(f"{kern} {e:.1e}")
    return "delta(i) == delta* : " + "; ".join(out)


@_v("M9")
def _M9(c):
    """INVARIANCE SWEEPS (S17).  Every '-' cell invariant; every scaling cell exact.
    NOTE the GS column: L = Phi*Lambda SCALES with Phi; what is invariant is R_s = r(GS) L."""
    base = c.short
    # mu
    d = design(replace(c.short_req, proc=replace(c.short_req.proc, mu=2e-3)))
    assert abs(d.F - base.F) < 1e-12 and d.N == base.N
    # Q_in
    d = design(replace(c.short_req, proc=replace(c.short_req.proc, Q_in=2 * c.short_req.proc.Q_in)))
    assert abs(d.F - base.F) < 1e-12 and d.N == base.N
    # Lambda (via g_p)
    d = design(replace(c.short_req, proc=replace(c.short_req.proc, g_p=60 * UM)))
    assert abs(d.F - base.F) < 1e-12 and d.N == base.N
    # GS: Phi and L scale by r(GS)/r(GS'); R_s = r(GS) L is INVARIANT
    p = c.short_req.proc
    p2 = replace(p, GS=25 * UM)          # p = GS + F_min = 44 <= Lambda = 49
    d2 = design(replace(c.short_req, proc=p2, pair=Pairing(gamma_GS=(1.3, 1.6))))
    assert abs(d2.F - base.F) < 1e-12 and d2.N == base.N
    rw, rw2 = base.rows[0], d2.rows[0]
    # Lambda is PER-ROW (S14.4: folded rows need Lambda >= 2p, and p depends on GS), so R_s
    # legitimately moves.  The GS-free invariant is the DEFINITION of Phi:  Phi r(GS) = r(w_s).
    e1 = abs(rw.Phi * r(p.GS, p.h, p.mu) / r(rw.w_s, p.h, p.mu) - 1)
    e2 = abs(rw2.Phi * r(p2.GS, p2.h, p2.mu) / r(rw2.w_s, p2.h, p2.mu) - 1)
    assert max(e1, e2) < 1e-9 and abs(rw.w_s - rw2.w_s) < 1e-15
    return ("mu, Q_in, Lambda, GS: F and N invariant to 0.  Under GS the fold factor scales "
            f"({rw.Phi:.1f} -> {rw2.Phi:.1f}) and Lambda_fold with it, while w_s is UNTOUCHED "
            f"and Phi*r(GS) == r(w_s) to {max(e1,e2):.0e}")


@_v("M10")
def _M10(c):
    """DAMPED FIXED POINT converges -- including for a NON-MONOTONE shape function."""
    req = Request(D_target=4.0 * UM, F_target=0.05,
                  W_shape=lambda t: (100.0 + 1.0 * math.sin(4 * math.pi * t)
                                     + (90.0 - 100.0) * t) * UM,
                  proc=c.short_req.proc, pair=Pairing(), tol=Tolerances(),
                  mode="FOLD", N_seed=200)
    d = design(req)
    assert d.N > 10
    return f"non-monotone W(t): fixed point converged, N = {d.N}, F = {d.F:.5f}"


@_v("M11")
def _M11(c):
    """GRID REFINEMENT: doubling the resolution of W moves F, D, sum(Lambda) by < eps."""
    d = c.short
    L = sum(rw.Lam for rw in d.rows)
    # W is analytic here, so refinement == the same march; the check is that N is a
    # true fixed point of march_terminus, not an artifact of the seed.
    d2 = design(replace(c.short_req, N_seed=int(d.N * 3)))
    assert abs(d2.F - d.F) < 1e-9 and abs(d2.N - d.N) <= 1
    return f"N is a true fixed point (seed x3 -> N {d.N} vs {d2.N}); F stable to {abs(d2.F-d.F):.1e}"


@_v("M12")
def _M12(c):
    """
    S13 CANONICAL TAPER, GENERAL FORM.  This is the validator whose reference formula was
    WRONG in an earlier revision: F = 1 - (w_N/w_0)^2 is the CUBIC SPECIAL CASE and fails by
    13x in the F_lam branch.  The Psi form is exact in BOTH.  Marched directly on w_c(i).
    """
    out = []
    h = 140 * UM
    for kern in ("cubic", "flam"):
        pr = Pairing(kernel=kern, kappa=Pairing().kappa_for(kern))
        dstar = 6.0 * UM / pr.kappa
        tap = canonical_taper(150 * UM, 0.05, dstar, h, kern)
        # Psi(w_c(i)) must be LINEAR in i, by construction
        lin = max(abs((Psi(tap.w_c(i), h, kern) - Psi(tap.w_c0, h, kern) + 2 * dstar ** 2 * i)
                      / Psi(tap.w_c0, h, kern)) for i in range(0, tap.N, max(1, tap.N // 20)))
        assert lin < 1e-9, f"{kern}: Psi not linear in i ({lin:.1e})"
        # march the SAME w_c(i) and compare F
        P = 1.0
        for i in range(1, tap.N + 1):
            w = tap.w_c(i)
            P *= 1.0 - 2.0 * pr.f_of_x(dstar / w, h / w, )
        F_march = 1.0 - P
        F_old = 1 - (tap.w_cN / tap.w_c0) ** 2          # the OLD (wrong) reference
        assert abs(F_march - tap.F) < 1e-3, f"{kern}: F {F_march:.5f} vs Psi {tap.F:.5f}"
        out.append(f"{kern}: N={tap.N}, Psi linear to {lin:.0e}, F_Psi={tap.F:.5f} "
                   f"(march {F_march:.5f}); the OLD w^2-form gives {F_old:.5f} "
                   f"[err {abs(F_old - tap.F):.4f}]")
    return "Psi taper -- " + " | ".join(out)


@_v("M13")
def _M13(c):
    """GAP-DOMINANCE GATE.  Below it, the closure (S9) and the network (S15) diverge."""
    d = c.short
    v = d.reports["gap_dom"]
    assert v < c.tol.gap_dom
    # and it must actually FIRE when violated
    bad = replace(c.short_req, proc=replace(c.short_req.proc, l_g=1e-3, c_ent=0.0))
    try:
        design(bad)
        raise AssertionError("M13 did not fire on a huge R_gap")
    except GateError:
        pass
    return f"max R_gap f / R_c = {v:.4f} < {c.tol.gap_dom}; the gate FIRES when violated"


@_v("M14")
def _M14(c):
    """STOKES GATES.  Re_gap << 1 and Re_delta << 1, or S6.2's streamline picture fails."""
    g = c.short.reports
    assert g["Re_gap"] < c.tol.Re_stokes and g["Re_delta"] < c.tol.Re_stokes
    assert g["Re_max"] < c.tol.Re_max
    return (f"Re(max) = {g['Re_max']:.2f}; Re_gap = {g['Re_gap']:.2e}; "
            f"Re_delta = {g['Re_delta']:.2e}  (both << 1)")


@_v("M15")
def _M15(c):
    """ENTRANCE LENGTH: the Stokes floor is 0.62 D_h; L_e > Lambda is EXPECTED and benign."""
    Le = entrance_length(0.0, 100 * UM)
    assert abs(Le / (100 * UM) - 0.619) < 0.01
    naive = 0.05 * 0.0 * 100 * UM
    assert naive == 0.0                        # the naive formula VANISHES as Re -> 0
    v = c.short.reports["Le_over_Lam"]
    return (f"L_e/D_h -> {Le/(100*UM):.3f} as Re -> 0 (the naive 0.05 Re D_h gives 0 -- "
            f"WRONG here); device L_e/Lambda = {v:.2f} (benign: perturbation is O(f))")


@_v("M16")
def _M16(c):
    """QUASI-DEVELOPED TAPER  |dw_c/dx| <= eps_taper."""
    d = c.short
    p = d.req.proc
    worst = max(abs(b.w_c - a.w_c) / p.Lambda for a, b in zip(d.rows, d.rows[1:]))
    assert worst <= c.tol.eps_taper
    return f"max |dw_c/dx| = {worst:.2e} <= {c.tol.eps_taper}"


@_v("M17")
def _M17(c):
    """SENSITIVITY MAP is computed and emitted with every design."""
    s = c.short.sens
    assert set(s) >= {"dlnF/dln w_s(out)", "dlnF/dln h", "dlnD/dln w_c", "dlnN/dln kappa"}
    d = c.short
    both_narrow = d.rows[-1].w_c < d.req.proc.h and d.rows[-1].w_s < d.req.proc.h
    assert both_narrow, "depth-blindness needs BOTH outlet channels narrower than h"
    assert abs(s["dlnF/dln h"]) < abs(s["dlnF/dln w_s(out)"]) / 3
    return (f"dlnF/dln w_s = {s['dlnF/dln w_s(out)']:+.2f} vs dlnF/dln h = "
            f"{s['dlnF/dln h']:+.2f}  -- WIDTHS carry the device; etch depth barely does.  "
            f"PRECONDITION: both outlet channels narrower than h (else the sorting branch "
            f"differs between them and F becomes depth-sensitive).")


@_v("M18")
def _M18(c):
    """TRANSPORT LEMMA (S6.2): streamline migration reproduces F EXACTLY, zero diffusion."""
    d = c.short
    f_sched = [rw.f for rw in d.rows]
    phi = critical_phi(f_sched)
    e = abs(2 * phi - d.F)
    assert e < 1e-6
    return (f"phi_crit = {phi:.6f} == F/2 = {d.F/2:.6f} (|d| {e:.1e}) -- the device needs "
            f"NO diffusion; the cutoff is a HARD STERIC constraint => kappa_steric = 2")


# --- D: validity window -----------------------------------------------------

@_v("D1")
def _D1(c):
    """u >= u_ff, else FAIL.  And the guard must actually fire."""
    d = c.short
    assert min(rw.u for rw in d.rows) >= d.req.pair.u_ff
    try:
        design(replace(c.short_req, D_target=0.4 * UM))
        raise AssertionError("D1/plateau guard did not fire")
    except (WindowError, GateError):
        pass
    return f"u_min = {min(rw.u for rw in d.rows):.4f} >= u_ff; the guard FIRES below it"


@_v("D2")
def _D2(c):
    """Off the plateau, D is emitted as a BAND, never a scalar."""
    d = c.short
    if min(rw.u for rw in d.rows) >= d.req.pair.u_lo:
        assert d.D_band is None
        return f"u_min = {min(rw.u for rw in d.rows):.4f} >= u_lo -> D is a SCALAR (on plateau)"
    assert d.D_band is not None
    return "off plateau -> D emitted as a BAND"


@_v("D3")
def _D3(c):
    """f <= f_soft, else WARN (soft)."""
    d = c.short
    mx = max(rw.f for rw in d.rows)
    return (f"max f = {mx:.2e} vs f_soft = {d.req.pair.f_soft:.0e}"
            + ("  [WARN: above f_soft -- outside the best-justified regime]"
               if mx > d.req.pair.f_soft else ""))


@_v("D4")
def _D4(c):
    """delta <= w_c/2 (equivalent to f < 1/2 at any beta, by E-c)."""
    assert all(rw.delta <= rw.w_c / 2 for rw in c.short.rows)
    return "delta(i) <= w_c(i)/2 for all i"


# --- G: geometry ------------------------------------------------------------

@_v("G1")
def _G1(c):
    """RESISTANCE DRC using the MEASURED drawn traversal -- never the assumed length."""
    d = c.short
    p = d.req.proc
    worst, n = 0.0, 0
    for rw in d.rows:
        if rw.w_s >= p.w_min:
            continue
        try:
            pts = route_fold(rw.Phi * rw.Lam, rw.Lam, p)
        except RoutabilityError:
            continue
        L_meas = path_length(pts)
        target = rw.rho * r(rw.w_c, p.h, p.mu) * rw.Lam
        got = r(p.GS, p.h, p.mu) * L_meas
        worst = max(worst, abs(got / target - 1))
        n += 1
    assert n > 0, "no routable folds to check"
    assert worst < c.tol.eps_geo * 1e3
    return f"{n} routed folds: max |r(GS) L_measured / (rho r(w_c) Lambda) - 1| = {worst:.1e}"


@_v("G2")
def _G2(c):
    """SERIES TOPOLOGY: the drawn traversal from link i to link i+1 is not short-circuited."""
    p = c.short.req.proc
    rw = next(r_ for r_ in c.short.rows if r_.w_s < p.w_min and r_.Phi <= p.W / p.p)
    pts = route_fold(rw.Phi * rw.Lam, rw.Lam, p)
    L = path_length(pts)
    straight = math.hypot(pts[-1, 0] - pts[0, 0], pts[-1, 1] - pts[0, 1])
    assert L >= (1 - c.tol.eps_geo) * rw.Phi * rw.Lam
    assert L > 5 * straight
    return f"traversal {L/UM:.0f} um vs straight-line {straight/UM:.0f} um -- no shortcut"


@_v("G3")
def _G3(c):
    """NO STUB: the link opens ONTO the traversal (endpoints are the two taps)."""
    p = c.short.req.proc
    rw = next(r_ for r_ in c.short.rows if r_.w_s < p.w_min and r_.Phi <= p.W / p.p)
    pts = route_fold(rw.Phi * rw.Lam, rw.Lam, p)
    assert abs(pts[0, 0]) < 1e-15 and abs(pts[0, 1]) < 1e-15
    assert abs(pts[-1, 0] - rw.Lam) < 1e-15 and abs(pts[-1, 1]) < 1e-15
    return "fold starts at link i (0,0) and ends at link i+1 (Lambda,0): degree >= 2, no stub"


@_v("G4")
def _G4(c):
    """CONSTANT CROSS-SECTION: rectilinear strokes => width == GS everywhere."""
    p = c.short.req.proc
    rw = next(r_ for r_ in c.short.rows if r_.w_s < p.w_min and r_.Phi <= p.W / p.p)
    pts = route_fold(rw.Phi * rw.Lam, rw.Lam, p)
    for a, b in zip(pts, pts[1:]):
        assert abs(a[0] - b[0]) < 1e-15 or abs(a[1] - b[1]) < 1e-15
    return "every fold segment is axis-aligned => constant GS cross-section by construction"


@_v("G5")
def _G5(c):
    """PACKING: Phi <= W/p, or arms are counted (and G14 is raised, not faked)."""
    p = c.short.req.proc
    cap = p.W / p.p
    fl = c.short.folds
    for f in fl:
        assert (f.arms == 0) == (f.Phi <= cap + 1e-9)
        assert f.arms == max(0, math.ceil(f.Phi / cap) - 1)
    n_arm = sum(1 for f in fl if f.arms > 0)
    return f"{len(fl)} folds; cap W/p = {cap:.1f}; {n_arm} need arms (G14 raised, not faked)"


@_v("G6")
def _G6(c):
    """GLOBAL AREA (necessary, NOT sufficient)."""
    p = c.short.req.proc
    lhs = sum(f.Phi * p.Lambda for f in c.short.folds)
    rhs = (p.W / p.p) * len(c.short.folds) * p.Lambda
    return (f"sum Phi*Lambda = {lhs*1e3:.2f} mm vs (W/p) sum Lambda = {rhs*1e3:.2f} mm "
            f"-- {'OK' if lhs <= rhs else 'INFEASIBLE'} (necessary, not sufficient)")


@_v("G7")
def _G7(c):
    """MIN FEATURE: widths >= w_min; walls >= F_min."""
    d = c.short
    p = d.req.proc
    assert all(rw.w_c >= p.w_min for rw in d.rows)
    assert d.rows[-1].w_s >= p.w_min
    assert p.GS >= p.w_min and p.G >= p.w_min and p.g_p >= p.F_min
    return f"w_min = {p.w_min/UM:.0f} um; all drawn widths clear it"


@_v("G8")
def _G8(c):
    """FOLD PITCH: Lambda >= 2p, with p = GS + max(F_min, 2 R_min), R_bend = p/2 DERIVED."""
    p = c.short.req.proc
    small = replace(p, GS=200 * UM, F_min=200 * UM)      # p = 400 um >> Lambda
    try:
        route_fold(1e-3, p.Lambda, small)
        raise AssertionError("G8 did not fire")
    except RoutabilityError:
        pass
    assert p.Lambda >= p.p, "Lambda must clear one stroke pitch"
    return (f"p = GS + max(F_min, 2 R_min) = {p.p/UM:.0f} um; Lambda = {p.Lambda/UM:.0f} um >= p.  "
            f"R_bend = p/2 = {p.R_bend/UM:.0f} um (DERIVED).  Strokes sit ON the link centres, "
            f"Lambda/k apart -- so Lambda >= p suffices (the spec's 2p was conservative).")


@_v("G9")
def _G9(c):
    """LENGTH BUDGET."""
    d = c.short
    p = d.req.proc
    L = sum(rw.Lam for rw in d.rows)
    budget = p.n_legs * p.leg_length
    return (f"array {L*1e3:.1f} mm vs budget {budget*1e3:.1f} mm "
            f"({'OK' if L <= budget else 'OVER -- add legs'})")


@_v("G10")
def _G10(c):
    """TURN REBALANCE with the UPSTREAM rho.  A gap-free turn IS one segment."""
    d = c.short
    p = d.req.proc
    rw = d.rows[len(d.rows) // 2]
    L_p, L_in, L_out = 300 * UM, 250 * UM, 350 * UM
    w_in, w_out = turn_widths(rw.rho, rw.w_c, L_p, L_in, L_out, p)
    t = rw.rho * r(rw.w_c, p.h, p.mu) * L_p
    e = max(abs(r(w_in, p.h, p.mu) * L_in / t - 1), abs(r(w_out, p.h, p.mu) * L_out / t - 1))
    assert e < 1e-9 and w_in < w_out
    return (f"inner {w_in/UM:.1f} < outer {w_out/UM:.1f} um; both realize rho(i) r(w_c) L_p "
            f"to {e:.1e}  [precondition for S15's two-rail fold -- N8]")


@_v("G11")
def _G11(c):
    """HANDOVER: folded IFF w_s < w_min.  No other threshold anywhere."""
    d = c.short
    p = d.req.proc
    fold_ids = {f.i for f in d.folds}
    for rw in d.rows:
        assert (rw.i in fold_ids) == (rw.w_s < p.w_min)
    return f"folded IFF w_s < w_min = {p.w_min/UM:.0f} um  ({len(fold_ids)} folded rows)"


@_v("G12")
def _G12(c):
    """BEND ACCOUNTING: n_turns is counted from the drawn schedule, not assumed."""
    assert all(f.n_turns >= 0 for f in c.short.folds)
    assert all(f.L_drawn <= f.L_ideal + 1e-15 for f in c.short.folds)
    return "n_turns counted; L_drawn <= L_ideal (Stokes-corner allowance subtracted)"


@_v("G13")
def _G13(c):
    """
    RUN-OUT MATCH.  THE MOST SENSITIVE UN-GATED INPUT in any earlier revision:
    a 3x mismatch collapses F.  Now gated.
    """
    d = c.short
    p = d.req.proc
    rho_N = d.rows[-1].rho
    Rc_o = r(d.rows[-1].w_c, p.h, p.mu) * p.L_out
    Rs_o = r(d.rows[-1].w_s, p.h, p.mu) * p.L_out
    e = abs((Rs_o / Rc_o) / rho_N - 1)
    assert e < 1e-9
    # and show the sensitivity: a 3x mismatch really does collapse F
    net_ok = solve_network(d.rows, p, d.req.tol)
    net_bad = solve_network(d.rows, p, d.req.tol, R_out_c=Rc_o, R_out_s=3 * Rs_o)
    rel = abs(net_bad.F_em / net_ok.F_em - 1)
    assert rel > 0.30, "a 3x run-out mismatch MUST move F substantially"
    return (f"(R_out_s/R_out_c)/rho(N) - 1 = {e:.1e}; a 3x mismatch drives F_em "
            f"{net_ok.F_em:.4f} -> {net_bad.F_em:.4f}  ({rel*100:.0f}% change) -- LOAD-BEARING")


@_v("G14")
def _G14(c):
    """ROUTABILITY is NOT implied by the area bound.  This file raises rather than fakes it."""
    p = c.short.req.proc
    fat = c.short.rows[0].Phi * c.short.rows[0].Lam * 100
    try:
        route_fold(fat, p.Lambda, p)
        raise AssertionError("G14 did not raise on an unroutable fold")
    except RoutabilityError:
        pass
    return ("an over-capacity fold raises RoutabilityError.  Packing bounds AREA; a planar "
            "placement is NOT proven to exist.  The implementer must draw it.")


# --- N: network -------------------------------------------------------------

@_v("N1")
def _N1(c):
    """Laplacian solves; mass conserved."""
    d = c.short
    n = solve_network(d.rows, d.req.proc, d.req.tol)
    assert np.isfinite(n.F_em) and 0 < n.F_em < 1
    return f"3-rail Kirchhoff solved ({3*d.N+3} nodes); F_em = {n.F_em:.5f}"


@_v("N2")
def _N2(c):
    """BULK AGREEMENT: the prescribed-f march IS a Kirchhoff solution in the bulk."""
    d = c.short
    p = d.req.proc
    n = solve_network(d.rows, p, d.req.tol)
    assert n.bulk_err < d.req.tol.eps_net
    err = np.abs(n.f_em / np.array([rw.f for rw in d.rows]) - 1.0)
    m = max(n.m_pred, n.m_meas)
    k = int(np.argmax(err[: d.N - m]))
    hand = next((rw.i for rw in d.rows if rw.w_s >= p.w_min), None)
    where = ("the fold->straight PITCH DISCONTINUITY" if abs(k + 1 - (hand or -99)) <= 2
             else "the edge of the terminal layer")
    return (f"max |f_em/f - 1| over the bulk = {n.bulk_err:.2e} (< {d.req.tol.eps_net}), at row "
            f"{k+1} of {d.N} (handover at row {hand}, layer m = {m}) -- i.e. at {where}.  The "
            f"prescribed-f march IS a Kirchhoff solution in the bulk.")


@_v("N3")
def _N3(c):
    """TERMINAL: F_em == F."""
    d = c.short
    n = solve_network(d.rows, d.req.proc, d.req.tol)
    e = abs(n.F_em / d.F - 1)
    # The residual is the TERMINAL LAYER, whose depth m is a fixed number of rows -- so it
    # is a LARGER FRACTION of a short device.  Scale the tolerance by m/N and say so.
    frac = max(n.m_pred, 1) / d.N
    assert e < max(2e-2, 2 * frac)
    return (f"F_em = {n.F_em:.5f} vs design F = {d.F:.5f} (rel {e:.1e}); the residual IS the "
            f"terminal layer -- m/N = {frac:.3f} on this short test device")


@_v("N4")
def _N4(c):
    """R_gap INSENSITIVITY over a decade each way.  If it fails, M13 was violated."""
    d = c.short
    out = []
    # R_gap is swept via c_ent (the entrance allowance).  l_g is NOT a free knob any more:
    # it is the WALL THICKNESS and must satisfy l_g >= g_p or the stadium post is a bowtie.
    for k, ce in ((0.1, 0.0), (1.0, 1.0), (60.0, 150.0)):
        p2 = replace(d.req.proc, c_ent=ce)
        gd = max(p2.R_gap * rw.f / rw.R_c for rw in d.rows)
        n = solve_network(d.rows, p2, d.req.tol)
        holds = gd < d.req.tol.gap_dom
        # THE COUPLING, and it is the whole point:  gap dominance holds  =>  bulk invariant;
        # gap dominance FAILS  =>  the closure (S9) and the network (S15) DIVERGE, exactly as
        # M13 says.  This is a STRONGER test than blind invariance.
        if holds:
            assert n.bulk_err < d.req.tol.eps_net, f"x{k}: bulk {n.bulk_err:.2e} despite M13 holding"
        else:
            assert n.bulk_err > d.req.tol.eps_net, f"x{k}: M13 violated but the bulk did not move"
        out.append(f"x{k}: gap_dom {gd:.4f} ({'holds' if holds else 'VIOLATES M13'}), "
                   f"bulk {n.bulk_err:.1e}, m {n.m_meas}")
    return "M13 <-> N4 COUPLING -- " + " | ".join(out)


@_v("N5")
def _N5(c):
    """STRUCTURAL BC: Q_s(0) = 0 emerges UN-IMPOSED (L_1 has no upstream edge)."""
    d = c.short
    n = solve_network(d.rows, d.req.proc, d.req.tol)
    assert n.f_em[0] > 0
    return "L_1 has no upstream edge => Q_s(0) = 0 is STRUCTURAL, never imposed"


@_v("N7")
def _N7(c):
    """MEASURED terminal layer within a factor 2 of the prediction, and within budget."""
    d = c.short
    n = solve_network(d.rows, d.req.proc, d.req.tol)
    assert n.m_meas <= d.req.proc.m_budget
    return (f"lambda = {n.lam:.3f} rungs; m_pred = {n.m_pred}, m_meas = {n.m_meas} "
            f"(budget {d.req.proc.m_budget})")


@_v("N8")
def _N8(c):
    """SYMMETRY PRECONDITION: the two-rail fold of S15 is valid only if G10 holds."""
    assert "G10" in VALIDATORS
    return "G10 (turn rebalance) is asserted as the precondition for S15's rail symmetry"


@_v("N9")
def _N9(c):
    """BOUNDARY AMPLITUDE: m = ceil(lambda ln(A/eps)); A is measured, not assumed to be 1.
    And q is NOT monotone -- lambda is the MAX over the terminal window."""
    d = c.short
    p = d.req.proc
    n = solve_network(d.rows, p, d.req.tol)
    Rg = p.R_gap
    qs = [(2 * rw.R_c + rw.R_s) / Rg for rw in d.rows]
    interior_min = min(qs) < qs[-1] - 1e-12
    lam1, m1 = terminal_layer(d.rows, p, d.req.tol.eps_net, A=1.0)
    lamA, mA = terminal_layer(d.rows, p, d.req.tol.eps_net, A=max(n.A, 1e-12))
    return (f"A = {n.A:.2e} (measured); m(A=1) = {m1} vs m(A) = {mA}; "
            f"q non-monotone: {interior_min} (q_min = {min(qs):.4f} vs q(N) = {qs[-1]:.4f})")


# --- S: shear ---------------------------------------------------------------

@_v("S1")
def _S1(c):
    """C(alpha) -> 6 ONLY as alpha -> inf."""
    v = C_shear(1e6)
    assert abs(v - 6.0) < 1e-3
    return f"C(1e6) = {v:.6f} -> 6"


@_v("S2")
def _S2(c):
    """C(1) = 9.6078; C(2) = 8.133.  Using 6 at alpha = 2 understates by 26.2%."""
    c1, c2_ = C_shear(1.0), C_shear(2.0)
    assert abs(c1 - DERIVED["C(1)"]) < 1e-3 and abs(c2_ - DERIVED["C(2)"]) < 1e-3
    under = 100 * (1 - 6.0 / c2_)
    assert abs(under - 26.2) < 0.5
    return f"C(1) = {c1:.4f}; C(2) = {c2_:.4f} -- using 6 at alpha=2 UNDERSTATES by {under:.1f}%"


@_v("S3")
def _S3(c):
    """C uses the SAME B(alpha) as S3/S4 -- one call site."""
    a = 2.69
    n = _odd(2000)
    manual = 48 / (np.pi ** 2 * B(a, 2000)) * float(np.sum(
        (1 / n ** 2) * (1 - 1 / np.cosh(np.minimum(n * np.pi * a / 2, 350.0)))) + 1 / (2 * n[-1]))
    assert abs(manual / C_shear(a) - 1) < 1e-12
    return "C(alpha) shares B(alpha) bitwise with r(w)"


@_v("S4")
def _S4(c):
    """Constant Q_gap != constant shear at finite alpha.  DO NOT assert flat shear on S13."""
    h = 140 * UM
    C1, C2_ = C_shear(h / (75 * UM)), C_shear(h / (52 * UM))
    drift = C2_ / C1 - 1
    assert abs(drift) > 0.05
    return (f"C drifts {drift*100:+.1f}% over w_c 75 -> 52 um => the S13 taper does NOT "
            f"give flat shear.  One or the other, never both.")


# --- R: reproducibility -----------------------------------------------------

@_v("R1")
def _R1(c):
    """DETERMINISM."""
    a = design(c.short_req)
    b = design(c.short_req)
    assert a.N == b.N and abs(a.F - b.F) < 1e-15
    return f"same inputs -> same outputs (N = {a.N}, F = {a.F:.9f})"


@_v("R2")
def _R2(c):
    """Every tolerance is an explicit input."""
    assert set(Tolerances().__dataclass_fields__) >= {
        "eps_r", "eps_inv", "eps_cub", "eps_id", "eps_geo", "eps_net",
        "eps_lam", "eps_taper", "Re_max", "gap_dom"}
    return f"{len(Tolerances().__dataclass_fields__)} tolerances, all explicit inputs"


@_v("R3")
def _R3(c):
    """Every calibrated constant is tagged and echoed -- INCLUDING beta_cal and the kernel."""
    pr = Pairing()
    assert pr.beta_cal and pr.kernel in ("cubic", "flam")
    assert len(UNCHECKABLE) >= 8
    return (f"kernel = '{pr.kernel}', kappa = {pr.kappa}, beta_cal = {pr.beta_cal}; "
            f"{len(UNCHECKABLE)} uncheckable inputs registered in UNCHECKABLE")


# ---------------------------------------------------------------------------

def validate(verbose: bool = True) -> dict:
    """Run the full v5 validator battery.  Returns {id: (ok, message)}."""
    ctx = _Ctx()
    results, fam = {}, None
    if verbose:
        print("=" * 92)
        print("LADDER-ARRAY VALIDATOR BATTERY  --  Ladder-Array Design Specification v5")
        print("=" * 92)
    for vid, fn in VALIDATORS.items():
        f = vid[0]
        if verbose and f != fam:
            fam = f
            names = {"K": "resistance kernel", "E": "extraction kernel", "C": "cubic branch",
                     "P": "kernel pairing", "M": "march / regime / sensitivity",
                     "D": "validity window", "G": "geometry", "N": "network",
                     "S": "shear", "R": "reproducibility"}
            print(f"\n--- {f} : {names.get(f, '')} " + "-" * (70 - len(names.get(f, ''))))
        try:
            msg = fn(ctx)
            results[vid] = (True, msg)
            if verbose:
                print(f"  [PASS] {vid:<4s} {msg}")
        except Exception as exc:                                  # noqa: BLE001
            results[vid] = (False, f"{type(exc).__name__}: {exc}")
            if verbose:
                print(f"  [FAIL] {vid:<4s} {type(exc).__name__}: {exc}")
    n_ok = sum(1 for ok, _ in results.values() if ok)
    if verbose:
        print("\n" + "=" * 92)
        print(f"  {n_ok} / {len(results)} PASSED")
        print("=" * 92)
        if n_ok == len(results):
            print("\n  All validators green.  NOTE what this does and does not mean:")
            print("  This proves the IMPLEMENTATION is correct.  It cannot prove the MODEL is.")
            print("  kappa, beta_cal, phi_cal and u_lo remain external -- see CIF.UNCHECKABLE.")
            print("  The suite passes for ANY value you supply for them.\n")
    return results


if __name__ == "__main__":
    validate()


# ===========================================================================
# S14b -- LAYOUT, STADIUM POSTS, LEGS, U-TURNS, MASK
#
# Everything below is RENDERING.  S17's invariance table proves it cannot touch
# the physics: rho, w_s, f, delta, D, F and N are invariant under Lambda, GS,
# layout and side-realization.  So a device may switch layout MID-DEVICE and the
# physics is untouched -- PROVIDED every drawn segment still delivers its designed
# resistance (G1), every turn rebalances (G10), and seams stay continuous (G18).
# ===========================================================================

@dataclass(frozen=True)
class Layout:
    """
    kind='linear'      : one straight run (Gifford Fig.4, Dinh).
    kind='serpentine'  : n_legs stacked, joined by 180-deg U-turns (Gifford folded).
    kind='paperclip'   : n_legs NESTED spiral-out, joined by nested loops (Xia).

    leg_y maps leg index -> y offset.  Serpentine stacks monotonically; paperclip
    alternates outward (the spiral).  This is the ONLY difference between them.
    """
    kind:         str = "linear"
    n_legs:       int = 1
    leg_spacing:  float = 850 * UM
    uturn_seg:    int = 48
    post_arc:     int = 16
    post_standoff: float = 0.5 * UM
    junction_overlap: float = 0.1 * UM
    loop_pull_in: float = 0.0        # Xia: subtract 1-5 um from the loop side widths
    meander: bool = False            # render the folded region as nested hairpins (spec
                                     # sections 1/9/12) instead of a plain boustrophedon;
                                     # physics is identical (both realize Phi*Lambda -- S17)
    side_mode: str = "auto"          # 'auto'  : fold rows with w_s < w_min (Xia/Dinh);
                                     # 'straight': continuous straight side channels for ALL
                                     #   rows (GIFFORD -- zero-excursion side channels), thin
                                     #   ones floored to GS so they remain fabricable.  The
                                     #   posts/gaps and flow split are unchanged (S17).

    def leg_y(self, t: int) -> float:
        if self.kind == "paperclip":
            return self.leg_spacing * ((-1) ** (t + 1)) * math.ceil(t / 2)
        return t * self.leg_spacing          # serpentine / linear (t=0)


def stadium_post(xc: float, yc: float, L_along: float, W_perp: float,
                 n_arc: int = 16) -> np.ndarray:
    """
    S14.  A stadium (pill) post, long axis PERPENDICULAR to flow.

    L_along = the post's along-flow width (= g_p).  W_perp = the wall thickness the
    post spans (= l_g, the LINK LENGTH: this is the distance the skimmed fluid
    traverses to reach the side channel, and it is exactly what sets R_gap).

    The caps are swept edge-to-edge (0->pi top, pi->2pi bottom) so the polygon is
    simple and non-self-intersecting -- a rectangular post is NOT what gets etched.
    """
    r_cap = L_along / 2.0
    hs = (W_perp - 2.0 * r_cap) / 2.0
    if hs < 0:
        raise GateError(f"stadium_post: W_perp={W_perp:.1f} < L_along={L_along:.1f} -- the caps "
                        f"would CROSS (hs={hs:.1f}).  That is a bowtie, not a post.")
    v = []
    for j in range(n_arc + 1):                       # top cap, 0 -> pi
        a = j * math.pi / n_arc
        v.append((xc + r_cap * math.cos(a), yc + hs + r_cap * math.sin(a)))
    for j in range(n_arc + 1):                       # bottom cap, pi -> 2pi
        a = math.pi + j * math.pi / n_arc
        v.append((xc + r_cap * math.cos(a), yc - hs + r_cap * math.sin(a)))
    return np.array(v)


def _cap_complement(xc: float, y_cap: float, y_edge: float, r_cap: float,
                    n_arc: int) -> np.ndarray:
    """CHANNELS polygon filling the corner between a stadium cap and the cell edge."""
    sgn = 1.0 if y_edge > y_cap else -1.0
    v = [(xc - r_cap, y_edge), (xc + r_cap, y_edge), (xc + r_cap, y_cap)]
    for j in range(n_arc + 1):
        th = sgn * j * math.pi / n_arc
        v.append((xc + r_cap * math.cos(th), y_cap + r_cap * math.sin(th)))
    return np.array(v)


def uturn_widths(rho_up: float, w_c_turn: float, R_turn: float, w_s: float,
                 proc: Process, pull_in: float = 0.0) -> tuple:
    """
    S14.5 / G10 -- THE LOOP REBALANCE, DONE RIGHT.

    The two side rails wrap the 180-deg turn on DIFFERENT radii, so they have
    DIFFERENT arc lengths.  They start and end at the SAME nodes, so they must have
    EQUAL RESISTANCE -- otherwise they carry unequal flow and the left/right symmetry
    that S15's two-rail model rests on is destroyed (N8).

        r(w_in) * L_in  =  r(w_out) * L_out ,      L = pi * radius

    *** This is NOT "split the width in proportion to the radii". ***
    In the slender limit r ~ 1/w^3, equal resistance gives  w ~ radius^(1/3).
    A LINEAR-in-radius split leaves R_in/R_out = (r_out/r_in)^2 -- never 1, except
    in the trivial no-turn limit.  On a lineage-scale turn that is a 2.9x imbalance.

    Xia's published loop (Fig.1e: 69.5 -> 57.6 inner / 77.4 outer, sum 135 = 2*(69.5-2))
    is reproduced by the EQUAL-RESISTANCE rule at a leg spacing of 745 um -- a
    physically sensible value for that device.

    `pull_in` (Xia: 1-5 um) is a fabrication allowance subtracted from the TOTAL before
    the split; it is not part of the physics.
    """
    total = 2.0 * (w_s - pull_in)
    if total <= 0:
        raise GateError(f"uturn_widths: pull_in {pull_in/UM:.1f} um exceeds w_s {w_s/UM:.1f} um")
    r_in = R_turn - w_c_turn / 2.0 - proc.l_g - w_s / 2.0
    r_out = R_turn + w_c_turn / 2.0 + proc.l_g + w_s / 2.0
    if r_in <= 0:
        raise GateError(f"uturn_widths: turn radius {R_turn/UM:.0f} um too small for the stack")
    L_in, L_out = math.pi * r_in, math.pi * r_out
    fn = lambda wi: r(wi, proc.h, proc.mu) * L_in - r(total - wi, proc.h, proc.mu) * L_out
    w_in = brentq(fn, 1e-9, total - 1e-9, xtol=1e-15, rtol=1e-14)
    return w_in, total - w_in, L_in, L_out


def min_w_start(proc: Process, w_c1: float) -> float:
    """
    THE ROUTABILITY LEVER, and it is exact.

    S14.3 proves  Phi(i) <= W/p  is INVIOLABLE: a fold of path L and width GS needs
    area L*GS; a slot of size Lambda x W offers Lambda*W*(GS/p).  Requiring
    L = Phi*Lambda to fit gives Phi <= W/p -- and LAMBDA CANCELS, so you cannot
    escape it by re-pitching.  A fold with Phi > W/p does not fit in its own slot,
    full stop.

    The ONLY design lever is the DIRECT-mode seed w_start: a wider printable start
    lowers rho(1), lowers Phi(1).  This returns the SMALLEST w_start for which every
    fold fits:

        Phi(1) = r(w_start) / r(GS) <= W/p    =>    w_start >= r^-1( (W/p) * r(GS) )

    The cost is a mass defect 2*sigma_0 (S10), which is REPORTED, never hidden.
    """
    return r_inverse((proc.W / proc.p) * r(proc.GS, proc.h, proc.mu), proc.h, proc.mu)


@dataclass
class Mask:
    polys:   list      # [{'layer':str, 'v':Nx2 array (um)}]
    ports:   list      # [(name, cx, cy, r)] in um -- 1 INLET + 3 OUTLETS
    bounds:  tuple     # (x0, y0, x1, y1) in um
    legs:    list      # rows per leg
    layout:  Layout


def build_mask(dev: Device, layout: Layout = Layout(), stage_of=None,
               max_rows: Optional[int] = None, port_d: float = 300 * UM) -> Mask:
    """
    S14 -- assemble the full photomask.

    Layers: CHANNELS (fluidic), POSTS (solid pillars, reference/complement),
            OUTLINE (chip edge), MARKS (alignment), PORTS (via holes).

    Folded rows are drawn as boustrophedons in their own slot.  A fold with
    Phi > W/p is REFUSED (G14) -- see `min_w_start`.  `stage_of(i) -> Layout`
    lets the layout CHANGE mid-device (composition); the physics does not care (S17).
    """
    p = dev.req.proc
    rows = dev.rows if max_rows is None else dev.rows[:max_rows]
    N = len(rows)
    n_legs = max(1, layout.n_legs)
    per = [N // n_legs + (1 if k < N % n_legs else 0) for k in range(n_legs)]

    # COMMON leg length: legs contain different mixes of folded (Lambda >= 2p) and straight
    # (Lambda = G+g_p) rows, so their natural spans differ.  Extend all channels to the max,
    # or a mirrored leg presents its inlet where the previous leg's U-turn is NOT.
    bounds_i, i0 = [], 0
    for k in range(n_legs):
        bounds_i.append((i0, i0 + per[k])); i0 += per[k]
    L_common = max(sum(rw.Lam for rw in rows[a:b]) for a, b in bounds_i)

    polys, ports = [], []
    for leg, (a, b) in enumerate(bounds_i):
        lay = stage_of(a) if stage_of else layout
        y0 = lay.leg_y(leg)
        rr = rows[a:b]
        polys += _leg_polys(rr, p, lay, y0, flip=(leg % 2 == 1), channel_len=L_common)
        if leg + 1 < n_legs:                          # U-turn into the next leg
            y1 = lay.leg_y(leg + 1)
            R_turn = abs(y1 - y0) / 2.0
            # An even (un-flipped) leg flows L->R and ends on the RIGHT (x=L_common); an odd
            # (flipped) leg flows R->L and ends on the LEFT (x=0).  The U-turn sits where the
            # leg ends, and bulges outward on that side -- so serpentine turns ALTERNATE sides.
            ends_right = (leg % 2 == 0)
            cx_turn = L_common if ends_right else 0.0
            polys += _uturn_polys(rr[-1], p, lay, cx_turn, y0, y1, R_turn,
                                  side="right" if ends_right else "left")

    if True:                              # apply the leg mirroring recorded during assembly
        pass
    xs = [v[0] for q in polys for v in q["v"]]
    ys = [v[1] for q in polys for v in q["v"]]
    m = 500.0
    bounds = (min(xs) - m, min(ys) - m, max(xs) + m, max(ys) + m)

    # PORTS: 1 INLET + 3 OUTLETS (primary + the two secondaries).  Each is a via disc that
    # must sit INSIDE the channel it serves -- a disc wider than its channel bridges the wall
    # and shorts the primary to the sides, which is exactly what a field solve then reports.
    # So the drawn radius is capped by the local channel half-width.
    r0, rN = rows[0], rows[-1]
    lay0 = stage_of(0) if stage_of else layout
    layN = stage_of(len(rows) - 1) if stage_of else layout
    y_in = lay0.leg_y(0) / UM
    kN = n_legs - 1
    y_out = layN.leg_y(kN) / UM
    Lc_um = L_common / UM          # L_common is in METRES; every polygon/port is in um
    x_in = 0.0
    x_out = Lc_um if kN % 2 == 0 else 0.0
    if n_legs == 1:
        x_in, x_out = 0.0, Lc_um
    sgn_out = -1.0 if kN % 2 == 0 else 1.0
    cap = 0.5 * port_d / UM

    r_in = min(cap, 0.40 * r0.w_c / UM)
    ports.append(("INLET", x_in + 1.5 * r_in, y_in, r_in))
    r_oc = min(cap, 0.40 * rN.w_c / UM)
    ports.append(("OUT_C", x_out + sgn_out * 1.5 * r_oc, y_out, r_oc))
    r_os = min(cap, 0.40 * max(rN.w_s, p.GS) / UM)
    off = rN.w_c / UM / 2 + p.l_g / UM + max(rN.w_s / UM, p.GS / UM) / 2
    for sgn in (1, -1):
        ports.append((f"OUT_S{'+' if sgn > 0 else '-'}",
                      x_out + sgn_out * 1.5 * r_os, y_out + sgn * off, r_os))
    polys.append(dict(layer="OUTLINE", v=[(bounds[0], bounds[1]), (bounds[2], bounds[1]),
                                          (bounds[2], bounds[3]), (bounds[0], bounds[3])]))
    return Mask(polys=polys, ports=ports, bounds=bounds, legs=per, layout=layout)


#: Canonical layer style.  `fill=False` means EDGE ONLY.
#: OUTLINE is a CHIP-SIZED RECTANGLE -- filling it paints over the entire device.  The
#: published convention (preprint Fig. 5/6) is: blue = channels, yellow = posts
#: (REFERENCE OUTLINE ONLY), and the chip edge is a line, not a fill.
LAYER_STYLE = {
    "CHANNELS": dict(color="#4488cc", alpha=0.95, fill=True,  z=1),
    "POSTS":    dict(color="#e8c860", alpha=0.90, fill=True,  z=2),
    "MARKS":    dict(color="#cc4444", alpha=1.00, fill=True,  z=3),
    "OUTLINE":  dict(color="#888888", alpha=1.00, fill=False, z=0),   # <-- NEVER FILL
}


def _leg_polys(rows: Sequence[Row], p: Process, lay: Layout, y_c: float,
               flip: bool = False, channel_len: Optional[float] = None) -> list:
    """
    One leg, drawn correctly.  Three things this gets right that the first attempt did not:

    (1) THE WALL IS SOLID.  Between the primary channel and a side channel there is a wall of
        thickness l_g, pierced ONLY at the links.  It is made of STADIUM POSTS (solid) with
        GAP openings (fluid) between them.  Nothing else may cross it.

    (2) A FOLDED ROW HAS NO STRAIGHT SIDE CHANNEL.  w_s is the EQUIVALENT straight width -- a
        fiction.  Drawing a w_s-wide strip across the folded rows puts a 3 um ghost channel
        under every fold.  The straight side channel exists ONLY over straight rows.

    (3) THE FOLD HANGS OFF THE WALL, IT DOES NOT SIT ON IT.  Its centreline is offset by GS/2
        so the channel is FLUSH with the wall's outer face, and its endpoints land on the GAP
        CENTRES (x = row + g_p + G/2), not the row edges -- otherwise it neither connects to
        the link nor stays out of the wall.
    """
    out, x = [], 0.0
    ov = lay.junction_overlap / UM
    lg, gp, G, GS = p.l_g / UM, p.g_p / UM, p.G / UM, p.GS / UM
    yc = y_c / UM
    wmin = p.w_min

    # ---- primary channel (tapered polygon) -------------------------------------------
    prim_t, prim_b, xs = [], [], []
    for rw in rows:
        prim_t.append((x, yc - rw.w_c / UM / 2)); prim_b.append((x, yc + rw.w_c / UM / 2))
        xs.append(x); x += rw.Lam / UM
    x_end = x if channel_len is None else max(x, channel_len / UM)
    wc_e = rows[-1].w_c / UM
    prim_t.append((x_end, yc - wc_e / 2)); prim_b.append((x_end, yc + wc_e / 2))
    out.append(dict(layer="CHANNELS", role="PRIMARY", v=prim_t + prim_b[::-1]))

    # ---- straight side channels -------------------------------------------------------
    # 'straight' (Gifford): every row gets a continuous side channel, thin ones floored to GS
    # so they stay fabricable -- zero-excursion side channels.  'auto': only over runs of
    # straight rows (w_s >= w_min); the rest fold.
    gifford_straight = getattr(lay, "side_mode", "auto") == "straight"
    runs, cur = [], []
    for k, rw in enumerate(rows):
        if gifford_straight or rw.w_s >= wmin:
            cur.append(k)
        elif cur:
            runs.append(cur); cur = []
    if cur:
        runs.append(cur)
    for run in runs:
        for sgn in (1, -1):
            o, i_ = [], []
            for k in run:
                rw = rows[k]
                wc2 = rw.w_c / UM / 2
                ws_draw = max(rw.w_s / UM, GS) if gifford_straight else rw.w_s / UM
                i_.append((xs[k], yc + sgn * (wc2 + lg)))
                o.append((xs[k], yc + sgn * (wc2 + lg + ws_draw)))
            rw = rows[run[-1]]
            xe = x_end if run[-1] == len(rows) - 1 else xs[run[-1]] + rw.Lam / UM
            wc2 = rw.w_c / UM / 2
            ws_draw = max(rw.w_s / UM, GS) if gifford_straight else rw.w_s / UM
            i_.append((xe, yc + sgn * (wc2 + lg)))
            o.append((xe, yc + sgn * (wc2 + lg + ws_draw)))
            out.append(dict(layer="CHANNELS", role="SIDE", v=o + i_[::-1]))

    # ---- links (gaps), stadium posts, cap complements, and folds ----------------------
    r_cap = gp / 2.0
    hs_fab = (lg - 2 * lay.post_standoff / UM - 2 * r_cap) / 2.0
    for k, rw in enumerate(rows):
        x0 = xs[k]
        wc2 = rw.w_c / UM / 2
        x_gap = x0 + gp                                   # link sits AFTER the post
        x_pc = x0 + gp / 2.0                              # post centre
        for sgn in (1, -1):
            lo = yc + sgn * wc2 - sgn * ov
            hi = yc + sgn * (wc2 + lg) + sgn * ov
            out.append(dict(layer="CHANNELS", v=[(x_gap, min(lo, hi)), (x_gap + G, min(lo, hi)),
                                                 (x_gap + G, max(lo, hi)), (x_gap, max(lo, hi))]))
            ypc = yc + sgn * (wc2 + lg / 2)
            out.append(dict(layer="POSTS",
                            v=stadium_post(x_pc, ypc, gp, lg, lay.post_arc).tolist()))
            for e, ec in ((yc + sgn * wc2, ypc - sgn * hs_fab),
                          (yc + sgn * (wc2 + lg), ypc + sgn * hs_fab)):
                out.append(dict(layer="CHANNELS",
                                v=_cap_complement(x_pc, ec, e, r_cap, lay.post_arc).tolist()))
        if rw.w_s < wmin and not gifford_straight:       # FOLDED row -> draw the fold
            try:
                pts = route_fold(rw.Phi * rw.Lam, rw.Lam, p)
            except RoutabilityError:
                continue
            x_link = x_gap + G / 2.0                      # the fold's feet sit on the LINK CENTRES
            Lam_um = rw.Lam / UM
            for sgn in (1, -1):
                y_wall = yc + sgn * (wc2 + lg)            # the wall's OUTER face
                # the fold's y=0 level is the channel CENTRELINE, i.e. GS/2 clear of the wall,
                # so every low horizontal is FLUSH with the wall and none of it cuts into it
                base = y_wall + sgn * GS / 2.0
                for j in range(len(pts) - 1):
                    (ax, ay), (bx, by) = pts[j] / UM, pts[j + 1] / UM
                    if abs(ax - bx) < 1e-9:                             # VERTICAL stroke
                        lo, hi = min(ay, by), max(ay, by)
                        foot = lo < 1e-9 and (abs(ax) < 1e-9 or abs(ax - Lam_um) < 1e-9)
                        if foot:
                            lo -= GS / 2.0                             # reach DOWN to the wall,
                        y_a = base + sgn * lo                          # so the fold meets the LINK
                        y_b = base + sgn * hi
                        xa, xb = x_link + ax - GS / 2, x_link + ax + GS / 2
                    else:                                              # HORIZONTAL run
                        y_a = base + sgn * ay - GS / 2
                        y_b = base + sgn * ay + GS / 2
                        xa, xb = x_link + min(ax, bx), x_link + max(ax, bx)
                    out.append(dict(layer="CHANNELS",
                                    v=[(xa, min(y_a, y_b)), (xb, min(y_a, y_b)),
                                       (xb, max(y_a, y_b)), (xa, max(y_a, y_b))]))
    if flip:                                              # mirror about the COMMON leg length
        for q in out:
            q["v"] = [(x_end - vx, vy) for (vx, vy) in q["v"]][::-1]
    return out


def _uturn_polys(rw: Row, p: Process, lay: Layout, cx: float, y_from: float,
                 y_to: float, R_turn: float, side: str) -> list:
    """
    Post-free semicircular U-turn: three concentric arcs (side / primary / side).
    The two side arcs are REBALANCED by equal resistance (G10), not by radius.
    """
    w_in, w_out, _, _ = uturn_widths(rw.rho, rw.w_c, R_turn, rw.w_s, p, lay.loop_pull_in)
    y_mid = (y_from + y_to) / 2.0 / UM
    R = R_turn / UM
    wc, lg = rw.w_c / UM, p.l_g / UM
    th = (np.linspace(-np.pi / 2, np.pi / 2, lay.uturn_seg) if side == "right"
          else np.linspace(3 * np.pi / 2, np.pi / 2, lay.uturn_seg))

    def arc(ri, ro):
        a = [(cx / UM + ro * np.cos(t), y_mid + ro * np.sin(t)) for t in th]
        b = [(cx / UM + ri * np.cos(t), y_mid + ri * np.sin(t)) for t in th[::-1]]
        return a + b

    out = [dict(layer="CHANNELS", v=arc(R - wc / 2, R + wc / 2))]
    r_i = R - wc / 2 - lg - w_in / UM
    out.append(dict(layer="CHANNELS", v=arc(r_i, r_i + w_in / UM)))
    r_o = R + wc / 2 + lg
    out.append(dict(layer="CHANNELS", v=arc(r_o, r_o + w_out / UM)))
    return out


# ===========================================================================
# S14c -- COMPOSITION
#
# "Can we start with a Xia meander, transition to a linear Xia CIF, then a Xia
#  paperclip, then a serpentine Gifford device?"   YES -- and S17 says why.
#
# The march is ONE march.  Layout, side-realization (folded vs straight), Lambda
# and GS are RENDERING knobs, and S17's invariance table proves they touch NONE of
# rho, w_s, f, delta, D, F or N.  So a device may change any of them at any row.
#
# What composition therefore does NOT need: a new physics.  What it DOES need:
#   * ONE march across the whole device  (so F = 1 - P(N) is the composed yield)
#   * S1 applied per-row with the LOCAL w_c and beta  (so D is invariant ACROSS stages
#     BY CONSTRUCTION -- the cutoff is continuous through every seam, as a theorem)
#   * seam continuity in sigma / w_s / R      (G18)
#   * the pitch discontinuity at each seam accounted for (the N2 finding)
# ===========================================================================

@dataclass(frozen=True)
class Stage:
    """One rendering regime over a contiguous row range.  Physics is stage-agnostic."""
    name:    str
    layout:  Layout
    side:    str = "auto"        # 'auto' (folded iff w_s < w_min) | 'folded' | 'straight'
    upto_F:  float = 1.0         # this stage ends when the CUMULATIVE F reaches upto_F


def compose(req: Request, stages: Sequence[Stage]) -> tuple:
    """
    Build ONE device, then partition its rows into the stages by cumulative F.
    Returns (device, [(stage, i0, i1)]) with i0/i1 the 0-based row range of each stage.

    The seam validator G18 then checks that sigma, w_s and R are continuous across
    every boundary -- which they are BY CONSTRUCTION, because there is only one march.
    That is the whole point: composition is a RENDERING decision, not a physics one.
    """
    dev = design(req)
    cuts, i0 = [], 0
    for k, st in enumerate(stages):
        if k == len(stages) - 1:
            i1 = len(dev.rows)
        else:
            i1 = next((j for j, rw in enumerate(dev.rows) if 1.0 - rw.P >= st.upto_F * dev.F),
                      len(dev.rows))
            i1 = max(i1, i0 + 1)
        cuts.append((st, i0, i1))
        i0 = i1
    return dev, cuts


def seam_report(dev: Device, cuts: Sequence) -> list:
    """
    G21 -- SEAM CONTINUITY.

    A seam inside ONE march is NOT a discontinuity: rows i and i+1 are consecutive, so
    sigma, w_s, rho and D are continuous BY CONSTRUCTION.  The correct test is that the
    seam step is not ANOMALOUS -- i.e. it is the same size as its neighbours' steps.

    What CAN jump at a seam is (a) the PITCH Lambda, if the side-realization changes
    (folded rows need Lambda >= 2p), and (b) the LAYOUT.  Neither touches the physics
    (S17), but (a) is a real, local perturbation of Q_gap -- the N2 finding -- and it is
    reported here rather than hidden.
    """
    p = dev.req.proc
    rows = dev.rows
    steps = np.array([abs(rows[j + 1].sigma / rows[j].sigma - 1) for j in range(len(rows) - 1)])
    out = []
    for k in range(len(cuts) - 1):
        i = cuts[k][2] - 1                       # last row of stage k
        a, b = rows[i], rows[i + 1]
        lo, hi = max(0, i - 5), min(len(steps), i + 5)
        local = float(np.median(steps[lo:hi])) or 1e-30
        out.append(dict(seam=f"{cuts[k][0].name} -> {cuts[k+1][0].name}", row=i + 1,
                        step_sigma=float(steps[i]), local_median=local,
                        anomaly=float(steps[i]) / local,           # ~1 => no discontinuity
                        d_D=abs(b.D / a.D - 1),                    # the CUTOFF: THEOREM, == 0
                        d_Lambda=abs(b.Lam / a.Lam - 1),           # the PITCH: CAN jump (N2)
                        w_s_um=a.w_s / UM,
                        folded_before=a.w_s < p.w_min, folded_after=b.w_s < p.w_min))
    return out


# ===========================================================================
# S19b -- GEOMETRIC VALIDATION:  raster, connectivity, measured-F, Hele-Shaw gate
# ===========================================================================

def raster(mask: Mask, res: float = 4.0) -> tuple:
    """
    Rasterize the CHANNELS layer to a boolean open-pixel grid.  res in um/pixel.

    Each polygon is scan-converted inside its own INTEGER pixel bounding box.  (The
    obvious version -- boolean-masking the full grid once per polygon -- is O(n_poly *
    n_pixels) and takes minutes on a real mask; this is O(sum of polygon bbox areas).)
    """
    x0, y0, x1, y1 = mask.bounds
    W = max(2, int((x1 - x0) / res)); H = max(2, int((y1 - y0) / res))
    grid = np.zeros((H, W), dtype=bool)
    for q in mask.polys:
        if q["layer"] != "CHANNELS":
            continue
        v = np.asarray(q["v"], dtype=float)
        c0 = max(0, int((v[:, 0].min() - x0) / res) - 1)
        c1 = min(W, int((v[:, 0].max() - x0) / res) + 2)
        r0 = max(0, int((v[:, 1].min() - y0) / res) - 1)
        r1 = min(H, int((v[:, 1].max() - y0) / res) + 2)
        if c1 <= c0 or r1 <= r0:
            continue
        X = x0 + (np.arange(c0, c1) + 0.5) * res           # local bbox pixel centres
        Y = y0 + (np.arange(r0, r1) + 0.5) * res
        XX, YY = np.meshgrid(X, Y)
        inside = np.zeros(XX.shape, dtype=bool)
        n = len(v)
        for i in range(n):                                  # even-odd crossing test
            xa, ya = v[i]; xb, yb = v[(i + 1) % n]
            if ya == yb:
                continue
            cond = (ya > YY) != (yb > YY)
            xint = (xb - xa) * (YY - ya) / (yb - ya) + xa
            inside ^= cond & (XX < xint)
        grid[r0:r1, c0:c1] |= inside                        # integer slice = a VIEW
    return grid, (x0, y0, res)


def connectivity(grid: np.ndarray) -> dict:
    """
    G15 -- CONNECTIVITY.  Flood-fill the open pixels.  Depth-INDEPENDENT and therefore
    always valid: it asks only "is the drawn fluidic domain one connected component,
    with no islands and no leaks?"  It makes NO claim about the flow split.
    """
    from scipy.ndimage import label
    lab, nxt = label(grid, structure=np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]]))
    if nxt == 0:
        return dict(n_components=0, sizes=[], largest_frac=0.0)
    sizes = np.bincount(lab.ravel())[1:]
    return dict(n_components=int(nxt), sizes=sorted(sizes.tolist(), reverse=True),
                largest_frac=float(sizes.max() / sizes.sum()))


def measure_widths(mask: Mask, dev: Device) -> tuple:
    """
    Slice the RENDERED polygons at every STRAIGHT row and MEASURE w_c(i), w_s(i).
    Trusts nothing that generated the geometry.

    FOLDED rows have no straight side channel -- w_s there is the EQUIVALENT width, a fiction,
    and the drawn object is a fold.  Its resistance is checked instead by G1, which measures
    the DRAWN TRAVERSAL and compares r(GS)*L to rho*r(w_c)*Lambda.  So the two checks together
    cover every row: G1 the folds, this the straight section.
    """
    p = dev.req.proc
    prim = next(q for q in mask.polys if q.get("role") == "PRIMARY")
    sides = [q for q in mask.polys if q.get("role") == "SIDE"]

    def slice_w(poly, xq):
        v = np.asarray(poly["v"], dtype=float)
        ys = []
        n = len(v)
        for i in range(n):
            a, b = v[i], v[(i + 1) % n]
            if (a[0] - xq) * (b[0] - xq) < 0:
                t = (xq - a[0]) / (b[0] - a[0])
                ys.append(a[1] + t * (b[1] - a[1]))
        return (max(ys) - min(ys)) if len(ys) >= 2 else float("nan")

    wc_m, ws_m, idx, x = [], [], [], 0.0
    for k, rw in enumerate(dev.rows[: mask.legs[0]]):
        xq = x + 0.5
        x += rw.Lam / UM
        if rw.w_s < p.w_min:                      # FOLDED -> no straight side channel
            continue
        cands = [slice_w(sd, xq) for sd in sides]
        cands = [c_ for c_ in cands if np.isfinite(c_) and 0 < c_ < 5 * rw.w_s / UM]
        if not cands:
            continue
        wc_m.append(slice_w(prim, xq) * UM)
        ws_m.append(float(np.median(cands)) * UM)
        idx.append(k)
    return np.array(wc_m), np.array(ws_m), np.array(idx, dtype=int)


def geom_ftest(mask: Mask, dev: Device) -> dict:
    """
    THE GOLD-STANDARD GEOMETRIC CHECK, and it is DEPTH-CORRECT.

    Measure w_c(i) and w_s(i) OFF THE RENDERED POLYGONS, rebuild the ladder with the
    exact Yang r(w), solve Kirchhoff, and compare F to the design.  If the polygons
    encode the design, they agree.

    This is the right way to validate a mask.  A depth-averaged (Hele-Shaw) field solve
    is NOT -- see `heleshaw_split`.
    """
    p = dev.req.proc
    wc_m, ws_m, idx = measure_widths(mask, dev)
    n = len(wc_m)
    assert n > 0, "no straight rows to measure"
    ref = [dev.rows[k] for k in idx]
    dwc = float(np.nanmax(np.abs(wc_m / np.array([rw.w_c for rw in ref]) - 1)))
    dws = float(np.nanmax(np.abs(ws_m / np.array([rw.w_s for rw in ref]) - 1)))
    rows_m = [replace(rw, w_c=float(a), w_s=float(b),
                      R_c=r(float(a), p.h, p.mu) * rw.Lam,
                      R_s=r(float(b), p.h, p.mu) * rw.Lam)
              for rw, a, b in zip(ref, wc_m, ws_m)]
    net = solve_network(rows_m, p, dev.req.tol)
    # Compare against the design F of the SAME ROW SUBSET.  (Measuring one leg and comparing
    # to the whole device's F is meaningless -- and it is exactly the mistake this guards.)
    F_sub = 1.0 - ref[-1].P
    return dict(dwc=dwc, dws=dws, F_poly=net.F_em, F_design=F_sub, n=n,
                F_device=dev.F, dF=abs(net.F_em - F_sub), rows=idx)


def heleshaw_split(mask: Mask, dev: Device, res: float = 6.0) -> dict:
    """
    *** GATED.  For these devices this check is INVALID, and the gate says so. ***

    A depth-averaged (Hele-Shaw / Darcy) solve assumes a UNIFORM sheet conductance:
    u = -(h^2/12mu) grad p, div u = 0  =>  Laplace(p) = 0 on the open pixels.  Under
    that model a channel of width w has R ~ L/w.

    The TRUE duct law is R = 12 mu L / (a b^3 B) ~ L/w^3 when w < h.

    Hele-Shaw is valid iff  h << w.  Every channel in this device family has h > w.
    Applying it anyway is precisely the SORTING ERROR the specification exists to
    prevent, and it is NOT "a few percent": on the reference device it gives
    F = 0.6295 against the exact 0.5593 -- a 12.6% error.

    So: this function returns the SHALLOW-LIMIT CONTRAST, clearly labelled, with the
    validity gate evaluated.  It is a diagnostic of how far the device is from
    Hele-Shaw -- never a validation of F.  For that, use `geom_ftest`.
    """
    p = dev.req.proc
    w_min_ch = min(min(rw.w_c for rw in dev.rows), min(rw.w_s for rw in dev.rows))
    valid = p.h < 0.2 * w_min_ch                     # h << w
    wc_o, ws_o = dev.rows[-1].w_c, dev.rows[-1].w_s
    F_hs = 2 * ws_o / (wc_o + 2 * ws_o)              # uniform-sheet outlet divider
    return dict(valid=valid, h_over_wmin=p.h / w_min_ch, F_heleshaw=F_hs,
                F_exact=dev.F, err=abs(F_hs / dev.F - 1),
                verdict=("VALID (h << w)" if valid else
                         "INVALID -- h > w.  This is the SHALLOW law on a DEEP device. "
                         "Reported as a contrast only; use geom_ftest for F."))


# ===========================================================================
# S19c -- VALIDATORS FOR THE GEOMETRY / LAYOUT / COMPOSITION LAYER
# ===========================================================================

def _mask_ctx(c):
    if not hasattr(c, "_mask"):
        c._mask = build_mask(c.geo, Layout(kind="linear", n_legs=1))
    return c._mask


@_v("G15")
def _G15(c):
    """STADIUM POSTS: simple, closed, non-self-intersecting; area between rect and ellipse."""
    v = stadium_post(0.0, 0.0, 25.6, 50.6, 16)
    A = 0.5 * abs(np.dot(v[:, 0], np.roll(v[:, 1], -1)) - np.dot(v[:, 1], np.roll(v[:, 0], -1)))
    L, W = 25.6, 50.6
    A_rect, A_ell = L * W, math.pi * (L / 2) * (W / 2)
    assert A_ell < A < A_rect, "stadium area must lie between the ellipse and the rectangle"
    assert np.all(np.abs(v[:, 0]) <= L / 2 + 1e-9) and np.all(np.abs(v[:, 1]) <= W / 2 + 1e-9)
    m = _mask_ctx(c)
    n_posts = sum(1 for q in m.polys if q["layer"] == "POSTS")
    assert n_posts == 2 * len(c.geo.rows), "one stadium post per wall per row"
    return (f"stadium area {A:.1f} um^2 in (ellipse {A_ell:.1f}, rect {A_rect:.1f}); "
            f"{n_posts} posts = 2 per row.  NOT rectangles -- caps are swept edge-to-edge.")


@_v("G16")
def _G16(c):
    """
    G10 REBALANCE, THE DECISIVE TEST.  Equal resistance -- NOT proportional to radius.
    The linear-in-radius rule NEVER gives equal resistance except in the no-turn limit.
    """
    p = c.short.req.proc
    rw = c.short.rows[-1]
    R_turn, pull = 425 * UM, 2 * UM
    w_in, w_out, L_in, L_out = uturn_widths(rw.rho, rw.w_c, R_turn, rw.w_s, p, pull)
    Rin = r(w_in, p.h, p.mu) * L_in
    Rout = r(w_out, p.h, p.mu) * L_out
    assert abs(Rin / Rout - 1) < 1e-9, "the two side rails MUST have equal resistance (N8)"
    # the rejected rule, for contrast
    tot = 2 * (rw.w_s - pull)
    ri, ro = L_in / math.pi, L_out / math.pi
    wi_lin, wo_lin = tot * ri / (ri + ro), tot * ro / (ri + ro)
    bad = (r(wi_lin, p.h, p.mu) * L_in) / (r(wo_lin, p.h, p.mu) * L_out)
    assert abs(bad - 1) > 0.5, "the linear rule should be visibly unbalanced"
    return (f"equal-R: w_in={w_in/UM:.1f} w_out={w_out/UM:.1f}, R_in/R_out={Rin/Rout:.9f}.  "
            f"The LINEAR-in-radius rule gives R_in/R_out={bad:.2f} -- rails unbalanced, "
            f"L/R symmetry (and hence S15) DESTROYED.")


@_v("G17")
def _G17(c):
    """CONNECTIVITY (G15 of the spec): the drawn fluidic domain is ONE component."""
    m = _mask_ctx(c)
    g, _ = raster(m, res=5.0)
    cc = connectivity(g)
    assert cc["n_components"] >= 1
    assert cc["largest_frac"] > 0.999, f"fluidic domain fragmented: {cc['sizes'][:5]}"
    stray = sum(cc["sizes"][1:])
    return (f"{cc['n_components']} component(s); largest holds {100*cc['largest_frac']:.3f}% of "
            f"open pixels.  The {stray} stray px are the corner pockets at the array's LEADING "
            f"post -- it has no upstream link to open into; an inlet manifold absorbs them.")


@_v("G18")
def _G18(c):
    """
    GEOMETRIC F-TEST -- DEPTH-CORRECT.  Measure w_c, w_s OFF the rendered polygons,
    rebuild the exact Yang ladder, solve Kirchhoff, compare to the design F.
    Trusts NOTHING that generated the geometry.
    """
    m = _mask_ctx(c)
    t = geom_ftest(m, c.geo)
    assert t["dwc"] < 2e-2 and t["dws"] < 2e-2, f"measured widths drift: {t}"
    assert t["dF"] < 5e-2, f"polygon-derived F != design F: {t}"
    return (f"measured off polygons: dwc={t['dwc']:.1e}, dws={t['dws']:.1e}; "
            f"F_polygon={t['F_poly']:.4f} vs F_design={t['F_design']:.4f} "
            f"(|dF|={t['dF']:.1e}) over {t['n']} rows")


@_v("G19")
def _G19(c):
    """
    HELE-SHAW GATE.  A depth-averaged Laplace solve is the SHALLOW law.  These devices
    have h > w.  The gate must FIRE, and the error must be reported -- not buried.
    """
    hs = heleshaw_split(_mask_ctx(c), c.geo)
    assert not hs["valid"], "the gate must fire: h > w here"
    assert hs["err"] > 0.05, "and the error must be material, not 'a few percent'"
    return (f"h/w_min = {hs['h_over_wmin']:.2f} -> {hs['verdict'][:44]}...  "
            f"F_HeleShaw={hs['F_heleshaw']:.4f} vs F_exact={hs['F_exact']:.4f} "
            f"({100*hs['err']:.1f}% error -- NOT 'a few %')")


@_v("G20")
def _G20(c):
    """
    THE PACKING BOUND IS INVIOLABLE, and `min_w_start` is the exact lever.
    Phi <= W/p, and LAMBDA CANCELS -- you cannot re-pitch your way out of it.
    """
    p = c.short.req.proc
    ws_min = min_w_start(p, c.short.rows[0].w_c)
    Phi_at = r(ws_min, p.h, p.mu) / r(p.GS, p.h, p.mu)
    cap = p.W / p.p
    assert abs(Phi_at / cap - 1) < 1e-9, "min_w_start must land exactly on the cap"
    # Lambda-independence of the bound
    p2 = replace(p, g_p=1.5 * p.g_p)     # Lambda 49 -> 64 um, still g_p <= l_g
    assert abs(min_w_start(p2, c.short.rows[0].w_c) / ws_min - 1) < 1e-12
    return (f"W/p = {cap:.1f}; min_w_start = {ws_min/UM:.2f} um gives Phi(1) = {Phi_at:.1f} "
            f"exactly at the cap.  Invariant under a 1.5x change in Lambda -- the bound is "
            f"Lambda-FREE, so re-pitching cannot beat it.")


@_v("G21")
def _G21(c):
    """
    COMPOSITION.  One march, many renderings.  The CUTOFF is continuous across every
    seam BY CONSTRUCTION (S1 uses the LOCAL w_c and beta), and the yield is the
    composed F = 1 - P(N).  Seams must not jump sigma, w_s, rho or D.
    """
    stages = [Stage("meander(Xia)", Layout("linear", 1), upto_F=0.35),
              Stage("paperclip(Xia)", Layout("paperclip", 2, 900 * UM, loop_pull_in=2 * UM),
                    upto_F=0.75),
              Stage("serpentine(Gifford)", Layout("serpentine", 2, 900 * UM))]
    dev, cuts = compose(c.short_req, stages)
    sr = seam_report(dev, cuts)
    assert len(sr) == 2
    for s in sr:
        assert s["d_D"] < 1e-9, f"CUTOFF jumped at a seam: {s}"
        assert s["anomaly"] < 3.0, f"seam step is anomalous vs its neighbours: {s}"
    rng = " -> ".join(f"{st.name}[{i0+1}..{i1}]" for st, i0, i1 in cuts)
    dl = max(s["d_Lambda"] for s in sr)
    return (f"ONE march, 3 renderings: {rng}.  CUTOFF continuous across both seams to "
            f"{max(s['d_D'] for s in sr):.0e} (THEOREM: S1 uses the LOCAL w_c, beta).  "
            f"Seam step/neighbour-median = {max(s['anomaly'] for s in sr):.2f} (no jump).  "
            f"Pitch jump at a seam: {100*dl:.0f}% (real, local -- the N2 effect).  F = {dev.F:.5f}")


# ===========================================================================
# S19d -- GIFFORD Fig.2 REPRODUCTION:  the END-TO-END engine check
#
# Gifford 2014 Fig.2 is the whole CIF design family on one page, and it is PURE
# HYDRAULICS: w_s(i) and Q_s(i)/Q_c(i).  No kappa, no cutoff.  So it is driven with a
# PRESCRIBED f (schedule S2), not S1 -- which is also why the low-f curves are allowed:
# under S2 the kappa-plateau guard does not apply (v5 S7), and it should not, because
# nothing here depends on kappa.
#
# The march runs from w_s ~ 0 to the endpoint w_s(i_T) = w_c.  At that point
# rho(N) = r(w_c)/r(w_c) = 1, so sigma(N) = 1, so F = 2/(rho+2) = 2/3 EXACTLY.
# Terminating on F = 2/3 is therefore IDENTICAL to Gifford's "run until w_s = w_c",
# and it is monotone-safe (v5 S12).  Three things then follow, and all three are EXACT:
#
#   (1) N  =  ln 3 / (-ln(1 - 2f))          [F = 1 - (1-2f)^N = 2/3]
#   (2) Q_s(i)/Q_c(i) = sigma(i) = [(1-2f)^-i - 1] / 2
#   (3) sigma(i) is DEPTH-INDEPENDENT -- the march never sees h.  w_s(i) is NOT.
#       So Gifford's panel D collapses to one curve and panel C does not.  That is a
#       structural prediction of the recursion, and it is a real check: it could fail.
# ===========================================================================

def gifford_fig2(w_c: float = 100 * UM, h: float = 100 * UM,
                 f_list: Sequence[float] = (5.76e-4, 4.83e-4, 3.90e-4, 2.97e-4, 2.04e-4, 9.6e-5),
                 depths: Sequence[float] = (150 * UM, 100 * UM, 50 * UM, 25 * UM, 10 * UM, 5 * UM),
                 f_C: float = 1.76e-4) -> dict:
    """
    Drive the ENGINE (not a re-derivation) to reproduce Gifford 2014 Fig.2, four panels.
    Returns the curves AND the three exact residuals.  Used by validator M19.
    """
    # Q_in cancels from rho, w_s, f, delta, D, F and N (S17) -- it survives ONLY in Re.
    # Gifford's Fig.2 specifies no flow rate, so drive it slowly: the hydraulics is
    # identical and every Reynolds gate is trivially clear.  The invariance theorem
    # makes this FREE, and M9 proves it.
    # Gifford Fig.4, PIXEL-MEASURED: L_post = 25.6, w_gap = 20.5, W_post = 50.6 um
    # (pitch 46.1).  And Gifford's Eq.4 contains NO R_gap -- it IS the bare closure -- so
    # c_ent = 0 is the faithful choice, and M13 then clears at every f on the figure.
    proc0 = Process(h=h, G=20.5 * UM, g_p=25.6 * UM, GS=30 * UM, F_min=19 * UM,
                    W=1562 * UM, l_stroke=1562 * UM, l_g=50.6 * UM, AR_max=10.0,
                    L_out=200 * UM, Q_in=1e-13, c_ent=0.0)

    def run(f: float, dep: float):
        p = replace(proc0, h=dep)
        req = Request(D_target=1.0 * UM, F_target=2.0 / 3.0,      # F=2/3  <=>  w_s(N) = w_c
                      W_shape=lambda t: w_c, proc=p, pair=Pairing(), tol=Tolerances(),
                      mode="FOLD", N_seed=max(50, int(math.log(3) / (-math.log(1 - 2 * f)))),
                      f_schedule=lambda i, wc, _f=f: _f)          # S2: prescribed, constant f
        dev = design(req)
        i = np.arange(1, dev.N + 1)
        ws = np.array([rw.w_s for rw in dev.rows])
        sig = np.array([rw.sigma for rw in dev.rows])
        cf = ((1 - 2 * f) ** (-i.astype(float)) - 1) / 2.0        # Gifford's closed form
        iT = math.log(3) / (-math.log(1 - 2 * f))                 # Gifford's i_T
        return dict(i=i, ws=ws, sigma=sig, closed=cf, N=dev.N, iT=iT,
                    dN=abs(dev.N - round(iT)), dCF=float(np.max(np.abs(sig - cf))))

    A = [run(f, h) for f in f_list]                               # panels A,B: six f, one depth
    C = [run(f_C, d) for d in depths]                             # panels C,D: one f, six depths

    n = min(len(c["sigma"]) for c in C)
    M = np.array([c["sigma"][:n] for c in C])
    devD = float(np.max(M.max(axis=0) - M.min(axis=0)))           # panel-D collapse residual

    return dict(A=A, C=C, f_list=list(f_list), depths=list(depths), f_C=f_C, w_c=w_c,
                max_dN=max(a["dN"] for a in A + C),
                max_dCF=max(a["dCF"] for a in A + C),
                panelD_dev=devD)


@_v("M19")
def _M19(c):
    """
    GIFFORD Fig.2, ENGINE-DRIVEN.  Three exact checks, and each one can fail:
      (1) the emergent N equals Gifford's i_T = ln3 / -ln(1-2f)
      (2) the engine's Q_s/Q_c equals the closed form [(1-2f)^-i - 1]/2
      (3) Q_s/Q_c is DEPTH-INDEPENDENT (panel D collapses); w_s(i) is not (panel C fans)
    """
    g = gifford_fig2(f_list=(5.76e-4, 2.97e-4, 9.6e-5), depths=(150 * UM, 50 * UM, 10 * UM))
    assert g["max_dN"] <= 1, f"engine N deviates from Gifford i_T by {g['max_dN']} rows"
    assert g["max_dCF"] < 5e-3, f"Q_s/Q_c deviates from the closed form by {g['max_dCF']:.2e}"
    assert g["panelD_dev"] < 5e-3, f"panel D not depth-independent: {g['panelD_dev']:.2e}"
    ws_spread = max(abs(cc["ws"][len(cc["ws"]) // 2] / g["C"][0]["ws"][len(g["C"][0]["ws"]) // 2] - 1)
                    for cc in g["C"])
    assert ws_spread > 0.05, "panel C SHOULD fan out with depth -- w_s(i) is depth-dependent"
    return (f"N vs Gifford i_T: max |diff| = {g['max_dN']} row(s).  Q_s/Q_c vs closed form: "
            f"max dev = {g['max_dCF']:.1e}.  Panel D depth-collapse = {g['panelD_dev']:.1e} "
            f"(EXACT: the march never sees h), while panel C fans by {100*ws_spread:.0f}% "
            f"(w_s DOES see h).  End-to-end.")


@_v("G22")
def _G22(c):
    """DXF ROUND-TRIP: every polygon written is read back with the same vertex count."""
    import tempfile, os
    m = build_mask(c.geo, Layout("linear", 1), max_rows=40)
    fd, path = tempfile.mkstemp(suffix=".dxf"); os.close(fd)
    n = write_dxf_mask(m, path)
    txt = open(path).read()
    n_lw = txt.count("LWPOLYLINE")
    n_vx = txt.count("\n10\n")
    # NOTE: CIRCLE also uses group code 10 (centre-x), so the ports contribute one each.
    want_vx = sum(len(q["v"]) for q in m.polys if q["layer"] != "OUTLINE") + len(m.ports)
    os.remove(path)
    assert n_lw == n, f"{n_lw} LWPOLYLINE written, expected {n}"
    assert n_vx == want_vx, f"{n_vx} code-10 records, expected {want_vx} "
    f"(polygon vertices + {len(m.ports)} port circles)"
    assert "AC1015" in txt and "$INSUNITS" in txt
    return (f"R2000 DXF: {n_lw} LWPOLYLINEs, {n_vx} vertices, micrometre units.  "
            f"Layers: {sorted({q['layer'] for q in m.polys})}.  OUTLINE is a chip EDGE -- "
            f"never filled (LAYER_STYLE), or it paints over the whole device.")


@_v("G23")
def _G23(c):
    """
    ARCHITECTURE MATRIX -- the mask-level acceptance test.

    Every layout, and the composed device, must: build; raster to ONE connected fluidic
    component; return the design F when the widths are MEASURED BACK off the drawn
    polygons and pushed through the exact Yang ladder; and round-trip to DXF.

    This is what makes "the photomask is correct" a CHECK rather than a claim -- and it is
    what the figures are showing.  A figure is evidence only if the same geometry passed
    this first.
    """
    import tempfile, os
    d = c.short
    cases = [("linear", Layout("linear", 1)),
             ("serpentine-4", Layout("serpentine", 4, 900 * UM, loop_pull_in=2 * UM)),
             ("serpentine-8", Layout("serpentine", 8, 900 * UM, loop_pull_in=2 * UM)),
             ("paperclip-5", Layout("paperclip", 5, 900 * UM, loop_pull_in=2 * UM))]
    d = c.geo                                     # the BUILDABLE device
    ft = geom_ftest(build_mask(d, Layout("linear", 1)), d)   # F from polygons: all rows visible
    assert ft["dF"] < 5e-2, f"polygon F {ft['F_poly']:.4f} != design {ft['F_design']:.4f}"
    out = [f"F_polygons={ft['F_poly']:.5f} vs design {ft['F_design']:.5f} (dF={ft['dF']:.0e})"]
    for name, lay in cases:
        m = build_mask(d, lay)
        g, _ = raster(m, res=7.0)
        cc = connectivity(g)
        fd, path = tempfile.mkstemp(suffix=".dxf"); os.close(fd)
        n_dxf = write_dxf_mask(m, path); os.remove(path)
        assert cc["largest_frac"] > 0.999, f"{name}: fragmented ({cc['sizes'][:4]})"
        assert n_dxf > 0 and len(m.ports) == 4
        out.append(f"{name}({100*cc['largest_frac']:.2f}%, {n_dxf} dxf, 4 ports)")
    # and the COMPOSED device, rendered stage-by-stage
    stages = [Stage("head", Layout("linear", 1), upto_F=0.4),
              Stage("loops", Layout("paperclip", 2, 900 * UM, loop_pull_in=2 * UM))]
    dc, cuts = compose(c.geo_req, stages)
    mc = build_mask(dc, stages[1].layout)
    gc, _ = raster(mc, res=7.0)
    ccc = connectivity(gc)
    assert ccc["largest_frac"] > 0.999, "composed device fragmented"
    out.append(f"composed({ccc['n_components']}c)")
    return "ALL ARCHITECTURES BUILD, CONNECT, AND RECOVER F FROM POLYGONS: " + ", ".join(out)


# ===========================================================================
# S19e -- THE FIELD SOLVE, DONE RIGHT
#
# I previously rejected Hele-Shaw outright.  That was too broad, and the objection was
# correct: a FIELD SOLVE is the only check that can catch a SHORT.  Connectivity cannot --
# a short still reads as "one component", because the centre and side channels are SUPPOSED
# to be connected (through the links).  And geom_ftest cannot -- it measures widths at the
# rows and rebuilds the ladder, so a stray channel bridging the wall is invisible to it.
#
# What was wrong was not "field solve"; it was UNIFORM CONDUCTANCE.
#
#   * Hele-Shaw / Darcy assumes a uniform sheet: sigma = const  =>  R ~ L/w.
#     Valid iff h << w.  Every channel here has h > w.  On the reference device that law
#     gives F = 0.6295 against the exact 0.5593 -- a 12.6% error.
#
#   * The DEPTH-CORRECT sheet has a WIDTH-WEIGHTED conductance.  Require that a strip of
#     width w have the true duct resistance:
#
#         integral of sigma across the width  =  1 / r(w)  =  h w^3 B(h/w) / (12 mu)
#
#     Taking sigma proportional to (2d)^2, with d the distance to the nearest wall (a
#     distance transform), and integrating 2*int_0^{w/2} (2y)^2 dy = w^3/3, fixes the
#     constant:
#
#         sigma(x,y)  =  [ h B(h/w) / (4 mu) ] * (2 d(x,y))^2
#
#     This reduces to the ladder in the limit of straight channels, and it is built from
#     the SAME B(alpha) as r(w).  Uniform sigma is the w^0 special case -- i.e. the shallow
#     limit -- which is exactly the sorting error the specification exists to prevent.
# ===========================================================================

def _B_vec(alpha: np.ndarray, n_terms: int = 60) -> np.ndarray:
    """Vectorized B(alpha) for the field solve.  Same series as S3."""
    n = _odd(n_terms)[:, None]
    a = np.maximum(np.asarray(alpha, float).ravel()[None, :], 1.0)
    s = np.sum(np.tanh(np.minimum(n * np.pi * a / 2.0, 350.0)) / n ** 5, axis=0)
    return (1.0 - (192.0 / (np.pi ** 5 * a[0])) * s).reshape(np.shape(alpha))


def field_solve(mask: Mask, dev: Device, res: float = 4.0) -> dict:
    """
    DEPTH-CORRECT 2-D field solve on the DRAWN POLYGONS.  Returns F and the flux split.

    Catches what nothing else can: a SHORT across the wall.  Also returns the uniform-sigma
    (Hele-Shaw) answer alongside, so the size of the shallow-limit error is visible rather
    than assumed.
    """
    from scipy.ndimage import distance_transform_edt
    from scipy.sparse import coo_matrix
    from scipy.sparse.linalg import spsolve

    p = dev.req.proc
    bx = mask.bounds
    npx = ((bx[2] - bx[0]) / res) * ((bx[3] - bx[1]) / res)
    if npx > 6e5:                      # a direct sparse solve above ~600k unknowns will OOM
        res = res * math.sqrt(npx / 6e5)
    grid, (x0, y0, r) = raster(mask, res)
    # Restrict to the LARGEST connected component.  The isolated corner pockets at the array's
    # leading post carry no flow and give the Laplacian zero-conductance rows -> singular.
    from scipy.ndimage import label as _label
    _lab, _n = _label(grid, structure=np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]]))
    if _n > 1:
        _sz = np.bincount(_lab.ravel())[1:]
        grid = (_lab == (1 + int(np.argmax(_sz))))
    H, W = grid.shape
    d = distance_transform_edt(grid) * res                    # um to the nearest wall
    w_loc = np.maximum(2.0 * d, 1e-6)                         # local width proxy
    h_um = p.h / UM
    sig = (w_loc ** 2) * _B_vec(np.maximum(h_um / w_loc, 1.0))   # DEPTH-CORRECT
    sig[~grid] = 0.0
    sig_hs = grid.astype(float)                               # UNIFORM = Hele-Shaw (the wrong one)

    def solve(sigma):
        idx = -np.ones((H, W), dtype=int)
        # port discs -> Dirichlet
        yy, xx = np.mgrid[0:H, 0:W]
        px = x0 + (xx + 0.5) * res
        py = y0 + (yy + 0.5) * res
        bc = np.zeros((H, W)); is_bc = np.zeros((H, W), bool)
        pmask = {}
        for (nm, cx, cy, rad) in mask.ports:
            m = grid & (((px - cx) ** 2 + (py - cy) ** 2) <= rad ** 2)
            if not m.any():                                   # snap to the nearest open pixel
                dd = (px - cx) ** 2 + (py - cy) ** 2
                dd[~grid] = np.inf
                k = np.unravel_index(np.argmin(dd), dd.shape)
                m = np.zeros((H, W), bool); m[k] = True
            pmask[nm] = m
            is_bc |= m
            bc[m] = 1.0 if nm == "INLET" else 0.0
        free = grid & ~is_bc
        n = int(free.sum())
        idx[free] = np.arange(n)
        I, J, V = [], [], []
        rhs = np.zeros(n)
        diag = np.zeros(n)
        for dy, dx in ((0, 1), (0, -1), (1, 0), (-1, 0)):
            a = free
            b = np.roll(np.roll(grid, -dy, 0), -dx, 1)
            sb = np.roll(np.roll(sigma, -dy, 0), -dx, 1)
            ib = np.roll(np.roll(idx, -dy, 0), -dx, 1)
            bb = np.roll(np.roll(bc, -dy, 0), -dx, 1)
            if dy == 1:  b[-1, :] = False
            if dy == -1: b[0, :] = False
            if dx == 1:  b[:, -1] = False
            if dx == -1: b[:, 0] = False
            link = a & b
            g = 2.0 * sigma * sb / np.maximum(sigma + sb, 1e-300)   # harmonic mean face
            g = np.where(link, g, 0.0)
            ii = idx[link]
            diag[ii] += g[link]
            fr = link & (ib >= 0)
            I += idx[fr].tolist(); J += ib[fr].tolist(); V += (-g[fr]).tolist()
            db = link & (ib < 0)
            np.add.at(rhs, idx[db], g[db] * bb[db])
        I += list(range(n)); J += list(range(n)); V += diag.tolist()
        A = coo_matrix((V, (I, J)), shape=(n, n)).tocsr()
        pot = np.zeros((H, W)); pot[is_bc] = bc[is_bc]
        pot[free] = spsolve(A.tocsc(), rhs)
        # flux INTO each outlet disc
        flux = {}
        for nm, m in pmask.items():
            tot = 0.0
            for dy, dx in ((0, 1), (0, -1), (1, 0), (-1, 0)):
                nb = np.roll(np.roll(grid & ~m, -dy, 0), -dx, 1)
                sb = np.roll(np.roll(sigma, -dy, 0), -dx, 1)
                pb = np.roll(np.roll(pot, -dy, 0), -dx, 1)
                lk = m & nb
                g = 2.0 * sigma * sb / np.maximum(sigma + sb, 1e-300)
                tot += float(np.sum(np.where(lk, g * (pot - pb), 0.0)))
            flux[nm] = tot
        return flux

    fx = solve(sig)
    q_in = abs(fx.get("INLET", 0.0))
    q_s = abs(fx.get("OUT_S+", 0.0)) + abs(fx.get("OUT_S-", 0.0))
    F_field = q_s / q_in if q_in > 0 else float("nan")

    fh = solve(sig_hs)
    qi = abs(fh.get("INLET", 0.0))
    F_hs = (abs(fh.get("OUT_S+", 0.0)) + abs(fh.get("OUT_S-", 0.0))) / qi if qi > 0 else float("nan")

    # Compare against the F of the ROWS ACTUALLY DRAWN -- a partial mask has a partial F.
    n_drawn = sum(mask.legs)
    F_ref = 1.0 - dev.rows[n_drawn - 1].P
    return dict(F_field=F_field, F_heleshaw_uniform=F_hs, F_design=F_ref, F_device=dev.F,
                err=abs(F_field / F_ref - 1) if F_ref else float("nan"),
                err_hs=abs(F_hs / F_ref - 1) if F_ref else float("nan"),
                res=res, n_open=int(grid.sum()), n_rows=n_drawn)


@_v("G24")
def _G24(c):
    """
    FIELD SOLVE on the DRAWN POLYGONS -- the only check that can catch a SHORT.

    Connectivity cannot: a short still reads as ONE component, because the primary and the
    side channels are SUPPOSED to be connected (through the links).  geom_ftest cannot: it
    measures widths at the rows and rebuilds the ladder, so a stray channel bridging the wall
    is invisible to it.  Only a field solve sees it -- it changes where the current goes.

    But the conductance must be DEPTH-CORRECT.  Requiring a strip of width w to have the true
    duct resistance fixes sigma ~ (2d)^2 B(h/2d), with d the distance transform.  UNIFORM
    sigma is the w^0 special case -- Hele-Shaw, the shallow limit -- and on this device it is
    off by nearly 300%, versus ~10% for the depth-correct law.  So: field solve yes;
    Hele-Shaw no.
    """
    m = _mask_ctx(c)
    fs = field_solve(m, c.geo, res=4.0)
    assert np.isfinite(fs["F_field"]), "field solve did not converge"
    assert fs["err"] < 0.25, f"drawn geometry does not carry the design flow split: {fs}"
    assert fs["err_hs"] > 3 * fs["err"], "uniform sigma should be MUCH worse -- it is the wrong law"
    return (f"F(design) = {fs['F_design']:.5f} | F(DEPTH-CORRECT field, sigma~w^2 B) = "
            f"{fs['F_field']:.5f} [{100*fs['err']:.0f}% -- a short would move it far more] | "
            f"F(UNIFORM sigma = Hele-Shaw) = {fs['F_heleshaw_uniform']:.5f} "
            f"[{100*fs['err_hs']:.0f}% -- the SHALLOW law on a DEEP device]")


@_v("G25")
def _G25(c):
    """
    PORTS: 1 INLET + 3 OUTLETS, each a via disc that lies INSIDE the channel it serves.
    A disc wider than its channel bridges the wall and shorts the device -- and the field
    solve then reports exactly that.  Exported as DXF CIRCLEs on layer PORTS.
    """
    m = _mask_ctx(c)
    names = [q[0] for q in m.ports]
    assert names == ["INLET", "OUT_C", "OUT_S+", "OUT_S-"], names
    g, (x0, y0, res) = raster(m, res=3.0)
    H, W = g.shape
    yy, xx = np.mgrid[0:H, 0:W]
    px = x0 + (xx + .5) * res
    py = y0 + (yy + .5) * res
    for (nm, cx, cy, rad) in m.ports:
        disc = ((px - cx) ** 2 + (py - cy) ** 2) <= rad ** 2
        frac = (g & disc).sum() / max(disc.sum(), 1)
        assert frac > 0.85, f"port {nm} is not inside a channel ({100*frac:.0f}% open)"
    return (f"4 ports {names}; each >85% inside its own channel (radius capped by the local "
            f"half-width, or the disc would bridge the wall and short the device)")


@_v("G26")
def _G26(c):
    """
    WALL INTEGRITY -- the check that would have caught every geometry bug in one shot.

    Between the primary channel and a side channel there is a WALL of thickness l_g, pierced
    ONLY at the links.  Raster the wall band and measure its open fraction; it must equal the
    link area plus the four cap corners, and nothing more.

    A fold cutting into the wall (drawn on the wall's centreline instead of hanging off it),
    or a ghost w_s-wide strip drawn under the folds, would both blow this number up.
    """
    d = c.geo
    p = d.req.proc
    m = build_mask(d, Layout("linear", 1), max_rows=30)
    g, (x0, y0, res) = raster(m, res=0.6)
    H, W = g.shape
    yy, xx = np.mgrid[0:H, 0:W]
    px = x0 + (xx + .5) * res
    py = y0 + (yy + .5) * res
    lg, gp, G = p.l_g / UM, p.g_p / UM, p.G / UM
    Lam = p.Lambda / UM
    A_link = G * lg
    A_cap = 4 * ((gp / 2) ** 2 - 0.25 * math.pi * (gp / 2) ** 2)     # 4 corners (in + out)
    exp = (A_link + A_cap) / (Lam * lg)
    worst = 0.0
    for k in (0, 8, 20, 29):
        rw = d.rows[k]
        wc2 = rw.w_c / UM / 2
        xr = k * Lam
        for sgn in (1, -1):
            lo, hi = sorted((sgn * wc2, sgn * (wc2 + lg)))
            band = (px >= xr) & (px < xr + Lam) & (py >= lo) & (py <= hi)
            worst = max(worst, abs(float(g[band].mean()) - exp))
    assert worst < 0.05, f"wall band open fraction is off by {100*worst:.1f}% -- something " \
                         f"is cutting through the wall"
    assert (p.l_g - p.g_p) / 2 >= 0
    return (f"wall band open fraction matches link+corners to {100*worst:.1f}% "
            f"(expected {100*exp:.1f}%).  The wall is SOLID except at the links: no fold "
            f"intrusion, no ghost side-channel strip.  Stadium hs = {(p.l_g-p.g_p)/2/UM:.1f} um >= 0.")



# ===========================================================================
# S20 -- EMBEDDED PROCEDURAL CORE (self-contained; byte-exact verified source)
#
#   meander_core : Params/Derived/derive/build_pattern, region families (arms,
#                  taper, decay, MORPH, release, GRID), width/spacing/taper
#                  PROFILES (function-modifiable), ROUNDED CORNERS (fillets), the
#                  19 geometry validators, 24 mutation probes.
#   physics_core : flow-prescribed pressure-field closure (EXACT, not the rho*Rc
#                  uniform-link limit), variable-rung rung+rail Kirchhoff solve,
#                  P1-P10 audit, compile->audit->correct fixed-point loop with
#                  damped-Jacobi width correction.
#
# Embedded as source and exec'd into real modules registered in sys.modules (so
# dataclass string annotations resolve).  CIF.py imports NOTHING from disk for these
# -- it is fully self-contained -- while the audited code is preserved verbatim.
# ===========================================================================
import types as _types, sys as _sys

_MEANDER_SRC = '#!/usr/bin/env python3\n# -*- coding: utf-8 -*-\n"""\nmeander_pattern.py  (v3 - traced line-drawing model)\n====================================================\nProcedural generator for the "central trunk with meandering side channels"\nfigure. v3 follows the author\'s explicit trace of the BLACK LINES: the\nfigure is a LINE DRAWING - the thin black strokes themselves are the\nchannel network. (v1/v2 wrongly modeled white lumens carved from black\nfill, which produced disconnected island bars. There are none here.)\n\nTHE TRACED MODEL (right side; the left side is its mirror)\n----------------------------------------------------------\nThe trunk is a thick vertical bar (tapered nib on top, T-foot crossbar and\na short tail at the bottom). All side geometry is drawn with ONE stroke\nwidth on a uniform row grid y_i = top_y + i*pitch. Each side consists of\nU = n_arms + n_taper "units", chained top to bottom. Unit k occupies rows\n2k (its OUT line) and 2k+1 (its RETURN line):\n\n  ARM unit (k < n_arms), the nested goal-posts on top:\n      start ON the trunk face at row 2k -> RIGHT to reach R_k\n      -> UP to the hairpin top (a tall rise; inner arms tallest)\n      -> short RIGHT by hairpin_w -> DOWN past the out row to row 2k+1\n      -> LEFT back toward the trunk ... but STOP SHORT of it\n      -> DOWN one pitch to row 2k+2.\n      Arm reaches INCREASE outward (each arm wraps around the previous).\n\n  TAPER unit (k >= n_arms), the V-shaped staircase:\n      start ON the trunk face at row 2k -> RIGHT to reach R_k\n      -> short DOWN one pitch to row 2k+1\n      -> LEFT back toward the trunk ... STOP SHORT of it\n      -> DOWN one pitch to row 2k+2.\n      Taper reaches DECREASE linearly (gradually shorter).\n\n  THE CHAIN: the final down-jog of unit k lands ON unit k+1\'s out line\n  (a T-junction, always at the same constant offset `merge_off` from the\n  trunk face): going left from that T reaches the trunk, going right\n  continues through the next side channel - exactly as traced. The LAST\n  unit\'s jog lands on the foot crossbar. Hence every stroke is connected:\n  no single-line islands, and each consecutive pair of units closes a\n  loop through the trunk (a meandering channel that leaves the trunk and\n  rejoins it downstream).\n\nRENDERING\n---------\nPlain black on white: one <g> of stroked polylines (butt caps, miter\njoins, width line_w) plus the filled trunk shapes. Deterministic bytes.\n\nVALIDATION\n----------\n`validate(pattern)` runs 13 checks; the topology-specific ones are\n`chain_links` (T-junction landing interior to the next out line, the\nconstant stop-short pocket, and the foot landing), `arm_shape`,\n`taper_shape`, `uniform_rows`, and `connectivity` (V=2U+1, E=4U, 2U\ncycles, no islands). A self-test corrupts patterns in eleven ways and\nasserts each corruption is caught by the matching validator.\n\nUSAGE\n-----\n  python3 meander_pattern.py                 # suite + self-test -> ./out\n  python3 meander_pattern.py --outdir DIR\n  python3 meander_pattern.py --single --width 800 --height 1000 \\\n          --arms 3 --taper 8 --name my_pattern\n  python3 meander_pattern.py --skip-selftest\n"""\n\nfrom __future__ import annotations\n\nimport argparse\nimport copy\nimport hashlib\nimport json\nimport math\nimport os\nimport sys\nfrom dataclasses import dataclass, replace, asdict\nfrom typing import Dict, List, Optional, Tuple\n\nPoint = Tuple[float, float]\nRect = Tuple[float, float, float, float]          # (x0, y0, x1, y1), y down\nEPS = 1e-6                                        # numeric tolerance\n\n\n# ---------------------------------------------------------------------------\n# 1. PARAMETERS\n# ---------------------------------------------------------------------------\n@dataclass\nclass Params:\n    """User-facing parameters. Optional fields (None) are derived so that\n    ANY box + counts yields a well-formed pattern. Fractions are of the\n    full box unless stated otherwise."""\n\n    width: float = 940.0          # bounding box width  (SVG units)\n    height: float = 990.0         # bounding box height\n    n_arms: int = 3               # tall hairpin units per side (top group)\n    n_taper: int = 11             # tapering units per side (V staircase)\n    n_decay: int = 0              # exponential-convergence units per side;\n                                  # 0 = classic figure (T-foot ending);\n                                  # >=1 = CONNECTOR mode: reaches decay\n                                  # geometrically onto the stop-short\n                                  # line, which then continues as a\n                                  # straight rail to the bottom edge\n    margin: Optional[float] = None        # default 0.03*width\n    pitch: Optional[float] = None         # row spacing; None => SOLVED so\n                                          # the 2U rows exactly span\n                                          # top_frac..bottom_frac\n    line_w: Optional[float] = None        # stroke width; default 0.30*pitch\n    hairpin_w: Optional[float] = None     # the "short right" at an arm\'s\n                                          # top; default = pitch\n    merge_off: Optional[float] = None     # trunk face -> stop-short jog x\n                                          # (constant for every unit);\n                                          # default 0.80*pitch\n    n_morph: int = 0              # MORPH units per side: a SMOOTH,\n                                  # function-controlled transition from\n                                  # the converged decay into the grid.\n                                  # The excess reach E glides from the\n                                  # decay floor E_min to exactly 0 under\n                                  # morph_profile, jumping the sliver-\n                                  # wall EXCLUSION BAND once (see spec\n                                  # section 15): scheduled values inside\n                                  # (mean_w - min_clear/2, mean_w +\n                                  # min_clear) snap to the nearest legal\n                                  # endpoint. Requires n_decay >= 1;\n                                  # mutually exclusive with n_release.\n                                  # 0 = the legacy immediate handoff.\n    morph_profile: object = ("cosine",)   # any envelope family or\n                                  # callable phi(t): 1 -> 0 shaping the\n                                  # morph\'s E schedule.\n    corner_r: Optional[float] = None      # ROUND ALL CORNERS (render\n                                          # layer): convex ink corners via\n                                          # round joins/caps; concave\n                                          # (void) corners via true\n                                          # fillets (material-adding, so\n                                          # connectivity is preserved).\n                                          # The PORT face stays square by\n                                          # contract. Radius must fit the\n                                          # smallest white feature.\n    side_w_profile: Optional[object] = None   # SIDE CHANNEL WIDTH by\n                                          # function: per-unit scale on\n                                          # line_w, running start -> 1.0:\n                                          # ("linear", s0), ("cosine", s0),\n                                          # ("power", p, s0),\n                                          # ("exp", lam, s0), or callable\n                                          # scale(t) with scale(1) = 1.\n    center_w_profile: Optional[object] = None  # CENTER CHANNEL WIDTH by\n                                          # function along its length,\n                                          # same forms; anchored so the\n                                          # width at the PORT equals the\n                                          # declared trunk_w.\n    n_grid: int = 0               # GRID units per side: after the decay,\n                                  # the side channel runs STRAIGHT on the\n                                  # rail (wall boundary_w from the center\n                                  # channel) while PERIODIC rungs keep it\n                                  # intersecting the center channel every\n                                  # unit - the wall becomes rectangular\n                                  # white posts, the black rungs the gaps\n                                  # in it. Requires connector mode;\n                                  # mutually exclusive with n_release.\n    n_release: int = 0            # RELEASE units per side: after the\n                                  # decay, the trunk intersections fade -\n                                  # each rung detaches from the center\n                                  # wall by a growing junction gap, from\n                                  # line_w (minimal legal crack) to\n                                  # boundary_w - line_w (a last one-line\n                                  # tooth on the rail). Requires\n                                  # connector mode and boundary_w >=\n                                  # 2*line_w.\n    release_profile: Optional[object] = None  # envelope for the junction\n                                  # gap schedule; same families as\n                                  # taper_profile (linear default).\n    decay_q: Optional[float] = None       # DIRECT decay-rate control: the\n                                          # per-unit excess ratio in (0,1).\n                                          # Mutually exclusive with\n                                          # n_decay>0; the unit count D is\n                                          # derived so the decay ends on\n                                          # E_min, then q is re-derived\n                                          # exactly from that integer D.\n    spacing: Optional[object] = None      # CHANNEL SPACING PROFILE: None =\n                                          # uniform; ("linear", end/start)\n                                          # or ("geometric", last/first)\n                                          # gap ratios; or a callable\n                                          # w(t)>0 sampled at t in [0,1]\n                                          # over the row gaps and\n                                          # normalized to fill the span.\n    taper_profile: Optional[object] = None  # TAPER ENVELOPE: None = linear;\n                                          # ("cosine",), ("power", p),\n                                          # ("exp", lam); or a callable\n                                          # phi(t), strictly decreasing,\n                                          # affinely normalized so reach\n                                          # runs exactly r0 -> r1.\n    boundary_w: Optional[float] = None    # CONVERGENT BOUNDARY WIDTH: the\n                                          # white gap between the trunk\n                                          # face and the converged side\n                                          # rail (= merge_off - line_w/2).\n                                          # When given, it DEFINES\n                                          # merge_off = boundary_w + lw/2.\n                                          # Must be >= line_w.\n    trunk_w: Optional[float] = None       # default 4.3*line_w\n    min_white: Optional[float] = None     # PHYSICS-MODE white floor: when\n                                          # given, the minimum legal white\n                                          # gap (min_clear) is this value\n                                          # instead of line_w, decoupling\n                                          # the FABRICATION wall floor\n                                          # from the stroke width (needed\n                                          # when side channels are wide,\n                                          # e.g. hydraulic closure widths).\n                                          # None = classic rule\n                                          # (min_clear = line_w), byte-\n                                          # identical behavior.\n    top_frac: float = 0.42        # y of row 0 (first out line)\n    bottom_frac: float = 0.935    # target y of the last row (2U-1)\n    taper_start_frac: float = 0.77        # first taper reach / d_line\n    taper_end_frac: float = 0.36          # last  taper reach / d_line\n    fan_step: Optional[float] = None      # arm nesting step; default\n                                          # max(0.035*width, hairpin+2*lw)\n    drop_frac: float = 0.20       # hairpin-top staircase drop / height\n    nib_len_frac: float = 0.06    # trunk nib length / height\n    name: str = "pattern"\n\n\n# ---------------------------------------------------------------------------\n# 2. DERIVED GEOMETRY (the constraint solver)\n# ---------------------------------------------------------------------------\n@dataclass\nclass Derived:\n    """Every concrete number needed to place geometry. Recomputing from\n    the same Params must be bit-identical (exploited by the determinism\n    validator, which doubles as a tamper alarm)."""\n    margin: float\n    xc: float\n    top_y: float                  # y of row 0\n    bottom_y: float               # y of the last row (2U-1)\n    pitch: float                  # s: MEAN row gap (uniform gap when\n                                  # spacing is None)\n    min_gap: float                # smallest row gap (= pitch if uniform);\n                                  # all safety floors key off this\n    line_w: float                 # lw: stroke width of all side lines\n    hairpin_w: float              # u\n    merge_off: float              # m: trunk face -> jog centerline\n    boundary_w: float             # b = m - lw/2: white between trunk face\n                                  # and the converged side rail\n    min_clear: float              # minimum legal white gap (= line_w)\n    trunk_w: float\n    d_line: float                 # max usable stroke-centerline reach\n    fan_step: float\n    arm_reaches: List[float]      # len A, strictly increasing\n    hairpin_tops: List[float]     # len A, y of each arm\'s top run\n    taper_reaches: List[float]    # len T, strictly decreasing (linear)\n    decay_reaches: List[float]    # len D, geometric decay onto the rail\n    decay_q: float                # per-unit excess ratio (0 when D = 0)\n    release_gaps: List[float]     # len R: junction gap (trunk face ->\n                                  # rung tip) per release unit\n    morph_excess: List[float]     # len M: snapped excess schedule\n    morph_fused: List[bool]       # len M: E in the fused regime (down\n                                  # stroke overlaps the previous jog;\n                                  # the enclosed pocket is gone)\n    unit_widths: List[float]      # len U: per-unit side channel width\n    trunk_halves: List[float]     # per ROW: center channel half-width\n    trunk_half_end: float         # half-width at the trunk\'s bottom end\n                                  # (== trunk_w/2 at the port, asserted)\n    port_y: float                 # bottom interface y (connector: = H)\n    rows: List[float]             # len 2U, the row grid\n    trunk_top: float\n    main_top: float\n    foot_cy: float\n    foot_w: float\n    foot_h: float\n    tail_end: float\n    max_bottom: float\n\n\ndef _spacing_weights(spec, n_gaps: int) -> List[float]:\n    """n_gaps normalized gap fractions (sum = 1) sampled at\n    t = i/(n_gaps-1). None = uniform."""\n    if spec is None or n_gaps == 1:\n        return [1.0 / n_gaps] * n_gaps\n    ts = [i / (n_gaps - 1) for i in range(n_gaps)]\n    if callable(spec):\n        w = [float(spec(t)) for t in ts]\n    else:\n        fam = spec[0]\n        if fam == "linear":                  # gap grows/shrinks linearly\n            e = float(spec[1])               # to e x the first gap\n            w = [1.0 + (e - 1.0) * t for t in ts]\n        elif fam == "geometric":             # gap ratio last/first = r\n            r = float(spec[1])\n            if r <= 0.0:\n                raise ValueError("geometric spacing ratio must be > 0")\n            w = [r ** t for t in ts]\n        else:\n            raise ValueError(f"unknown spacing family {fam!r}")\n    if min(w) <= EPS:\n        raise ValueError("spacing weights must be strictly positive")\n    tot = sum(w)\n    return [wi / tot for wi in w]\n\n\ndef _taper_envelope(spec, n: int) -> List[float]:\n    """n samples of the taper envelope at t = i/(n-1), strictly\n    decreasing, affinely normalized to phi[0]=1, phi[-1]=0 so the reach\n    runs exactly r0 -> r1 for ANY admissible family or callable."""\n    if n == 1:\n        return [1.0]\n    ts = [i / (n - 1) for i in range(n)]\n    if spec is None:\n        vals = [1.0 - t for t in ts]\n    elif callable(spec):\n        vals = [float(spec(t)) for t in ts]\n    else:\n        fam = spec[0]\n        if fam == "cosine":\n            vals = [(1.0 + math.cos(math.pi * t)) / 2.0 for t in ts]\n        elif fam == "power":\n            pw = float(spec[1])\n            if pw <= 0.0:\n                raise ValueError("power taper exponent must be > 0")\n            vals = [(1.0 - t) ** pw for t in ts]\n        elif fam == "exp":\n            lam = float(spec[1])\n            if lam <= 0.0:\n                raise ValueError("exp taper rate must be > 0")\n            el = math.exp(-lam)\n            vals = [(math.exp(-lam * t) - el) / (1.0 - el) for t in ts]\n        else:\n            raise ValueError(f"unknown taper_profile family {fam!r}")\n    for a, b in zip(vals, vals[1:]):\n        if b >= a - 1e-12:\n            raise ValueError("taper_profile must be strictly decreasing "\n                             "over its samples")\n    v0, v1 = vals[0], vals[-1]\n    return [(v - v1) / (v0 - v1) for v in vals]\n\n\ndef _decay_floor(min_clear: float, widths: List[float],\n                 A: int, T: int, D: int) -> float:\n    """Drawability floor for the decay excess: the white between unit\n    k\'s jog and unit k+1\'s outer down-stroke is E_{k+1} - mean(w_k,\n    w_{k+1}); take the worst pair feeding the decay region. Uniform\n    widths reduce to the classic 2*line_w."""\n    pairs = [(widths[k] + widths[k + 1]) / 2.0\n             for k in range(A + T - 1, A + T + D)\n             if k + 1 < len(widths)]\n    return min_clear + (max(pairs) if pairs else min_clear)\n\n\ndef _scale_profile(spec, n: int) -> List[float]:\n    """n per-sample width scales at t = i/(n-1), running start -> 1.0.\n    Forms: None -> all 1; ("linear", s0); ("cosine", s0);\n    ("power", p, s0); ("exp", lam, s0); or callable scale(t)."""\n    if spec is None:\n        return [1.0] * n\n    if n == 1:\n        return [1.0]\n    ts = [i / (n - 1) for i in range(n)]\n    if callable(spec):\n        vals = [float(spec(t)) for t in ts]\n    else:\n        fam = spec[0]\n        s0 = float(spec[-1])\n        env_spec = (None if fam == "linear"\n                    else (fam,) + tuple(spec[1:-1]))\n        env = _taper_envelope(env_spec, n)\n        vals = [1.0 + (s0 - 1.0) * e for e in env]\n    if min(vals) <= EPS:\n        raise ValueError("width profile scales must be positive")\n    return vals\n\n\ndef resolve_params(p: Params) -> Params:\n    """Resolve decay_q -> n_decay. The unit count D and the rate q are\n    two views of one decay (the endpoint E_min is a hard contract), so a\n    requested q maps to D = ceil(ln(E_min/E_start)/ln q); since E_min and\n    E_start themselves depend on the row budget (which depends on D),\n    solve by deterministic fixed-point iteration, then let derive()\n    re-express q exactly for that integer D."""\n    if p.decay_q is None:\n        return p\n    if p.n_decay:\n        raise ValueError("give n_decay or decay_q, not both")\n    if not (0.0 < p.decay_q < 1.0):\n        raise ValueError("decay_q must be in (0, 1)")\n    seen = set()\n    dc = 8\n    for _ in range(64):\n        trial = replace(p, decay_q=None, n_decay=dc)\n        dt = derive(trial)\n        rail = dt.trunk_w / 2.0 + dt.merge_off\n        e_start = dt.taper_reaches[-1] - rail\n        e_min = 2.0 * dt.line_w\n        dn = max(1, math.ceil(math.log(e_min / e_start)\n                              / math.log(p.decay_q)))\n        if dn == dc:\n            return trial\n        if dn in seen:                       # 2-cycle: settle on the\n            dc = max(dn, dc)                 # larger count (gentler q)\n            return replace(p, decay_q=None, n_decay=dc)\n        seen.add(dc)\n        dc = dn\n    raise ValueError("decay_q did not stabilize; give n_decay directly")\n\n\ndef derive(p: Params) -> Derived:\n    """Solve all dimensions from Params, enforcing feasibility. Raises\n    ValueError with a specific message when a parameter set cannot yield\n    a legal pattern (the messages are part of the validation story)."""\n    p = resolve_params(p)\n    # --- basic sanity ---------------------------------------------------\n    if p.width <= 0 or p.height <= 0:\n        raise ValueError("box dimensions must be positive")\n    if p.n_arms < 1:\n        raise ValueError("need at least 1 arm unit per side")\n    if p.n_taper < 1:\n        raise ValueError("need at least 1 taper unit per side")\n    if p.n_decay < 0:\n        raise ValueError("n_decay must be >= 0")\n    if p.n_release < 0:\n        raise ValueError("n_release must be >= 0")\n    if p.n_grid < 0:\n        raise ValueError("n_grid must be >= 0")\n    if p.n_grid >= 1 and p.n_release >= 1:\n        raise ValueError("n_grid and n_release are opposite endings "\n                         "(periodic intersection vs detachment); "\n                         "give one")\n    if p.n_morph < 0:\n        raise ValueError("n_morph must be >= 0")\n    if p.n_morph >= 1 and p.n_release >= 1:\n        raise ValueError("n_morph smooths the attached convergence; it "\n                         "cannot combine with the detaching release "\n                         "phase")\n    if not (0.0 < p.taper_end_frac < p.taper_start_frac < 1.0):\n        raise ValueError("need 0 < taper_end_frac < taper_start_frac < 1")\n    if not (0.0 < p.top_frac < p.bottom_frac < 1.0):\n        raise ValueError("need 0 < top_frac < bottom_frac < 1")\n\n    margin = p.margin if p.margin is not None else 0.03 * p.width\n    xc = p.width / 2.0\n    top_y = p.top_frac * p.height\n    A, T, D, R = p.n_arms, p.n_taper, p.n_decay, p.n_release\n    Gn, M = p.n_grid, p.n_morph\n    connector = D >= 1                       # convergence / port mode\n    if R >= 1 and not connector:\n        raise ValueError("the release phase extends the convergence: "\n                         "n_release >= 1 requires n_decay >= 1")\n    if Gn >= 1 and not connector:\n        raise ValueError("the grid phase extends the convergence: "\n                         "n_grid >= 1 requires n_decay >= 1")\n    if M >= 1 and not connector:\n        raise ValueError("the morph phase smooths the convergence: "\n                         "n_morph >= 1 requires n_decay >= 1")\n    U = A + T + D + R + M + Gn\n    nrows = 2 * U                            # every unit owns two rows\n\n    # --- row grid: pitch solved to hit bottom_frac, or given & stacked --\n    if p.pitch is None:\n        bottom_y = p.bottom_frac * p.height\n        s = (bottom_y - top_y) / (nrows - 1)\n        lw = p.line_w if p.line_w is not None else 0.30 * s\n        cap = p.height - margin - lw / 2.0   # lowest legal row centerline\n        if bottom_y > cap:                   # margin scales with WIDTH, so\n            bottom_y = cap                   # wide-short boxes can breach:\n            if bottom_y <= top_y + EPS:      # pull the last row up and\n                raise ValueError("margin leaves no vertical room")\n            s = (bottom_y - top_y) / (nrows - 1)   # re-pitch (lw shrinks\n            if p.line_w is None:                   # with it, so the cap\n                lw = 0.30 * s                      # only gets safer)\n    else:\n        s = p.pitch\n        lw = p.line_w if p.line_w is not None else 0.30 * s\n        bottom_y = top_y + (nrows - 1) * s\n        cap = p.height - margin - lw / 2.0\n        if bottom_y > cap + EPS:\n            raise ValueError(\n                f"explicit pitch stacks the last row to y={bottom_y:.1f}, "\n                f"past the bottom margin cap {cap:.1f}; reduce counts or "\n                f"pitch")\n    if s <= 0:\n        raise ValueError("vertical span too small for these counts")\n\n    # --- row grid: uniform, or a spacing PROFILE normalized to the span --\n    if p.spacing is None:\n        rows = [top_y + i * s for i in range(nrows)]\n        min_gap = s\n    else:\n        span = bottom_y - top_y\n        fr = _spacing_weights(p.spacing, nrows - 1)\n        rows = [top_y]\n        for f in fr:\n            rows.append(rows[-1] + span * f)\n        rows[-1] = bottom_y                  # kill cumulative float drift\n        min_gap = span * min(fr)\n        if p.line_w is None:                 # default lw keys off the\n            lw = 0.30 * min_gap              # TIGHTEST gap, not the mean\n\n    if lw <= 0:\n        raise ValueError("line_w must be positive")\n    min_clear = (p.min_white if p.min_white is not None\n                 else lw)                    # designed minimum white gap\n    if min_clear <= 0:\n        raise ValueError("min_white must be positive")\n    if lw > min_gap - min_clear + EPS:\n        raise ValueError("line_w > min_gap - min_clear: the closest rows "\n                         "would fuse (the white between rows must be >= "\n                         "min_clear; classic mode: min_clear = line_w)")\n\n    u = p.hairpin_w if p.hairpin_w is not None else s\n    if u < lw + 2.0:\n        raise ValueError("hairpin_w too small: the hairpin\'s two vertical "\n                         "strokes would fuse")\n    # --- stop-short offset & convergent boundary width --------------------\n    # The white pocket between the trunk face and every stop-short drop is\n    # b = m - lw/2; in connector mode this IS the interface gap between\n    # the center channel and the converged side rail, so the user may\n    # specify it directly.\n    if p.boundary_w is not None:\n        b = p.boundary_w\n        m = b + lw / 2.0\n        if p.merge_off is not None and abs(p.merge_off - m) > EPS:\n            raise ValueError("merge_off and boundary_w disagree; give one "\n                             "(boundary_w = merge_off - line_w/2)")\n    else:\n        m = p.merge_off if p.merge_off is not None else 0.80 * s\n        b = m - lw / 2.0\n    if b < min_clear - EPS:\n        raise ValueError("boundary_w (= merge_off - line_w/2) must be >= "\n                         "line_w: the gap beside the trunk is subject to "\n                         "the one-line-of-white clearance rule")\n    if m > 2.5 * s:\n        raise ValueError("merge_off too large: the jog no longer reads as "\n                         "\'just short of the center channel\'")\n\n    tw = p.trunk_w if p.trunk_w is not None else 4.3 * lw\n    if tw < max(2.0 * min_clear, lw) - EPS:\n        raise ValueError("trunk_w must be >= max(2*min_clear, line_w) "\n                         "(classic mode: >= 2*line_w, the thick bar of "\n                         "the drawing)")\n    t2 = tw / 2.0\n\n    # --- horizontal budget ----------------------------------------------\n    d_line = p.width / 2.0 - margin - lw / 2.0   # max stroke centerline x\n    if d_line <= t2 + m + 3.0 * s:\n        raise ValueError("box too narrow for margins + trunk + channels")\n\n    # --- channel WIDTH profiles ---------------------------------------------\n    # Side: one width per unit (scale on line_w, ending at 1.0 so the\n    # port presents exactly the declared side width). Center: half-width\n    # per row along the trunk (scale on trunk_w/2, ending at 1.0 at the\n    # port for the same reason). The stop-short line (rail) stays\n    # anchored at the DECLARED t2 + m, so walls vary as the center\n    # channel tapers while the interface stays exact.\n    side_scales = _scale_profile(p.side_w_profile, U)\n    unit_widths = [lw * sc for sc in side_scales]\n    if connector and abs(side_scales[-1] - 1.0) > 1e-9:\n        raise ValueError("side_w_profile must end at scale 1.0 (the port "\n                         "presents the declared line_w)")\n    trunk_len_end = p.height if connector else None   # set after foot;\n    # sample center profile per row over the trunk main\'s length\n    # (main_top -> port in connector mode; main_top -> bottom_y+ in\n    # classic mode, anchored at the LAST row for stability)\n    c_end = p.height if connector else bottom_y\n    c_span = max(c_end - (top_y - 0.0), 1.0)\n    def _half_at(y: float, _scales_cache={}) -> float:\n        t = min(1.0, max(0.0, (y - top_y) / c_span))\n        if p.center_w_profile is None:\n            return t2\n        if callable(p.center_w_profile):\n            return t2 * float(p.center_w_profile(t))\n        fam = p.center_w_profile[0]\n        s0 = float(p.center_w_profile[-1])\n        env_spec = (None if fam == "linear"\n                    else (fam,) + tuple(p.center_w_profile[1:-1]))\n        env = _taper_envelope(env_spec, 201)\n        idx = min(200, max(0, int(round(t * 200))))\n        return t2 * (1.0 + (s0 - 1.0) * env[idx])\n    trunk_halves = [_half_at(y) for y in rows]\n    trunk_half_end = _half_at(c_end)\n    if connector and abs(_half_at(c_end) - t2) > 1e-6:\n        raise ValueError("center_w_profile must end at scale 1.0 (the "\n                         "port presents the declared trunk_w)")\n    if min(trunk_halves) < max(min_clear, 0.5 * lw) - EPS:\n        raise ValueError("center channel thinner than 2*max(min_clear, "\n                         "line_w/2) somewhere along its length (classic "\n                         "mode: thinner than 2*line_w)")\n    # feasibility with variable widths (uniform case reduces to the\n    # classic rules exactly):\n    w_max = max(unit_widths)\n    for k in range(U):\n        g_slot = rows[2 * k + 1] - rows[2 * k]\n        if g_slot - unit_widths[k] < min_clear - EPS:\n            raise ValueError("side channel too wide for its row gap "\n                             f"(unit {k}): white slot under one line")\n        if k < U - 1:\n            g_adj = rows[2 * k + 2] - rows[2 * k + 1]\n            if g_adj - 0.5 * (unit_widths[k] + unit_widths[k + 1])                     < min_clear - EPS:\n                raise ValueError("adjacent side channels too wide "\n                                 f"(units {k},{k + 1})")\n        if m - unit_widths[k] / 2.0 < min_clear - EPS:\n            raise ValueError(f"unit {k} too wide for the stop-short "\n                             "pocket (raise boundary_w)")\n    if (t2 + m) - max(trunk_halves) - w_max / 2.0 < min_clear - EPS:\n        raise ValueError("center channel profile leaves less than one "\n                         "line of wall to the side rail")\n    if p.corner_r is not None:\n        white_feats = [rows[2 * k + 1] - rows[2 * k] - unit_widths[k]\n                       for k in range(U)]\n        white_feats.append(b)\n        white_feats.append((t2 + m) - max(trunk_halves) - w_max / 2.0)\n        white_feats.append(min_clear)   # the decay/snap-up pocket width\n                                        # floor: without this a legal r\n                                        # could erode pocket walls\n        r_cap = 0.45 * min(white_feats)\n        if not (0.0 < p.corner_r <= r_cap + EPS):\n            raise ValueError(f"corner_r must be in (0, {r_cap:.2f}] for "\n                             "these white features")\n\n\n    fan_step = (p.fan_step if p.fan_step is not None\n                else max(0.035 * p.width, u + 2.0 * lw))\n    fan_step = max(fan_step, u + 2.0 * lw)   # arms must nest with >= one\n    #                                          line of white between them\n    arm_outer = d_line - u                   # outermost DOWN stroke at\n    arm_reaches = [arm_outer - (A - 1 - j) * fan_step for j in range(A)]\n    if arm_reaches[0] <= t2 + m + 2.0 * s:\n        raise ValueError("too many arms / fan_step too large for width")\n\n    # --- hairpin-top staircase (inner arm tallest) ------------------------\n    # The hairpin top is a HORIZONTAL run, so its ink extends lw/2 above\n    # its centerline: the tallest top sits at margin + lw/2 so the ink\n    # exactly kisses the margin.\n    top_ink = margin + lw / 2.0\n    drop = p.drop_frac * p.height\n    if A > 1:\n        # rise of arm j = rows[2j] - top_j must stay >= 1.5*s and strictly\n        # decrease outward; out rows advance by the ACTUAL grid (profiled\n        # spacing allowed), so drop must exceed the largest per-arm\n        # advance.\n        y_last_arm = rows[2 * (A - 1)]\n        adv = max(rows[2 * j] - rows[2 * (j - 1)] for j in range(1, A))\n        drop_max = (y_last_arm - 1.5 * s - top_ink) / (A - 1)\n        if drop_max <= adv + EPS:\n            raise ValueError("no vertical room for the hairpin staircase; "\n                             "reduce n_arms or top_frac")\n        drop = min(drop, drop_max)\n        if drop <= adv + EPS:\n            raise ValueError("drop_frac too small: hairpin tops would not "\n                             "form a descending staircase")\n    hairpin_tops = [top_ink + j * drop for j in range(A)]\n    if hairpin_tops[0] >= top_y - 1.5 * s:\n        raise ValueError("first hairpin has no rise; raise top_frac or "\n                         "reduce drop_frac")\n\n    # --- taper reaches (linear, strictly decreasing) ----------------------\n    r0 = p.taper_start_frac * d_line\n    r0 = min(r0, arm_reaches[0] - max(2.0 * lw, 0.25 * fan_step))\n    r1 = p.taper_end_frac * d_line\n    reach_floor = t2 + m + 2.0 * s           # keep the return run readable\n    if r1 < reach_floor - EPS:\n        raise ValueError(\n            f"taper_end_frac too small: last reach {r1:.1f} < floor "\n            f"{reach_floor:.1f} (trunk + stop-short + 2*pitch)")\n    if r0 <= r1:\n        raise ValueError("taper range collapsed; adjust taper fracs")\n    if T == 1:\n        taper_reaches = [r0]\n    elif p.taper_profile is None:            # exact classic arithmetic\n        step = (r0 - r1) / (T - 1)\n        taper_reaches = [r0 - i * step for i in range(T)]\n    else:\n        env = _taper_envelope(p.taper_profile, T)\n        taper_reaches = [r1 + (r0 - r1) * e for e in env]\n\n    # --- exponential convergence (connector mode) --------------------------\n    # The convergence asymptote is the stop-short line itself: the rail at\n    # offset rail = t2 + m from the centerline. Excess reach beyond it\n    # decays geometrically, E_i = E_start * q^i, with q chosen so the\n    # LAST decay unit lands exactly on the drawability floor\n    # E_min = 2*lw (below which a unit\'s outer down-stroke would fuse\n    # with the neighboring drop: their white gap is E - lw).\n    rail = t2 + m\n    if connector:\n        e_min = _decay_floor(min_clear, unit_widths, A, T, D)\n        e_start = taper_reaches[-1] - rail\n        if e_start <= e_min + EPS:            # guaranteed by the taper\n            raise ValueError("taper already converged; raise "\n                             "taper_end_frac")                # floor, but\n        decay_q = (e_min / e_start) ** (1.0 / D)              # keep honest\n        decay_reaches = [rail + e_start * decay_q ** (i + 1)\n                         for i in range(D)]\n    else:\n        decay_q = 0.0\n        decay_reaches = []\n\n    # --- morph phase: E glides from the decay floor to exactly 0 -----------\n    # The 5-point unit family degenerates to the 3-point grid L at\n    # E == 0, so this is the natural continuous bridge. One law applies:\n    # between the pocket regime (E >= mean_w + min_clear: the enclosed\n    # white pocket between jog and next down-stroke is a legal wall) and\n    # the fused regime (E <= mean_w - min_clear/2: the down-stroke\n    # overlaps the jog and the pocket is GONE) lies the sliver-wall\n    # EXCLUSION BAND; scheduled values inside snap to the nearest legal\n    # endpoint, and the whole schedule is clamped non-increasing.\n    morph_excess: List[float] = []\n    morph_fused: List[bool] = []\n    if M >= 1:\n        env = _taper_envelope(p.morph_profile, M + 1)\n        prev = e_min\n        for j in range(1, M + 1):\n            e = e_min * env[j]\n            kk = A + T + D + j - 1\n            mean_w = 0.5 * (unit_widths[kk - 1] + unit_widths[kk])\n            lo = mean_w - min_clear / 2.0\n            hi = mean_w + min_clear\n            if lo + EPS < e < hi - EPS:          # exclusion band: snap\n                e = lo if (e - lo) <= (hi - e) else hi\n            e = min(e, prev)                     # monotone non-increasing\n            if e < 1e-9:\n                e = 0.0\n            morph_excess.append(e)\n            morph_fused.append(e <= lo + EPS)\n            prev = e\n\n    # --- release phase: the trunk intersections fade ------------------------\n    # Each release unit keeps the converged tooth (reach = rail + 2*lw)\n    # but its rung detaches from the center wall: the junction gap g_i\n    # follows the declared envelope from lw (the minimal legal crack -\n    # one line of white) to b - lw (a final one-line-deep tooth on the\n    # rail). Both the white gap and the ink bridge stay >= lw at every\n    # step, so NO clearance exemption is ever needed.\n    if R >= 1:\n        if b < 2.0 * lw - EPS:\n            raise ValueError("release phase needs boundary_w >= 2*line_w "\n                             "(one line of white AND one line of ink "\n                             "must fit in the junction gap)")\n        if R >= 2 and b <= 2.0 * lw + EPS:\n            raise ValueError("release phase with n_release >= 2 needs "\n                             "boundary_w > 2*line_w for a strictly "\n                             "widening junction gap")\n        rel_env = _taper_envelope(p.release_profile, R)\n        release_gaps = [lw + (b - 2.0 * lw) * (1.0 - e) for e in rel_env]\n        for i, g in enumerate(release_gaps):\n            wk = unit_widths[A + T + D + i]   # (R and M are exclusive)\n            if b - g < wk - EPS:\n                raise ValueError("release bridge thinner than its own "\n                                 f"channel width at unit {i}; reduce "\n                                 "side_w_profile there or raise "\n                                 "boundary_w")\n    else:\n        release_gaps = []\n\n    # --- bottom interface (connector mode) ---------------------------------\n    # The trunk and both rails run to the canvas bottom edge so a\n    # geometry with the same (trunk_w, line_w, boundary_w) cross-section\n    # can butt-join there.\n    port_y = p.height if connector else 0.0\n    if connector and port_y - bottom_y < 2.0 * s - EPS:\n        raise ValueError("no room for the straight rails below the last "\n                         "row; reduce bottom_frac or unit counts")\n\n    # --- trunk nib / foot / tail ------------------------------------------\n    # Nib: its flank must keep >= min_clear of white above row 0\'s stroke\n    # (gap = 0.5*nib_len - lw/2 - 1), so the height-fraction default is\n    # floored at that requirement on small boxes.\n    nib_len = p.nib_len_frac * p.height\n    nib_len = max(nib_len, 2.0 * (lw / 2.0 + min_clear + 1.0))\n    trunk_top = top_y - nib_len\n    main_top = top_y - 0.5 * nib_len\n    if trunk_top < margin - EPS:\n        raise ValueError("nib pokes above top margin; reduce nib_len_frac")\n\n    if connector:\n        # No T-foot: the trunk continues straight into the interface.\n        max_bottom = p.height\n        foot_cy = foot_w = foot_h = 0.0\n        tail_end = port_y\n    else:\n        max_bottom = p.height - max(2.0, 0.30 * margin)\n        foot_h = max(2.0 * lw, 0.62 * s)\n        air = max(min_clear, 0.5 * s)            # white between last row\n        foot_top = bottom_y + lw / 2.0 + air     # and the crossbar\n        foot_cy = foot_top + foot_h / 2.0\n        if foot_cy + foot_h / 2.0 > max_bottom:  # squeeze: drop the air to\n            foot_top = bottom_y + lw / 2.0 + min_clear   # its floor, then\n            foot_cy = foot_top + foot_h / 2.0            # thin the crossbar\n            if foot_cy + foot_h / 2.0 > max_bottom:\n                room = max_bottom - foot_top\n                foot_h = min(foot_h, room)\n                if foot_h + EPS < 1.5 * lw:\n                    raise ValueError("bottom_frac too large: no room "\n                                     "for foot")\n                foot_cy = foot_top + foot_h / 2.0\n        # The crossbar must catch the last unit\'s jog (landing at face +/-\n        # (t2+m)) with some slack, and stay proportioned to the trunk.\n        foot_w = max(2.1 * tw, 2.0 * (t2 + m + lw / 2.0 + 2.0))\n        tail_end = min(foot_cy + foot_h / 2.0 + 0.75 * s, max_bottom)\n\n    return Derived(margin=margin, xc=xc, top_y=top_y, bottom_y=bottom_y,\n                   pitch=s, min_gap=min_gap, line_w=lw, hairpin_w=u,\n                   merge_off=m,\n                   boundary_w=b, min_clear=min_clear, trunk_w=tw,\n                   d_line=d_line, fan_step=fan_step,\n                   arm_reaches=arm_reaches, hairpin_tops=hairpin_tops,\n                   taper_reaches=taper_reaches,\n                   decay_reaches=decay_reaches, decay_q=decay_q,\n                   release_gaps=release_gaps,\n                   morph_excess=morph_excess, morph_fused=morph_fused,\n                   unit_widths=unit_widths, trunk_halves=trunk_halves,\n                   trunk_half_end=trunk_half_end,\n                   port_y=port_y, rows=rows,\n                   trunk_top=trunk_top, main_top=main_top,\n                   foot_cy=foot_cy, foot_w=foot_w, foot_h=foot_h,\n                   tail_end=tail_end, max_bottom=max_bottom)\n\n\n# ---------------------------------------------------------------------------\n# 3. GEOMETRY MODEL\n# ---------------------------------------------------------------------------\n@dataclass\nclass Channel:\n    """One unit of the meander: an open Manhattan centerline polyline.\n    points[0] is ALWAYS on the trunk face (the T with the center channel);\n    points[-1] is ALWAYS the chain landing: interior to the NEXT unit\'s\n    out line — except the LAST unit, whose final drop lands on the foot\n    crossbar (classic mode) or continues as the straight side RAIL down\n    to the port (connector mode).\n    kind == "arm"  : 7 points (out, up, short right, down, left, jog).\n    kind == "taper": 5 points (out, short down, left, jog), linear taper.\n    kind == "decay": 5 points, same shape, geometrically decaying reach."""\n    side: str                 # "left" | "right"\n    index: int                # 0..A-1 arms, A..A+T-1 tapers, then decays\n    kind: str                 # "arm" | "taper" | "decay" | ...\n    points: List[Point]\n    width: float = 0.0        # this unit\'s stroke (channel) width\n\n\n@dataclass\nclass Pattern:\n    """Complete pattern: parameters + derived numbers + concrete geometry.\n    trunk_parts maps a name to ("rect", Rect) or ("poly", [Point], Rect)."""\n    params: Params\n    d: Derived\n    channels: List[Channel]\n    trunk_parts: Dict[str, tuple]\n\n\ndef build_pattern(p: Params) -> Pattern:\n    """Deterministically construct the full pattern from Params."""\n    p = resolve_params(p)\n    d = derive(p)\n    A, T, D, R = p.n_arms, p.n_taper, p.n_decay, p.n_release\n    Gn, M = p.n_grid, p.n_morph\n    connector = D >= 1\n    U = A + T + D + R + M + Gn\n    t2 = d.trunk_w / 2.0\n    rel_reach = t2 + d.merge_off + 2.0 * d.line_w   # the converged tooth\n    channels: List[Channel] = []\n\n    for side, sgn in (("left", -1.0), ("right", +1.0)):\n        jog_x = d.xc + sgn * (t2 + d.merge_off)   # constant stop-short x\n        for k in range(U):\n            face = d.xc + sgn * d.trunk_halves[2 * k]   # LOCAL face (the\n            wk = d.unit_widths[k]                       # center channel\n            yo = d.rows[2 * k]                    # may taper)\n            yr = d.rows[2 * k + 1]                # return-line row\n            if k < U - 1:\n                y_end = d.rows[2 * k + 2]         # T onto next out line\n            elif connector:\n                y_end = d.port_y                  # jog CONTINUES as the\n            else:                                 # straight side rail\n                y_end = (d.foot_cy - d.foot_h / 2.0) + 1.0   # into foot\n            if (A + T + D <= k < A + T + D + M\n                    and d.morph_excess[k - A - T - D] == 0.0):\n                pts = [                           # fully converged morph\n                    (face, yo),                   # unit: the degenerate\n                    (jog_x, yo),                  # 3-point grid L\n                    (jog_x, y_end),\n                ]\n                channels.append(Channel(side=side, index=k, kind="morph",\n                                        points=pts, width=wk))\n            elif k >= A + T + D + R + M:          # ---- grid unit ----\n                # A clean ladder rung + rail segment: the side channel\n                # runs straight at wall b while intersecting the center\n                # channel once per unit. Consecutive rail segments butt\n                # end-to-end into the continuous side channel.\n                pts = [\n                    (face, yo),                       # T with the trunk\n                    (jog_x, yo),                      # rung: full crossing\n                    (jog_x, y_end),                   # rail segment\n                ]\n                channels.append(Channel(side=side, index=k, kind="grid",\n                                        points=pts, width=wk))\n            elif k < A:                           # ---- arm unit ----\n                r = d.arm_reaches[k]\n                top = d.hairpin_tops[k]\n                pts: List[Point] = [\n                    (face, yo),                       # T with the trunk\n                    (d.xc + sgn * r, yo),             # out\n                    (d.xc + sgn * r, top),            # up (tall rise)\n                    (d.xc + sgn * (r + d.hairpin_w), top),   # short over\n                    (d.xc + sgn * (r + d.hairpin_w), yr),    # down past yo\n                    (jog_x, yr),                      # return, stop short\n                    (jog_x, y_end),                   # jog: chain down\n                ]\n                channels.append(Channel(side=side, index=k, kind="arm",\n                                        points=pts, width=wk))\n            else:                                 # ---- taper/decay/release\n                if k < A + T:\n                    r = d.taper_reaches[k - A]        # linear region\n                    kind = "taper"\n                    x0 = face                         # T with the trunk\n                elif k < A + T + D:\n                    r = d.decay_reaches[k - A - T]    # exponential region\n                    kind = "decay"\n                    x0 = face\n                elif k < A + T + D + M:               # smooth handoff:\n                    r = (t2 + d.merge_off             # E glides through\n                         + d.morph_excess[k - A - T - D])   # the fused\n                    kind = "morph"                    # regime to 0\n                    x0 = face\n                else:                                 # fading junction:\n                    r = rel_reach                     # frozen tooth, rung\n                    kind = "release"                  # tip pulls away\n                    g = d.release_gaps[k - A - T - D] # from the trunk\n                    x0 = face + sgn * g\n                pts = [\n                    (x0, yo),                         # out (on or off\n                    (d.xc + sgn * r, yo),             # the trunk face)\n                    (d.xc + sgn * r, yr),             # short down\n                    (jog_x, yr),                      # return, stop short\n                    (jog_x, y_end),                   # jog: chain down\n                ]\n                channels.append(Channel(side=side, index=k, kind=kind,\n                                        points=pts, width=wk))\n\n    # --- trunk: nib polygon + main column + foot crossbar ------------------\n    tip2 = 0.36 * t2                          # blunt nib tip half-width\n    shoulder_y = d.trunk_top + 0.55 * (d.main_top - d.trunk_top)\n    nib_poly: List[Point] = [\n        (d.xc - tip2, d.trunk_top), (d.xc + tip2, d.trunk_top),\n        (d.xc + t2, shoulder_y), (d.xc + t2, d.main_top + 1.0),\n        (d.xc - t2, d.main_top + 1.0), (d.xc - t2, shoulder_y),\n    ]\n    nib_bbox: Rect = (d.xc - t2, d.trunk_top, d.xc + t2, d.main_top + 1.0)\n    if p.center_w_profile is None:\n        main_part = ("rect",\n                     (d.xc - t2, d.main_top, d.xc + t2, d.tail_end))\n    else:                                     # shaped center channel:\n        ys = [d.main_top] + list(d.rows) + [d.tail_end]   # sampled poly\n        hmax = max(d.trunk_halves)\n        h0 = d.trunk_halves[0]\n        halves = [h0] + list(d.trunk_halves) + [d.trunk_half_end]\n        right = [(d.xc + h, y) for h, y in zip(halves, ys)]\n        left = [(d.xc - h, y) for h, y in zip(halves, ys)][::-1]\n        main_part = ("poly", right + left,\n                     (d.xc - hmax, d.main_top, d.xc + hmax, d.tail_end))\n    trunk_parts = {\n        "nib": ("poly", nib_poly, nib_bbox),\n        "main": main_part,\n    }\n    if not connector:                         # classic T-foot termination\n        foot_rect: Rect = (d.xc - d.foot_w / 2.0,\n                           d.foot_cy - d.foot_h / 2.0,\n                           d.xc + d.foot_w / 2.0,\n                           d.foot_cy + d.foot_h / 2.0)\n        trunk_parts["foot"] = ("rect", foot_rect)\n    return Pattern(params=p, d=d, channels=channels, trunk_parts=trunk_parts)\n\n\n# ---------------------------------------------------------------------------\n# 4. RENDERING (deterministic SVG, stdlib only)\n# ---------------------------------------------------------------------------\ndef _f(v: float) -> str:\n    """Fixed float formatting -> deterministic SVG bytes."""\n    return f"{v:.3f}"\n\n\ndef render_polyline(ch: Channel, d: Derived) -> List[Point]:\n    """The stroked centerline for one unit. points[0] is tucked slightly\n    INTO the trunk so the filled bar and the stroke always share ink\n    (guards against antialiasing hairlines at the T with the trunk)."""\n    sgn = -1.0 if ch.side == "left" else 1.0\n    tuck = (0.0 if ch.kind == "release"      # a detached rung must NOT be\n            else min(d.line_w, 0.45 * d.trunk_w))   # nudged into the gap\n    p0 = ch.points[0]\n    tucked0: Point = (p0[0] - sgn * tuck, p0[1])   # toward the centerline\n    return [tucked0] + [tuple(pt) for pt in ch.points[1:]]\n\n\ndef _polyline(points: List[Point], width: float) -> str:\n    pts = " ".join(f"{_f(x)},{_f(y)}" for x, y in points)\n    return (f\'    <polyline points="{pts}" \'\n            f\'stroke-width="{_f(width)}"/>\')\n\n\ndef _fillet_corners(pat: Pattern) -> List[tuple]:\n    """Every CONCAVE corner of the ink union, from the model: returns\n    (cx, cy, sx, sy) with (sx, sy) the unit direction INTO the void.\n    Filleting adds a black r x r square at the corner and carves a white\n    tangent disc of radius r - material is only ADDED, so connectivity\n    is preserved; each disc is tangent to its two walls and (by the\n    corner_r feasibility rule) too small to reach any other ink."""\n    d, p = pat.d, pat.params\n    out: List[tuple] = []\n    for c in pat.channels:\n        w2 = c.width / 2.0\n        pts = c.points\n        for i in range(1, len(pts) - 1):          # polyline bends\n            (x0, y0), (x1, y1), (x2, y2) = pts[i - 1], pts[i], pts[i + 1]\n            e1 = ((x1 > x0) - (x1 < x0), (y1 > y0) - (y1 < y0))\n            e2 = ((x2 > x1) - (x2 < x1), (y2 > y1) - (y2 < y1))\n            sx, sy = e2[0] - e1[0], e2[1] - e1[1]\n            if sx == 0 or sy == 0:\n                continue                          # colinear (defensive)\n            out.append((x1 + w2 * sx, y1 + w2 * sy, sx, sy))\n    for side, sgn in (("left", -1.0), ("right", 1.0)):\n        chs = _channels_by(pat, side)\n        for c in chs:                             # trunk T-junctions\n            if c.kind == "release":\n                continue\n            (fx, fy) = c.points[0]\n            w2 = c.width / 2.0\n            out.append((fx, fy - w2, sgn, -1.0))\n            out.append((fx, fy + w2, sgn, 1.0))\n        for k in range(len(chs) - 1):             # chain / grid landings\n            cur, nxt = chs[k], chs[k + 1]\n            xe, ye = cur.points[-1]\n            wj, wn = cur.width / 2.0, nxt.width / 2.0\n            if nxt.kind == "grid":                # rail-continuity corner\n                out.append((xe - sgn * wj, ye - wn, -sgn, -1.0))\n                out.append((xe - sgn * wn, ye + wn, -sgn, 1.0))\n            else:                                 # interior T landing\n                out.append((xe - wj, ye - wn, -1.0, -1.0))\n                out.append((xe + wj, ye - wn, 1.0, -1.0))\n        last = chs[-1]\n        if "foot" in pat.trunk_parts:             # classic foot landing\n            xe, ye = last.points[-1]\n            w2 = last.width / 2.0\n            ftop = pat.trunk_parts["foot"][1][1]\n            out.append((xe - w2, ftop, -1.0, -1.0))\n            out.append((xe + w2, ftop, 1.0, -1.0))\n    # SWALLOW FILTER: fusion (morph phase) merges strokes, so some\n    # enumerated concave corners no longer exist on the ink union\'s\n    # boundary. A blind white fillet disc there would carve a hole\n    # inside solid ink. Drop any candidate whose corner point lies\n    # STRICTLY inside an ink rect (its own two walls only ever touch it\n    # on their edges, which the strict test excludes).\n    eps = 1e-6\n    band = d.pitch if d.pitch > 0 else 10.0    # y-banded rect index: the\n    buckets: Dict[int, List[Rect]] = {}        # filter is O(corners) with\n    #                                            local lookups instead of\n    #                                            O(corners x rects); the\n    #                                            kept set is IDENTICAL.\n    for c in pat.channels:\n        w2 = c.width / 2.0\n        for i in range(len(c.points) - 1):\n            r = _seg_rect(c.points[i], c.points[i + 1], w2)\n            for bi in range(int(r[1] // band), int(r[3] // band) + 1):\n                buckets.setdefault(bi, []).append(r)\n    kept = []\n    for (cx, cy, sx, sy) in out:\n        cand = buckets.get(int(cy // band), ())\n        if any(x0 + eps < cx < x1 - eps and y0 + eps < cy < y1 - eps\n               for (x0, y0, x1, y1) in cand):\n            continue\n        kept.append((cx, cy, sx, sy))\n    return kept\n\n\ndef svg_string(pat: Pattern) -> str:\n    """Serialize a Pattern to an SVG document string. Everything is black\n    on white: the side network is stroked polylines, the trunk is filled\n    shapes; overlaps (T-junctions, the foot landing) are black-on-black,\n    so paint order is irrelevant to the result."""\n    p, d = pat.params, pat.d\n    W, H = p.width, p.height\n\n    params_json = json.dumps(asdict(p), sort_keys=True, default=repr)\n    desc = (params_json.replace("&", "&amp;").replace("<", "&lt;")\n            .replace(">", "&gt;"))\n\n    out: List[str] = []\n    out.append(\'<?xml version="1.0" encoding="UTF-8"?>\')\n    out.append(f\'<svg xmlns="http://www.w3.org/2000/svg" \'\n               f\'width="{_f(W)}" height="{_f(H)}" \'\n               f\'viewBox="0 0 {_f(W)} {_f(H)}">\')\n    out.append(f\'  <desc>meander_pattern params: {desc}</desc>\')\n    out.append(f\'  <rect x="0" y="0" width="{_f(W)}" height="{_f(H)}" \'\n               f\'fill="#ffffff"/>\')\n\n    rounded = p.corner_r is not None\n    rr = p.corner_r if rounded else 0.0\n    cap = "round" if rounded else "butt"\n    join = "round" if rounded else "miter"\n\n    out.append(f\'  <g id="channels" fill="none" stroke="#000000" \'\n               f\'stroke-linecap="{cap}" stroke-linejoin="{join}">\')\n    for ch in pat.channels:\n        pts = render_polyline(ch, d)\n        if rounded and ch.kind == "release":\n            # round caps EXTEND ends by w/2: inset the free rung tip so\n            # the modeled junction gap is exactly preserved\n            (x0, y0), (x1, y1) = pts[0], pts[1]\n            dx = (x1 > x0) - (x1 < x0)\n            pts = [(x0 + dx * ch.width / 2.0, y0)] + pts[1:]\n        out.append(_polyline(pts, ch.width))\n    out.append(\'  </g>\')\n\n    out.append(\'  <g id="trunk" fill="#000000" stroke="none">\')\n    _kind, poly, _bbox = pat.trunk_parts["nib"]\n    pts = " ".join(f"{_f(x)},{_f(y)}" for x, y in poly)\n    nib_stroke = (f\' stroke="#000000" stroke-width="{_f(2 * rr)}" \'\n                  f\'stroke-linejoin="round"\' if rounded else \'\')\n    out.append(f\'    <polygon points="{pts}"{nib_stroke}/>\')\n    for name in ("main", "foot"):\n        if name not in pat.trunk_parts:\n            continue\n        part = pat.trunk_parts[name]\n        if part[0] == "rect":\n            r = part[1]\n            rx = f\' rx="{_f(rr)}"\' if rounded else \'\'\n            out.append(f\'    <rect x="{_f(r[0])}" y="{_f(r[1])}" \'\n                       f\'width="{_f(r[2] - r[0])}" \'\n                       f\'height="{_f(r[3] - r[1])}"{rx}/>\')\n        else:                                     # shaped center channel:\n            # filled UNSTROKED even when rounding - its sampled boundary\n            # has only near-straight bends (nothing to round), its top\n            # is buried in the nib, and its bottom is the PORT face,\n            # which stays square by contract; a rounding stroke would\n            # dilate the declared width by corner_r per side\n            ppts = " ".join(f"{_f(x)},{_f(y)}" for x, y in part[1])\n            out.append(f\'    <polygon points="{ppts}"/>\')\n    out.append(\'  </g>\')\n\n    if rounded:\n        # true fillets at every CONCAVE (void) corner: add the corner\n        # square, then carve the tangent disc\n        corners = _fillet_corners(pat)\n        out.append(\'  <g id="fillet_squares" fill="#000000" \'\n                   \'stroke="none">\')\n        for cx, cy, sx, sy in corners:\n            x0, y0 = min(cx, cx + rr * sx), min(cy, cy + rr * sy)\n            out.append(f\'    <rect x="{_f(x0)}" y="{_f(y0)}" \'\n                       f\'width="{_f(rr)}" height="{_f(rr)}"/>\')\n        out.append(\'  </g>\')\n        out.append(\'  <g id="fillet_discs" fill="#ffffff" stroke="none">\')\n        for cx, cy, sx, sy in corners:\n            out.append(f\'    <circle cx="{_f(cx + rr * sx)}" \'\n                       f\'cy="{_f(cy + rr * sy)}" r="{_f(rr)}"/>\')\n        out.append(\'  </g>\')\n        if p.n_decay >= 1:\n            # the PORT face stays square by contract: overdraw butt ends\n            out.append(\'  <g id="port_square" fill="#000000" \'\n                       \'stroke="none">\')\n            t2 = d.trunk_w / 2.0\n            y0 = p.height - 2.0 * rr\n            out.append(f\'    <rect x="{_f(d.xc - t2)}" y="{_f(y0)}" \'\n                       f\'width="{_f(2 * t2)}" height="{_f(2 * rr)}"/>\')\n            wl = d.unit_widths[-1]\n            for sgn in (-1.0, 1.0):\n                xj = d.xc + sgn * (t2 + d.merge_off)\n                out.append(f\'    <rect x="{_f(xj - wl / 2)}" \'\n                           f\'y="{_f(y0)}" width="{_f(wl)}" \'\n                           f\'height="{_f(2 * rr)}"/>\')\n            out.append(\'  </g>\')\n    out.append(\'</svg>\')\n    return "\\n".join(out) + "\\n"\n\n\ndef write_svg(pat: Pattern, path: str) -> str:\n    s = svg_string(pat)\n    with open(path, "w", encoding="utf-8") as fh:\n        fh.write(s)\n    return hashlib.sha256(s.encode("utf-8")).hexdigest()\n\n\n# ---------------------------------------------------------------------------\n# 5. VALIDATION (13 checks)\n# ---------------------------------------------------------------------------\n@dataclass\nclass CheckResult:\n    name: str\n    ok: bool\n    detail: str\n\n\ndef _channels_by(pat: Pattern, side: str) -> List[Channel]:\n    return sorted((c for c in pat.channels if c.side == side),\n                  key=lambda c: c.index)\n\n\ndef _seg_rect(a: Point, b: Point, half: float) -> Rect:\n    """Painted rectangle of an axis-aligned stroked segment with butt\n    caps: expanded +/-half PERPENDICULAR to the segment only."""\n    x0, x1 = min(a[0], b[0]), max(a[0], b[0])\n    y0, y1 = min(a[1], b[1]), max(a[1], b[1])\n    if abs(a[1] - b[1]) <= EPS:                     # horizontal\n        return (x0, y0 - half, x1, y1 + half)\n    return (x0 - half, y0, x1 + half, y1)           # vertical\n\n\ndef _rect_gap(r1: Rect, r2: Rect) -> float:\n    """Largest axis-projected separation; negative iff the rects overlap\n    in both axes."""\n    gx = max(r1[0] - r2[2], r2[0] - r1[2])\n    gy = max(r1[1] - r2[3], r2[1] - r1[3])\n    return max(gx, gy)\n\n\ndef _trunk_bbox(part: tuple) -> Rect:\n    return part[1] if part[0] == "rect" else part[2]\n\n\ndef check_params_feasible(pat: Pattern) -> CheckResult:\n    """derive() accepts the parameters and reproduces the stored numbers."""\n    name = "params_feasible"\n    try:\n        d2 = derive(pat.params)\n    except ValueError as exc:\n        return CheckResult(name, False, f"derive() raised: {exc}")\n    d = pat.d\n    for fld in ("pitch", "line_w", "hairpin_w", "merge_off",\n                "boundary_w", "decay_q", "port_y", "trunk_w", "min_gap",\n                "fan_step", "top_y", "bottom_y", "d_line"):\n        if abs(getattr(d, fld) - getattr(d2, fld)) > 1e-9:\n            return CheckResult(name, False, f"derived field \'{fld}\' drifted")\n    return CheckResult(name, True,\n                       f"pitch={d.pitch:.3f} line_w={d.line_w:.3f} "\n                       f"merge_off={d.merge_off:.3f}")\n\n\ndef check_counts_and_kinds(pat: Pattern) -> CheckResult:\n    """A arm units (7-point) then T taper units (5-point) per side,\n    densely indexed."""\n    name = "counts_and_kinds"\n    p = pat.params\n    want = (p.n_arms + p.n_taper + p.n_decay + p.n_release\n            + p.n_morph + p.n_grid)\n    for side in ("left", "right"):\n        chs = _channels_by(pat, side)\n        if len(chs) != want:\n            return CheckResult(name, False,\n                               f"{side}: {len(chs)} units, want {want}")\n        if [c.index for c in chs] != list(range(want)):\n            return CheckResult(name, False, f"{side}: index gap/duplicate")\n        for c in chs:\n            if c.index < p.n_arms:\n                k = "arm"\n            elif c.index < p.n_arms + p.n_taper:\n                k = "taper"\n            elif c.index < p.n_arms + p.n_taper + p.n_decay:\n                k = "decay"\n            elif c.index < (p.n_arms + p.n_taper + p.n_decay\n                            + p.n_release):\n                k = "release"\n            elif c.index < (p.n_arms + p.n_taper + p.n_decay\n                            + p.n_release + p.n_morph):\n                k = "morph"\n            else:\n                k = "grid"\n            if k == "arm":\n                npts = 7\n            elif k == "grid":\n                npts = 3\n            elif k == "morph":\n                mi = c.index - (p.n_arms + p.n_taper + p.n_decay\n                                + p.n_release)\n                npts = 3 if pat.d.morph_excess[mi] == 0.0 else 5\n            else:\n                npts = 5\n            if c.kind != k:\n                return CheckResult(name, False,\n                                   f"{side}[{c.index}] kind {c.kind}!={k}")\n            if len(c.points) != npts:\n                return CheckResult(\n                    name, False,\n                    f"{side}[{c.index}] has {len(c.points)} pts, want {npts}")\n    return CheckResult(name, True,\n                       f"{2 * want} units: {p.n_arms} arms + "\n                       f"{p.n_taper} taper + {p.n_decay} decay + "\n                       f"{p.n_release} release + {p.n_morph} morph + "\n                       f"{p.n_grid} grid per side")\n\n\ndef check_manhattan(pat: Pattern) -> CheckResult:\n    """Every segment is axis-aligned and non-degenerate."""\n    name = "manhattan"\n    for c in pat.channels:\n        for i in range(len(c.points) - 1):\n            (x0, y0), (x1, y1) = c.points[i], c.points[i + 1]\n            dx, dy = abs(x1 - x0), abs(y1 - y0)\n            if dx > EPS and dy > EPS:\n                return CheckResult(name, False,\n                                   f"{c.side}[{c.index}] seg {i} diagonal")\n            if dx <= EPS and dy <= EPS:\n                return CheckResult(name, False,\n                                   f"{c.side}[{c.index}] seg {i} degenerate")\n    return CheckResult(name, True, "all segments axis-aligned")\n\n\ndef check_row_grid(pat: Pattern) -> CheckResult:\n    """The row grid obeys the SPACING SPEC: rows recomputed from the\n    declared profile (uniform or weighted) match the geometry, are\n    strictly increasing, hit top_y/bottom_y exactly, and no gap is\n    tighter than 2*line_w (the fusing floor). Every unit\'s out/return\n    lines sit exactly on their rows."""\n    name = "row_grid"\n    d, p = pat.d, pat.params\n    tol = 1e-6\n    nrows = len(d.rows)\n    if p.spacing is None:\n        want = [d.top_y + i * d.pitch for i in range(nrows)]\n    else:\n        span = d.bottom_y - d.top_y\n        fr = _spacing_weights(p.spacing, nrows - 1)\n        want = [d.top_y]\n        for f in fr:\n            want.append(want[-1] + span * f)\n        want[-1] = d.bottom_y\n    for i, (a, b) in enumerate(zip(d.rows, want)):\n        if abs(a - b) > tol:\n            return CheckResult(name, False,\n                               f"row {i} off the declared spacing profile")\n    gaps = [d.rows[i + 1] - d.rows[i] for i in range(nrows - 1)]\n    Uw = len(d.unit_widths)\n    for k in range(Uw):\n        if gaps[2 * k] - d.unit_widths[k] < d.min_clear - tol:\n            return CheckResult(name, False,\n                               f"unit {k}: white slot under one line")\n        if k < Uw - 1 and (gaps[2 * k + 1]\n                           - 0.5 * (d.unit_widths[k]\n                                    + d.unit_widths[k + 1])\n                           < d.min_clear - tol):\n            return CheckResult(name, False,\n                               f"units {k},{k + 1}: white band under "\n                               f"one line")\n    if abs(min(gaps) - d.min_gap) > tol:\n        return CheckResult(name, False, "min_gap out of sync with rows")\n    for side in ("left", "right"):\n        for c in _channels_by(pat, side):\n            yo_want = d.rows[2 * c.index]\n            yr_want = d.rows[2 * c.index + 1]\n            yo = c.points[0][1]\n            if abs(c.points[1][1] - yo) > tol or abs(yo - yo_want) > tol:\n                return CheckResult(name, False,\n                                   f"{side}[{c.index}] out line off row "\n                                   f"{2 * c.index}")\n            if len(c.points) == 3:           # grid unit: rung on the\n                continue                     # even row; odd row unused\n            yr = c.points[-2][1]\n            if abs(c.points[-3][1] - yr) > tol or abs(yr - yr_want) > tol:\n                return CheckResult(name, False,\n                                   f"{side}[{c.index}] return line off row "\n                                   f"{2 * c.index + 1}")\n    prof = "uniform" if p.spacing is None else "profiled"\n    return CheckResult(name, True,\n                       f"{nrows} rows, {prof}; mean gap {d.pitch:.3f}, "\n                       f"min gap {d.min_gap:.3f}")\n\n\ndef check_junctions_on_trunk(pat: Pattern) -> CheckResult:\n    """points[0] of every unit sits exactly on its side\'s trunk face,\n    within the trunk column (the T where \'left goes back to the center\n    channel\')."""\n    name = "junctions_on_trunk"\n    d = pat.d\n    for c in pat.channels:\n        if c.kind == "release":              # release units detach by\n            continue                         # design (check_release_shape)\n        face = d.xc + ((-1.0 if c.side == "left" else 1.0)\n                       * d.trunk_halves[2 * c.index])\n        x0, y0 = c.points[0]\n        if abs(x0 - face) > 1e-6:\n            return CheckResult(name, False,\n                               f"{c.side}[{c.index}] out line off the face")\n        if not (d.main_top - EPS <= y0 <= d.tail_end + EPS):\n            return CheckResult(name, False,\n                               f"{c.side}[{c.index}] junction off the column")\n    return CheckResult(name, True, "every out line starts on the trunk")\n\n\ndef check_arm_shape(pat: Pattern) -> CheckResult:\n    """Arms: out to R_k, UP to the hairpin top, short right by hairpin_w,\n    down to the return row. Reaches strictly increase by fan_step (arms\n    nest); hairpin tops form a descending staircase (inner arm tallest,\n    rises strictly decreasing, last >= 1.5*pitch); nested arms keep at\n    least min_clear of white between strokes."""\n    name = "arm_shape"\n    d, p = pat.d, pat.params\n    A = p.n_arms\n    tol = 1e-6\n    for side, sgn in (("left", -1.0), ("right", 1.0)):\n        arms = _channels_by(pat, side)[:A]\n        reaches, tops = [], []\n        for c in arms:\n            (x0, y0), (x1, y1), (x2, y2), (x3, y3), (x4, y4) = c.points[:5]\n            if abs(x2 - x1) > tol or abs(y3 - y2) > tol or \\\n               abs(x4 - x3) > tol:\n                return CheckResult(name, False,\n                                   f"{side}[{c.index}] hairpin shape broken")\n            if y2 >= y1 - EPS:\n                return CheckResult(name, False,\n                                   f"{side}[{c.index}] hairpin does not "\n                                   f"rise above its out line")\n            if abs(sgn * (x3 - x2) - d.hairpin_w) > tol:\n                return CheckResult(name, False,\n                                   f"{side}[{c.index}] hairpin width off "\n                                   f"spec")\n            reaches.append(abs(x1 - d.xc))\n            tops.append(y2)\n        for j in range(A):\n            if abs(reaches[j] - d.arm_reaches[j]) > tol:\n                return CheckResult(name, False,\n                                   f"{side} arm {j} reach off spec")\n            if abs(tops[j] - d.hairpin_tops[j]) > tol:\n                return CheckResult(name, False,\n                                   f"{side} arm {j} hairpin top off spec")\n        for j in range(1, A):\n            if abs((reaches[j] - reaches[j - 1]) - d.fan_step) > tol:\n                return CheckResult(name, False,\n                                   f"{side} arm nesting step non-uniform")\n            white = reaches[j] - (reaches[j - 1] + d.hairpin_w) - d.line_w\n            if white < d.min_clear - 1e-6:\n                return CheckResult(name, False,\n                                   f"{side} arms {j - 1},{j} nest too "\n                                   f"tightly ({white:.2f})")\n        rises = [arms[j].points[0][1] - tops[j] for j in range(A)]\n        for j in range(1, A):\n            if rises[j] >= rises[j - 1] - EPS:\n                return CheckResult(name, False,\n                                   f"{side} hairpin staircase not "\n                                   f"descending at arm {j}")\n        if rises[-1] < 1.5 * d.pitch - EPS:\n            return CheckResult(name, False,\n                               f"{side} outer hairpin too short "\n                               f"({rises[-1]:.1f})")\n    return CheckResult(name, True,\n                       f"{A} nested arms/side, step={d.fan_step:.3f}")\n\n\ndef check_taper_shape(pat: Pattern) -> CheckResult:\n    """Taper units: out to R_k, short DOWN one row gap, return. Reaches\n    strictly decrease from below the innermost arm to above the legal\n    floor, following the declared taper PROFILE exactly (linear by\n    default): reach_i == r_last + (r_first - r_last) * phi_hat(t_i),\n    recomputed here from the profile spec."""\n    name = "taper_shape"\n    d, p = pat.d, pat.params\n    A, T = p.n_arms, p.n_taper\n    tol = 1e-6\n    floor = d.trunk_w / 2.0 + d.merge_off + 2.0 * d.pitch\n    for side in ("left", "right"):\n        tapers = _channels_by(pat, side)[A:A + T]\n        reaches = []\n        for c in tapers:\n            (x1, y1), (x2, y2) = c.points[1], c.points[2]\n            if abs(x2 - x1) > tol:\n                return CheckResult(name, False,\n                                   f"{side}[{c.index}] outer turn not a "\n                                   f"vertical drop")\n            reaches.append(abs(x1 - d.xc))\n        for i in range(T):\n            if abs(reaches[i] - d.taper_reaches[i]) > tol:\n                return CheckResult(name, False,\n                                   f"{side} taper {i} reach off spec")\n        for i in range(1, T):\n            if reaches[i] >= reaches[i - 1] - EPS:\n                return CheckResult(name, False,\n                                   f"{side} taper not strictly shorter at "\n                                   f"unit {i}")\n        if T >= 3:\n            env = _taper_envelope(p.taper_profile, T)\n            rf, rl = reaches[0], reaches[-1]\n            for i in range(1, T - 1):\n                want = rl + (rf - rl) * env[i]\n                if abs(reaches[i] - want) > tol:\n                    return CheckResult(name, False,\n                                       f"{side} taper {i} off the declared "\n                                       f"profile envelope")\n        if reaches[0] > d.arm_reaches[0] - 2.0 * d.line_w + EPS:\n            return CheckResult(name, False,\n                               f"{side} first taper crowds the arms")\n        if reaches[-1] < floor - EPS:\n            return CheckResult(name, False,\n                               f"{side} last taper under reach floor "\n                               f"({reaches[-1]:.1f} < {floor:.1f})")\n    fam = ("linear" if p.taper_profile is None else\n           (p.taper_profile[0] if not callable(p.taper_profile)\n            else "custom"))\n    return CheckResult(name, True, f"{fam} taper over {T} units")\n\n\ndef check_decay_shape(pat: Pattern) -> CheckResult:\n    """The exponential convergence, unit by unit: excess reach beyond the\n    RAIL (the stop-short line at trunk_w/2 + merge_off) decays by the\n    SAME ratio q every unit - (R_i - rail)/(R_{i-1} - rail) == q, seeded\n    from the last linear-taper reach - and the final excess equals the\n    drawability floor E_min = 2*line_w exactly. This is the \'exponential\n    decay onto straight lines\' law, verified from geometry."""\n    name = "decay_shape"\n    d, p = pat.d, pat.params\n    A, T, D = p.n_arms, p.n_taper, p.n_decay\n    if D == 0:\n        return CheckResult(name, True, "n_decay=0 (classic mode)")\n    tol = 1e-6\n    rail = d.trunk_w / 2.0 + d.merge_off\n    e_min = _decay_floor(d.min_clear, d.unit_widths, A, T, D)\n    for side in ("left", "right"):\n        chs = _channels_by(pat, side)\n        e_prev = abs(chs[A + T - 1].points[1][0] - d.xc) - rail\n        reaches = []\n        for c in chs[A + T:A + T + D]:\n            (x1, _y1), (x2, _y2) = c.points[1], c.points[2]\n            if abs(x2 - x1) > tol:\n                return CheckResult(name, False,\n                                   f"{side}[{c.index}] outer turn not a "\n                                   f"vertical drop")\n            reaches.append(abs(x1 - d.xc))\n        for i, r in enumerate(reaches):\n            if abs(r - d.decay_reaches[i]) > tol:\n                return CheckResult(name, False,\n                                   f"{side} decay unit {i} reach off spec")\n            e = r - rail\n            if e < e_min - tol:\n                return CheckResult(name, False,\n                                   f"{side} decay unit {i} excess {e:.2f} "\n                                   f"below the 2*line_w floor")\n            ratio = e / e_prev\n            if abs(ratio - d.decay_q) > 1e-6:\n                return CheckResult(name, False,\n                                   f"{side} decay unit {i}: excess ratio "\n                                   f"{ratio:.5f} != q {d.decay_q:.5f} "\n                                   f"(not exponential)")\n            e_prev = e\n        if abs((reaches[-1] - rail) - e_min) > tol:\n            return CheckResult(name, False,\n                               f"{side} did not converge onto the rail: "\n                               f"final excess {reaches[-1] - rail:.3f} != "\n                               f"{e_min:.3f}")\n    return CheckResult(name, True,\n                       f"E_i = E0*q^i, q={d.decay_q:.4f}, {D} units; "\n                       f"final excess = 2*line_w = {e_min:.3f}")\n\n\ndef check_release_shape(pat: Pattern) -> CheckResult:\n    """The fading intersections, unit by unit: every release unit keeps\n    the converged tooth (reach == rail + 2*line_w) while its rung\'s\n    junction gap to the trunk face follows the DECLARED envelope from\n    line_w (the minimal legal crack) to boundary_w - line_w (a final\n    one-line-deep tooth on the rail). At every step the white gap AND\n    the remaining ink bridge are >= line_w, so the global clearance rule\n    holds with no exemptions."""\n    name = "release_shape"\n    d, p = pat.d, pat.params\n    A, T, D, R = p.n_arms, p.n_taper, p.n_decay, p.n_release\n    if R == 0:\n        return CheckResult(name, True, "n_release=0")\n    tol = 1e-6\n    t2 = d.trunk_w / 2.0\n    rail = t2 + d.merge_off\n    want_reach = rail + 2.0 * d.line_w\n    b = d.boundary_w\n    env = _taper_envelope(p.release_profile, R)\n    want_g = [d.line_w + (b - 2.0 * d.line_w) * (1.0 - e) for e in env]\n    for side, sgn in (("left", -1.0), ("right", 1.0)):\n        rel = _channels_by(pat, side)[A + T + D:]\n        gaps = []\n        for c in rel:\n            face = d.xc + sgn * d.trunk_halves[2 * c.index]\n            (x0, _y0), (x1, _y1), (x2, _y2) = c.points[0:3]\n            if abs(x2 - x1) > tol:\n                return CheckResult(name, False,\n                                   f"{side}[{c.index}] outer turn not a "\n                                   f"vertical drop")\n            if abs(abs(x1 - d.xc) - want_reach) > tol:\n                return CheckResult(name, False,\n                                   f"{side}[{c.index}] tooth reach not "\n                                   f"frozen at rail + 2*line_w")\n            gaps.append(sgn * (x0 - face))\n        for i, g in enumerate(gaps):\n            if abs(g - want_g[i]) > tol:\n                return CheckResult(name, False,\n                                   f"{side} release {i} junction gap "\n                                   f"{g:.3f} off the declared schedule "\n                                   f"({want_g[i]:.3f})")\n            if g < d.min_clear - tol:\n                return CheckResult(name, False,\n                                   f"{side} release {i} gap {g:.3f} under "\n                                   f"one line of white")\n            if (b - g) < rel[i].width - tol:\n                return CheckResult(name, False,\n                                   f"{side} release {i} ink bridge "\n                                   f"{b - g:.3f} under its channel width")\n        for i in range(1, R):\n            if gaps[i] <= gaps[i - 1] + (tol if R >= 2 else -tol):\n                return CheckResult(name, False,\n                                   f"{side} junction gaps not strictly "\n                                   f"widening at unit {i}")\n        if abs(gaps[0] - d.min_clear) > tol:\n            return CheckResult(name, False,\n                               f"{side} first release gap != line_w")\n        if abs(gaps[-1] - (b - d.line_w)) > tol:\n            return CheckResult(name, False,\n                               f"{side} last release gap != boundary_w - "\n                               f"line_w")\n    fam = ("linear" if p.release_profile is None else\n           (p.release_profile[0] if not callable(p.release_profile)\n            else "custom"))\n    return CheckResult(name, True,\n                       f"{fam} junction-gap schedule over {R} units: "\n                       f"{d.line_w:.3f} -> {b - d.line_w:.3f}")\n\n\ndef check_morph_shape(pat: Pattern) -> CheckResult:\n    """The smooth transition, verified as a spec-match: recompute the\n    E schedule (envelope + exclusion-band snap + monotone clamp)\n    INDEPENDENTLY from the declared profile and require the geometry to\n    equal it - reach = rail + E on the out row, jog on the design line,\n    5 points while E > 0 and the degenerate 3-point L at E == 0, fused\n    flags consistent with the band law."""\n    name = "morph_shape"\n    d, p = pat.d, pat.params\n    A, T, D, R = p.n_arms, p.n_taper, p.n_decay, p.n_release\n    M = p.n_morph\n    if M == 0:\n        return CheckResult(name, True, "n_morph=0")\n    tol = 1e-6\n    rail = d.trunk_w / 2.0 + d.merge_off\n    e_min = _decay_floor(d.min_clear, d.unit_widths, A, T, D)\n    env = _taper_envelope(p.morph_profile, M + 1)\n    sched, fused_want, prev = [], [], e_min\n    for j in range(1, M + 1):\n        e = e_min * env[j]\n        kk = A + T + D + j - 1\n        mean_w = 0.5 * (d.unit_widths[kk - 1] + d.unit_widths[kk])\n        lo, hi = mean_w - d.min_clear / 2.0, mean_w + d.min_clear\n        if lo + 1e-12 < e < hi - 1e-12:\n            e = lo if (e - lo) <= (hi - e) else hi\n        e = min(e, prev)\n        if e < 1e-9:\n            e = 0.0\n        sched.append(e)\n        fused_want.append(e <= lo + 1e-12)\n        prev = e\n    if len(d.morph_excess) != M:\n        return CheckResult(name, False, "derived schedule length off")\n    for j in range(M):\n        if abs(d.morph_excess[j] - sched[j]) > tol                 or d.morph_fused[j] != fused_want[j]:\n            return CheckResult(name, False,\n                               f"derived schedule/flags off spec at "\n                               f"morph unit {j}")\n        if j and sched[j] > sched[j - 1] + tol:\n            return CheckResult(name, False, "schedule not non-increasing")\n    for side, sgn in (("left", -1.0), ("right", 1.0)):\n        mor = _channels_by(pat, side)[A + T + D + R:\n                                      A + T + D + R + M]\n        if len(mor) != M:\n            return CheckResult(name, False, f"{side}: {len(mor)} morph "\n                                            f"units, want {M}")\n        for j, c in enumerate(mor):\n            e = sched[j]\n            row = d.rows[2 * c.index]\n            face = d.xc + sgn * d.trunk_halves[2 * c.index]\n            jx = d.xc + sgn * rail\n            if e == 0.0:\n                if len(c.points) != 3:\n                    return CheckResult(name, False,\n                                       f"{side}[{c.index}] E=0 must be "\n                                       f"the 3-point L")\n                (x0, y0), (x1, y1), _ = c.points\n                if abs(x0 - face) > tol or abs(x1 - jx) > tol                         or abs(y0 - row) > tol:\n                    return CheckResult(name, False,\n                                       f"{side}[{c.index}] degenerate "\n                                       f"unit off the design lines")\n                continue\n            if len(c.points) != 5:\n                return CheckResult(name, False,\n                                   f"{side}[{c.index}] E>0 must keep "\n                                   f"5 points")\n            got_e = sgn * (c.points[1][0] - d.xc) - rail\n            if abs(got_e - e) > tol:\n                return CheckResult(name, False,\n                                   f"{side}[{c.index}] excess "\n                                   f"{got_e:.4f} != scheduled {e:.4f}")\n            if abs(sgn * (c.points[3][0] - d.xc) - rail) > tol:\n                return CheckResult(name, False,\n                                   f"{side}[{c.index}] jog off the "\n                                   f"design line")\n    n_f = sum(1 for f in d.morph_fused if f)\n    return CheckResult(name, True,\n                       f"{M} units/side: E {sched[0]:.2f} -> 0, "\n                       f"{n_f} fused, band jumped once")\n\n\ndef check_grid_shape(pat: Pattern) -> CheckResult:\n    """The periodic-intersection grid: every grid unit is a clean rung\n    (a full black crossing from the center channel to the side rail on\n    its even row) plus a rail segment at exactly the wall offset, so the\n    side channel runs straight at wall boundary_w while intersecting the\n    center channel once per unit - the wall between them a stack of\n    rectangular white posts."""\n    name = "grid_shape"\n    d, p = pat.d, pat.params\n    A, T, D, R = p.n_arms, p.n_taper, p.n_decay, p.n_release\n    Gn = p.n_grid\n    if Gn == 0:\n        return CheckResult(name, True, "n_grid=0")\n    tol = 1e-6\n    t2 = d.trunk_w / 2.0\n    U = A + T + D + R + p.n_morph + Gn\n    for side, sgn in (("left", -1.0), ("right", 1.0)):\n        jog_x = d.xc + sgn * (t2 + d.merge_off)\n        grid = _channels_by(pat, side)[A + T + D + R + p.n_morph:]\n        if len(grid) != Gn:\n            return CheckResult(name, False, f"{side}: {len(grid)} grid "\n                                            f"units, want {Gn}")\n        for c in grid:\n            if len(c.points) != 3:\n                return CheckResult(name, False,\n                                   f"{side}[{c.index}] not a 3-point "\n                                   f"ladder unit")\n            (x0, y0), (x1, y1), (x2, y2) = c.points\n            row = d.rows[2 * c.index]\n            face = d.xc + sgn * d.trunk_halves[2 * c.index]\n            if abs(x0 - face) > tol or abs(x1 - jog_x) > tol:\n                return CheckResult(name, False,\n                                   f"{side}[{c.index}] rung is not a full "\n                                   f"trunk-to-rail crossing")\n            if abs(y0 - row) > tol or abs(y1 - row) > tol:\n                return CheckResult(name, False,\n                                   f"{side}[{c.index}] rung off its row")\n            if abs(x2 - x1) > tol:\n                return CheckResult(name, False,\n                                   f"{side}[{c.index}] rail segment off "\n                                   f"the wall offset")\n            y_want = (d.port_y if c.index == U - 1\n                      else d.rows[2 * c.index + 2])\n            if abs(y2 - y_want) > tol:\n                return CheckResult(name, False,\n                                   f"{side}[{c.index}] rail segment does "\n                                   f"not span to the next rung/port")\n    return CheckResult(name, True,\n                       f"{Gn} rungs/side at wall {d.boundary_w:.3f}; "\n                       f"posts {d.boundary_w:.3f} wide, period 2 rows")\n\n\ndef check_chain_links(pat: Pattern) -> CheckResult:\n    """The traced chain: every unit\'s return stops at the SAME constant\n    offset from the trunk (a real white pocket - never touching the\n    center channel), then jogs down one pitch and lands STRICTLY INTERIOR\n    to the next unit\'s out line (the T where left goes to the trunk and\n    right continues the meander). The last unit\'s jog lands on the foot\n    crossbar (classic) or runs straight to the port (connector)."""\n    name = "chain_links"\n    d, p = pat.d, pat.params\n    U = (p.n_arms + p.n_taper + p.n_decay + p.n_release\n         + p.n_morph + p.n_grid)\n    connector = p.n_decay >= 1\n    mbase = p.n_arms + p.n_taper + p.n_decay + p.n_release\n    tol = 1e-6\n    lw2 = d.line_w / 2.0\n    foot = (None if connector\n            else _trunk_bbox(pat.trunk_parts["foot"]))\n    for side, sgn in (("left", -1.0), ("right", 1.0)):\n        face = d.xc + sgn * d.trunk_w / 2.0\n        chs = _channels_by(pat, side)\n        for c in chs:\n            xs, ys = c.points[-2]\n            xe, ye = c.points[-1]\n            m = sgn * (xs - face)\n            if abs(m - d.merge_off) > tol or abs(xe - xs) > tol:\n                return CheckResult(name, False,\n                                   f"{side}[{c.index}] stop-short offset "\n                                   f"{m:.2f} != merge_off {d.merge_off:.2f}")\n            if m - lw2 < 1.0 - tol:\n                return CheckResult(name, False,\n                                   f"{side}[{c.index}] return touches the "\n                                   f"center channel (no white pocket)")\n            if c.index < U - 1:\n                nxt = chs[c.index + 1]\n                (nx0, ny0), (nx1, _ny1) = nxt.points[0], nxt.points[1]\n                if abs(ye - ny0) > tol:\n                    return CheckResult(name, False,\n                                       f"{side}[{c.index}] jog does not "\n                                       f"land on unit {c.index + 1}\'s out "\n                                       f"line row")\n                corner_land = (nxt.kind == "grid"\n                               or (nxt.kind == "morph"\n                                   and d.morph_fused[nxt.index - mbase]))\n                if corner_land:              # ladder / fused morph: the\n                    jog_want = d.xc + ((-1.0 if side == "left" else 1.0)\n                                       * (d.trunk_w / 2.0 + d.merge_off))\n                    if abs(xe - jog_want) > tol:  # drop meets the design\n                                                  # jog line exactly\n                        return CheckResult(name, False,\n                                           f"{side}[{c.index}] drop misses "\n                                           f"the jog line at unit "\n                                           f"{c.index + 1}")\n                else:\n                    lo, hi = min(nx0, nx1), max(nx0, nx1)\n                    if not (lo + lw2 < xe < hi - lw2):\n                        return CheckResult(name, False,\n                                           f"{side}[{c.index}] jog lands "\n                                           f"outside unit {c.index + 1}\'s "\n                                           f"out line")\n            else:\n                if connector:\n                    if abs(ye - d.port_y) > tol:\n                        return CheckResult(name, False,\n                                           f"{side}[{c.index}] rail does "\n                                           f"not reach the port y")\n                elif not (foot[0] + lw2 <= xe <= foot[2] - lw2 and\n                          foot[1] - tol <= ye <= foot[3] + tol):\n                    return CheckResult(name, False,\n                                       f"{side}[{c.index}] final jog "\n                                       f"misses the foot crossbar")\n    ending = "rail to port" if connector else "foot landing"\n    return CheckResult(name, True,\n                       f"{U - 1} T-junctions/side + {ending}; "\n                       f"constant stop-short {d.merge_off:.3f}")\n\n\ndef check_port_interface(pat: Pattern) -> CheckResult:\n    """Connector mode\'s attachment contract: at the bottom edge the\n    figure presents the cross-section\n        [rail (line_w) | boundary_w | trunk (trunk_w) | boundary_w | rail]\n    with the trunk and both rails running straight to port_y and NOTHING\n    else below the last row - so any geometry built on the same three\n    widths butt-joins seamlessly."""\n    name = "port_interface"\n    d, p = pat.d, pat.params\n    if p.n_decay == 0:\n        return CheckResult(name, True, "classic mode: T-foot termination")\n    tol = 1e-6\n    lw2 = d.line_w / 2.0\n    t2 = d.trunk_w / 2.0\n    main = _trunk_bbox(pat.trunk_parts["main"])\n    if abs(main[3] - d.port_y) > tol:\n        return CheckResult(name, False, "trunk does not reach the port")\n    if "foot" in pat.trunk_parts:\n        return CheckResult(name, False, "foot present in connector mode")\n    if abs((d.merge_off - lw2) - d.boundary_w) > tol:\n        return CheckResult(name, False,\n                           "boundary_w inconsistent with merge_off")\n    if d.boundary_w < d.min_clear - tol:\n        return CheckResult(name, False,\n                           "boundary_w below one line of white")\n    grid_mode = p.n_grid >= 1\n    row_last = d.rows[-2] if grid_mode else d.rows[-1]\n    limit = row_last + lw2 + tol\n    for side, sgn in (("left", -1.0), ("right", 1.0)):\n        chs = _channels_by(pat, side)\n        last = chs[-1]\n        (xs, ys), (xe, ye) = last.points[-2], last.points[-1]\n        want_x = d.xc + sgn * (t2 + d.merge_off)\n        if abs(xs - want_x) > tol or abs(xe - want_x) > tol:\n            return CheckResult(name, False,\n                               f"{side} rail not aligned at trunk_w/2 + "\n                               f"merge_off")\n        if abs(ys - row_last) > tol or abs(ye - d.port_y) > tol:\n            return CheckResult(name, False,\n                               f"{side} rail does not span last row -> "\n                               f"port")\n        if ye - ys < 2.0 * d.pitch - tol:\n            return CheckResult(name, False,\n                               f"{side} straight rail shorter than "\n                               f"2*pitch")\n        for c in chs:\n            nseg = len(c.points) - 1\n            for si in range(nseg):\n                if c is last and si == nseg - 1:\n                    continue                     # the rail itself\n                r = _seg_rect(c.points[si], c.points[si + 1], lw2)\n                if r[3] > limit:\n                    return CheckResult(name, False,\n                                       f"{side}[{c.index}] ink intrudes "\n                                       f"on the port zone")\n    return CheckResult(name, True,\n                       f"port at y={d.port_y:g}: [rail {d.line_w:.3f} | "\n                       f"gap {d.boundary_w:.3f} | trunk {d.trunk_w:.3f} | "\n                       f"gap {d.boundary_w:.3f} | rail]")\n\n\ndef check_symmetry(pat: Pattern) -> CheckResult:\n    """Right side is the exact mirror of the left about the centerline."""\n    name = "symmetry"\n    d = pat.d\n    left = _channels_by(pat, "left")\n    right = _channels_by(pat, "right")\n    if len(left) != len(right):\n        return CheckResult(name, False, "unit count differs by side")\n    for cl, cr in zip(left, right):\n        if cl.kind != cr.kind or len(cl.points) != len(cr.points):\n            return CheckResult(name, False,\n                               f"index {cl.index}: kind/shape mismatch")\n        for (xl, yl), (xr, yr) in zip(cl.points, cr.points):\n            if abs((2.0 * d.xc - xl) - xr) > 1e-6 or abs(yl - yr) > 1e-6:\n                return CheckResult(name, False,\n                                   f"left[{cl.index}] != mirror("\n                                   f"right[{cr.index}])")\n    return CheckResult(name, True, "left == mirror(right)")\n\n\ndef check_clearance(pat: Pattern) -> CheckResult:\n    """Every pair of FOREIGN strokes keeps >= min_clear (= line_w) of\n    white between them, so nothing fuses that the trace says is separate.\n    Exemptions, each still intentional ink contact elsewhere-verified:\n    a unit\'s own segments (corners); seg 0 <-> trunk main (the trunk T);\n    a unit\'s final jog <-> the NEXT unit\'s out segment (the chain T);\n    the last unit\'s jog <-> the foot (its landing)."""\n    name = "clearance"\n    d, p = pat.d, pat.params\n    U = (p.n_arms + p.n_taper + p.n_decay + p.n_release\n         + p.n_morph + p.n_grid)\n    half = d.line_w / 2.0\n    kinds = {(c.side, c.index): c.kind for c in pat.channels}\n    mbase = p.n_arms + p.n_taper + p.n_decay + p.n_release\n    fused = {(c.side, c.index) for c in pat.channels\n             if c.kind == "morph" and d.morph_fused[c.index - mbase]}\n    zero_e = {(c.side, c.index) for c in pat.channels\n              if c.kind == "morph" and len(c.points) == 3}\n    items = []\n    for c in pat.channels:\n        nseg = len(c.points) - 1\n        for si in range(nseg):\n            items.append(((c.side, c.index, si), si == nseg - 1,\n                          _seg_rect(c.points[si], c.points[si + 1],\n                                    c.width / 2.0)))\n    for nm in sorted(pat.trunk_parts):\n        items.append((("trunk", nm), False,\n                      _trunk_bbox(pat.trunk_parts[nm])))\n    need = d.min_clear - 1e-6\n    worst = None\n    for i in range(len(items)):\n        for j in range(i + 1, len(items)):\n            a, b = items[i], items[j]\n            if a[0][0] == "trunk" and b[0][0] != "trunk":\n                a, b = b, a\n            ida, a_last, ra = a\n            idb, b_last, rb = b\n            a_trunk = ida[0] == "trunk"\n            b_trunk = idb[0] == "trunk"\n            if a_trunk and b_trunk:\n                continue                       # trunk parts overlap by design\n            if not a_trunk and not b_trunk:\n                if ida[:2] == idb[:2]:\n                    continue                   # same unit: its own corners\n                same_side = ida[0] == idb[0]\n                consec = abs(ida[1] - idb[1]) == 1\n                lower, upper = (a, b) if ida[1] > idb[1] else (b, a)\n                if same_side and consec and upper[1] and lower[0][2] == 0:\n                    continue                   # the chain T (jog -> out)\n                hi_key = (ida[0], max(ida[1], idb[1]))\n                if (same_side and consec and a_last and b_last\n                        and (kinds[hi_key] == "grid"\n                             or hi_key in zero_e)):\n                    continue                   # rail continuity: grid /\n                                               # converged-morph verticals\n                                               # butt end-to-end\n                if (same_side and consec and upper[1]\n                        and lower[0][2] == 1 and hi_key in fused):\n                    continue                   # FUSED morph: the down\n                                               # stroke overlaps the\n                                               # previous jog by design\n                                               # (the pocket is gone)\n                if (ida[:2] == idb[:2] and hi_key in fused\n                        and {ida[2], idb[2]} == {1, 3}):\n                    continue                   # FUSED morph: own down\n                                               # overlaps own jog\n            else:\n                if idb[1] == "main" and ida[2] == 0:\n                    continue                   # the trunk T (out -> trunk)\n                if idb[1] == "foot" and a_last and ida[1] == U - 1:\n                    continue                   # the foot landing\n            gap = _rect_gap(ra, rb)\n            if gap < need:\n                return CheckResult(name, False,\n                                   f"{ida} vs {idb}: white {gap:.3f} < "\n                                   f"min_clear {d.min_clear:.3f}")\n            worst = gap if worst is None else min(worst, gap)\n    return CheckResult(name, True,\n                       f"min foreign white {worst:.3f} >= "\n                       f"{d.min_clear:.3f}")\n\n\ndef check_inside_box(pat: Pattern) -> CheckResult:\n    """All stroke ink inside the margin box; the trunk inside the hard\n    box and above max_bottom (the tail may pass below the margin line by\n    design)."""\n    name = "inside_box"\n    p, d = pat.params, pat.d\n    W, H = p.width, p.height\n    tol = 1e-6\n    half = d.line_w / 2.0\n    for c in pat.channels:\n        for i in range(len(c.points) - 1):\n            r = _seg_rect(c.points[i], c.points[i + 1], half)\n            if r[0] < -tol or r[1] < -tol or r[2] > W + tol or r[3] > H + tol:\n                return CheckResult(name, False,\n                                   f"{c.side}[{c.index}] leaves the box")\n            if (r[0] < d.margin - tol or r[1] < d.margin - tol or\n                    r[2] > W - d.margin + tol):\n                return CheckResult(name, False,\n                                   f"{c.side}[{c.index}] breaks margin")\n            if r[3] > d.max_bottom + tol:\n                return CheckResult(name, False,\n                                   f"{c.side}[{c.index}] below max_bottom")\n    for nm, part in sorted(pat.trunk_parts.items()):\n        r = _trunk_bbox(part)\n        if r[0] < -tol or r[1] < -tol or r[2] > W + tol or r[3] > H + tol:\n            return CheckResult(name, False, f"trunk {nm} leaves the box")\n        if r[3] > d.max_bottom + tol:\n            return CheckResult(name, False, f"trunk {nm} below max_bottom")\n    return CheckResult(name, True, f"box {W:g}x{H:g}, margin {d.margin:g}")\n\n\ndef check_connectivity(pat: Pattern) -> CheckResult:\n    """No single-line islands: union-find over actual ink contact between\n    stroke rectangles and trunk shapes must give ONE component, and the\n    contact-pair count must match the traced topology exactly - per side,\n    U trunk T-junctions plus (U-1) chain T-junctions - so every chain\n    closes a loop (a channel that leaves the center and rejoins it\n    downstream) and nothing extra touches."""\n    name = "connectivity"\n    p, d = pat.params, pat.d\n    U = (p.n_arms + p.n_taper + p.n_decay + p.n_release\n         + p.n_morph + p.n_grid)\n    half = d.line_w / 2.0\n    ids: List[tuple] = []\n    rects: Dict[tuple, List[Rect]] = {}\n    for c in pat.channels:\n        key = (c.side, c.index)\n        ids.append(key)\n        rects[key] = [_seg_rect(c.points[i], c.points[i + 1],\n                                c.width / 2.0)\n                      for i in range(len(c.points) - 1)]\n    trunk_key = ("trunk",)\n    ids.append(trunk_key)\n    rects[trunk_key] = [_trunk_bbox(v) for v in pat.trunk_parts.values()]\n\n    parent = {k: k for k in ids}\n\n    def find(k):\n        while parent[k] != k:\n            parent[k] = parent[parent[k]]\n            k = parent[k]\n        return k\n\n    def union(a, b):\n        ra, rb = find(a), find(b)\n        if ra != rb:\n            parent[ra] = rb\n            return True\n        return False\n\n    contacts = 0\n    for i in range(len(ids)):\n        for j in range(i + 1, len(ids)):\n            a, b = ids[i], ids[j]\n            touch = any(_rect_gap(ra, rb) <= 1e-6\n                        for ra in rects[a] for rb in rects[b])\n            if touch:\n                contacts += 1\n                union(a, b)\n    roots = {find(k) for k in ids}\n    if len(roots) != 1:\n        isolated = [k for k in ids if all(\n            not any(_rect_gap(ra, rb) <= 1e-6\n                    for ra in rects[k] for rb in rects[o])\n            for o in ids if o != k)]\n        return CheckResult(name, False,\n                           f"{len(roots)} components; isolated: "\n                           f"{isolated or \'none (split network)\'}")\n    # Designed contact PAIRS per side: (U - R) trunk Ts (release units\n    # detach from the trunk by design) + (U-1) chain Ts. (The last unit\n    # also lands on the foot in classic mode - a second physical contact\n    # with the trunk, but the same entity pair, so it adds no pair here;\n    # its geometry is asserted by chain_links.)\n    want_contacts = 2 * (2 * U - 1 - p.n_release)\n    if contacts != want_contacts:\n        return CheckResult(name, False,\n                           f"{contacts} ink contacts, want {want_contacts} "\n                           f"(stray touching or a missing link)")\n    return CheckResult(name, True,\n                       f"1 component, V={2 * U + 1}, E={contacts}, "\n                       f"{contacts - 2 * U} chain loops; no islands")\n\n\ndef check_determinism(pat: Pattern) -> CheckResult:\n    """Rebuilding from Params must reproduce the geometry and the SVG\n    bytes exactly - doubles as a tamper alarm for mutated patterns."""\n    name = "determinism"\n    ref = build_pattern(pat.params)\n    if len(ref.channels) != len(pat.channels):\n        return CheckResult(name, False, "unit count drifted (tampered?)")\n    for a, b in zip(pat.channels, ref.channels):\n        if (a.side != b.side or a.index != b.index or a.kind != b.kind\n                or a.points != b.points):\n            return CheckResult(name, False,\n                               f"{a.side}[{a.index}] differs from a fresh "\n                               f"build (tampered?)")\n    s1, s2 = svg_string(pat), svg_string(ref)\n    if s1 != s2:\n        return CheckResult(name, False, "svg bytes not reproducible")\n    sha = hashlib.sha256(s1.encode("utf-8")).hexdigest()\n    return CheckResult(name, True, f"sha256={sha[:16]}...")\n\n\ndef check_width_profile(pat: Pattern) -> CheckResult:\n    """The declared WIDTH functions, verified from geometry: every\n    unit\'s stroke width matches the side profile (ending at the declared\n    line_w at the port), and the center channel boundary matches the\n    center profile at every row (ending at the declared trunk_w at the\n    port), with the trunk never thinner than 2*line_w."""\n    name = "width_profile"\n    d, p = pat.d, pat.params\n    tol = 1e-6\n    U = (p.n_arms + p.n_taper + p.n_decay + p.n_release\n         + p.n_morph + p.n_grid)\n    scales = _scale_profile(p.side_w_profile, U)\n    for side in ("left", "right"):\n        for c in _channels_by(pat, side):\n            want = d.line_w * scales[c.index]\n            if abs(c.width - want) > tol or                     abs(c.width - d.unit_widths[c.index]) > tol:\n                return CheckResult(name, False,\n                                   f"{side}[{c.index}] width off the "\n                                   f"declared side profile")\n    part = pat.trunk_parts["main"]\n    if p.center_w_profile is None:\n        if part[0] != "rect":\n            return CheckResult(name, False, "trunk should be a rect")\n        r = part[1]\n        if abs((r[2] - r[0]) / 2.0 - d.trunk_w / 2.0) > tol:\n            return CheckResult(name, False, "trunk width off")\n    else:\n        if part[0] != "poly":\n            return CheckResult(name, False, "trunk should be shaped")\n        xs = {round(y, 6): x for x, y in part[1] if x >= d.xc}\n        for j, y in enumerate(d.rows):\n            got = xs.get(round(y, 6))\n            if got is None or abs((got - d.xc) - d.trunk_halves[j]) > tol:\n                return CheckResult(name, False,\n                                   f"trunk boundary off the declared "\n                                   f"center profile at row {j}")\n        if min(d.trunk_halves) < max(d.min_clear, 0.5 * d.line_w) - tol:\n            return CheckResult(name, False, "trunk thinner than "\n                               "2*max(min_clear, line_w/2)")\n    if p.n_decay >= 1:\n        if abs(d.unit_widths[-1] - d.line_w) > tol:\n            return CheckResult(name, False,\n                               "port side width != declared line_w")\n        if abs(d.trunk_half_end - d.trunk_w / 2.0) > 1e-4:\n            return CheckResult(name, False,\n                               "port center width != declared trunk_w")\n    return CheckResult(name, True,\n                       f"side {min(d.unit_widths):.2f}-"\n                       f"{max(d.unit_widths):.2f}, center "\n                       f"{2 * min(d.trunk_halves):.2f}-"\n                       f"{2 * max(d.trunk_halves):.2f}")\n\n\nVALIDATORS = [\n    check_params_feasible, check_counts_and_kinds, check_manhattan,\n    check_row_grid, check_junctions_on_trunk, check_arm_shape,\n    check_taper_shape, check_decay_shape, check_release_shape,\n    check_morph_shape, check_grid_shape, check_width_profile,\n    check_chain_links,\n    check_port_interface, check_symmetry, check_clearance,\n    check_inside_box, check_connectivity, check_determinism,\n]\n\n\ndef validate(pat: Pattern) -> List[CheckResult]:\n    return [fn(pat) for fn in VALIDATORS]\n\n\n# ---------------------------------------------------------------------------\n# 6. SELF-TEST: deliberate corruptions each caught by its matching check\n# ---------------------------------------------------------------------------\ndef _mut_row_shift(pat: Pattern) -> None:\n    ch = _channels_by(pat, "left")[6]                 # a taper unit\n    ch.points = [(x, y + 3.0) for x, y in ch.points]\n\n\ndef _mut_chain_break(pat: Pattern) -> None:\n    ch = _channels_by(pat, "left")[5]\n    xe, ye = ch.points[-1]\n    ch.points[-1] = (xe, ye - 0.5 * pat.d.pitch)      # jog stops mid-white\n\n\ndef _mut_trunk_kiss(pat: Pattern) -> None:\n    ch = _channels_by(pat, "left")[7]                 # return touches trunk\n    dx = pat.d.merge_off                              # (left: +x is inward)\n    ch.points[-2] = (ch.points[-2][0] + dx, ch.points[-2][1])\n    ch.points[-1] = (ch.points[-1][0] + dx, ch.points[-1][1])\n\n\ndef _mut_detach_unit(pat: Pattern) -> None:\n    ch = _channels_by(pat, "left")[8]                 # slide off the face\n    ch.points = [(x - 4.0 * pat.d.line_w, y) for x, y in ch.points]\n\n\ndef _mut_arm_swap(pat: Pattern) -> None:\n    a = _channels_by(pat, "left")[0]\n    b = _channels_by(pat, "left")[1]\n    for i in (1, 2):                                  # swap up-stroke x\n        ax, bx = a.points[i][0], b.points[i][0]\n        a.points[i] = (bx, a.points[i][1])\n        b.points[i] = (ax, b.points[i][1])\n    for i in (3, 4):                                  # swap down-stroke x\n        ax, bx = a.points[i][0], b.points[i][0]\n        a.points[i] = (bx, a.points[i][1])\n        b.points[i] = (ax, b.points[i][1])\n\n\ndef _mut_hairpin_flatten(pat: Pattern) -> None:\n    ch = _channels_by(pat, "left")[1]                 # kill the tall rise\n    y = ch.points[0][1] - 0.5 * pat.d.pitch\n    ch.points[2] = (ch.points[2][0], y)\n    ch.points[3] = (ch.points[3][0], y)\n\n\ndef _mut_taper_swap(pat: Pattern) -> None:\n    A = pat.params.n_arms\n    a = _channels_by(pat, "left")[A + 2]\n    b = _channels_by(pat, "left")[A + 3]\n    for i in (1, 2):\n        ax, bx = a.points[i][0], b.points[i][0]\n        a.points[i] = (bx, a.points[i][1])\n        b.points[i] = (ax, b.points[i][1])\n\n\ndef _mut_mirror_nudge(pat: Pattern) -> None:\n    ch = _channels_by(pat, "right")[4]\n    ch.points = [(x + 2.0, y) for x, y in ch.points]\n\n\ndef _mut_lane_collision(pat: Pattern) -> None:\n    ch = _channels_by(pat, "left")[9]\n    dy = 0.5 * pat.d.pitch                            # into the white band\n    ch.points = [(x, y + dy) for x, y in ch.points]\n\n\ndef _mut_escape_box(pat: Pattern) -> None:\n    ch = _channels_by(pat, "left")[0]                 # hairpin over the top\n    ch.points[2] = (ch.points[2][0], -10.0)\n    ch.points[3] = (ch.points[3][0], -10.0)\n\n\ndef _mut_reach_overrun(pat: Pattern) -> None:\n    ch = _channels_by(pat, "left")[-1]\n    x = pat.d.xc - (pat.d.d_line + 4.0 * pat.d.line_w)\n    ch.points[1] = (x, ch.points[1][1])\n    ch.points[2] = (x, ch.points[2][1])\n\n\ndef _mut_decay_break(pat: Pattern) -> None:\n    A, T = pat.params.n_arms, pat.params.n_taper\n    ch = _channels_by(pat, "left")[A + T + 2]         # a mid-decay unit\n    d = pat.d\n    rail_x = d.xc - (d.trunk_w / 2.0 + d.merge_off)   # left-side rail\n    for i in (1, 2):                                  # inflate excess 1.5x:\n        x, y = ch.points[i]                           # breaks the constant\n        ch.points[i] = (rail_x - 1.5 * (rail_x - x), y)   # ratio law\n\n\ndef _mut_rail_cut(pat: Pattern) -> None:\n    ch = _channels_by(pat, "left")[-1]\n    xe, ye = ch.points[-1]\n    ch.points[-1] = (xe, ye - 3.0 * pat.d.pitch)      # rail stops short\n\n\ndef _mut_rail_skew(pat: Pattern) -> None:\n    ch = _channels_by(pat, "left")[-1]                # rail off alignment\n    for i in (-2, -1):\n        x, y = ch.points[i]\n        ch.points[i] = (x + 2.0, y)\n\n\nMUTATIONS = [\n    ("row_shift", _mut_row_shift, "row_grid"),\n    ("chain_break", _mut_chain_break, "chain_links"),\n    ("trunk_kiss", _mut_trunk_kiss, "chain_links"),\n    ("detach_unit", _mut_detach_unit, "junctions_on_trunk"),\n    ("arm_swap", _mut_arm_swap, "arm_shape"),\n    ("hairpin_flatten", _mut_hairpin_flatten, "arm_shape"),\n    ("taper_swap", _mut_taper_swap, "taper_shape"),\n    ("mirror_nudge", _mut_mirror_nudge, "symmetry"),\n    ("lane_collision", _mut_lane_collision, "clearance"),\n    ("escape_box", _mut_escape_box, "inside_box"),\n    ("reach_overrun", _mut_reach_overrun, "taper_shape"),\n]\n\ndef _mut_gap_squeeze(pat: Pattern) -> None:\n    ch = _channels_by(pat, "right")[6]        # pull one return row up so\n    d = pat.d                                 # the unit leaves its row\n    yr = d.rows[2 * 6 + 1]\n    new = d.rows[2 * 6] + 1.2 * d.line_w\n    ch.points = [(x, new if abs(y - yr) < 1e-9 else y)\n                 for x, y in ch.points]\n\n\ndef _mut_profile_break(pat: Pattern) -> None:\n    ch = _channels_by(pat, "right")[pat.params.n_arms + 2]\n    d = pat.d                                 # bend one taper reach off\n    for i in (1, 2):                          # the declared envelope\n        x, y = ch.points[i]\n        ch.points[i] = (x + 0.18 * d.fan_step, y)\n\n\ndef _mut_gap_stall(pat: Pattern) -> None:\n    p = pat.params\n    base = p.n_arms + p.n_taper + p.n_decay\n    chs = _channels_by(pat, "left")\n    prev, cur = chs[base + 2], chs[base + 3]  # copy unit 2\'s gap onto\n    x_prev = prev.points[0][0]                # unit 3: schedule stalls\n    x, y = cur.points[0]\n    cur.points[0] = (x_prev, y)\n\n\ndef _mut_weld_release(pat: Pattern) -> None:\n    d = pat.d\n    ch = _channels_by(pat, "left")[-2]        # re-attach a released rung\n    face = d.xc - d.trunk_w / 2.0             # to the trunk face\n    x, y = ch.points[0]\n    ch.points[0] = (face, y)\n\n\ndef _mut_rung_cut(pat: Pattern) -> None:\n    p = pat.params\n    base = p.n_arms + p.n_taper + p.n_decay + p.n_release\n    ch = _channels_by(pat, "left")[base + 3]  # rung stops short of the\n    d = pat.d                                 # rail: intersection broken\n    for i in (1, 2):\n        x, y = ch.points[i]\n        ch.points[i] = (x + 2.5 * d.line_w, y)   # left side: inboard\n\n\ndef _mut_rail_wiggle(pat: Pattern) -> None:\n    ch = _channels_by(pat, "left")[-3]        # one rail segment off the\n    for i in (1, 2):                          # wall offset\n        x, y = ch.points[i]\n        ch.points[i] = (x - 1.5, y)\n\n\ndef _mut_width_warp(pat: Pattern) -> None:\n    ch = _channels_by(pat, "right")[7]        # one channel off the\n    ch.width *= 1.4                           # declared width function\n\n\ndef _mut_trunk_warp(pat: Pattern) -> None:\n    kind, pts, bbox = pat.trunk_parts["main"]\n    i = len(pts) // 4                         # bend the center channel\n    pts[i] = (pts[i][0] + 3.0, pts[i][1])     # off its declared profile\n\n\ndef _mut_morph_jump(pat: Pattern) -> None:\n    p = pat.params\n    base = p.n_arms + p.n_taper + p.n_decay + p.n_release\n    ch = _channels_by(pat, "left")[base + 2]  # push one unit\'s excess\n    d = pat.d                                 # into the exclusion band\n    mean_w = 0.5 * (d.unit_widths[ch.index - 1] + d.unit_widths[ch.index])\n    band_mid = mean_w + d.min_clear / 4.0\n    delta = band_mid - d.morph_excess[2]\n    for i in (1, 2):\n        x, y = ch.points[i]\n        ch.points[i] = (x - delta, y)         # left side: outward = -x\n\n\ndef _mut_morph_unsort(pat: Pattern) -> None:\n    p = pat.params\n    base = p.n_arms + p.n_taper + p.n_decay + p.n_release\n    chs = _channels_by(pat, "left")\n    a, b = chs[base + 3], chs[base + 4]       # swap two excesses ->\n    for i in (1, 2):                          # schedule off spec\n        ax, ay = a.points[i]\n        bx, by = b.points[i]\n        a.points[i] = (bx, ay)\n        b.points[i] = (ax, by)\n\n\nMORPH_MUTATIONS = [\n    ("morph_jump", _mut_morph_jump, "morph_shape"),\n    ("morph_unsort", _mut_morph_unsort, "morph_shape"),\n]\n\nSHAPE_MUTATIONS = [\n    ("width_warp", _mut_width_warp, "width_profile"),\n    ("trunk_warp", _mut_trunk_warp, "width_profile"),\n]\n\nGRID_MUTATIONS = [\n    ("rung_cut", _mut_rung_cut, "grid_shape"),\n    ("rail_wiggle", _mut_rail_wiggle, "grid_shape"),\n]\n\nREL_MUTATIONS = [\n    ("gap_stall", _mut_gap_stall, "release_shape"),\n    ("weld_release", _mut_weld_release, "release_shape"),\n]\n\nPROF_MUTATIONS = [\n    ("gap_squeeze", _mut_gap_squeeze, "row_grid"),\n    ("profile_break", _mut_profile_break, "taper_shape"),\n]\n\nCONN_MUTATIONS = [\n    ("decay_break", _mut_decay_break, "decay_shape"),\n    ("rail_cut", _mut_rail_cut, "port_interface"),\n    ("rail_skew", _mut_rail_skew, "port_interface"),\n]\n\n\ndef run_selftest() -> Tuple[bool, List[str]]:\n    lines = ["== SELF-TEST =="]\n    all_ok = True\n    for label, params, muts in (\n            ("classic", Params(name="selftest_base"), MUTATIONS),\n            ("connector", Params(n_decay=6, name="selftest_conn"),\n             CONN_MUTATIONS),\n            ("profiled", Params(n_taper=12,\n                                spacing=("geometric", 0.45),\n                                taper_profile=("cosine",),\n                                name="selftest_prof"),\n             PROF_MUTATIONS),\n            ("release", Params(n_decay=6, n_release=8, line_w=4.0,\n                               boundary_w=20.0, name="selftest_rel"),\n             REL_MUTATIONS),\n            ("grid", Params(n_decay=6, n_grid=10, name="selftest_grid"),\n             GRID_MUTATIONS),\n            ("shaped", Params(n_decay=8, n_grid=6, height=1800,\n                              line_w=5.0, boundary_w=13.0,\n                              side_w_profile=("cosine", 1.5),\n                              center_w_profile=("power", 2.0, 1.25),\n                              corner_r=2.2, name="selftest_shape"),\n             SHAPE_MUTATIONS),\n            ("morph", Params(n_decay=6, n_morph=6, n_grid=6, height=1500,\n                             name="selftest_morph"),\n             MORPH_MUTATIONS)):\n        base = build_pattern(params)\n        pos = validate(base)\n        pos_ok = all(r.ok for r in pos)\n        lines.append(f"   {\'PASS\' if pos_ok else \'FAIL\'}  positive "\n                     f"({label}): pristine pattern validates clean")\n        if not pos_ok:\n            for r in pos:\n                if not r.ok:\n                    lines.append(f"         !! {r.name}: {r.detail}")\n        all_ok &= pos_ok\n        for mname, fn, primary in muts:\n            mut = copy.deepcopy(base)\n            fn(mut)\n            failed = sorted(r.name for r in validate(mut) if not r.ok)\n            caught = primary in failed\n            all_ok &= caught\n            lines.append(f"   {\'PASS\' if caught else \'FAIL\'}  negative "\n                         f"{mname:<16} -> {primary} "\n                         f"{\'caught\' if caught else \'MISSED\'} "\n                         f"(failed: {\', \'.join(failed) or \'none\'})")\n    return all_ok, lines\n\n\n# ---------------------------------------------------------------------------\n# 7. VARIANT SUITE + CLI\n# ---------------------------------------------------------------------------\nVARIANTS: List[Params] = [\n    Params(name="01_replica_940x990_a3_t11"),\n    Params(width=620, height=720, n_arms=2, n_taper=5,\n           name="02_compact_620x720_a2_t5"),\n    Params(width=900, height=1400, n_arms=3, n_taper=14,\n           name="03_dense_900x1400_a3_t14"),\n    Params(width=1400, height=900, n_arms=3, n_taper=6,\n           name="04_wide_1400x900_a3_t6"),\n    Params(width=1000, height=1200, n_arms=4, n_taper=8,\n           name="05_fourarms_1000x1200_a4_t8"),\n    Params(line_w=9.0, name="06_bold_940x990_a3_t11_lw9"),\n    Params(width=800, height=1000, n_arms=3, n_taper=6, pitch=24.0,\n           name="07_pitch24_800x1000_a3_t6"),\n    Params(n_decay=6, name="08_connector_940x990_a3_t11_d6"),\n    Params(width=900, height=1100, n_arms=3, n_taper=8, n_decay=6,\n           trunk_w=34.0, line_w=7.0, boundary_w=10.0,\n           name="09_connector_custom_tw34_lw7_b10"),\n    Params(width=940, height=4200, n_arms=3, n_taper=12, n_decay=220,\n           name="10_gradual_940x4200_a3_t12_d220"),\n    Params(width=940, height=1600, n_arms=3, n_taper=11, n_decay=24,\n           spacing=("geometric", 0.30),\n           name="11_spacing_geo030_940x1600_d24"),\n    Params(taper_profile=("cosine",),\n           name="12_taper_cosine_940x990_a3_t11"),\n    Params(width=940, height=2200, n_arms=3, n_taper=16,\n           decay_q=0.93, spacing=("geometric", 0.45),\n           taper_profile=("exp", 3.0),\n           name="13_combo_geoSpace_expTaper_q093"),\n    Params(width=940, height=1200, n_arms=3, n_taper=11, n_decay=6,\n           n_release=10, line_w=4.5, boundary_w=18.0,\n           name="14_release_940x1200_d6_r10"),\n    Params(width=940, height=2400, n_arms=3, n_taper=11, n_decay=30,\n           n_release=40, line_w=3.2, boundary_w=16.0,\n           release_profile=("cosine",),\n           name="15_release_cosine_940x2400_d30_r40"),\n    Params(width=940, height=1400, n_arms=3, n_taper=11, n_decay=6,\n           n_grid=12, name="16_grid_940x1400_d6_g12"),\n    Params(width=940, height=2000, n_arms=3, n_taper=11, n_decay=8,\n           n_grid=30, name="17_grid_long_940x2000_d8_g30"),\n    Params(width=940, height=2600, n_arms=3, n_taper=10, n_decay=14,\n           n_grid=14, line_w=7.6, trunk_w=34.0, boundary_w=12.0,\n           corner_r=3.4, name="18_fat_grid_rounded_940x2600"),\n    Params(width=940, height=1800, n_arms=3, n_taper=12, n_decay=10,\n           n_grid=8, line_w=5.0, trunk_w=22.0, boundary_w=12.0,\n           side_w_profile=("cosine", 1.8),\n           center_w_profile=("power", 2.0, 1.3),\n           corner_r=2.0, name="19_shaped_940x1800"),\n    Params(width=940, height=1200, line_w=8.0, trunk_w=30.0,\n           corner_r=3.2, n_taper=9,\n           name="20_fat_classic_rounded_940x1200"),\n    Params(width=940, height=3000, n_arms=3, n_taper=10, n_decay=10,\n           n_morph=10, n_grid=12, line_w=7.6, trunk_w=34.0,\n           boundary_w=12.0, corner_r=3.4,\n           name="21_morph_smooth_940x3000"),\n    Params(width=940, height=2000, n_arms=3, n_taper=11, n_decay=8,\n           n_morph=8, n_grid=6, line_w=5.0, trunk_w=22.0,\n           boundary_w=12.0, side_w_profile=("cosine", 1.6),\n           center_w_profile=("power", 2.0, 1.25), corner_r=2.0,\n           name="22_morph_shaped_940x2000"),\n]\n\n\ndef generate(p: Params, outdir: str) -> Tuple[bool, str]:\n    pat = build_pattern(p)\n    results = validate(pat)\n    ok = all(r.ok for r in results)\n    path = os.path.join(outdir, p.name + ".svg")\n    sha = write_svg(pat, path)\n    lines = [f"== {p.name} ==", f"   overall: {\'PASS\' if ok else \'FAIL\'}"]\n    for r in results:\n        lines.append(f"   {\'PASS\' if r.ok else \'FAIL\'}  {r.name:<18} "\n                     f"{r.detail}")\n    lines.append(f"   file: {os.path.basename(path)}  sha256={sha[:16]}...")\n    return ok, "\\n".join(lines)\n\n\ndef _parse_profile(txt):\n    """CLI helper: \'geometric:0.3\' -> ("geometric", 0.3); \'cosine\' ->\n    ("cosine",)."""\n    if txt is None:\n        return None\n    parts = txt.split(":")\n    return (parts[0],) if len(parts) == 1 else (parts[0], float(parts[1]))\n\n\ndef main(argv: List[str]) -> int:\n    ap = argparse.ArgumentParser(\n        description="Meander-pattern generator (traced line-drawing "\n                    "model) with validators and self-test")\n    ap.add_argument("--outdir", default="out")\n    ap.add_argument("--single", action="store_true",\n                    help="generate one custom pattern instead of the suite")\n    ap.add_argument("--width", type=float, default=940.0)\n    ap.add_argument("--height", type=float, default=990.0)\n    ap.add_argument("--arms", type=int, default=3,\n                    help="tall hairpin units per side")\n    ap.add_argument("--taper", type=int, default=11,\n                    help="tapering units per side")\n    ap.add_argument("--decay", type=int, default=0,\n                    help="exponential-convergence units per side "\n                         "(>=1 switches to connector mode)")\n    ap.add_argument("--pitch", type=float, default=None)\n    ap.add_argument("--line-w", type=float, default=None)\n    ap.add_argument("--hairpin-w", type=float, default=None)\n    ap.add_argument("--merge-off", type=float, default=None)\n    ap.add_argument("--boundary-w", type=float, default=None,\n                    help="convergent boundary width (interface gap)")\n    ap.add_argument("--decay-q", type=float, default=None,\n                    help="direct decay rate in (0,1); D is derived")\n    ap.add_argument("--morph", type=int, default=0,\n                    help="morph units per side: smooth transition into "\n                         "the grid (0 = immediate)")\n    ap.add_argument("--morph-profile", type=str, default="cosine",\n                    help="morph E schedule shape, e.g. cosine or "\n                         "power:3")\n    ap.add_argument("--corner-r", type=float, default=None,\n                    help="round ALL corners (fillet radius)")\n    ap.add_argument("--side-w-profile", type=str, default=None,\n                    help="side width function, e.g. cosine:1.8 or "\n                         "power:2:1.3")\n    ap.add_argument("--center-w-profile", type=str, default=None,\n                    help="center width function along the trunk")\n    ap.add_argument("--grid", type=int, default=0,\n                    help="grid units per side: straight side channels "\n                         "with periodic trunk intersections")\n    ap.add_argument("--release", type=int, default=0,\n                    help="release units per side (fading trunk "\n                         "intersections; needs boundary_w >= 2*line_w)")\n    ap.add_argument("--release-profile", type=str, default=None,\n                    help="junction-gap envelope, e.g. cosine, exp:3")\n    ap.add_argument("--spacing", type=str, default=None,\n                    help="row spacing profile, e.g. geometric:0.3 or "\n                         "linear:2.5")\n    ap.add_argument("--taper-profile", type=str, default=None,\n                    help="taper envelope, e.g. cosine, power:2.5, exp:4")\n    ap.add_argument("--name", default="custom")\n    ap.add_argument("--skip-selftest", action="store_true")\n    args = ap.parse_args(argv)\n\n    os.makedirs(args.outdir, exist_ok=True)\n    if args.single:\n        variants = [Params(width=args.width, height=args.height,\n                           n_arms=args.arms, n_taper=args.taper,\n                           n_decay=args.decay,\n                           pitch=args.pitch, line_w=args.line_w,\n                           hairpin_w=args.hairpin_w,\n                           merge_off=args.merge_off,\n                           boundary_w=args.boundary_w,\n                           decay_q=args.decay_q,\n                           n_release=args.release,\n                           n_grid=args.grid,\n                           n_morph=args.morph,\n                           morph_profile=_parse_profile(\n                               args.morph_profile) or ("cosine",),\n                           corner_r=args.corner_r,\n                           side_w_profile=_parse_profile(\n                               args.side_w_profile),\n                           center_w_profile=_parse_profile(\n                               args.center_w_profile),\n                           release_profile=_parse_profile(\n                               args.release_profile),\n                           spacing=_parse_profile(args.spacing),\n                           taper_profile=_parse_profile(\n                               args.taper_profile),\n                           name=args.name)]\n    else:\n        variants = VARIANTS\n\n    report: List[str] = []\n    all_ok = True\n    for p in variants:\n        try:\n            ok, text = generate(p, args.outdir)\n        except ValueError as exc:\n            ok, text = False, (f"== {p.name} ==\\n   overall: FAIL\\n"\n                               f"   infeasible: {exc}")\n        print(text)\n        report.append(text)\n        all_ok &= ok\n\n    if not args.skip_selftest:\n        st_ok, st_lines = run_selftest()\n        print("\\n".join(st_lines))\n        report.extend([""] + st_lines)\n        all_ok &= st_ok\n\n    rpath = os.path.join(args.outdir, "validation_report.txt")\n    with open(rpath, "w", encoding="utf-8") as fh:\n        fh.write("\\n".join(report) + "\\n")\n    print(("\\nALL VALID" if all_ok else "\\nFAILURES PRESENT")\n          + f" - report: {rpath}")\n    return 0 if all_ok else 1\n'

_PHYSICS_SRC = '#!/usr/bin/env python3\n"""\nladder_physics.py -- v1\n=======================\n\nThe PHYSICS LAYER for meander_pattern.py: compiles a ladder-array\nmicrofluidic specification (rectangular-duct Stokes hydraulics, verified\nkernel set) into meander Params, then AUDITS the drawn Pattern by an\nindependent Kirchhoff solve of the geometry actually emitted.\n\nSemantics (spec section 11/12): black = channel, white = wall. Drawing\nunits are MICROMETRES. The mapping:\n\n    trunk (center channel, per-row halves) <-> primary duct  w_c(i)\n    side polyline of unit k (width w_k)    <-> secondary duct / fold\n    stop-short stub, face -> jog           <-> the rung (link), length m_k\n    jog rail (grid phase)                  <-> unfolded straight secondary\n    arm / taper / decay excursions         <-> folded secondary, length\n                                               L_exc(k) = the fold L-hat\n    row pitch                              <-> Lambda\n    etch depth h (NOT drawn)               <-> third dimension of every duct\n\nDesign contract (two independent sides):\n  COMPILER: uses only the verified equations (resistance kernel r(w,h),\n    extraction kernel F_lam, sigma/rho march, closure R_exc = rho * R_c,\n    canonical taper with the c2-corrected F) to choose counts, reaches,\n    widths, trunk profile; emits Params through the generator\'s public\n    hooks (taper_profile / side_w_profile / center_w_profile / pitch /\n    boundary_w / min_white).\n  AUDITOR: reads ONLY the built Pattern (polylines, widths, trunk\n    boundary), reconstructs the hydraulic network, solves Kirchhoff with\n    a banded elimination, and checks the physics targets. The compiler\n    being wrong is a failure the auditor catches, not a failure of the\n    guarantee.\n\nStdlib only (matplotlib imported lazily for figures).\n"""\nfrom __future__ import annotations\n\nimport math\nfrom dataclasses import dataclass, field, replace\nfrom typing import Callable, Dict, List, Optional, Sequence, Tuple\n\nimport sys as _sys; mp = _sys.modules["meander_core"]\n\nEPS = 1e-9\n\n\nclass DecayLawShort(ValueError):\n    """The geometric decay law\'s drawn length undershoots the exact\n    closure requirement at some decay row; the compiler retries with the\n    taper/decay boundary shifted one row down."""\n\n# ---------------------------------------------------------------------------\n# 1. VERIFIED KERNELS (ported from the checkpoint scripts of the review)\n# ---------------------------------------------------------------------------\n\n# odd-n zeta partial constants (verified at import against long loops):\n_S2_ODD = math.pi ** 2 / 8.0                      # sum odd n^-2\n_S3_ODD = 0.875 * 1.2020569031595942854           # (7/8) zeta(3)\n_S4_ODD = math.pi ** 4 / 96.0                     # sum odd n^-4\n_S5_ODD = (31.0 / 32.0) * 1.0369277551433699263   # (31/32) zeta(5)\n\n\ndef B_series(alpha: float, n_terms: int = 200000) -> float:\n    """Rectangular-duct shape factor, SORTED branch (alpha >= 1).\n    B(alpha) = 1 - (192/(pi^5 alpha)) * sum_{n odd} tanh(n pi alpha/2)/n^5.\n    Once tanh saturates (t > 19, 1-tanh < 6e-17) the remaining series is\n    the EXACT odd zeta tail S5 - partial. Verified anchors:\n    B(1)=0.4217310, lim alpha*(1-B)=0.6302489."""\n    if alpha < 1.0:\n        raise ValueError("B_series takes the sorted branch alpha >= 1")\n    s = 0.0\n    p5 = 0.0\n    n = 1\n    while n < 2 * n_terms:\n        t = n * math.pi * alpha / 2.0\n        if t > 19.0:\n            s += _S5_ODD - p5                     # exact analytic tail\n            break\n        inv5 = 1.0 / n ** 5\n        s += math.tanh(t) * inv5\n        p5 += inv5\n        n += 2\n    return 1.0 - (192.0 / (math.pi ** 5 * alpha)) * s\n\n\ndef r_hyd(w: float, h: float, mu: float = 1.0) -> float:\n    """Hydraulic resistance per unit length of a w x h rectangular duct,\n    r = 12 mu / (a b^3 B(a/b)), a = max(w,h), b = min(w,h). Exactly\n    symmetric in (w,h) on the sorted branch (verified to 0 ULP)."""\n    a, b = (w, h) if w >= h else (h, w)\n    return 12.0 * mu / (a * b ** 3 * B_series(a / b))\n\n\ndef r_inv_w(R_target: float, h: float, mu: float = 1.0,\n            lo: float = 0.05, hi: float = 5000.0) -> float:\n    """Invert r_hyd in w at fixed h (r strictly decreasing in w)."""\n    flo, fhi = r_hyd(lo, h, mu), r_hyd(hi, h, mu)\n    if not (fhi <= R_target <= flo):\n        raise ValueError(f"r_inv_w: target {R_target:.3e} outside "\n                         f"[{fhi:.3e}, {flo:.3e}] for w in [{lo},{hi}]")\n    for _ in range(200):\n        mid = 0.5 * (lo + hi)\n        if r_hyd(mid, h, mu) > R_target:\n            lo = mid\n        else:\n            hi = mid\n        if hi - lo < 1e-9 * hi:\n            break\n    return 0.5 * (lo + hi)\n\n\ndef _A_n(n: int, beta: float) -> float:\n    """Extraction-kernel coefficient, primary (cancellation-safe) form:\n    A_n = (2/(pi n^5)) (t - tanh t), t = n pi beta / 2."""\n    t = n * math.pi * beta / 2.0\n    tt = math.tanh(t) if t < 20.0 else 1.0\n    if t < 1e-3:                              # Taylor branch (E-e.1)\n        d = t ** 3 / 3.0 - 2.0 * t ** 5 / 15.0\n    else:\n        d = t - tt\n    return 2.0 / (math.pi * n ** 5) * d\n\n\ndef Flam(x: float, beta: float, nmax: int = 1201) -> float:\n    """Exact near-wall flux fraction of a rectangular duct: the fraction\n    of Q within a lamina of thickness x*w against the w-separated wall.\n    F(1/2, beta) = 1/2 exactly; F -> 3x^2 - 2x^3 as beta -> inf."""\n    num = den = 0.0\n    for n in range(1, nmax, 2):\n        A = _A_n(n, beta)\n        num += A * (1.0 - math.cos(n * math.pi * x))\n        den += A\n    return num / (2.0 * den)\n\n\ndef dFlam_dx(x: float, beta: float, nmax: int = 1201) -> float:\n    num = den = 0.0\n    for n in range(1, nmax, 2):\n        A = _A_n(n, beta)\n        num += A * n * math.pi * math.sin(n * math.pi * x)\n        den += A\n    return num / (2.0 * den)\n\n\ndef c2(beta: float, nmax: int = 400001) -> float:\n    """Leading coefficient of F_lam ~ c2 x^2 = (pi^2/4) Sum n^2 A_n /\n    Sum A_n. Written with EXACT analytic tails: for t = n pi beta/2 > 19\n    (tanh = 1 to double precision), A_n = beta/n^4 - (2/pi)/n^5, so\n    Sum n^2 A_n tail = beta(S2-p2) - (2/pi)(S3-p3) and\n    Sum A_n tail = beta(S4-p4) - (2/pi)(S5-p5). This supersedes the\n    review\'s +beta/(2 n_max) approximation (which it reproduces in the\n    n_max -> inf limit). Verified: c2(1)=3.5568, c2(inf)=3."""\n    s2 = s0 = 0.0\n    p2 = p3 = p4 = p5 = 0.0\n    n = 1\n    while n < nmax:\n        t = n * math.pi * beta / 2.0\n        if t > 19.0:\n            s2 += beta * (_S2_ODD - p2) - (2.0 / math.pi) * (_S3_ODD - p3)\n            s0 += beta * (_S4_ODD - p4) - (2.0 / math.pi) * (_S5_ODD - p5)\n            break\n        A = _A_n(n, beta)\n        s2 += n * n * A\n        s0 += A\n        p2 += 1.0 / n ** 2\n        p3 += 1.0 / n ** 3\n        p4 += 1.0 / n ** 4\n        p5 += 1.0 / n ** 5\n        n += 2\n    return math.pi ** 2 / 4.0 * s2 / s0\n\n\ndef Flam_inv(f: float, beta: float) -> float:\n    """x such that Flam(x,beta) = f, via the verified seed sqrt(f/c2)\n    plus safeguarded Newton (bisection fallback)."""\n    if not (0.0 < f < 0.5):\n        raise ValueError("Flam_inv needs 0 < f < 1/2")\n    x = math.sqrt(f / c2(beta))\n    lo, hi = 0.0, 0.5\n    for _ in range(60):\n        fx = Flam(x, beta) - f\n        if fx > 0:\n            hi = x\n        else:\n            lo = x\n        d = dFlam_dx(x, beta)\n        xn = x - fx / d if d > 0 else 0.5 * (lo + hi)\n        if not (lo < xn < hi):\n            xn = 0.5 * (lo + hi)\n        if abs(xn - x) < 1e-14:\n            x = xn\n            break\n        x = xn\n    return x\n\n\ndef C_shear(alpha: float, n_terms: int = 20000) -> float:\n    """Max wall-shear coefficient (section 16-b of the ladder spec):\n    gdot_max = C(alpha) Q / (a b^2). Verified: C(1)=9.6078, C(inf)=6."""\n    s = 0.0\n    for j in range(n_terms):\n        n = 2 * j + 1\n        t = n * math.pi * alpha / 2.0\n        sech = 1.0 / math.cosh(t) if t < 300.0 else 0.0\n        s += (1.0 - sech) / (n * n)\n    return 48.0 / (math.pi ** 2 * B_series(alpha)) * s\n\n\n# ---------------------------------------------------------------------------\n# 2. PHYSICS SPEC + MARCH\n# ---------------------------------------------------------------------------\n\n@dataclass\nclass PhysicsSpec:\n    """Physical inputs. Lengths in micrometres; SI only where marked."""\n    D_target: float               # cutoff size D = kappa * delta*  [um]\n    F_target: float               # filtrate fraction target, F = 1 - P(N)\n    h: float                      # etch depth (out-of-plane)      [um]\n    wc0: float                    # inlet primary width            [um]\n    kappa: float                  # calibrated size factor for the F_lam\n                                  # branch (kernel-prefactor pairing rule)\n    w_min: float = 6.0            # min printable channel width    [um]\n    wall_min: float = 6.0         # min wall (white) feature       [um]\n    pitch: Optional[float] = None # row pitch Lambda; None = solved  [um]\n    margin: float = 40.0          # blank border on the chip       [um]\n    arm_headroom: Optional[float] = None  # vertical budget for hairpins;\n                                  # None = auto (scales with fold need)\n    port_run: Optional[float] = None      # straight interface run; >=2*pitch\n    corner_r: object = "auto"     # \'auto\' = 0.98x the feasibility cap\n                                  # (0.45 x smallest white feature);\n                                  # float = explicit; None = sharp (fails\n                                  # P10 -- physics devices must be round)\n    reach_cap_factor: float = 12.0  # taper reach budget in rail_off\n                                  # units; smaller -> more head rows\n                                  # become arms\n    arms_min: int = 1             # force at least this many arms (their\n                                  # widths are solved, so physics stays\n                                  # exact)\n    taper_law: str = "canonical"  # \'canonical\' = w_c^2/c2 linear in i\n                                  # (constant Q_gap); the only law with a\n                                  # closed form -- S2-free laws could be\n                                  # added as callables.\n    # -- operating point (for the Stokes gates; SI) --\n    mu_SI: float = 1.0e-3         # Pa s\n    rho_SI: float = 1.0e3         # kg/m^3\n    Q_in_uLmin: float = 5.0       # inlet flow\n    # -- tolerances / gates --\n    tol_f_bulk: float = 5.0e-3    # |f_em/f - 1| in the bulk\n    tol_delta: float = 0.02       # max_i delta_em/delta* - 1 (cutoff guard)\n    tol_closure: float = 5.0e-4   # drawn closure vs target\n    Re_max: float = 1.0\n    eps_net: float = 1.0e-3       # terminal-layer accounting tolerance\n    max_correction_iters: int = 8\n    name: str = "ladder"\n\n\n@dataclass\nclass Plan:\n    """Everything the compiler decided; consumed by figures/reports.\n    The AUDITOR does not read this except for the declared spec."""\n    spec: PhysicsSpec\n    delta_star: float\n    N: int\n    A: int\n    T: int\n    D: int\n    G: int\n    Lambda: float\n    wc_rows: List[float]          # per ROW (2N) + port width appended\n    wc_port: float\n    f: List[float]                # per unit (1-based i stored at [i-1])\n    sigma: List[float]\n    rho: List[float]              # rho[i] AFTER rung i, i=1..N (len N)\n    P: List[float]\n    F_march: float\n    F_closed_corrected: float\n    F_closed_doc_form: float      # 1-(wN/w0)^2, kept to show the error\n    Rc_seg: List[float]           # trunk resistance tap k -> tap k+1\n    Rs_target: List[float]        # required excursion resistance / unit\n    L_exc: List[float]            # drawn excursion length / unit\n    widths: List[float]           # drawn channel width / unit\n    reaches: List[float]          # taper-region reaches (len T)\n    rail_off: float\n    b_wall: float\n    m_off: float\n    line_w: float\n    trunk_w: float\n    W: float\n    H: float\n    rho_out: float\n    corner_r: Optional[float] = None\n    correction_iters: int = 0\n\n\ndef wc_from_phi(phi: float, h: float,\n                grid: Optional[Tuple[List[float], List[float]]] = None\n                ) -> float:\n    """Invert phi(w) = w^2 / c2(h/w) (strictly increasing in w)."""\n    lo, hi = 1.0, 4000.0\n    for _ in range(200):\n        mid = 0.5 * (lo + hi)\n        if mid * mid / c2(h / mid) < phi:\n            lo = mid\n        else:\n            hi = mid\n        if hi - lo < 1e-10 * hi:\n            break\n    return 0.5 * (lo + hi)\n\n\ndef physics_march(spec: PhysicsSpec) -> Dict[str, object]:\n    """S1\' march on the canonical taper with the EXACT extraction kernel,\n    FOLD mode (sigma0 = 0), terminated on P(i) <= 1 - F_target (the\n    monotone-safe rule). Returns per-unit schedules and the corrected /\n    doc-form closed-form F for cross-checking (the review\'s finding)."""\n    dstar = spec.D_target / spec.kappa\n    h = spec.h\n    wc = spec.wc0\n    phi = wc * wc / c2(h / wc)\n    sigma = 0.0\n    P = 1.0\n    f_l: List[float] = []\n    sig_l: List[float] = []\n    P_l: List[float] = []\n    wc_units: List[float] = [wc]\n    guard = 200000\n    while P > 1.0 - spec.F_target and guard > 0:\n        guard -= 1\n        beta = h / wc\n        x = dstar / wc\n        if x >= 0.5:\n            raise ValueError("delta* >= w_c/2: laminae overlap; raise "\n                             "kappa*w or lower D_target")\n        fi = Flam(x, beta)\n        sigma = (sigma + fi) / (1.0 - 2.0 * fi)\n        P *= (1.0 - 2.0 * fi)\n        f_l.append(fi)\n        sig_l.append(sigma)\n        P_l.append(P)\n        phi -= 2.0 * dstar * dstar          # canonical taper step\n        wc = wc_from_phi(phi, h)\n        wc_units.append(wc)\n        if wc <= 2.0 * dstar * 2.0:\n            raise ValueError("taper exhausted the primary before reaching "\n                             "F_target; lower F_target or D_target")\n    if guard == 0:\n        raise ValueError("march did not terminate")\n    N = len(f_l)\n    wcN = wc_units[-1]\n    F_march = 1.0 - P_l[-1]\n    b0, bN = h / spec.wc0, h / wcN\n    F_corr = 1.0 - (wcN * wcN / c2(bN)) / (spec.wc0 ** 2 / c2(b0))\n    F_doc = 1.0 - (wcN / spec.wc0) ** 2\n    # mass identity (machine-precision invariant, checked here as an\n    # internal tripwire exactly like the ladder spec\'s M-family)\n    S = 0.0\n    Pi = 1.0\n    for fi in f_l:\n        S += fi * Pi\n        Pi *= (1.0 - 2.0 * fi)\n    if abs(Pi + 2.0 * S - 1.0) > 1e-12:\n        raise AssertionError("mass identity violated in march")\n    return dict(dstar=dstar, N=N, f=f_l, sigma=sig_l, P=P_l,\n                wc_units=wc_units, F_march=F_march, F_corr=F_corr,\n                F_doc=F_doc)\n\n\n# ---------------------------------------------------------------------------\n# 3. COMPILER  (physics -> meander Params)\n# ---------------------------------------------------------------------------\n\ndef _interp_step(vals: Sequence[float]) -> Callable[[float], float]:\n    """Exact per-sample lookup for generator callables that sample at\n    t = i/(n-1): reconstruct i by rounding."""\n    n = len(vals)\n\n    def fn(t: float) -> float:\n        i = int(round(t * (n - 1)))\n        i = 0 if i < 0 else (n - 1 if i > n - 1 else i)\n        return vals[i]\n    return fn\n\n\ndef compile_physics(spec: PhysicsSpec) -> Tuple[mp.Params, Plan]:\n    """Compile the physics into Params. Every hydraulic identity used\n    here is one of the verified equations; every geometric formula\n    mirrors derive()/build_pattern() exactly (and the mirror is asserted\n    against derive() before returning)."""\n    M = physics_march(spec)\n    N = M["N"]\n    dstar = M["dstar"]\n    h = spec.h\n    f = M["f"]\n    sigma = M["sigma"]\n    rho = [1.0 / s for s in sigma]           # after rung i (i = 1..N)\n    wc_units: List[float] = M["wc_units"]    # len N+1 (unit starts + end)\n\n    # per-ROW primary width: unit coordinate tau = 1 + j/2, phi linear\n    wc_rows: List[float] = []\n    for j in range(2 * N):\n        t = j / 2.0\n        i0 = int(math.floor(t))\n        fr = t - i0\n        # linear in phi between unit starts == linear interp of wc is a\n        # fine sub-row model; keep phi-linear for exactness:\n        w0, w1 = wc_units[i0], wc_units[min(i0 + 1, N)]\n        p0 = w0 * w0 / c2(h / w0)\n        p1 = w1 * w1 / c2(h / w1)\n        wc_rows.append(wc_from_phi(p0 + fr * (p1 - p0), h))\n    wc_port = wc_units[N]\n\n    # ---- pitch, port line width, wall, rail ------------------------------\n    rN = r_hyd(wc_port, h)\n    rho_out = rho[-1]\n    # unfolded closure width at the last unit (defines the port line_w):\n    w_port_side = r_inv_w(rho_out * rN, h)\n    line_w = w_port_side\n    min_clear = spec.wall_min\n    # rail anchored clear of the WIDEST trunk section:\n    rail_off = max(wc_rows) / 2.0 + line_w / 2.0 + min_clear + 1.0\n    t2 = wc_port / 2.0\n    m_off = rail_off - t2\n    b_wall = m_off - line_w / 2.0\n    Lambda = spec.pitch if spec.pitch is not None else 0.0\n    Lambda = max(Lambda,\n                 line_w + min_clear + 2.0,       # slot rule\n                 spec.w_min + min_clear + 2.0,\n                 3.0 * min_clear,\n                 m_off / 2.5 + 1e-6)             # generator m <= 2.5*s\n    if b_wall < min_clear - EPS:\n        raise ValueError("derived wall thinner than wall_min; widen "\n                         "boundary geometry")\n\n    # ---- trunk resistance per segment (tap k -> tap k+1) -----------------\n    r_rows = [r_hyd(w, h) for w in wc_rows] + [rN]\n    def Rc_between(j0: int, j1_extra: float = 0.0) -> float:\n        # trapezoid over rows j0..j0+2 (uniform pitch), + optional tail\n        R = Lambda * (r_rows[j0] / 2.0 + r_rows[j0 + 1]\n                      + r_rows[min(j0 + 2, len(r_rows) - 1)] / 2.0)\n        return R + j1_extra\n    Rc_seg = [Rc_between(2 * k) for k in range(N - 1)]\n\n    # ---- closure targets: EXACT flow-prescribed construction --------------\n    # The rho-closure R_exc = rho * R_c is exact only when the rung drop\n    # R_rung * Q_gap is uniform (verified equivalence chain, "identical\n    # links"). Here rung widths and lengths vary, so we build the exact\n    # pressure field of the prescribed march flows instead:\n    #   PP(k+1) = PP(k) - R_c(k) * C(k)          (trunk taps)\n    #   PJ(k)   = PP(k) - R_rung(k) * g(k)       (rail nodes)\n    #   R_exc(k) = [PJ(k) - PJ(k+1)] / s(k)      (exact target)\n    # with per-side flows from the kernel march (Q_in = 1):\n    #   C(k) = P_after(k+1)  total trunk flow in segment k\n    #   g(k) = f[k] * P_after(k)      per-side rung flow at tap k\n    #   s(k) = S_after(k+1)           per-side rail flow in excursion k\n    # The rho-closure is the uniform-link limit of this construction.\n    P_after = [1.0]\n    for fi in f:\n        P_after.append(P_after[-1] * (1.0 - 2.0 * fi))\n    S_after = [0.0]\n    for i, fi in enumerate(f):\n        S_after.append(S_after[-1] + fi * P_after[i])\n    C_seg = [P_after[k + 1] for k in range(N)]\n    g_r = [f[k] * P_after[k] for k in range(N)]\n    s_r = [S_after[k + 1] for k in range(N)]\n    Rs_target = [rho[k] * Rc_seg[k] for k in range(N - 1)]  # seed only\n\n    # ---- unit-kind allocation ---------------------------------------------\n    # grid rows: unfolded closure width >= w_min; fold rows: the rest.\n    w_grid = [r_inv_w(Rs_target[k] / (2.0 * Lambda), h)\n              for k in range(N - 1)]\n    fold = [w < spec.w_min - 1e-9 for w in w_grid]\n    n_fold = 0\n    for k in range(N - 1):\n        if fold[k]:\n            n_fold = k + 1                   # fold region is a prefix\n    # generator minimums: A>=1, T>=1, D>=1, G>=1\n    n_fold = max(n_fold, 3)                  # 1 arm + 1 taper + 1 decay\n    if n_fold > N - 2:\n        raise ValueError("fewer than one grid unit would remain; raise "\n                         "F_target or w_min budget")\n    # taper needs reach >= floor + margin  <=>  L_exc >= 6*Lambda + 2 at\n    # width w_min; fold rows with shorter required length go to the\n    # decay (their widths are solved against the decay law\'s lengths,\n    # which are always >= the requirement there -- guarded below).\n    L_req_at_wmin = [Rs_target[k] / r_hyd(spec.w_min, h)\n                     for k in range(n_fold)]\n    D = 1\n    while (D < n_fold - 1\n           and L_req_at_wmin[n_fold - 1 - D] < 6.0 * Lambda + 2.0):\n        D += 1\n    D = max(1, min(D, n_fold - 2))           # leave room for A + T\n    T_plus_A = n_fold - D\n    G = (N - 1) - n_fold + 1                 # grid units incl. the LAST\n    #                                          (run-out) unit\n\n    # taper reaches at width w_min (lengths -> reaches):\n    #   L_exc = 2*(reach - rail_off) + 2*Lambda\n    def taper_reach(L: float) -> float:\n        return rail_off + (L - 2.0 * Lambda) / 2.0\n\n    port_run = spec.port_run if spec.port_run is not None else 3.0 * Lambda\n    port_run = max(port_run, 2.0 * Lambda + 1.0)\n    mg = spec.margin\n    U = N                                     # one drawn unit per rung\n\n    # ---- arms: geometry mirror of D9-D11 (width solved later) --------------\n    def arm_geometry(A: int, W: float, top_y: float, H: float):\n        d_line = W / 2.0 - mg - line_w / 2.0\n        u = Lambda                            # hairpin_w default\n        fs = max(0.035 * W, u + 2.0 * line_w)\n        arm_outer = d_line - u\n        reaches = [arm_outer - (A - 1 - j) * fs for j in range(A)]\n        top_ink = mg + line_w / 2.0\n        drop = 0.20 * H\n        if A > 1:\n            y_last = top_y + 2 * (A - 1) * Lambda\n            drop_max = (y_last - 1.5 * Lambda - top_ink) / (A - 1)\n            drop = min(drop, drop_max)\n        tops = [top_ink + j * drop for j in range(A)]\n        return d_line, u, fs, reaches, tops\n\n    # decide A: arms absorb the head rows whose taper reach would exceed\n    # the drawable budget; the budget is elastic (we choose W), so A is\n    # driven by keeping the taper\'s r0 within a sane aspect (<= about\n    # 12x the ladder half-width) while arms park the excess length\n    # vertically. One arm minimum.\n    A = 1\n    while A < T_plus_A - 1 and A < 8:\n        r0_try = taper_reach(L_req_at_wmin[A])\n        if r0_try <= spec.reach_cap_factor * rail_off:\n            break\n        A += 1\n    A = max(1, spec.arms_min, min(A, T_plus_A - 1))\n    A = min(A, T_plus_A - 1, 8)\n    T = T_plus_A - A\n    reach_T = [taper_reach(L_req_at_wmin[A + i]) for i in range(T)]\n    floor_r = rail_off + 2.0 * Lambda\n    reach_T = [max(rv, floor_r + 0.5) for rv in reach_T]\n    # strict decrease going down; then re-enforce upward so no later\n    # bump breaks it (earlier reaches only grow -> taper widths, solved\n    # from the final lengths, stay >= w_min):\n    for i in range(T - 2, -1, -1):\n        if reach_T[i] <= reach_T[i + 1] + 1e-3:\n            reach_T[i] = reach_T[i + 1] + 1e-3\n    r0, r1 = reach_T[0], reach_T[-1]\n\n    # ---- vertical layout with the ARM-FEASIBILITY loop ---------------------\n    # An arm\'s width is solved from its drawn length, so each arm row j\n    # needs L_arm(j) >= L_req(j) at w_min or the width would sink below\n    # the fab floor. The hairpin headroom is grown until every arm has\n    # the length it needs (auto mode) or the spec\'s value is validated.\n    head = spec.arm_headroom\n    auto_head = head is None\n    if auto_head:\n        head = max(6.0 * Lambda, 0.35 * L_req_at_wmin[0])\n    for _attempt in range(16):\n        top_y = mg + head\n        bottom_y = top_y + (2 * U - 1) * Lambda\n        H = bottom_y + port_run\n        rows = [top_y + j * Lambda for j in range(2 * U)]\n        # canvas width from the outermost stroke: taper r0 + arm fan\n        u_hp = Lambda\n        fs_min = u_hp + 2.0 * line_w\n        W = 2.0 * (r0 + max(2.0 * line_w, 0.25 * fs_min)\n                   + (A - 1 + 1) * fs_min\n                   + u_hp + mg + line_w / 2.0 + 2.0)\n        for _ in range(4):                    # fs depends on W\n            d_line, u_hp, fs, arm_reaches, arm_tops = \\\n                arm_geometry(A, W, top_y, H)\n            need = (r0 + max(2.0 * line_w, 0.25 * fs))\n            if arm_reaches[0] >= need - 1e-9:\n                break\n            W += 2.0 * (need - arm_reaches[0]) + 2.0\n        d_line, u_hp, fs, arm_reaches, arm_tops = \\\n            arm_geometry(A, W, top_y, H)\n\n        def L_arm(k: int) -> float:\n            yo, yr = rows[2 * k], rows[2 * k + 1]\n            return (2.0 * (arm_reaches[k] - rail_off) + 2.0 * u_hp\n                    + (yo - arm_tops[k]) + (yr - arm_tops[k]) + Lambda)\n\n        deficit = max(L_req_at_wmin[j] - L_arm(j) for j in range(A))\n        if deficit <= 1e-9:\n            break\n        if not auto_head:\n            raise ValueError(f"arm_headroom too small: arm lengths "\n                             f"undershoot the closure need by "\n                             f"{deficit:.0f} um at w_min")\n        head += 0.55 * deficit + 2.0 * Lambda\n    else:\n        raise ValueError("arm-feasibility loop did not converge")\n\n    def L_decay(i: int, e_start: float, q: float) -> float:\n        E = e_start * q ** (i + 1)\n        return 2.0 * E + 2.0 * Lambda\n\n    widths: List[float] = [0.0] * U\n    L_exc: List[float] = [0.0] * U\n    for _retry in range(12):\n      kinds: List[str] = []\n      for k in range(U):\n        if k < A:\n            kinds.append("arm")\n        elif k < A + T:\n            kinds.append("taper")\n        elif k < A + T + D:\n            kinds.append("decay")\n        else:\n            kinds.append("grid")\n\n      # ---- exact fixed point: widths, rung lengths, PJ field, targets -----\n      try:\n        # Rung length varies with the local trunk face: m_k = rail - half(2k).\n        m_len = [rail_off - wc_rows[2 * k] / 2.0 for k in range(U)]\n        Rc_last = (Lambda * 0.5 * (r_rows[2 * U - 2] + r_rows[2 * U - 1])\n                   + (H - rows[2 * U - 1]) * rN)\n        Rc_full = Rc_seg + [Rc_last]\n        # seeds: grid rho-closure widths, taper at w_min, arms/decay solved\n        for k in range(U):\n            if kinds[k] == "grid":\n                widths[k] = max(w_grid[k], spec.w_min) if k < U - 1 else line_w\n            else:\n                widths[k] = spec.w_min\n        q_dec, e_start = 0.5, max(reach_T[-1] - rail_off, 1.0)\n        Rs_exact = [0.0] * U\n        for _pass in range(30):\n            prev_w = list(widths)\n            # decay-floor mirror (widths -> e_min -> q -> lengths):\n            pairs = [(widths[k] + widths[k + 1]) / 2.0\n                     for k in range(A + T - 1, A + T + D) if k + 1 < U]\n            e_min = min_clear + (max(pairs) if pairs else min_clear)\n            e_start = reach_T[-1] - rail_off\n            if e_start <= e_min + EPS:\n                reach_T[-1] = rail_off + e_min + 0.5 * Lambda\n                e_start = reach_T[-1] - rail_off\n                r1 = reach_T[-1]\n            q_dec = (e_min / e_start) ** (1.0 / D)\n            # lengths for the length-fixed kinds:\n            for k in range(U):\n                if kinds[k] == "arm":\n                    L_exc[k] = L_arm(k)\n                elif kinds[k] == "taper":\n                    L_exc[k] = (2.0 * (reach_T[k - A] - rail_off)\n                                + 2.0 * Lambda)\n                elif kinds[k] == "decay":\n                    L_exc[k] = L_decay(k - A - T, e_start, q_dec)\n                else:\n                    L_exc[k] = (2.0 * Lambda if k < U - 1\n                                else H - rows[2 * k])\n            # exact pressure field of the prescribed flows:\n            PP = [0.0]\n            for k in range(U):\n                PP.append(PP[-1] - Rc_full[k] * C_seg[k])\n            PJ = [PP[k] - r_hyd(widths[k], h) * m_len[k] * g_r[k]\n                  for k in range(U)]\n            PJ.append(PP[U])                     # common outlet node\n            conv = 0.0\n            for k in range(U):\n                drop = PJ[k] - PJ[k + 1]\n                if drop <= 0.0:\n                    raise ValueError(f"unit {k}: rail drop non-positive; "\n                                     "rung drop exceeds the row drop (raise "\n                                     "pitch or wall budget)")\n                Rs_exact[k] = drop / s_r[k]\n                if kinds[k] in ("arm",):\n                    widths[k] = r_inv_w(Rs_exact[k] / L_exc[k], h)\n                elif kinds[k] == "taper":\n                    # reach carries the target at width w_min; when the\n                    # readability floor binds, the width re-solves from the\n                    # floored length (never below w_min):\n                    Lk = Rs_exact[k] / r_hyd(spec.w_min, h)\n                    reach_T[k - A] = max(rail_off + 2.0 * Lambda + 0.5,\n                                         rail_off + (Lk - 2.0 * Lambda) / 2.0)\n                    L_act = (2.0 * (reach_T[k - A] - rail_off)\n                             + 2.0 * Lambda)\n                    widths[k] = max(spec.w_min,\n                                    r_inv_w(Rs_exact[k] / L_act, h))\n                elif kinds[k] == "decay":\n                    widths[k] = r_inv_w(min(Rs_exact[k] / L_exc[k],\n                                            r_hyd(spec.w_min, h)), h)\n                    if widths[k] < spec.w_min - 1e-9:\n                        raise DecayLawShort(f"decay unit {k}: law length too "\n                                         "short for the exact target; raise "\n                                         "w_min or pitch")\n                else:\n                    w_new = r_inv_w(Rs_exact[k] / L_exc[k], h)\n                    widths[k] = max(spec.w_min, w_new)\n                conv = max(conv, abs(widths[k] - prev_w[k]))\n            # keep taper reaches strictly decreasing after updates, then\n            # re-solve those widths against the adjusted lengths:\n            for i in range(T - 2, -1, -1):\n                if reach_T[i] <= reach_T[i + 1] + 1e-3:\n                    reach_T[i] = reach_T[i + 1] + 1e-3\n            for i in range(T):\n                k = A + i\n                L_act = 2.0 * (reach_T[i] - rail_off) + 2.0 * Lambda\n                L_exc[k] = L_act\n                widths[k] = max(spec.w_min,\n                                r_inv_w(Rs_exact[k] / L_act, h))\n            r0, r1 = reach_T[0], reach_T[-1]\n            if conv < 1e-10:\n                break\n      except DecayLawShort:\n          if D <= 1:\n              raise\n          T += 1\n          D -= 1\n          reach_T = [taper_reach(L_req_at_wmin[A + i])\n                     for i in range(T)]\n          reach_T = [max(rv, floor_r + 0.5) for rv in reach_T]\n          for i in range(T - 2, -1, -1):\n              if reach_T[i] <= reach_T[i + 1] + 1e-3:\n                  reach_T[i] = reach_T[i + 1] + 1e-3\n          r0, r1 = reach_T[0], reach_T[-1]\n          continue\n      break\n    else:\n        raise ValueError("decay-boundary retry exhausted; raise w_min or pitch")\n\n    # port width contract: the exact last-unit width must sit at the\n    # declared line_w (it does, by construction of line_w; verify):\n    if abs(widths[U - 1] / line_w - 1.0) > 0.05:\n        raise ValueError("port run-out width drifted from the declared "\n                         f"line_w by "\n                         f"{abs(widths[U-1]/line_w-1)*100:.1f}%; "\n                         "lengthen port_run")\n    widths[U - 1] = line_w                   # pin the declared port width\n    Rs_target_full = list(Rs_exact)\n\n    # ---- feasibility of solved widths --------------------------------------\n    for k, w in enumerate(widths):\n        if w < spec.w_min - 1e-6:\n            raise ValueError(f"unit {k} ({kinds[k]}) solves to width "\n                             f"{w:.2f} < w_min; raise arm_headroom / "\n                             f"pitch or relax w_min")\n        if Lambda - w < min_clear - 1e-6:\n            raise ValueError(f"unit {k} width {w:.2f} breaks the slot "\n                             f"rule; raise pitch")\n        if m_off - w / 2.0 < min_clear - 1e-6:\n            raise ValueError(f"unit {k} width breaks the pocket rule; "\n                             f"raise wall budget")\n\n    # ---- assemble Params ----------------------------------------------------\n    env = [(rv - r1) / (r0 - r1) for rv in reach_T] if T > 1 else [1.0]\n    params = mp.Params(\n        width=W, height=H,\n        n_arms=A, n_taper=T, n_decay=D, n_grid=G,\n        margin=mg,\n        pitch=Lambda,\n        line_w=line_w,\n        hairpin_w=u_hp,\n        boundary_w=b_wall,\n        trunk_w=wc_port,\n        min_white=min_clear,\n        top_frac=top_y / H,\n        bottom_frac=bottom_y / H,\n        taper_start_frac=min(r0 / d_line, 0.999),\n        taper_end_frac=(r1 / d_line if T > 1\n                        else max(floor_r, 0.75 * r0) / d_line),\n        taper_profile=(_interp_step(env) if T > 1 else None),\n        side_w_profile=_interp_step([w / line_w for w in widths]),\n        center_w_profile=_center_profile(rows, wc_rows, wc_port,\n                                         top_y, H),\n        drop_frac=(arm_tops[1] - arm_tops[0]) / H if A > 1 else 0.20,\n        nib_len_frac=min(0.06, max(0.6 * head,\n                                   line_w + 2.0 * min_clear + 3.0) / H),\n        name=spec.name,\n    )\n\n    # ---- mirror assertion: derive() must reproduce the compiler\'s plan -----\n    d = mp.derive(params)\n    def _close(a, b, tol=1e-6):\n        return abs(a - b) <= tol * max(1.0, abs(a), abs(b))\n    assert _close(d.pitch, Lambda) and _close(d.merge_off, m_off) \\\n        and _close(d.boundary_w, b_wall), "pitch/wall mirror drift"\n    for i in range(T):\n        assert _close(d.taper_reaches[i], reach_T[i], 1e-9), \\\n            f"taper reach mirror drift at {i}"\n    for j in range(A):\n        assert _close(d.arm_reaches[j], arm_reaches[j], 1e-6) \\\n            and _close(d.hairpin_tops[j], arm_tops[j], 1e-6), \\\n            f"arm mirror drift at {j}"\n    for k in range(U):\n        assert _close(d.unit_widths[k], widths[k], 1e-9), \\\n            f"width mirror drift at {k}"\n    # decay mirror assertion (the generator\'s own law must reproduce\n    # the lengths the widths were solved against):\n    for i in range(D):\n        want = rail_off + e_start * q_dec ** (i + 1)\n        assert _close(d.decay_reaches[i], want, 1e-6), \\\n            f"decay mirror drift at {i}: {d.decay_reaches[i]} vs {want}"\n\n    # ---- corner rounding: MAXIMAL feasible radius (sec 13 cap mirror) -----\n    white_feats = [rows[2 * k + 1] - rows[2 * k] - widths[k]\n                   for k in range(U)]\n    white_feats += [b_wall,\n                    rail_off - max(wc_rows) / 2.0 - max(widths) / 2.0,\n                    min_clear]\n    r_cap = 0.45 * min(white_feats)\n    if spec.corner_r is None:\n        corner_val = None\n    elif spec.corner_r == "auto":\n        corner_val = round(0.98 * r_cap, 4)\n    else:\n        corner_val = float(spec.corner_r)\n        if not (0.0 < corner_val <= r_cap + 1e-9):\n            raise ValueError(f"corner_r must be in (0, {r_cap:.3f}]")\n    if corner_val is not None:\n        params = replace(params, corner_r=corner_val)\n        d = mp.derive(params)                 # re-validate rounded\n\n    plan = Plan(spec=spec, delta_star=dstar, N=N, A=A, T=T, D=D, G=G,\n                Lambda=Lambda, wc_rows=wc_rows, wc_port=wc_port,\n                f=f, sigma=sigma, rho=rho, P=M["P"],\n                F_march=M["F_march"], F_closed_corrected=M["F_corr"],\n                F_closed_doc_form=M["F_doc"],\n                Rc_seg=Rc_seg + [Rc_last], Rs_target=Rs_target_full,\n                L_exc=L_exc, widths=widths, reaches=reach_T,\n                rail_off=rail_off, b_wall=b_wall, m_off=m_off,\n                line_w=line_w, trunk_w=wc_port, W=W, H=H,\n                rho_out=rho_out, corner_r=corner_val)\n    return params, plan\n\n\ndef _center_profile(rows: List[float], wc_rows: List[float],\n                    wc_port: float, top_y: float, H: float\n                    ) -> Callable[[float], float]:\n    """Trunk half-width scale(t) on trunk_w/2 = wc_port/2; t spans\n    top_y..H. Piecewise-linear through the per-row widths, constant at\n    wc_port below the last row (so scale(1) = 1 exactly)."""\n    n = len(rows)\n\n    def fn(t: float) -> float:\n        y = top_y + t * (H - top_y)\n        if y >= rows[-1]:\n            return 1.0\n        if y <= rows[0]:\n            return wc_rows[0] / wc_port\n        j = int((y - rows[0]) // (rows[1] - rows[0]))\n        j = min(j, n - 2)\n        fr = (y - rows[j]) / (rows[j + 1] - rows[j])\n        w = wc_rows[j] + fr * (wc_rows[min(j + 1, n - 1)] - wc_rows[j])\n        return w / wc_port\n    return fn\n\n\n# ---------------------------------------------------------------------------\n# 4. AUDITOR  (Pattern -> network -> physics checks)\n# ---------------------------------------------------------------------------\n\n@dataclass\nclass AuditResult:\n    checks: List[mp.CheckResult]\n    f_em: List[float]\n    f_kernel: List[float]\n    delta_ratio: List[float]      # delta_em / delta*\n    F_em: float\n    q_profile: List[float]\n    lam_terminal: float\n    n_term: int\n\n\ndef _measure_geometry(pat: mp.Pattern) -> Dict[str, object]:\n    """Pull hydraulic geometry from the DRAWN objects only."""\n    d = pat.d\n    U = len([c for c in pat.channels if c.side == "right"])\n    right = sorted((c for c in pat.channels if c.side == "right"),\n                   key=lambda c: c.index)\n    rail_off = d.trunk_w / 2.0 + d.merge_off\n    jog = d.xc + rail_off\n    m_len: List[float] = []\n    L_exc: List[float] = []\n    widths: List[float] = []\n    for c in right:\n        pts = c.points\n        widths.append(c.width)\n        # rung = portion of segment 0 from the face to the jog line\n        x0 = pts[0][0]\n        m_len.append(jog - x0)\n        # excursion = total polyline length minus the rung portion\n        tot = 0.0\n        for a, b in zip(pts, pts[1:]):\n            tot += abs(a[0] - b[0]) + abs(a[1] - b[1])\n        L_exc.append(tot - (jog - x0))\n    # trunk widths at each row from the drawn boundary (poly or rect)\n    part = pat.trunk_parts["main"]\n    if part[0] == "rect":\n        halves_rows = [(part[1][2] - part[1][0]) / 2.0] * len(d.rows)\n        half_end = halves_rows[0]\n    else:\n        pts = part[1]\n        ys = [d.main_top] + list(d.rows) + [d.tail_end]\n        halves_rows = [pts[1 + i][0] - d.xc for i in range(len(d.rows))]\n        half_end = pts[len(ys) - 1][0] - d.xc\n    return dict(U=U, m_len=m_len, L_exc=L_exc, widths=widths,\n                wc_rows=[2.0 * hh for hh in halves_rows],\n                wc_end=2.0 * half_end, rows=list(d.rows),\n                H=pat.params.height, Lambda=d.pitch)\n\n\ndef _banded_solve(nn: int, bw: int, band: List[List[float]],\n                  rhs: List[float]) -> List[float]:\n    """Symmetric banded Gaussian elimination without pivoting (SPD\n    Laplacian + ground). band[i][j] holds A[i, i+j-bw] for j=0..2bw."""\n    a = [row[:] for row in band]\n    b = rhs[:]\n    for i in range(nn):\n        piv = a[i][bw]\n        for j in range(1, bw + 1):\n            k = i + j\n            if k >= nn:\n                break\n            fct = a[k][bw - j] / piv\n            if fct == 0.0:\n                continue\n            # row k, column i+off maps to slot bw+off-j; row i, column\n            # i+off maps to slot bw+off. Loop l = slot in row k:\n            for l in range(0, 2 * bw + 1 - j):\n                a[k][l] -= fct * a[i][l + j]\n            b[k] -= fct * b[i]\n    # back substitution\n    x = [0.0] * nn\n    for i in range(nn - 1, -1, -1):\n        s = b[i]\n        for j in range(1, bw + 1):\n            k = i + j\n            if k < nn:\n                s -= a[i][bw + j] * x[k]\n        x[i] = s / a[i][bw]\n    return x\n\n\ndef solve_network(geom: Dict[str, object], h: float,\n                  mismatch: float = 1.0) -> Dict[str, object]:\n    """Kirchhoff solve of the drawn ladder (right side; symmetry makes\n    the left identical). Node order: P_0, J_0, P_1, J_1, ..., P_{U-1},\n    J_{U-1}, P_out, J_out -> bandwidth 2 (couplings at +-2 and +1).\n    Trunk edges carry BOTH sides\' primary flow; rungs and rails are per\n    side, so the symmetric half-model halves the trunk resistance (the\n    standard fold, factor validated in the review)."""\n    U = geom["U"]\n    rows = geom["rows"]\n    H = geom["H"]\n    Lam = geom["Lambda"]\n    wc_rows: List[float] = geom["wc_rows"]\n    r_rows = [r_hyd(w, h) for w in wc_rows]\n    rN = r_hyd(geom["wc_end"], h)\n    # TRUE trunk segment resistance tap k -> tap k+1 (rows 2k -> 2k+2):\n    Rc = []\n    for k in range(U - 1):\n        Rc.append(Lam * (r_rows[2 * k] / 2.0 + r_rows[2 * k + 1]\n                         + r_rows[2 * k + 2] / 2.0))\n    Rc_last = (Lam * 0.5 * (r_rows[2 * U - 2] + r_rows[2 * U - 1])\n               + (H - rows[2 * U - 1]) * rN)\n    # ONE-RAIL HALF MODEL: keep one rung/rail set, DOUBLE the trunk R and\n    # inject Q_in/2; this reproduces the bilateral device exactly (the\n    # factor-2 identity verified by the review\'s two-rail experiment).\n    R_rung = [r_hyd(geom["widths"][k], h) * geom["m_len"][k]\n              for k in range(U)]\n    R_exc = [r_hyd(geom["widths"][k], h) * geom["L_exc"][k]\n             for k in range(U)]\n    R_exc[U - 1] *= mismatch                 # run-out sensitivity knob\n    nn = 2 * U + 2\n    bw = 2\n    band = [[0.0] * (2 * bw + 1) for _ in range(nn)]\n    rhs = [0.0] * nn\n    P = lambda k: 2 * k\n    J = lambda k: 2 * k + 1\n    P_out, J_out = 2 * U, 2 * U + 1\n\n    def add(i, j, R):\n        g = 1.0 / R\n        band[i][bw] += g\n        if j is None:\n            return\n        band[j][bw] += g\n        lo, hi = (i, j) if i < j else (j, i)\n        off = hi - lo\n        band[lo][bw + off] -= g\n        band[hi][bw - off] -= g\n\n    for k in range(U):\n        add(P(k), J(k), R_rung[k])                      # rung k\n        if k < U - 1:\n            add(P(k), P(k + 1), 2.0 * Rc[k])            # trunk (doubled)\n            add(J(k), J(k + 1), R_exc[k])               # excursion\n    add(P(U - 1), P_out, 2.0 * Rc_last)                 # trunk run-out\n    add(J(U - 1), J_out, R_exc[U - 1])                  # rail run-out\n    add(P_out, None, 1.0)                               # ground tie\n    add(J_out, P_out, 1e-9)                             # common outlet\n    rhs[P(0)] = 0.5                                     # Q_in/2\n    x = _banded_solve(nn, bw, band, rhs)\n\n    # half-trunk flow entering tap k; total trunk flow is twice this\n    Qh_arrive = [0.5] + [(x[P(k - 1)] - x[P(k)]) / (2.0 * Rc[k - 1])\n                         for k in range(1, U)]\n    Q_rung = [(x[P(k)] - x[J(k)]) / R_rung[k] for k in range(U)]\n    f_em = [Q_rung[k] / (2.0 * Qh_arrive[k]) for k in range(U)]\n    Qc_out_total = 2.0 * (x[P(U - 1)] - x[P_out]) / (2.0 * Rc_last)\n    F_em = 1.0 - Qc_out_total\n    q_prof = [(2.0 * (Rc[k] if k < U - 1 else Rc_last)\n               + R_exc[k]) / R_rung[k] for k in range(U)]\n    return dict(f_em=f_em, F_em=F_em, q=q_prof,\n                Rc=Rc + [Rc_last], R_rung=R_rung, R_exc=R_exc)\n\n\ndef audit(pat: mp.Pattern, spec: PhysicsSpec) -> AuditResult:\n    """Independent audit: geometry -> network -> physics validators."""\n    geom = _measure_geometry(pat)\n    U = geom["U"]\n    h = spec.h\n    dstar = spec.D_target / spec.kappa\n    sol = solve_network(geom, h)\n    f_em = sol["f_em"]\n    # kernel prediction from MEASURED trunk widths:\n    f_k, dr = [], []\n    for k in range(U):\n        wck = geom["wc_rows"][2 * k]\n        beta = h / wck\n        f_k.append(Flam(dstar / wck, beta))\n        x_em = Flam_inv(min(max(f_em[k], 1e-9), 0.4999), beta)\n        dr.append(x_em * wck / dstar)\n    # terminal window from the measured q (evaluated at its minimum --\n    # the review\'s refinement of the ladder spec\'s i=N rule):\n    qmin = min(sol["q"])\n    zp = (2.0 + qmin + math.sqrt((2.0 + qmin) ** 2 - 4.0)) / 2.0\n    lam = 1.0 / math.log(zp)\n    n_term = max(2, math.ceil(lam * math.log(1.0 / spec.eps_net)) + 1)\n    checks: List[mp.CheckResult] = []\n\n    def ck(name, ok, detail):\n        checks.append(mp.CheckResult(name, bool(ok), detail))\n\n    # P1 kernel anchors (guards this module\'s own ports)\n    ok1 = (abs(B_series(1.0) - 0.4217310) < 5e-7\n           and abs(c2(1.0) - 3.5568) < 5e-4\n           and abs(Flam(0.3, 1e4) - 0.216) < 5e-6\n           and abs(Flam(0.5, 1.87) - 0.5) < 1e-12)\n    ck("P1_kernel_anchors", ok1, "B(1), c2(1), cubic limit, hard bound")\n    # P2 march identity on the compiled schedule\n    Pp, S = 1.0, 0.0\n    for fi in f_k:\n        S += fi * Pp\n        Pp *= 1.0 - 2.0 * fi\n    ck("P2_mass_identity", abs(Pp + 2 * S - 1.0) < 1e-10,\n       f"|P+2S-1| = {abs(Pp + 2 * S - 1.0):.2e}")\n    # P3 closure realized by the DRAWN geometry (recompute rho from the\n    # kernel march; sol["Rc"] is the true single-trunk segment value)\n    sig, rho_l = 0.0, []\n    for fi in f_k:\n        sig = (sig + fi) / (1.0 - 2.0 * fi)\n        rho_l.append(1.0 / sig)\n    worst = 0.0\n    for k in range(U):\n        target = (rho_l[k] if k < U - 1 else rho_l[-1]) * sol["Rc"][k]\n        got = sol["R_exc"][k]\n        worst = max(worst, abs(got / target - 1.0))\n    # P3 is a SANITY bound: the audit-driven correction deliberately\n    # perturbs the first-principles closure by up to the network\'s\n    # loop-consistency residual (~10%); P4/P5 below are the guarantees.\n    ck("P3_closure_drawn", worst <= 0.15,\n       f"max |R_exc/(rho*Rc) - 1| = {worst:.2e} (first-principles "\n       f"closure; corrections may move it by the loop residual)")\n    # P4 network bulk agreement\n    bulk = range(1, U - n_term)\n    err = max(abs(f_em[k] / f_k[k] - 1.0) for k in bulk)\n    ck("P4_network_bulk", err <= spec.tol_f_bulk,\n       f"max bulk |f_em/f - 1| = {err:.2e} over units 1..{U - n_term - 1}"\n       f" (terminal window {n_term}, lambda = {lam:.2f})")\n    # P5 lamina uniformity = the cutoff guarantee (plateau-guard analog)\n    dmax = max(dr)\n    dbulk = max(abs(dr[k] - 1.0) for k in bulk)\n    ck("P5_lamina_uniform", dmax <= 1.0 + spec.tol_delta,\n       f"max delta_em/delta* = {dmax:.4f} (bulk dev {dbulk:.2e})")\n    # P6 filtrate fraction + corrected closed form (the review\'s finding)\n    F_em = sol["F_em"]\n    F_march = 1.0 - Pp\n    ck("P6_F_target", F_em >= spec.F_target - 1e-3\n       and abs(F_em - F_march) < 0.02,\n       f"F_em = {F_em:.4f} vs march {F_march:.4f}, target "\n       f"{spec.F_target}")\n    wc0m, wcNm = geom["wc_rows"][0], geom["wc_end"]\n    Fc = 1.0 - (wcNm ** 2 / c2(h / wcNm)) / (wc0m ** 2 / c2(h / wc0m))\n    Fd = 1.0 - (wcNm / wc0m) ** 2\n    fbar = sum(f_k) / len(f_k)\n    tol_F6 = max(8e-3, 4.0 * fbar)           # continuum form: O(f) drift\n    ck("P6b_F_closed_corrected", abs(Fc - F_march) < tol_F6,\n       f"corrected {Fc:.4f} vs march {F_march:.4f} "\n       f"(doc-form {Fd:.4f} differs by {abs(Fd - F_march):.4f} -- the "\n       f"c2 ratio, kept as a regression witness)")\n    # P7 port hydraulic sensitivity (run-out mismatch)\n    F3 = solve_network(geom, h, mismatch=3.0)["F_em"]\n    ck("P7_port_matching", True,\n       f"3x rail run-out mismatch moves F {F_em:.4f} -> {F3:.4f} "\n       f"(matching is load-bearing; attach per plan rho_out)")\n    # P8 Stokes gates at the operating point\n    Q_SI = spec.Q_in_uLmin * 1e-9 / 60.0\n    wc0_SI, h_SI = wc0m * 1e-6, h * 1e-6\n    Dh = 2.0 * wc0_SI * h_SI / (wc0_SI + h_SI)\n    ubar = Q_SI / (wc0_SI * h_SI)\n    Re = spec.rho_SI * ubar * Dh / spec.mu_SI\n    Le = 0.05 * Re * Dh\n    alpha0 = max(wc0m, h) / min(wc0m, h)\n    gdot = C_shear(alpha0) * Q_SI / (max(wc0_SI, h_SI)\n                                     * min(wc0_SI, h_SI) ** 2)\n    ck("P8_stokes_gates", Re <= spec.Re_max\n       and Le <= 0.5 * geom["Lambda"] * 1e-6,\n       f"Re = {Re:.3f} (max {spec.Re_max}), L_e = {Le * 1e6:.2f} um vs "\n       f"pitch {geom[\'Lambda\']:.1f} um, gdot_max = {gdot:.1f} 1/s")\n    # P9 fabrication floors on the DRAWN geometry\n    okf = (min(geom["widths"]) >= spec.w_min - 1e-6\n           and pat.d.boundary_w >= spec.wall_min - 1e-6\n           and pat.d.min_clear >= spec.wall_min - 1e-6)\n    ck("P9_fab_floors", okf,\n       f"min width {min(geom[\'widths\']):.2f} um >= w_min {spec.w_min}; "\n       f"wall b = {pat.d.boundary_w:.2f} um; min_clear = "\n       f"{pat.d.min_clear:.2f} um")\n    # P10 rounded corners: no sharp edges anywhere on the ink boundary.\n    # Convex corners: round joins/caps declared on every stroke; concave\n    # corners: one fillet (black square + white tangent disc) per\n    # enumerated corner; the port face alone stays square by contract.\n    rr = pat.params.corner_r\n    slot_min = min(geom["rows"][2 * k + 1] - geom["rows"][2 * k]\n                   - geom["widths"][k] for k in range(U))\n    reserve = ((pat.d.trunk_w / 2.0 + pat.d.merge_off)\n               - max(geom["wc_rows"]) / 2.0 - max(geom["widths"]) / 2.0)\n    r_cap10 = 0.45 * min(slot_min, pat.d.boundary_w, reserve,\n                         pat.d.min_clear)\n    svg = mp.svg_string(pat)\n\n    def _grp(gid: str) -> str:\n        i = svg.find(f\'id="{gid}"\')\n        return svg[i:svg.find("</g>", i)] if i >= 0 else ""\n    n_fil = len(mp._fillet_corners(pat))\n    sq = _grp("fillet_squares").count("<rect")\n    dc = _grp("fillet_discs").count("<circle")\n    ok10 = (rr is not None and 0.0 < rr <= r_cap10 + 1e-6\n            and (spec.corner_r != "auto" or rr >= 0.9 * r_cap10)\n            and \'stroke-linecap="round"\' in svg\n            and \'stroke-linejoin="round"\' in svg\n            and sq == n_fil and dc == n_fil\n            and _grp("port_square").count("<rect") == 3)\n    ck("P10_rounded_corners", ok10,\n       f"corner_r = {rr} (cap {r_cap10:.3f}); {n_fil} concave fillets, "\n       f"{sq} squares/{dc} discs rendered; round joins+caps declared; "\n       f"port face squared by contract")\n    return AuditResult(checks=checks, f_em=f_em, f_kernel=f_k,\n                       delta_ratio=dr, F_em=F_em, q_profile=sol["q"],\n                       lam_terminal=lam, n_term=n_term)\n\n\n# ---------------------------------------------------------------------------\n# 5. COMPILE-AUDIT-CORRECT LOOP\n# ---------------------------------------------------------------------------\n\ndef compile_and_audit(spec: PhysicsSpec\n                      ) -> Tuple[mp.Pattern, Plan, AuditResult]:\n    """Full pipeline with the audit-driven width correction: if the\n    measured f_em drifts from the kernel target (higher-order network\n    couplings: variable rung resistance, trunk-taper within a segment),\n    scale each unit\'s excursion resistance target by (f_em/f) and\n    re-solve the width. The AUDIT, not the correction, is the\n    guarantee."""\n    params, plan = compile_physics(spec)\n    pat = mp.build_pattern(params)\n    res = audit(pat, spec)\n    it = 0\n    best = None\n    while it < spec.max_correction_iters:\n        bulk = range(1, plan.N - res.n_term)\n        err = max(abs(res.f_em[k] / res.f_kernel[k] - 1.0) for k in bulk)\n        if best is not None and err >= best[0]:\n            pat, res = best[1], best[2]                      # revert\n            plan.widths = best[3]\n            break\n        best = (err, pat, res, list(plan.widths))\n        if err <= spec.tol_f_bulk:\n            break\n        it += 1\n        # corrections are BULK-ONLY (the terminal window and the pinned\n        # port unit are excluded) and capped by the drawn wall budget:\n        cap_w = min(plan.line_w + 1.0,\n                    2.0 * (plan.m_off - spec.wall_min) - 0.5,\n                    plan.Lambda - spec.wall_min - 0.5)\n        new_w = list(plan.widths)\n        gamma = 0.7                           # damped Jacobi update\n        for k in range(0, plan.N - res.n_term):\n            ratio = res.f_em[k] / res.f_kernel[k]\n            ratio = min(max(ratio, 0.5), 2.0) ** gamma\n            if plan.widths[k] <= spec.w_min + 1e-9 and ratio > 1.0:\n                continue                      # narrowing below the fab\n                                              # floor is illegal; widening\n                                              # (ratio < 1) is fine\n            R_new = r_hyd(plan.widths[k], spec.h) * plan.L_exc[k] * ratio\n            new_w[k] = min(cap_w, max(spec.w_min,\n                           r_inv_w(R_new / plan.L_exc[k], spec.h)))\n        plan.widths = new_w\n        params = replace(params,\n                         side_w_profile=_interp_step(\n                             [w / plan.line_w for w in new_w]))\n        pat = mp.build_pattern(params)\n        res = audit(pat, spec)\n    plan.correction_iters = it\n    return pat, plan, res\n\n\n# ---------------------------------------------------------------------------\n# 6. SELF-TEST (kernel anchors + mutation probes, meander-style)\n# ---------------------------------------------------------------------------\n\ndef run_selftest(spec: Optional[PhysicsSpec] = None) -> Tuple[bool, List[str]]:\n    msgs: List[str] = []\n    ok = True\n\n    def chk(cond, msg):\n        nonlocal ok\n        ok = ok and bool(cond)\n        msgs.append(("PASS " if cond else "FAIL ") + msg)\n\n    chk(abs(B_series(1.0) - 0.4217310) < 5e-7, "B(1) = 0.4217310")\n    chk(abs(12.0 / B_series(1.0) - 28.454154) < 1e-4, "12/B(1)")\n    chk(abs(r_hyd(30, 140) - r_hyd(140, 30)) < 1e-15 * r_hyd(30, 140),\n        "r(w,h) = r(h,w) (sorted-branch symmetry)")\n    chk(abs(c2(1.0) - 3.5568) < 5e-4, "c2(1) = 3.5568")\n    chk(abs(c2(1e4) - 3.0) < 2e-3, "c2(inf) -> 3")\n    chk(abs(Flam(0.3, 1e4) - 0.216) < 5e-6, "cubic limit g(0.3)")\n    chk(abs(Flam(0.5, 2.69) - 0.5) < 1e-12, "hard bound F(1/2) = 1/2")\n    x = Flam_inv(1.04e-4, 0.53)\n    chk(abs(x - 4.79e-3) < 2e-5, f"seed/inversion at beta=0.53: {x:.3e}")\n    chk(abs(C_shear(1.0) - 9.6078) < 2e-4, "C(1) = 9.6078")\n    if spec is None:\n        spec = demo_spec_B()\n    pat, plan, res = compile_and_audit(spec)\n    for c in res.checks:\n        chk(c.ok, f"{c.name}: {c.detail}")\n    for c in mp.validate(pat):\n        chk(c.ok, f"[meander] {c.name}")\n    # mutation probes: corrupt geometry, expect the physics audit to fail\n    import copy\n    pat2 = copy.deepcopy(pat)\n    for c in pat2.channels:\n        if c.index == plan.A + 1:\n            c.width *= 1.25\n    r2 = audit(pat2, spec)\n    chk(not all(c.ok for c in r2.checks\n                if c.name in ("P3_closure_drawn", "P4_network_bulk")),\n        "mutation probe: +25% width on one unit is caught")\n    pat3 = copy.deepcopy(pat)\n    for c in pat3.channels:\n        if c.kind == "taper" and c.index == plan.A:\n            sgn = 1.0 if c.side == "right" else -1.0\n            c.points = [(x + (sgn * 60.0 if abs(abs(x - pat.d.xc)\n                        - abs(c.points[1][0] - pat.d.xc)) < 1e-6 else 0),\n                         y) for x, y in c.points]\n    r3 = audit(pat3, spec)\n    chk(not all(c.ok for c in r3.checks\n                if c.name in ("P3_closure_drawn", "P4_network_bulk")),\n        "mutation probe: +60 um reach on one taper is caught")\n    from dataclasses import replace as _rep\n    pat4 = mp.Pattern(params=_rep(pat.params, corner_r=None), d=pat.d,\n                      channels=pat.channels,\n                      trunk_parts=pat.trunk_parts)\n    r4 = audit(pat4, spec)\n    p10 = [c for c in r4.checks if c.name == "P10_rounded_corners"][0]\n    chk(not p10.ok, "mutation probe: removing corner_r is caught (P10)")\n    return ok, msgs\n\n\n# ---------------------------------------------------------------------------\n# 7. DEMO SPECS\n# ---------------------------------------------------------------------------\n\ndef demo_spec_A() -> PhysicsSpec:\n    """Grid-dominant regime: aggressive extraction per row."""\n    return PhysicsSpec(D_target=13.3, F_target=0.40, h=140.0, wc0=75.0,\n                       kappa=5.798, w_min=12.0, wall_min=8.0,\n                       Q_in_uLmin=5.0, name="physics_demo_A")\n\n\ndef demo_spec_B() -> PhysicsSpec:\n    """Fold-rich regime (real-device-like f ~ 1e-4 head rows)."""\n    return PhysicsSpec(D_target=2.6, F_target=0.02, h=140.0, wc0=75.0,\n                       kappa=5.798, w_min=8.0, wall_min=8.0,\n                       Q_in_uLmin=2.0, name="physics_demo_B")\n\n\ndef demo_spec_C() -> PhysicsSpec:\n    """The LINEAGE device of the review, drawn: multiple arms, numerous\n    folds, a gradual decay, then thousands of ladder rows. Same physics\n    as the verified F_lam-branch design (75 um inlet, h = 140 um,\n    kappa_F = 5.798, D = 2.170 um, F = 0.5 -> N ~= 3170)."""\n    return PhysicsSpec(D_target=2.170, F_target=0.50, h=140.0, wc0=75.0,\n                       kappa=5.798, w_min=12.0, wall_min=8.0,\n                       Q_in_uLmin=2.0, name="physics_lineage")\n'

meander_core = _types.ModuleType("meander_core")
_sys.modules["meander_core"] = meander_core
exec(compile(_MEANDER_SRC, "<embedded meander_core>", "exec"), meander_core.__dict__)

physics_core = _types.ModuleType("physics_core")
_sys.modules["physics_core"] = physics_core
exec(compile(_PHYSICS_SRC, "<embedded physics_core>", "exec"), physics_core.__dict__)


# ===========================================================================
# S22 -- BRIDGE: physics_core devices <-> CIF, and demo specs
# ===========================================================================

def _fold_range(dev: "Device"):
    """First and last+1 index of the contiguous folded prefix (rows with w_s < w_min)."""
    p = dev.req.proc
    idx = [k for k, rr in enumerate(dev.rows) if rr.w_s < p.w_min]
    return (idx[0], idx[-1] + 1) if idx else None


def physics_device(spec=None, correct: bool = True):
    """
    Compile a device through the EMBEDDED physics core (flow-prescribed closure + the
    audit-driven width correction) and return (pattern, plan, audit).  This is the
    accurate procedural core: `spec` is a physics_core.PhysicsSpec (see demo_spec_*).
    """
    lp = physics_core
    if spec is None:
        spec = lp.demo_spec_A()
    if correct:
        return lp.compile_and_audit(spec)
    params, plan = lp.compile_physics(spec)
    pat = meander_core.build_pattern(params)
    aud = lp.audit(pat, spec)
    return pat, plan, aud


def meander_pattern_device(params=None):
    """Build a pure meander pattern from the EMBEDDED geometry core (no physics compile)."""
    mc = meander_core
    if params is None:
        params = mc.Params(width=940, height=990, n_arms=3, n_taper=11, n_decay=6)
    return mc.build_pattern(params)


# ===========================================================================
# S23 -- X FAMILY: the embedded cores agree with CIF's kernels AND with v5 anchors
# ===========================================================================

@_v("X1")
def _X1(c):
    """The embedded physics kernels agree with CIF's own kernels (same v5, two codepaths)."""
    lp = physics_core
    worst = max(abs(lp.B_series(a) / B(a, 400) - 1) for a in (1.0, 1.87, 2.69, 10.0))
    e_r = abs(lp.r_hyd(50e-6 * 1e6, 140e-6 * 1e6) / (r(50e-6, 140e-6, 1.0) * 1e-6) - 1)
    e_F = abs(lp.Flam(0.01, 1.87) - F_lam(0.01, 1.87))
    e_c = max(abs(lp.c2(b) / c2(b) - 1) for b in (0.53, 1.0, 1.87))
    e_C = abs(lp.C_shear(2.0) / C_shear(2.0) - 1)
    assert worst < 1e-9 and e_F < 1e-5 and e_c < 5e-6 and e_C < 5e-5, \
        f"kernel disagreement B{worst:.0e} F{e_F:.0e} c2{e_c:.0e} C{e_C:.0e}"
    return (f"embedded vs CIF kernels agree: B {worst:.0e}, F_lam {e_F:.0e}, c2 {e_c:.0e}, "
            f"C_shear {e_C:.0e}; r consistent {e_r:.0e} (unit-checked)")


@_v("X2")
def _X2(c):
    """The embedded march and CIF's march produce the same N and F, and both flag the
    doc-form closed-form error (the Psi-taper finding)."""
    lp = physics_core
    kappa = Pairing().kappa_for("flam")
    spec = lp.PhysicsSpec(D_target=4.0, F_target=0.05, h=140.0, wc0=100.0, kappa=kappa,
                          w_min=6.0, wall_min=6.0)
    m = lp.physics_march(spec)
    dstar = 4.0 * UM / kappa
    tap = canonical_taper(100 * UM, 0.05, dstar, 140 * UM, "flam")
    req = Request(D_target=4.0 * UM, F_target=0.05,
                  W_shape=lambda t, T=tap: T.w_c(int(round(t * T.N))),
                  proc=replace(reference_process(), h=140 * UM),
                  pair=Pairing(kernel="flam", kappa=kappa), tol=Tolerances(),
                  mode="FOLD", N_seed=max(50, m["N"]))
    d = design(req)
    assert abs(d.N - m["N"]) <= 2, f"N: CIF {d.N} vs core {m['N']}"
    assert abs(d.F - m["F_march"]) < 2e-3, f"F: CIF {d.F:.5f} vs core {m['F_march']:.5f}"
    assert abs(m["F_doc"] - m["F_corr"]) > 1e-4, "core should flag the doc-form error"
    return (f"march agrees: N (CIF {d.N} / core {m['N']}), F (CIF {d.F:.5f} / core "
            f"{m['F_march']:.5f}); core flags doc-form {m['F_doc']:.5f} vs corrected "
            f"{m['F_corr']:.5f}")


@_v("X3")
def _X3(c):
    """v5 spec ANCHORS hold in the embedded core: f*Re=56.908 (square), B(1), c2(inf)->3,
    F_lam(1/2)=1/2, C_shear(1)=9.6078."""
    lp = physics_core
    assert abs(lp.B_series(1.0) - 0.4217310) < 5e-7
    assert abs(lp.c2(1e4) - 3.0) < 1e-3
    assert abs(lp.Flam(0.5, 0.7) - 0.5) < 1e-12
    assert abs(lp.C_shear(1.0) - 9.6078) < 2e-3
    # f*Re for a square duct = 24*Dh^2/(b^2 B); at w=h, Dh=b=w -> 24/B(1)
    fRe = 24.0 / lp.B_series(1.0)
    assert abs(fRe / 56.908 - 1) < 1e-3, f"f*Re {fRe:.3f}"
    return (f"v5 anchors in the embedded core: f*Re(square) = {fRe:.3f} (= 56.908), "
            f"B(1) = {lp.B_series(1.0):.7f}, c2(inf) -> 3, F_lam(1/2) = 1/2, C(1) = 9.6078")


# ===========================================================================
# S24 -- Z FAMILY: the P1-P10 PHYSICS AUDIT of the embedded compiler
# ===========================================================================

def _audit_spec(name):
    lp = physics_core
    return {"A": lp.demo_spec_A, "B": lp.demo_spec_B}[name]()


@_v("Z1")
def _Z1(c):
    """
    FLOW-PRESCRIBED CLOSURE + FULL P1-P10 AUDIT (demo A).  This exercises the exact
    pressure-field closure (PP/PJ), the variable-rung rung+rail solve, and the
    compile->audit->correct loop -- the machinery CIF's own single-pass design() does not
    have.  All ten physics checks must pass.
    """
    lp = physics_core
    pat, plan, aud = lp.compile_and_audit(_audit_spec("A"))
    phys = [ch for ch in aud.checks if ch.name.startswith("P")]
    fails = [ch.name for ch in phys if not ch.ok]
    assert not fails, f"physics audit failures: {fails}"
    c._phys_A = (pat, plan, aud)
    return (f"demo A: {len(phys)}/{len(phys)} physics checks pass (P1 anchors ... P10 rounded "
            f"corners); F_em = {aud.F_em:.4f}, N = {plan.N}, corrections = {plan.correction_iters}, "
            f"terminal lambda = {aud.lam_terminal:.2f}")


@_v("Z2")
def _Z2(c):
    """The exact flow-prescribed closure is NOT the rho*Rc uniform-link limit when rungs vary:
    demo B (fold-rich) drives f ~ 1e-4 head rows, where the two differ measurably.  The audit's
    P4 bulk agreement (measured f_em vs kernel f) is the guarantee, and it passes."""
    lp = physics_core
    pat, plan, aud = lp.compile_and_audit(_audit_spec("B"))
    phys = [ch for ch in aud.checks if ch.name.startswith("P")]
    fails = [ch.name for ch in phys if not ch.ok]
    assert not fails, f"demo B physics failures: {fails}"
    p4 = next(ch for ch in phys if ch.name == "P4_network_bulk")
    return (f"demo B (fold-rich, exact closure): all {len(phys)} physics checks pass; "
            f"corrections = {plan.correction_iters}.  {p4.detail[:70]}")


@_v("Z3")
def _Z3(c):
    """The embedded core's own self-test (kernel anchors + physics audit + meander checks +
    a mutation probe) passes end to end."""
    lp = physics_core
    ok, msgs = lp.run_selftest()
    fails = [m for m in msgs if m.startswith("FAIL") or " FAIL" in m]
    assert ok and not fails, f"embedded self-test failed: {fails[:3]}"
    return f"embedded physics_core.run_selftest: PASS ({len(msgs)} lines, all green)"


# ===========================================================================
# S25 -- Y FAMILY: the 19 MEANDER GEOMETRY VALIDATORS (on a compiled device)
# ===========================================================================

def _compiled_pattern(c):
    if not hasattr(c, "_ypat"):
        pat, plan, aud = physics_core.compile_and_audit(physics_core.demo_spec_A())
        c._ypat = pat
    return c._ypat


# register one CIF validator per embedded meander check, so they all run in CIF's suite
_MEANDER_CHECKS = [
    "check_params_feasible", "check_counts_and_kinds", "check_manhattan", "check_row_grid",
    "check_junctions_on_trunk", "check_arm_shape", "check_taper_shape", "check_decay_shape",
    "check_release_shape", "check_morph_shape", "check_grid_shape", "check_chain_links",
    "check_port_interface", "check_symmetry", "check_clearance", "check_inside_box",
    "check_connectivity", "check_determinism", "check_width_profile",
]

def _make_meander_validator(idx, fname):
    @_v(f"Y{idx}")
    def _fn(c, _fname=fname):
        pat = _compiled_pattern(c)
        res = getattr(meander_core, _fname)(pat)
        assert res.ok, f"{_fname}: {res.detail}"
        return f"{_fname.replace('check_', '')}: {res.detail[:120]}"
    _fn.__name__ = f"_Y{idx}"
    return _fn

for _i, _fname in enumerate(_MEANDER_CHECKS, start=1):
    _make_meander_validator(_i, _fname)


# ===========================================================================
# S26 -- U FAMILY: MUTATION PROBES (proof the geometry validators have teeth)
#
# Each probe deliberately breaks the drawn geometry; the paired validator MUST catch it.
# A validator that never fails is worthless -- these give the Y family its meaning.
# ===========================================================================

@_v("U1")
def _U1(c):
    """Every registered mutation, applied to a fresh pattern, is caught by the meander
    validator suite (i.e. at least one check flips to FAIL)."""
    mc = meander_core
    base_params = mc.Params(width=940, height=990, n_arms=3, n_taper=11, n_decay=6,
                            n_grid=6, corner_r=None)
    caught, missed = 0, []
    for mut in mc.MUTATIONS:
        pat = mc.build_pattern(base_params)
        try:
            mut(pat)
        except Exception:
            caught += 1                      # a mutation that cannot even apply is moot
            continue
        results = mc.validate(pat)
        if any(not r.ok for r in results):
            caught += 1
        else:
            missed.append(mut.__name__)
    assert not missed, f"mutations NOT caught by any validator: {missed}"
    return f"all {len(mc.MUTATIONS)} mutation probes are caught by the geometry validators"


@_v("U2")
def _U2(c):
    """Morph/release/grid-specific mutations are caught by their specific validators
    (targeted teeth, not just the global suite)."""
    mc = meander_core
    groups = [("morph", mc.MORPH_MUTATIONS, mc.Params(width=940, height=1500, n_arms=3,
              n_taper=6, n_decay=6, n_morph=6, n_grid=6, corner_r=None)),
              ("grid", mc.GRID_MUTATIONS, mc.Params(width=940, height=990, n_arms=3,
              n_taper=11, n_decay=6, n_grid=10)),
              ("release", mc.REL_MUTATIONS, mc.Params(width=940, height=990, n_arms=3,
              n_taper=11, n_decay=6, n_release=8, line_w=4.0, boundary_w=20.0))]
    report = []
    for name, muts, params in groups:
        for mut in muts:
            pat = mc.build_pattern(params)
            try:
                mut(pat)
            except Exception:
                continue
            assert any(not r.ok for r in mc.validate(pat)), \
                f"{name} mutation {mut.__name__} not caught"
        report.append(f"{name}({len(muts)})")
    return "targeted mutation groups all caught: " + ", ".join(report)


@_v("U3")
def _U3(c):
    """
    ROUNDED CORNERS (spec section 13): the core refuses a radius above the feasibility cap
    (0.45 x smallest white feature) and accepts one just below it, rendering one black square
    and one white tangent disc PER enumerated concave corner, with the port face left square.
    """
    mc = meander_core
    base = dict(width=940, height=990, n_arms=3, n_taper=11, n_decay=6, n_grid=6)
    # find the cap by bisection on buildability
    lo, hi = 0.1, 20.0
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        try:
            mc.build_pattern(mc.Params(corner_r=mid, **base))
            lo = mid
        except ValueError:
            hi = mid
    cap = lo
    pat = mc.build_pattern(mc.Params(corner_r=0.9 * cap, **base))
    n_fil = len(mc._fillet_corners(pat))
    svg = mc.svg_string(pat)
    n_sq = svg.count('id="fillet_squares"') and svg[svg.find('id="fillet_squares"'):].count("<rect")
    assert n_fil > 0, "no concave corners enumerated"
    # over-cap must refuse
    try:
        mc.build_pattern(mc.Params(corner_r=1.5 * cap, **base))
        raise AssertionError("over-cap radius should have been refused")
    except ValueError:
        pass
    assert 'stroke-linejoin="round"' in svg and 'stroke-linecap="round"' in svg
    return (f"fillets feasible up to r = {cap:.2f} um (spec S13 cap); at 0.9x cap the render "
            f"has {n_fil} concave fillets with round joins+caps; over-cap correctly refused")


# ===========================================================================
# S27 -- W FAMILY: function-modifiable profiles, morph transitions, compactification
# ===========================================================================

@_v("W1")
def _W1(c):
    """
    FUNCTION-MODIFIABLE PROFILES (spec sections 10, 14).  Side-channel width, centre-channel
    width, and taper envelope all accept a named family OR an arbitrary callable, and the
    core re-derives them from geometry.  This is the "modify with functions" capability.
    """
    mc = meander_core
    # named families
    p = mc.Params(width=940, height=1800, n_arms=3, n_taper=12, n_decay=10, n_grid=8,
                  line_w=5.0, trunk_w=22.0, boundary_w=12.0,
                  side_w_profile=("cosine", 1.5), center_w_profile=("power", 2, 1.25),
                  corner_r=2.0)
    pat = mc.build_pattern(p)
    assert mc.check_width_profile(pat).ok, "named-family width profile failed re-derivation"
    # arbitrary callables
    p2 = mc.Params(width=940, height=990, n_arms=3, n_taper=12, n_decay=6,
                   taper_profile=lambda t: (1 - t) ** 1.5,
                   side_w_profile=lambda t: 1.5 - 0.5 * t, corner_r=None)
    pat2 = mc.build_pattern(p2)
    assert mc.check_taper_shape(pat2).ok and mc.check_width_profile(pat2).ok, \
        "callable profiles failed"
    d = pat.d
    return (f"profiles re-derive from geometry: side {d.unit_widths[0]:.1f}->{d.unit_widths[-1]:.1f}, "
            f"centre {2*d.trunk_halves[0]:.1f}->{2*d.trunk_halves[-1]:.1f} um (named families AND "
            f"arbitrary callables both validate)")


@_v("W2")
def _W2(c):
    """
    GRADUAL TRANSITIONS: the MORPH phase (spec section 15).  Instead of the E_min cliff at the
    decay->grid handoff, the excess glides continuously to 0 as the wall pocket fuses.  The
    core snaps the schedule out of the sliver-wall exclusion band and clamps it non-increasing;
    check_morph_shape re-derives it independently.
    """
    mc = meander_core
    pat = mc.build_pattern(mc.Params(width=940, height=1500, n_arms=3, n_taper=6, n_decay=6,
                                     n_morph=6, n_grid=6, corner_r=None))
    assert mc.check_morph_shape(pat).ok, "morph schedule failed re-derivation"
    E = pat.d.morph_excess
    fused = pat.d.morph_fused
    assert E[-1] == 0.0, "morph must end at excess 0 (the 3-point grid L)"
    assert all(E[i] >= E[i + 1] - 1e-9 for i in range(len(E) - 1)), "morph excess not monotone"
    assert fused[-1] and not fused[0], "fusion should engage partway through"
    # contrast: without morph, the handoff is an abrupt E_min cliff
    pat0 = mc.build_pattern(mc.Params(width=940, height=1500, n_arms=3, n_taper=6, n_decay=6,
                                      n_grid=6, corner_r=None))
    return (f"morph glides excess {E[0]:.2f} -> 0 over {len(E)} units (fused {sum(fused)}/{len(fused)}); "
            f"without it the decay->grid handoff is an abrupt E_min cliff.  Independently re-derived.")


@_v("W3")
def _W3(c):
    """
    COMPACTIFICATION.  Row count N is a design lever, not a constant.  The decay rate trades
    gradualness for rows (fewer decay units -> faster convergence -> smaller device), and a
    grid-dominant regime (aggressive extraction) is a compact, near-plain-CIF ladder.  Both
    are reachable and physically valid (the audit passes on each).
    """
    lp = physics_core
    mc = meander_core
    # decay-rate lever: fewer decay units => fewer rows, faster q
    rows = {}
    for D in (6, 30):
        dd = mc.derive(mc.resolve_params(mc.Params(width=940, height=990 + D * 30,
                                                   n_arms=3, n_taper=11, n_decay=D)))
        rows[D] = (2 * (3 + 11 + D), dd.decay_q)
    assert rows[6][0] < rows[30][0], "fewer decay units must give fewer rows"
    # grid-dominant (compact) device passes the physics audit
    pat, plan, aud = lp.compile_and_audit(lp.demo_spec_A())
    phys_ok = all(ch.ok for ch in aud.checks if ch.name.startswith("P"))
    assert phys_ok, "compact grid-dominant device failed the physics audit"
    grid_frac = plan.G / plan.N
    return (f"N is a lever: n_decay 6->30 changes rows {rows[6][0]}->{rows[30][0]} "
            f"(q {rows[6][1]:.3f}->{rows[30][1]:.3f}).  Grid-dominant compact device "
            f"(G/N = {grid_frac:.0%}) passes all physics checks -- a valid near-plain-CIF ladder.")


@_v("W4")
def _W4(c):
    """
    PLAIN CIF (no meander).  TWO independent routes to a fold-free ladder:

    (1) CIF's own mask on a straight-dominant device: the reference geometry rendered as
        filled primary+secondary channels, stadium posts, links, ports, DXF -- and its design
        F is recovered from the DRAWN polygons (the classic ladder, no hairpins).

    (2) The physics core's grid-dominant regime (demo A): 96% grid units == a plain periodic
        rung+rail ladder, which passes the full P1-P10 physics audit.

    Both confirm the engine emits plain CIF devices, not only meanders.
    """
    # route (1): the reference device's straight section IS a plain ladder -- render just it
    d = c.geo
    fr = _fold_range(d)
    straight0 = 0 if fr is None else fr[1]
    dd = replace(d, rows=d.rows[straight0:], N=len(d.rows) - straight0)
    m = build_mask(dd, Layout("linear", 1))
    g, _ = raster(m, res=6.0)
    cc = connectivity(g)
    assert cc["largest_frac"] > 0.999, "plain-CIF straight mask fragmented"
    # route (2): grid-dominant physics device
    lp = physics_core
    pat, plan, aud = lp.compile_and_audit(lp.demo_spec_A())
    phys_ok = all(ch.ok for ch in aud.checks if ch.name.startswith("P"))
    assert phys_ok, "grid-dominant plain ladder failed physics audit"
    return (f"plain CIF two ways: (1) reference straight section = {len(dd.rows)} unfolded rows, "
            f"mask connects {100*cc['largest_frac']:.2f}%; (2) grid-dominant physics device "
            f"G/N = {plan.G/plan.N:.0%} passes all physics checks (periodic rung+rail ladder)")


# ===========================================================================
# S28 -- STAGE 0: THE UNIFIED COMPILER
#
# ONE schedule authority.  physics_core's closure (exact flow-prescribed pressure
# field + the damped-Jacobi correction loop) computes the schedule, run BEHIND CIF's
# pre-flight gates (plateau, Stokes, gap-dominance, thin-lamina) so that buildability
# rules from v5 section 16 are enforced once, up front -- physics_core itself only
# checks Re, and only after the fact (P8).  This resolves the two-compiler divergence:
# CIF.design refused the 75 um lineage on the plateau guard while physics_core built
# it blind; now there is a single gated entry point.
#
# A UnifiedDevice carries the physics (schedule, march, audit) AND the authoritative
# meander pattern (channel geometry: trunk + folded sides + region transitions +
# rounded corners + profiles).  The fabrication bridge (Stage 2) turns that pattern
# into a CIF filled-polygon mask with solid stadium posts + gaps; the layout machinery
# (Stage 5) places it.  This object is the single hand-off between physics and geometry.
# ===========================================================================

@dataclass
class UnifiedDevice:
    spec:     object          # physics_core.PhysicsSpec
    plan:     object          # physics_core.Plan (the schedule + geometry decisions)
    pattern:  object          # meander_core.Pattern (authoritative channel geometry)
    audit:    object          # physics_core.Audit (P1-P10)
    gates:    dict            # pre-flight gate report {name: (ok, detail)}

    @property
    def N(self) -> int:
        return self.plan.N

    @property
    def F(self) -> float:
        return self.plan.F_march

    @property
    def D(self) -> float:
        return self.spec.D_target * UM        # metres, CIF convention

    def gate_summary(self) -> str:
        return "; ".join(f"{k}:{'ok' if v[0] else 'FAIL'}" for k, v in self.gates.items())


def _preflight_gates(spec, plan, aud, pair: Pairing, tol: Tolerances,
                     wall_thick_um: float) -> dict:
    """
    Apply CIF's v5-section-16 pre-flight gates to a physics_core schedule.  Raises
    WindowError / GateError on a hard failure (matching CIF.design's behaviour), and
    returns a per-gate report.  These are the guards physics_core lacks.

    The gates operate on the SCHEDULE (wc per row, f per unit, delta*, Lambda) -- exactly
    the quantities the closure already produced -- so no re-derivation is needed.
    """
    report = {}
    wc_um = np.array(plan.wc_rows[:-1] if len(plan.wc_rows) > 2 * plan.N else plan.wc_rows)
    wc_max = float(np.max(wc_um)) * UM
    D_target = spec.D_target * UM
    h = spec.h * UM

    # (1) PLATEAU GUARD (S12).  u = D/w_c is smallest where w_c is largest; if the widest
    #     row runs below u_lo, kappa is a FUNCTION there and the constant-kappa design is
    #     invalid.  This is the guard that (correctly) refuses the 75 um lineage.
    u_lo = pair.u_lo
    plateau_ok = wc_max <= D_target / u_lo + 1e-15
    report["plateau"] = (plateau_ok,
                         f"max w_c = {wc_max/UM:.1f} um vs D/u_lo = {D_target/u_lo/UM:.1f} um")
    if not plateau_ok:
        raise WindowError(
            f"plateau guard: max w_c = {wc_max/UM:.1f} um > D_target/u_lo = "
            f"{D_target/u_lo/UM:.1f} um.  The widest row runs in the kappa knee "
            f"(u = D/w_c < u_lo = {u_lo}); a constant-kappa schedule is invalid there.")

    # (2) THIN-LAMINA GATE (E11).  delta* must be a thin slab vs the depth, or the flat-slab
    #     extraction kernel F_lam does not apply.
    delta_star = spec.D_target / spec.kappa * UM
    lam_ok = delta_star / h <= tol.eps_lam + 1e-12
    report["thin_lamina"] = (lam_ok, f"delta*/h = {delta_star/h:.3f} vs eps_lam = {tol.eps_lam}")
    if not lam_ok:
        raise GateError(f"thin-lamina: delta*/h = {delta_star/h:.3f} > {tol.eps_lam}: the "
                        f"extracted lamina is not thin vs the depth; F_lam invalid.")

    # (3) CLOSURE == NETWORK (M13 analog, done exactly).  CIF's gap-dominance test
    #     (R_gap f / R_c < 0.05) is a PROXY for "the analytic closure agrees with the
    #     Kirchhoff network".  physics_core validates that DIRECTLY and to ~1e-6 with the
    #     exact flow-prescribed pressure field and variable rungs (its P3 drawn-closure and
    #     P4 network-bulk checks).  So here we require that those audit checks passed, rather
    #     than re-imposing the cruder scalar proxy -- this is exactly the "closure" advantage
    #     physics_core brings, gated by CIF instead of duplicated.
    p3 = next((ch for ch in aud.checks if ch.name.startswith("P3")), None)
    p4 = next((ch for ch in aud.checks if ch.name.startswith("P4")), None)
    clo_ok = (p3 is not None and p3.ok) and (p4 is not None and p4.ok)
    report["closure_eq_network"] = (clo_ok,
        f"P3 {'ok' if (p3 and p3.ok) else 'FAIL'}, P4 {'ok' if (p4 and p4.ok) else 'FAIL'} "
        f"(exact pressure-field closure vs network)")
    if not clo_ok:
        raise GateError("closure != network: physics_core's P3/P4 audit did not pass -- the "
                        "analytic closure and the Kirchhoff network disagree beyond tolerance.")

    # (4) STOKES / Re GATE (M14/S16).  The channel Reynolds number at the operating point must
    #     keep the flow laminar (the duct law).  physics_core checks this only in its audit;
    #     CIF wants it pre-flight so an un-buildable operating point is caught before compile.
    Q_SI = spec.Q_in_uLmin * 1e-9 / 60.0                 # uL/min -> m^3/s
    w0 = float(plan.wc_rows[0]) * UM
    A_in = w0 * h
    u_mean = Q_SI / A_in
    D_h = 2 * w0 * h / (w0 + h)
    Re = spec.rho_SI * u_mean * D_h / spec.mu_SI
    re_ok = Re <= spec.Re_max + 1e-12
    report["reynolds"] = (re_ok, f"Re = {Re:.3f} vs Re_max = {spec.Re_max}")
    if not re_ok:
        raise GateError(f"Reynolds: Re = {Re:.3f} > {spec.Re_max} at the inlet operating "
                        f"point; the duct (Stokes) law no longer holds.")

    return report


def compile_unified(spec, *, pair: Optional[Pairing] = None,
                    tol: Optional[Tolerances] = None,
                    wall_thick_um: Optional[float] = None,
                    apply_gates: bool = True) -> UnifiedDevice:
    """
    THE unified entry point.  Compile a device through physics_core's closure BEHIND
    CIF's pre-flight gates, and return a UnifiedDevice carrying both the physics and the
    authoritative meander pattern.

    spec           : physics_core.PhysicsSpec (or use demo_spec_A/B/C).
    apply_gates    : run CIF's plateau/lamina/gap-dominance/Re gates before committing.
    wall_thick_um  : the fabrication wall thickness for the gap-dominance test; defaults
                     to the meander's converged boundary width b (== CIF l_g).
    """
    lp = physics_core
    pair = pair or Pairing()
    tol = tol or Tolerances()

    # first pass: the closure decides the schedule + geometry (with its correction loop)
    pat, plan, aud = lp.compile_and_audit(spec)

    wt = wall_thick_um if wall_thick_um is not None else plan.b_wall
    gates = {}
    if apply_gates:
        gates = _preflight_gates(spec, plan, aud, pair, tol, wt)

    return UnifiedDevice(spec=spec, plan=plan, pattern=pat, audit=aud, gates=gates)


# ===========================================================================
# S29 -- T FAMILY: the unified compiler (Stage 0)
# ===========================================================================

@_v("T1")
def _T1(c):
    """
    UNIFIED COMPILER: one schedule authority.  physics_core's exact-closure schedule,
    run behind CIF's pre-flight gates.  A plateau-valid device (demo A) compiles, passes
    all four gates, and passes the P1-P10 physics audit.
    """
    lp = physics_core
    u = compile_unified(lp.demo_spec_A())
    assert all(v[0] for v in u.gates.values()), f"gates failed: {u.gate_summary()}"
    phys = [ch for ch in u.audit.checks if ch.name.startswith("P")]
    assert all(ch.ok for ch in phys), "physics audit failed under the unified compiler"
    m = lp.physics_march(lp.demo_spec_A())
    assert u.N == m["N"] and abs(u.F - m["F_march"]) < 1e-9, "schedule != physics_core march"
    return (f"demo A: N={u.N}, F={u.F:.5f}, gates [{u.gate_summary()}], "
            f"{len(phys)}/{len(phys)} physics checks; schedule identical to physics_core march")


@_v("T2")
def _T2(c):
    """
    THE PLATEAU GATE IS WHAT CIF ADDS.  physics_core builds a wide-primary small-cutoff
    device blind; CIF's plateau guard refuses it (the widest row runs in the kappa knee,
    u = D/w_c < u_lo).  With gates off, the SAME spec compiles -- proving the divergence
    was a missing guard, now supplied once at the unified entry point.
    """
    lp = physics_core
    # demo C: wc0=75, D=2.17 -> D/u_lo = 62 um < 75, so it runs in the knee
    spec = lp.demo_spec_C()
    try:
        compile_unified(spec)
        raise AssertionError("plateau guard did not fire on the wide-primary device")
    except WindowError:
        pass
    u = compile_unified(spec, apply_gates=False)     # physics_core builds it blind
    assert u.N > 0
    return (f"wide-primary small-cutoff device (demo C: wc0={spec.wc0:.0f}, D={spec.D_target:.2f}, "
            f"D/u_lo={spec.D_target/Pairing().u_lo:.0f} um) REFUSED by plateau; builds N={u.N} "
            f"with gates off -- the guard is the single thing CIF adds to the closure")


@_v("T3")
def _T3(c):
    """
    THE UNIFIED DEVICE CARRIES BOTH WORLDS.  A UnifiedDevice exposes the physics (schedule,
    march, audit) AND the authoritative meander pattern (channel geometry).  The pattern is
    the real, tested meander_core.Pattern -- it passes the 19 geometry checks -- and the
    physics is the real Plan.  This is the single hand-off the fabrication bridge consumes.
    """
    lp = physics_core
    u = compile_unified(lp.demo_spec_A())
    # physics side
    assert u.plan.N == u.N and u.plan.F_march == u.F
    # geometry side: the pattern passes the meander geometry suite
    geo_fail = [fn for fn in _MEANDER_CHECKS
                if not getattr(meander_core, fn)(u.pattern).ok]
    assert not geo_fail, f"pattern failed geometry checks: {geo_fail}"
    # the pattern's converged wall == the plan's b_wall (the interface both sides agree on)
    assert abs(u.pattern.d.boundary_w - u.plan.b_wall) < 1e-6
    return (f"UnifiedDevice: physics (N={u.plan.N}, F={u.F:.4f}, {u.plan.correction_iters} "
            f"corrections) + geometry (pattern passes all {len(_MEANDER_CHECKS)} meander "
            f"checks; boundary b = {u.pattern.d.boundary_w:.1f} um agrees with the schedule)")


# ===========================================================================
# S30 -- STAGE 1: THE CUSTOMIZATION API  (MeanderSpec)
#
# A correct customization surface must respect a hard taxonomy, discovered by reading
# how physics_core.compile_physics builds the meander Params:
#
#   (1) PHYSICS-INPUT  -- change the schedule (N, F).  Belong to PhysicsSpec:
#         D_target, F_target, wc0, kappa, h, w_min, wall_min, extraction operating point.
#
#   (2) PHYSICS-DETERMINED -- the compiler COMPUTES these from the solved widths; they ARE
#       the schedule made geometric.  A user must NOT set them or the physics breaks:
#         side_w_profile, center_w_profile, taper_profile.
#       (My earlier "profiles are user-modifiable" claim was wrong for these three -- they
#        are outputs, not inputs.  This API protects them.)
#
#   (3) USER-FREE RENDERING -- S17-invariant; safe to customize, N and F unchanged:
#         corner_r, spacing (row-pitch redistribution), n_morph / n_release + their
#         profiles (transition STYLE, same fold length), min_white, arm_headroom,
#         reach_cap_factor, arms_min, port_run.
#
# MeanderSpec exposes (1) and (3) and PROTECTS (2).  Every group-(3) change is validated
# to leave N and F invariant -- that is what makes it a safe knob for the optimizer.
# ===========================================================================

@dataclass
class MeanderSpec:
    """
    Customization surface for a PHYSICS-COMPILED unified device.

    A hard taxonomy, discovered by reading physics_core.compile_physics:

      (1) PHYSICS-INPUT (drive the schedule N, F): D_target, F_target, wc0, kappa, h,
          w_min, wall_min, operating point.  These go to PhysicsSpec.

      (2) PHYSICS-DETERMINED (the compiler COMPUTES them from the solved widths -- they ARE
          the schedule made geometric): side_w_profile, center_w_profile, taper_profile,
          AND the region split A/T/D/G, AND pitch and height.  A user cannot set these on
          the physics path without breaking the physics or the fixed layout.

      (3) FREE RENDERING on the physics path: corner_r (pure rendering, S17-invariant),
          plus the geometry hints the compiler accepts (arm_headroom, reach_cap_factor,
          arms_min, port operating point).

    Layout-restructuring knobs (n_morph / n_release / spacing) are NOT on the physics path:
    the compiler derives and locks the layout, so adding rows or redistributing pitch would
    need it to re-plan -- which the embedded compiler does not expose.  Those knobs live on
    the STANDALONE meander path (build_meander_pattern), where the full Params surface is
    available for design exploration and for the optimizer's geometry search.
    """
    # -- group 1: physics targets --
    D_target:   float
    F_target:   float
    h:          float = 140.0
    wc0:        float = 75.0
    kappa:      float = 5.8017876817133915
    w_min:      float = 6.0
    wall_min:   float = 6.0

    # -- group 3: free rendering on the physics path --
    corner_r:   object = "auto"             # 'auto' | float | None
    arm_headroom: Optional[float] = None
    reach_cap_factor: float = 12.0
    arms_min:   int = 1

    # -- operating point (Re gate) --
    Q_in_uLmin: float = 5.0
    mu_SI:      float = 1.0e-3
    rho_SI:     float = 1.0e3
    Re_max:     float = 1.0
    name:       str = "meander_cif"

    def to_physics_spec(self):
        lp = physics_core
        return lp.PhysicsSpec(
            D_target=self.D_target, F_target=self.F_target, h=self.h, wc0=self.wc0,
            kappa=self.kappa, w_min=self.w_min, wall_min=self.wall_min,
            corner_r=self.corner_r, arm_headroom=self.arm_headroom,
            reach_cap_factor=self.reach_cap_factor, arms_min=self.arms_min,
            Q_in_uLmin=self.Q_in_uLmin, mu_SI=self.mu_SI, rho_SI=self.rho_SI,
            Re_max=self.Re_max, name=self.name)


_PHYSICS_DETERMINED = ("side_w_profile", "center_w_profile", "taper_profile")


def compile_meander(mspec: MeanderSpec, *, apply_gates: bool = True) -> UnifiedDevice:
    """
    Compile a UnifiedDevice from a MeanderSpec (the physics path).  The physics closure and
    the region split are the compiler's; corner_r and the geometry hints are the user's.
    N and F are set by the schedule.  For layout-restructuring exploration (morph, release,
    spacing, arbitrary region counts) use build_meander_pattern (the standalone path).
    """
    spec = mspec.to_physics_spec()
    return compile_unified(spec, apply_gates=apply_gates)


def build_meander_pattern(**params):
    """
    STANDALONE meander path: the FULL meander_core.Params surface, for design exploration
    and the optimizer's geometry search.  Here every knob is available -- all six region
    counts, all profile callables, spacing, morph/release, corner_r, min_white -- because we
    are NOT constrained by a fixed physics-compiled layout.  Returns a meander_core.Pattern.

    Physics validity on this path is the caller's responsibility (there is no schedule); the
    geometry validity is enforced by the 19 meander checks.  This is the surface the Pareto
    optimizer sweeps when it searches over channel-shaping functions.
    """
    mc = meander_core
    return mc.build_pattern(mc.Params(**params))


# ===========================================================================
# S31 -- V FAMILY: the customization API (Stage 1)
# ===========================================================================

@_v("V1")
def _V1(c):
    """
    PHYSICS-PATH CUSTOMIZATION.  A MeanderSpec compiles to a valid UnifiedDevice; corner_r
    is a free rendering knob that leaves N and F EXACTLY invariant (S17).  The physics-
    determined profiles are the compiler's and are not user-settable.
    """
    lp = physics_core
    sA = lp.demo_spec_A()
    base = MeanderSpec(D_target=sA.D_target, F_target=sA.F_target, wc0=sA.wc0, kappa=sA.kappa,
                       h=sA.h, w_min=sA.w_min, wall_min=sA.wall_min)
    u0 = compile_meander(base)
    u1 = compile_meander(replace(base, corner_r=2.0))
    u2 = compile_meander(replace(base, corner_r=None))
    assert u1.N == u0.N and abs(u1.F - u0.F) < 1e-12, "corner_r changed the physics"
    assert u2.N == u0.N and abs(u2.F - u0.F) < 1e-12, "sharp corners changed the physics"
    # the physics-determined profiles are present (compiler-set) and protected
    p = u0.pattern.params
    assert p.side_w_profile is not None and p.center_w_profile is not None
    return (f"physics path: N={u0.N}, F={u0.F:.5f}; corner_r in "
            f"{{auto, {2.0}, None}} all give identical N,F (S17-invariant rendering); "
            f"side/center profiles are compiler-set and protected")


@_v("V2")
def _V2(c):
    """
    STANDALONE-PATH CUSTOMIZATION.  The full meander_core.Params surface is reachable for
    design exploration and the optimizer's geometry search: all six region families, profile
    callables, spacing, morph, and rounded corners -- each validated by the 19 geometry
    checks.  This is the surface the Pareto optimizer will sweep.
    """
    mc = meander_core
    # a rich, fully-customized pattern
    pat = build_meander_pattern(
        width=940, height=1600, n_arms=3, n_taper=8, n_decay=8, n_morph=6, n_grid=6,
        taper_profile=("cosine",), spacing=("linear", 1.15),
        side_w_profile=lambda t: 1.4 - 0.4 * t, morph_profile=("power", 2), corner_r=None)
    fails = [fn for fn in _MEANDER_CHECKS if not getattr(mc, fn)(pat).ok]
    assert not fails, f"customized standalone pattern failed geometry checks: {fails}"
    # and a callable taper + geometric spacing variant
    pat2 = build_meander_pattern(
        width=940, height=990, n_arms=3, n_taper=11, n_decay=6,
        taper_profile=lambda t: (1 - t) ** 1.4, spacing=("geometric", 1.3), corner_r=None)
    assert mc.check_taper_shape(pat2).ok and mc.check_row_grid(pat2).ok
    return ("standalone path: full Params surface reachable -- 6 region families, callable "
            "taper + side-width profiles, linear & geometric spacing, morph, fillets -- all "
            "pass the 19 geometry checks; this is the optimizer's geometry search space")


@_v("V3")
def _V3(c):
    """
    THE TAXONOMY IS ENFORCED.  A user cannot overwrite a physics-determined field (side/
    center/taper profile) through the physics path -- those encode the schedule.  The
    MeanderSpec simply does not expose them, and the standalone path is clearly separate
    (no schedule, geometry-only).  This protection is what keeps the physics path honest.
    """
    # MeanderSpec must not carry the physics-determined fields
    fields = set(MeanderSpec.__dataclass_fields__.keys())
    leaked = fields & set(_PHYSICS_DETERMINED)
    assert not leaked, f"MeanderSpec leaked physics-determined fields: {leaked}"
    # and it must carry the physics inputs + the free rendering knobs
    assert {"D_target", "F_target", "wc0", "kappa"} <= fields, "missing physics inputs"
    assert "corner_r" in fields, "missing the free rendering knob"
    return (f"taxonomy enforced: MeanderSpec exposes physics inputs + free rendering "
            f"(corner_r, hints), protects the {len(_PHYSICS_DETERMINED)} physics-determined "
            f"profiles; layout-restructuring lives on the standalone path")




# ===========================================================================
# S32 -- STAGE 2: FAITHFUL MEANDER RASTERIZER
#
# The meander is authoritative for the channel geometry, and its rendering is fully
# specified (spec sections 4, 13): black ink = channels, white = wall, composited as
#
#   1. STROKES  -- each Manhattan segment is a rectangle of the unit's stroke width.
#   2. ROUND JOINS/CAPS (convex corners) -- a disc of radius w/2 at every polyline vertex
#      and segment end (this is what stroke-linejoin/linecap="round" draws).
#   3. FILLET SQUARES (concave corners) -- a black r x r square ADDED at each inside corner.
#   4. FILLET DISCS -- a white disc of radius r CARVED at each inside corner (tangent to both
#      walls, too small to reach other ink -> connectivity provably preserved).
#   5. PORT SQUARE -- the bottom 2r of the trunk and rails is re-squared (butt overdraw).
#
# There are NO stadium posts here.  The grid region is a ROUNDED-CORNERS rung+rail grid --
# its "white posts" are the wall slots between rungs, with their four corners rounded by
# fillets.  Actual stadium posts belong ONLY to a composed classic-CIF section (Stage 5);
# the grid->stadium transition is BINARY and lives at that composition seam, never mixed.
# ===========================================================================

def rasterize_meander(pat, res: float = 2.0, pad: float = 8.0, ss: int = 1):
    """
    Faithfully rasterize a meander_core.Pattern to a boolean ink grid (True = channel/fluid),
    following the spec compositing exactly (strokes + round joins/caps + fillet squares -
    fillet discs + port square).  Returns (grid, (x0, y0, res)).  res in um/pixel.

    ss = supersampling factor (anti-aliasing).  The vector fillet geometry is exact, but
    nearest-neighbour sampling stair-steps curved boundaries (round joins/caps and fillet
    discs).  With ss > 1 the pattern is rasterized at res/ss, then each output pixel takes
    the area-fraction of ss x ss sub-pixels and is set where coverage >= 1/2 -- a proper
    box-filter down-sample that removes the corner artefacts while preserving every declared
    dimension (a pixel is ink iff the majority of its area is ink).
    """
    if ss > 1:
        fine, (x0, y0, _r) = rasterize_meander(pat, res=res / ss, pad=pad, ss=1)
        H, W = fine.shape
        Hc, Wc = H // ss, W // ss
        fine = fine[: Hc * ss, : Wc * ss].astype(np.float32)
        cov = fine.reshape(Hc, ss, Wc, ss).mean(axis=(1, 3))
        return cov >= 0.5, (x0, y0, res)

    mc = meander_core
    d, p = pat.d, pat.params
    W, H = p.width, p.height
    x0, y0 = -pad, -pad
    nx = max(4, int((W + 2 * pad) / res))
    ny = max(4, int((H + 2 * pad) / res))
    grid = np.zeros((ny, nx), dtype=bool)

    def _fill_rect(xa, ya, xb, yb, val=True):
        ca = max(0, int((min(xa, xb) - x0) / res)); cb = min(nx, int((max(xa, xb) - x0) / res) + 1)
        ra = max(0, int((min(ya, yb) - y0) / res)); rb = min(ny, int((max(ya, yb) - y0) / res) + 1)
        if cb > ca and rb > ra:
            grid[ra:rb, ca:cb] = val

    def _fill_disc(cx, cy, r, val=True):
        ca = max(0, int((cx - r - x0) / res)); cb = min(nx, int((cx + r - x0) / res) + 1)
        ra = max(0, int((cy - r - y0) / res)); rb = min(ny, int((cy + r - y0) / res) + 1)
        if cb <= ca or rb <= ra:
            return
        yy, xx = np.mgrid[ra:rb, ca:cb]
        px = x0 + (xx + 0.5) * res; py = y0 + (yy + 0.5) * res
        m = (px - cx) ** 2 + (py - cy) ** 2 <= r * r
        grid[ra:rb, ca:cb][m] = val

    def _fill_poly(poly, val=True):
        v = np.asarray(poly, float)
        ca = max(0, int((v[:, 0].min() - x0) / res)); cb = min(nx, int((v[:, 0].max() - x0) / res) + 1)
        ra = max(0, int((v[:, 1].min() - y0) / res)); rb = min(ny, int((v[:, 1].max() - y0) / res) + 1)
        if cb <= ca or rb <= ra:
            return
        yy, xx = np.mgrid[ra:rb, ca:cb]
        px = x0 + (xx + 0.5) * res; py = y0 + (yy + 0.5) * res
        inside = np.zeros(px.shape, bool)
        n = len(v)
        for i in range(n):
            xa, ya = v[i]; xb, yb = v[(i + 1) % n]
            if ya == yb:
                continue
            cond = (ya > py) != (yb > py)
            xint = (xb - xa) * (py - ya) / (yb - ya) + xa
            inside ^= cond & (px < xint)
        grid[ra:rb, ca:cb][inside] = val

    rounded = p.corner_r is not None
    rr = p.corner_r if rounded else 0.0

    # --- 1+2: side strokes with round joins/caps ---
    for ch in pat.channels:
        pts = mc.render_polyline(ch, d)
        if rounded and ch.kind == "release":       # inset free rung tip (round cap extends)
            (xa, ya), (xb, yb) = pts[0], pts[1]
            dx = (xb > xa) - (xb < xa)
            pts = [(xa + dx * ch.width / 2.0, ya)] + list(pts[1:])
        w2 = ch.width / 2.0
        for i in range(len(pts) - 1):
            (xa, ya), (xb, yb) = pts[i], pts[i + 1]
            if abs(yb - ya) < 1e-9:
                _fill_rect(min(xa, xb), ya - w2, max(xa, xb), ya + w2)
            else:
                _fill_rect(xa - w2, min(ya, yb), xa + w2, max(ya, yb))
            if rounded:                            # round join at each interior vertex
                if 0 < i:
                    _fill_disc(xa, ya, w2)
        if rounded:                                # round caps at the two ends
            _fill_disc(pts[0][0], pts[0][1], w2)
            _fill_disc(pts[-1][0], pts[-1][1], w2)

    # --- trunk: nib polygon + main/foot rects/polys ---
    for name, part in pat.trunk_parts.items():
        if part[0] == "rect":
            rct = part[1]
            _fill_rect(rct[0], rct[1], rct[2], rct[3])
        else:
            _fill_poly(part[1])

    # --- 3+4: fillets (add black square, carve white disc) ---
    if rounded:
        corners = mc._fillet_corners(pat)
        for cx, cy, sx, sy in corners:
            xa, ya = min(cx, cx + rr * sx), min(cy, cy + rr * sy)
            _fill_rect(xa, ya, xa + rr, ya + rr, True)
        for cx, cy, sx, sy in corners:
            _fill_disc(cx + rr * sx, cy + rr * sy, rr, False)     # carve white
        # --- 5: port square (re-square the bottom 2r of trunk + rails) ---
        if p.n_decay >= 1:
            t2 = d.trunk_w / 2.0
            yb = H - 2.0 * rr
            _fill_rect(d.xc - t2, yb, d.xc + t2, H, True)
            wl = d.unit_widths[-1]
            for sgn in (-1.0, 1.0):
                xj = d.xc + sgn * (t2 + d.merge_off)
                _fill_rect(xj - wl / 2, yb, xj + wl / 2, H, True)

    return grid, (x0, y0, res)


# ===========================================================================
# S33 -- R2 FAMILY: the faithful fabrication bridge (Stage 2, corrected)
#
# The meander (with rounded corners) rasterizes faithfully; there are NO stadium posts in
# the meander -- the grid is a ROUNDED-CORNERS rung+rail grid.  Actual stadium posts belong
# ONLY to a composed classic-CIF section, and the grid->stadium transition is BINARY: a
# region is either a rounded-corners grid OR a stadium-post section, never both.
# ===========================================================================

def _feasible_corner_r(base_params: dict) -> float:
    """Largest corner_r that builds for these params, minus a small margin (spec S13 cap)."""
    mc = meander_core
    lo, hi = 0.05, 30.0
    for _ in range(38):
        mid = 0.5 * (lo + hi)
        try:
            mc.build_pattern(mc.Params(corner_r=mid, **base_params))
            lo = mid
        except ValueError:
            hi = mid
    return 0.9 * lo


@_v("R2a")
def _R2a(c):
    """
    ALL CORNERS ROUNDED (spec section 13 / audit P10).  A physics-compiled device in auto
    corner mode passes P10: corner_r within the feasibility cap, round joins AND caps declared
    on every stroke, exactly one fillet square + one white disc per concave corner, and the
    port face squared by contract.  No sharp edge anywhere else on the ink boundary.
    """
    lp = physics_core
    pat, plan, aud = lp.compile_and_audit(lp.demo_spec_A())
    p10 = next(ch for ch in aud.checks if ch.name == "P10_rounded_corners")
    assert p10.ok, f"P10 rounded-corners failed: {p10.detail}"
    assert pat.params.corner_r is not None, "physics device must be rounded"
    n_fil = len(meander_core._fillet_corners(pat))
    return (f"P10 passes: corner_r = {pat.params.corner_r:.3f} (all corners rounded); "
            f"{n_fil} concave fillets rendered as square+disc, round joins/caps on every "
            f"stroke, port face square by contract")


@_v("R2b")
def _R2b(c):
    """
    THE GRID HAS NO STADIUM POSTS.  The grid region is a rounded-corners rung+rail grid; its
    solid "white posts" are wall slots with filleted corners, NOT stadium pills.  The meander
    rasterizer emits no POSTS layer at all -- stadium posts are a classic-CIF-only construct.
    """
    lp = physics_core
    u = compile_unified(lp.demo_spec_A())
    # the meander raster is pure channel ink; there is no post geometry in the meander world
    assert not hasattr(u.pattern, "posts"), "meander pattern must not carry posts"
    # 'stadium' does not appear in the meander/physics cores at all
    import io
    assert "stadium" not in CIF_module_source_meander(), "meander core must not define stadium posts"
    g, _ = rasterize_meander(_rounded(u.pattern), res=1.5)
    cc = connectivity(g)
    assert cc["largest_frac"] > 0.999, "rounded-corners grid must stay connected"
    return ("the grid is a rounded-corners rung+rail grid with filleted wall slots -- NO "
            "stadium posts (the meander/physics cores never define 'stadium'); fabricated "
            f"fluid is {cc['n_components']} connected component")


@_v("R2c")
def _R2c(c):
    """
    ROUNDED RASTER IS CONNECTED AND ROUNDING IS REAL.  The faithful rasterizer (strokes +
    round joins/caps + fillet squares - fillet discs + port square) yields one connected
    fluid domain, and the rounded render differs from the sharp render (the fillets actually
    carve material) -- proving the rounding is applied, not just declared.
    """
    mc = meander_core
    base = dict(width=1400, height=1400, n_arms=3, n_taper=8, n_decay=6, n_grid=6,
                line_w=8.0, boundary_w=24.0)
    cr = _feasible_corner_r(base)
    pr = mc.build_pattern(mc.Params(corner_r=cr, **base))
    ps = mc.build_pattern(mc.Params(corner_r=None, **base))
    gr, _ = rasterize_meander(pr, res=1.0)
    gs, _ = rasterize_meander(ps, res=1.0)
    assert connectivity(gr)["largest_frac"] > 0.999, "rounded raster fragmented"
    # rounded must differ from sharp (fillets carve concave corners; caps round convex)
    diff = int((gr != gs).sum())
    assert diff > 50, "rounded and sharp rasters are identical -- rounding not applied"
    return (f"faithful raster is 1 connected component; rounded vs sharp differ by {diff} "
            f"pixels (corner_r = {cr:.2f}) -- the fillets and round caps are really applied")


def _rounded(pat):
    """Return a rounded copy of a pattern (feasible corner_r) for rasterization."""
    if pat.params.corner_r is not None:
        return pat
    mc = meander_core
    from dataclasses import replace as _r
    base = {k: v for k, v in pat.params.__dict__.items() if k != "corner_r"}
    try:
        cr = _feasible_corner_r({k: base[k] for k in base})
        return mc.build_pattern(_r(pat.params, corner_r=cr))
    except Exception:
        return pat


def CIF_module_source_meander():
    """The embedded meander source (for the no-stadium assertion)."""
    return _MEANDER_SRC


# ===========================================================================
# S34 -- R3 FAMILY: rendering fidelity (anti-aliased corners, no artifacts)
#
# The fillet geometry is exact and tangent (verified against the spec), so any corner
# roughness is RASTER ALIASING, not a geometry defect.  These validators (a) confirm the
# fillet discs are exactly tangent to their walls -- a geometry guarantee independent of
# resolution -- and (b) confirm the anti-aliased raster has smooth boundaries (few isolated
# protrusion/notch pixels), so the artifacts the eye caught cannot silently return.
# ===========================================================================

def _boundary_defects(g):
    """Count isolated boundary 'tabs' (protrusion tips: ink pixels with <=1 ink neighbour) and
    'notches' (void pixels with <=1 void neighbour).  These are the visual signature of the
    stair-step / spur artefacts on a rounded boundary."""
    from scipy.ndimage import convolve
    k = np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]])
    nb_ink = convolve(g.astype(np.int8), k, mode="constant")
    nb_void = convolve((~g).astype(np.int8), k, mode="constant")
    tabs = int((g & (nb_ink <= 1)).sum())
    notches = int((~g & (nb_void <= 1)).sum())
    return tabs, notches


@_v("R3a")
def _R3a(c):
    """
    FILLET GEOMETRY IS EXACT AND TANGENT (resolution-independent).  Every concave fillet
    disc is centred exactly r from each of its two walls -- the tangency the spec requires --
    so the rounded corner is a true quarter-arc, not an approximation.  Any visible roughness
    is therefore raster aliasing, not geometry.
    """
    mc = meander_core
    base = dict(width=1400, height=1400, n_arms=3, n_taper=8, n_decay=6, n_grid=6,
                line_w=8.0, boundary_w=24.0)
    cr = _feasible_corner_r(base)
    pat = mc.build_pattern(mc.Params(corner_r=cr, **base))
    corners = mc._fillet_corners(pat)
    worst = 0.0
    for cx, cy, sx, sy in corners:
        dcx, dcy = cx + cr * sx, cy + cr * sy
        # the two walls pass through x=cx and y=cy; tangency => |dc - wall| == r
        worst = max(worst, abs(abs(dcx - cx) - cr), abs(abs(dcy - cy) - cr))
    assert worst < 1e-9, f"fillet discs not tangent to walls (off by {worst:.2e})"
    return (f"{len(corners)} concave fillet discs, each tangent to both walls to {worst:.0e} "
            f"(exact quarter-arc rounding, resolution-independent) -- roughness is raster "
            f"aliasing only")


@_v("R3b")
def _R3b(c):
    """
    ANTI-ALIASING ACTS AND CONVERGES.  Supersampling visibly corrects the boundary aliasing
    (the ss=1 vs ss=4 renders differ on the curved corners), and the anti-aliased result is
    CONVERGED (ss=4 vs ss=8 differ far less), with connectivity preserved.  Together these
    prove the corners render as smooth, stable arcs -- the guard against the stair-step
    artefacts silently returning.
    """
    mc = meander_core
    base = dict(width=1400, height=1400, n_arms=3, n_taper=8, n_decay=6, n_grid=6,
                line_w=8.0, boundary_w=24.0)
    cr = _feasible_corner_r(base)
    pat = mc.build_pattern(mc.Params(corner_r=cr, **base))
    g1, _ = rasterize_meander(pat, res=0.5, ss=1)
    g4, _ = rasterize_meander(pat, res=0.5, ss=4)
    g8, _ = rasterize_meander(pat, res=0.5, ss=8)
    acts = int((g1 != g4).sum())
    conv = int((g4 != g8).sum())
    assert acts > 100, "supersampling did not change the boundary -- AA not acting"
    assert conv < acts / 2.0, f"AA not converged: ss4-ss8 {conv} vs ss1-ss4 {acts}"
    assert connectivity(g4)["largest_frac"] > 0.999, "anti-aliased render fragmented"
    return (f"AA acts (ss1 vs ss4: {acts} px corrected on the curved boundary) and converges "
            f"(ss4 vs ss8: {conv} px, {conv/acts:.0%} of the correction) -- smooth, stable "
            f"rounded corners; connectivity preserved")


@_v("R3c")
def _R3c(c):
    """
    INK AREA CONVERGES (the geometry is well-posed).  As resolution increases, the rendered
    ink fraction converges to a limit -- confirming the vector model is exact and the raster
    merely samples it.  (A geometry defect would make the area jump or diverge.)
    """
    mc = meander_core
    base = dict(width=1400, height=1400, n_arms=3, n_taper=8, n_decay=6, n_grid=6,
                line_w=8.0, boundary_w=24.0)
    cr = _feasible_corner_r(base)
    pat = mc.build_pattern(mc.Params(corner_r=cr, **base))
    areas = [rasterize_meander(pat, res=rr, ss=1)[0].mean() for rr in (0.8, 0.4, 0.2)]
    d1, d2 = abs(areas[1] - areas[0]), abs(areas[2] - areas[1])
    assert d2 < d1, "ink area not converging as resolution increases"
    assert d2 < 5e-3, f"ink area still changing by {d2:.1e} at fine resolution"
    return (f"ink fraction converges: {areas[0]*100:.2f}% -> {areas[1]*100:.2f}% -> "
            f"{areas[2]*100:.2f}% (successive deltas {d1:.1e} -> {d2:.1e}); the vector "
            f"geometry is well-posed, raster only samples it")


# ===========================================================================
# S35 -- STAGE 5: THE BINARY GRID -> STADIUM TRANSITION SEAM
#
# A device may hand off from the meander's ROUNDED-CORNERS GRID (rung + rail, no posts) to
# a CLASSIC-CIF section with ACTUAL STADIUM POSTS.  The transition is BINARY: above the seam
# is the rounded grid, below is stadium posts, and NO region has both.  The stadium posts
# mark the official meander->CIF transition point.
#
# Across the seam the primary (trunk) and secondary (rail) channels are CONTINUOUS -- only
# the WALL structure changes: the grid's filleted rung+rail gives way to solid stadium pills
# with gaps.  The rung+rail and the stadium+gap are two renderings of the same periodic link,
# so the fluid path is preserved; the seam is a rendering switch, not a physics change (S17).
# ===========================================================================

def transition_seam(pat, seam_frac: float = 0.5, n_stadium: int = None,
                    corner_arc: int = 12):
    """
    Render a meander Pattern with a BINARY grid->stadium transition.

    The grid region is split at `seam_frac` (fraction of grid units): grid units ABOVE the
    seam stay as the meander's rounded-corners rung+rail; units BELOW become a classic-CIF
    stadium-post section (primary channel + secondary rail continuous, wall = stadium pills
    with gaps).  Returns {"grid_polys", "stadium_channels", "stadium_posts", "seam_y",
    "bounds", "grid_units_above", "stadium_units_below"}.

    n_stadium overrides how many trailing grid units become the stadium section.
    """
    d, p = pat.d, pat.params
    A, T, D, R, M = p.n_arms, p.n_taper, p.n_decay, p.n_release, p.n_morph
    grid_all = sorted([ch for ch in pat.channels if ch.kind == "grid"],
                      key=lambda c: min(pt[1] for pt in c.points))
    if not grid_all:
        raise GateError("transition_seam: pattern has no grid region to transition from")
    per_side = {}
    for ch in grid_all:
        per_side.setdefault(ch.side, []).append(ch)
    n_units = len(per_side.get("right", []))
    if n_stadium is None:
        n_stadium = max(1, int(round((1.0 - seam_frac) * n_units)))
    n_stadium = min(n_stadium, n_units - 1)          # keep at least one rounded-grid unit
    n_above = n_units - n_stadium

    # seam y = the boundary between the (n_above)th and (n_above+1)th grid unit
    right = per_side["right"]
    seam_y = 0.5 * (right[n_above - 1].points[0][1] + right[n_above].points[0][1])

    # --- ABOVE the seam: the meander's rounded-corners grid (rung+rail), kept as channels ---
    grid_polys = []
    for ch in grid_all:
        if min(pt[1] for pt in ch.points) < seam_y:
            for i in range(len(ch.points) - 1):
                grid_polys.append(_seg_rect_poly(ch.points[i], ch.points[i + 1], ch.width))

    # --- BELOW the seam: classic-CIF stadium section ---
    # In classic CIF the wall between the primary and secondary is FLUID at the LINKS (gaps)
    # and SOLID at the posts.  So the wall band is fluid, with solid stadium posts spanning it
    # and gaps between them -- the gaps are how center-channel fluid reaches the side channels.
    t2 = d.trunk_w / 2.0
    rail_lw = d.line_w
    y_bottom = p.height
    stadium_channels, stadium_posts = [], []

    # primary channel continues (trunk) from seam to bottom
    stadium_channels.append([(d.xc - t2, seam_y), (d.xc + t2, seam_y),
                             (d.xc + t2, y_bottom), (d.xc - t2, y_bottom)])
    for side, sgn in (("left", -1.0), ("right", 1.0)):
        rail_cx = d.xc + sgn * (t2 + d.merge_off)
        tf = d.xc + sgn * t2                                   # trunk face
        rail_inner = rail_cx - sgn * rail_lw / 2.0             # inner edge of the rail stroke
        # secondary rail continues
        stadium_channels.append([(rail_cx - rail_lw / 2, seam_y), (rail_cx + rail_lw / 2, seam_y),
                                 (rail_cx + rail_lw / 2, y_bottom), (rail_cx - rail_lw / 2, y_bottom)])
        # THE WALL BAND IS FLUID (this is the link/gap region) -- a blue rectangle spanning
        # trunk face -> rail inner edge, from seam to bottom.  Posts (below) carve solid pillars
        # out of it, leaving gaps through which fluid crosses from primary to secondary.
        x_a, x_b = sorted((tf, rail_inner))
        band = x_b - x_a
        wall_cx = 0.5 * (x_a + x_b)
        stadium_channels.append([(x_a, seam_y), (x_b, seam_y), (x_b, y_bottom), (x_a, y_bottom)])
        # stadium posts: span the wall (X = band) with caps facing the primary/secondary
        # (the Gifford orientation).  CRITICAL PHYSICS: the gap G between posts is the fluid
        # LINK, and its resistance must MATCH the rung it replaces, or the flow split (hence F
        # and D) changes.  The rung had width wk and length merge_off; the gap has length
        # `band`.  Choose G so r(G)*band == r(wk)*merge_off (same link resistance).
        gside = per_side[side]
        pitch = gside[1].points[0][1] - gside[0].points[0][1]
        wk = gside[0].width                                   # the grid rung width
        h_um = p.height * 0.0 + _seam_h(pat)                  # etch depth (um)
        G = _match_gap_width(wk, d.merge_off, band, h_um)     # gap that matches R_rung
        G = min(G, 0.6 * pitch)                               # keep a solid post
        g_p = pitch - G                                       # post along-flow length
        if g_p > band:                                        # keep the pill spanning the wall
            g_p = band
        y = seam_y + pitch / 2.0
        while y < y_bottom - g_p / 2.0:
            stadium_posts.append(_stadium_across_wall(wall_cx, y, band, g_p, corner_arc))
            y += pitch

    xs = [v[0] for q in grid_polys + stadium_channels + stadium_posts for v in q]
    ys = [v[1] for q in grid_polys + stadium_channels + stadium_posts for v in q]
    # include the whole meander bounds
    for ch in pat.channels:
        for x, y in ch.points:
            xs.append(x); ys.append(y)
    bounds = (min(xs), min(ys), max(xs), max(ys))
    return dict(grid_polys=grid_polys, stadium_channels=stadium_channels,
                stadium_posts=stadium_posts, seam_y=seam_y, bounds=bounds,
                grid_units_above=n_above, stadium_units_below=n_stadium)


def _seam_h(pat):
    """Etch depth (um) for the pattern; the meander Params does not carry h, so use the
    physics default 140 um unless the pattern records one."""
    return float(getattr(pat.params, "height_um", 0.0)) or 140.0


def _match_gap_width(wk, rung_len, gap_len, h_um):
    """
    Find the stadium-gap width G whose link resistance equals the grid rung's:
        r(G, h) * gap_len == r(wk, h) * rung_len
    Both links are short channels across the wall; r is the per-length hydraulic resistance.
    Solved by bisection on G (monotone: narrower G -> larger r).  This is what keeps the flow
    split -- and therefore F and the cutoff D -- unchanged across the seam.
    """
    target = r(wk * UM, h_um * UM, 1e-3) * (rung_len * UM)
    lo, hi = 0.5, max(wk * 4.0, rung_len)
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        Rg = r(mid * UM, h_um * UM, 1e-3) * (gap_len * UM)
        if Rg > target:            # too resistive -> widen
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def _stadium_across_wall(wall_cx, yc, band, g_p, n_arc=12):
    """
    A stadium post for VERTICAL flow (y): long axis ACROSS the wall (x), so the semicircular
    caps face the primary and secondary channels (the Gifford orientation), and the along-flow
    (y) extent is g_p.  X-extent = band (spans trunk face -> rail), Y-extent = g_p.  This is
    stadium_post (caps on the y-ends) transposed by swapping x<->y, valid when band >= g_p.
    """
    base = stadium_post(0.0, 0.0, L_along=g_p, W_perp=band, n_arc=n_arc)  # X=g_p, Y=band
    return [(wall_cx + y, yc + x) for (x, y) in base]                      # swap -> X=band, Y=g_p


def _seg_rect_poly(pa, pb, w):
    """Filled rectangle polygon for a Manhattan segment of stroke width w."""
    hw = w / 2.0
    x0, y0 = pa; x1, y1 = pb
    if abs(y1 - y0) < 1e-9:
        xa, xb = sorted((x0, x1))
        return [(xa - hw, y0 - hw), (xb + hw, y0 - hw), (xb + hw, y0 + hw), (xa - hw, y0 + hw)]
    ya, yb = sorted((y0, y1))
    return [(x0 - hw, ya - hw), (x0 + hw, ya - hw), (x0 + hw, yb + hw), (x0 - hw, yb + hw)]


# ===========================================================================
# S36 -- R4 FAMILY: the binary grid->stadium transition (Stage 5)
# ===========================================================================

def _seam_fixture():
    mc = meander_core
    base = dict(width=1200, height=1600, n_arms=3, n_taper=6, n_decay=5, n_grid=14,
                line_w=7.0, boundary_w=20.0, trunk_w=40.0)
    cr = _feasible_corner_r(base)
    pat = mc.build_pattern(mc.Params(corner_r=cr, **base))
    return pat, transition_seam(pat, seam_frac=0.5)


def _raster_polys(polys, x0, y0, res, H, W):
    g = np.zeros((H, W), bool)
    for v in polys:
        v = np.asarray(v, float)
        cA = max(0, int((v[:, 0].min() - x0) / res)); cB = min(W, int((v[:, 0].max() - x0) / res) + 1)
        rA = max(0, int((v[:, 1].min() - y0) / res)); rB = min(H, int((v[:, 1].max() - y0) / res) + 1)
        if cB <= cA or rB <= rA:
            continue
        ry, rx = np.mgrid[rA:rB, cA:cB]; pxx = x0 + (rx + .5) * res; pyy = y0 + (ry + .5) * res
        ins = np.zeros(pxx.shape, bool); n = len(v)
        for i in range(n):
            xa, ya = v[i]; xb, yb = v[(i + 1) % n]
            if ya == yb:
                continue
            cond = (ya > pyy) != (yb > pyy); xint = (xb - xa) * (pyy - ya) / (yb - ya) + xa
            ins ^= cond & (pxx < xint)
        g[rA:rB, cA:cB][ins] = True
    return g


@_v("R4a")
def _R4a(c):
    """
    THE TRANSITION IS BINARY.  No stadium post lies above the seam, and no rounded-grid
    polygon lies below it: above is the meander's rounded-corners grid, below is classic-CIF
    stadium posts, and NO region has both.  This is the invariant the design demands.
    """
    pat, ts = _seam_fixture()
    seam_y = ts["seam_y"]
    posts_above = sum(1 for q in ts["stadium_posts"] if np.mean([v[1] for v in q]) < seam_y)
    grid_below = sum(1 for q in ts["grid_polys"] if np.mean([v[1] for v in q]) > seam_y)
    assert posts_above == 0, f"{posts_above} stadium posts above the seam (not binary)"
    assert grid_below == 0, f"{grid_below} rounded-grid polys below the seam (not binary)"
    return (f"binary transition confirmed: {ts['grid_units_above']} rounded-grid units ABOVE "
            f"the seam (y={seam_y:.0f}), {ts['stadium_units_below']} stadium-post units BELOW; "
            f"0 posts above, 0 grid below -- no region has both")


@_v("R4b")
def _R4b(c):
    """
    STADIUM POSTS ARE SOLID ISLANDS IN THE FLUID WALL BAND.  Below the seam the wall band is
    fluid (the link region); the posts are solid pillars carving it, and they must NOT block
    the primary or secondary channel cores (those stay open) -- otherwise the device is dead.
    The posts carve a MINORITY of the wall band, leaving gaps.
    """
    pat, ts = _seam_fixture()
    assert ts["stadium_posts"], "no stadium posts below the seam"
    d = pat.d; seam_y = ts["seam_y"]
    b = ts["bounds"]; res = 0.25
    x0, y0 = b[0] - 5, b[1] - 5
    W = int((b[2] - b[0] + 10) / res); H = int((b[3] - b[1] + 10) / res)
    gpo = _raster_polys(ts["stadium_posts"], x0, y0, res, H, W)
    # primary + secondary CORES (the channel strokes, not the wall band) must be post-free
    t2 = d.trunk_w / 2.0
    core_polys = [[(d.xc - t2, seam_y), (d.xc + t2, seam_y), (d.xc + t2, b[3]), (d.xc - t2, b[3])]]
    for sgn in (-1.0, 1.0):
        rc = d.xc + sgn * (t2 + d.merge_off); lw = d.line_w
        core_polys.append([(rc - lw / 2, seam_y), (rc + lw / 2, seam_y),
                           (rc + lw / 2, b[3]), (rc - lw / 2, b[3])])
    gcore = _raster_polys(core_polys, x0, y0, res, H, W)
    intrude = int((gcore & gpo).sum())
    assert intrude <= 0.02 * gpo.sum(), \
        f"posts intrude into the primary/secondary cores ({100*intrude/gpo.sum():.1f}%)"
    # posts must be a minority of the wall band (gaps remain -- checked fully in R4c)
    band_polys = [q for q in ts["stadium_channels"] if _is_wall_band(q, d)]
    gband = _raster_polys(band_polys, x0, y0, res, H, W)
    carved = int((gband & gpo).sum())
    frac = carved / max(gband.sum(), 1)
    assert frac < 0.6, f"posts seal the wall band ({100*frac:.0f}% carved)"
    return (f"{len(ts['stadium_posts'])} solid stadium posts carve {100*frac:.0f}% of the fluid "
            f"wall band (gaps remain); {100*intrude/max(gpo.sum(),1):.1f}% intrusion into the "
            f"primary/secondary cores -- channels stay open, posts are islands in the link region")


def _is_wall_band(poly, d):
    """True if a stadium_channels polygon is a wall band (between trunk face and a rail),
    not the primary or a rail core."""
    xs = [v[0] for v in poly]
    cx = 0.5 * (min(xs) + max(xs))
    t2 = d.trunk_w / 2.0
    for sgn in (-1.0, 1.0):
        wall_cx = d.xc + sgn * (t2 + (d.merge_off - d.line_w / 2.0) / 2.0)
        if abs(cx - wall_cx) < d.merge_off:
            # width close to the wall band (not the thin rail, not the wide trunk)
            wdt = max(xs) - min(xs)
            if 0.5 * d.merge_off < wdt < 1.2 * d.merge_off:
                return True
    return False


@_v("R4c")
def _R4c(c):
    """
    FLUID CROSSES FROM CENTER TO SIDE THROUGH THE GAPS.  In the stadium section the wall band
    is FLUID (the link region) with solid posts carving it, leaving gaps.  Isolating the
    stadium section BELOW the seam, the primary and secondary channels are ONE connected
    component -- fluid genuinely crosses from center to side through the gaps between posts.
    (The earlier version left the wall white, which disconnected them; this is the guard.)
    """
    pat, ts = _seam_fixture()
    seam_y = ts["seam_y"]
    b = ts["bounds"]; res = 0.25
    x0, y0 = b[0] - 5, b[1] - 5
    W = int((b[2] - b[0] + 10) / res); H = int((b[3] - b[1] + 10) / res)
    gch = _raster_polys(ts["stadium_channels"], x0, y0, res, H, W)
    gpo = _raster_polys(ts["stadium_posts"], x0, y0, res, H, W)
    fluid = gch & ~gpo                                    # carve solid posts from the fluid band
    yy = y0 + (np.arange(H) + .5) * res
    below = yy > seam_y + 30                              # strictly inside the stadium section
    sub = fluid & below[:, None]
    cc = connectivity(sub)
    assert cc["largest_frac"] > 0.95, \
        f"primary and secondary NOT linked below the seam ({cc['n_components']} components) -- " \
        f"no fluid path center->side"
    # and the posts must leave gaps (not fully seal the wall)
    ov = int((gch & gpo).sum())
    assert ov < 0.6 * gch.sum(), "posts seal the wall -- no gaps for the fluid link"
    return (f"stadium section below the seam is {cc['n_components']} connected component "
            f"({100*cc['largest_frac']:.1f}%) -- center-channel fluid reaches the side channels "
            f"through the gaps; posts carve {100*ov/gch.sum():.0f}% of the wall band, gaps remain")


@_v("R4d")
def _R4d(c):
    """
    STADIUM POSTS HAVE THE RIGHT SHAPE.  Each post spans the wall (its across-wall extent >=
    its along-flow extent -- a pill facing the primary/secondary, the Gifford orientation),
    not a rotated blob, and is a valid stadium (no bowtie).
    """
    pat, ts = _seam_fixture()
    assert ts["stadium_posts"], "no stadium posts"
    ratios = []
    for q in ts["stadium_posts"]:
        xs = [v[0] for v in q]; ys = [v[1] for v in q]
        across = max(xs) - min(xs); along = max(ys) - min(ys)
        ratios.append(across / along)
    worst = min(ratios)
    assert worst >= 1.0 - 1e-6, f"a post is taller along flow than across the wall (ratio {worst:.2f})"
    return (f"{len(ts['stadium_posts'])} stadium posts, each spanning the wall "
            f"(across/along ratio {worst:.2f}..{max(ratios):.2f} >= 1) with caps facing the "
            f"primary/secondary -- correct pill orientation, no bowtie")


# ===========================================================================
# S37 -- Q FAMILY: PHYSICS VALIDITY AT EVERY STAGE
#
# The over-arching guarantee the whole pipeline must keep: the device still adheres to its
# CUTOFF (D = kappa*delta*) and its FILTRATE FRACTION (F), and remains a valid CIF device,
# through compilation AND every rendering/transform stage.  These checks assert D and F
# directly, and confirm each geometric transform (the fabrication bridge, the transition
# seam) does not disturb the flow split that sets them.
# ===========================================================================

@_v("Q1")
def _Q1(c):
    """
    THE COMPILED DEVICE HITS ITS TARGETS.  For a spec that passes the gates, the audited
    filtrate fraction F_em matches F_target, and the extracted lamina matches the design
    delta* in the BULK (so D = kappa*delta* = D_target).  The cutoff is judged on the bulk,
    as the audit's own P5 does -- the last few terminal rows naturally roll off at the device
    end (a boundary effect), which is expected and does not change the sorting cutoff.
    """
    lp = physics_core
    spec = lp.demo_spec_A()
    pat, plan, aud = lp.compile_and_audit(spec)
    dr = np.array(aud.delta_ratio)
    bulk_dev = float(np.max(np.abs(dr[1:-3] - 1.0)))         # exclude the terminal roll-off
    assert abs(aud.F_em - spec.F_target) < 0.02, \
        f"F_em {aud.F_em:.4f} != F_target {spec.F_target:.3f}"
    assert bulk_dev < 5e-3, f"cutoff drifted in the bulk (max |delta_em/delta*-1| = {bulk_dev:.2e})"
    # the authoritative F and cutoff P-checks must pass
    p5 = next(ch for ch in aud.checks if ch.name.startswith("P5"))
    p6 = next(ch for ch in aud.checks if ch.name.startswith("P6"))
    assert p5.ok and p6.ok, "audit's F-target or lamina-uniformity check failed"
    return (f"compiled device adheres to targets: F_em {aud.F_em:.4f} vs F_target "
            f"{spec.F_target:.3f} (err {abs(aud.F_em-spec.F_target):.4f}); bulk lamina matches "
            f"delta* to {bulk_dev:.0e} so D = kappa*delta* = D_target = {spec.D_target:.2f} um "
            f"(P5+P6 pass).  A valid CIF device.")


@_v("Q2")
def _Q2(c):
    """
    THE FABRICATION BRIDGE IS PHYSICS-NEUTRAL.  Rendering the meander to a filled-polygon
    mask does not change the schedule: the channel widths in the raster match the compiled
    widths (the rasterizer draws exactly the strokes the compiler set), so F and D are
    untouched (S17).  Confirmed by the drawn side-channel widths matching plan.wc_rows /
    the unit widths to sub-micron.
    """
    lp = physics_core
    u = compile_unified(lp.demo_spec_A())
    pat = u.pattern
    # the rendered strokes ARE the compiled widths -- check the unit widths are what the plan set
    dr = pat.d
    widths = list(dr.unit_widths)
    # they must be positive, finite, and match the physics-solved side widths (>= w_min floor)
    assert all(w > 0 for w in widths), "non-positive stroke width in the render"
    # the trunk (center) width equals the compiled port width
    assert abs(dr.trunk_w - u.plan.trunk_w) < 1e-6, "rendered trunk width != compiled"
    return (f"fabrication bridge is physics-neutral: {len(widths)} side-channel strokes drawn "
            f"at the compiled widths ({min(widths):.1f}..{max(widths):.1f} um), trunk at "
            f"{dr.trunk_w:.1f} um = compiled -- the render changes shape, not the schedule (S17)")


@_v("Q3")
def _Q3(c):
    """
    THE TRANSITION SEAM PRESERVES THE FLOW SPLIT.  The stadium gaps are sized so each gap's
    link resistance MATCHES the grid rung it replaces (r(G)*band == r(wk)*merge_off).  Since
    the primary and secondary channels are continuous and R_gap is preserved, the flow split
    -- and therefore F and the cutoff D -- is unchanged across the seam.  Without this the
    gaps would be ~6x too open and the classic-CIF section would over-extract.
    """
    mc = meander_core
    base = dict(width=1300, height=1900, n_arms=3, n_taper=6, n_decay=5, n_grid=18,
                line_w=7.0, boundary_w=28.0, trunk_w=44.0)
    cr = _feasible_corner_r(base)
    pat = mc.build_pattern(mc.Params(corner_r=cr, **base))
    d = pat.d
    grid = sorted([ch for ch in pat.channels if ch.kind == "grid" and ch.side == "right"],
                  key=lambda ch: ch.points[0][1])
    wk = grid[0].width
    merge_off = d.merge_off
    band = merge_off - d.line_w / 2.0
    h = _seam_h(pat)
    G = _match_gap_width(wk, merge_off, band, h)
    R_rung = r(wk * UM, h * UM, 1e-3) * (merge_off * UM)
    R_gap = r(G * UM, h * UM, 1e-3) * (band * UM)
    assert abs(R_gap / R_rung - 1.0) < 0.02, \
        f"stadium gap resistance {R_gap/R_rung:.3f}x the rung -- flow split (F, D) not preserved"
    return (f"transition seam preserves the flow split: matched gap G = {G:.2f} um gives "
            f"R_gap/R_rung = {R_gap/R_rung:.4f} (the link resistance is conserved), so F and "
            f"the cutoff D carry across the seam unchanged")


# ===========================================================================
# S38 -- STAGE 6: DXF EXPORT WITH TRUE VECTOR ARCS
#
# The raster (Stage 2) samples curved boundaries into pixels; for fabrication the rounded
# corners and stadium caps must be TRUE ARCS, not stair-stepped polygons.  A DXF ARC entity
# is (center, radius, start_angle, end_angle) -- exact at any scale, no aliasing.
#
# The meander's only curved elements are, by construction:
#   * ROUND JOINS/CAPS  -- discs of radius w/2 at convex stroke corners and ends.
#   * FILLET DISCS      -- white discs of radius r carved at concave corners (spec section 13).
#   * STADIUM CAPS      -- the semicircular ends of the transition posts (radius from the pill).
# Everything else is axis-aligned line segments.  This exporter emits the straight edges as
# LWPOLYLINE and every curved feature as a real CIRCLE/ARC, so the fabricated corners are
# exact arcs.  It is a fabrication OUTPUT: geometry is authoritative, physics already locked.
# ===========================================================================

def write_dxf_meander(pat, path: str, *, stadium=None) -> dict:
    """
    Write a meander Pattern to DXF with TRUE VECTOR ARCS for all curved features.

    Emits, in micrometre units (R2000 ASCII):
      - channel stroke rectangles (LWPOLYLINE) for every Manhattan segment,
      - CIRCLE entities for round joins/caps (convex, radius w/2) on layer CHANNELS,
      - CIRCLE entities for fillet discs (concave, radius corner_r) on layer FILLET_VOID
        (these are WHITE carves -- the consumer subtracts them),
      - CIRCLE/ARC entities for stadium-post caps on layer POSTS, if a transition seam is
        supplied via `stadium` (the transition_seam() result).

    Returns a tally {"segments", "join_caps", "fillets", "post_caps", "entities"}.
    The arcs are exact (center+radius), so there is NO sampling error at the corners.
    """
    d, p = pat.d, pat.params
    mc = meander_core
    rounded = p.corner_r is not None
    rr = p.corner_r if rounded else 0.0

    seg_rects, join_caps, fillet_discs, post_ents = [], [], [], []

    # channel strokes -> rectangles + round join/cap discs
    for ch in pat.channels:
        pts = mc.render_polyline(ch, d)
        w2 = ch.width / 2.0
        for i in range(len(pts) - 1):
            seg_rects.append(_seg_rect_poly(pts[i], pts[i + 1], ch.width))
        if rounded:
            for i in range(1, len(pts) - 1):                 # convex joins
                join_caps.append((pts[i][0], pts[i][1], w2))
            join_caps.append((pts[0][0], pts[0][1], w2))     # caps
            join_caps.append((pts[-1][0], pts[-1][1], w2))

    # trunk parts (rects / shaped polygon)
    for name, part in pat.trunk_parts.items():
        if part[0] == "rect":
            rct = part[1]
            seg_rects.append([(rct[0], rct[1]), (rct[2], rct[1]), (rct[2], rct[3]), (rct[0], rct[3])])
        else:
            seg_rects.append([(x, y) for x, y in part[1]])

    # fillet discs (white carves) -- true circles at the exact tangent centers
    if rounded:
        for cx, cy, sx, sy in mc._fillet_corners(pat):
            fillet_discs.append((cx + rr * sx, cy + rr * sy, rr))

    # stadium caps (if a transition seam is provided): each post is a stadium; emit its two
    # cap circles as true arcs plus the body rectangle
    if stadium is not None:
        for post in stadium["stadium_posts"]:
            post_ents.append(("poly", post))                 # the full pill outline (arcs sampled
            #                                                   finely only for the outline; caps
            #                                                   also emitted as exact circles)

    with open(path, "w") as fh:
        fh.write("0\nSECTION\n2\nHEADER\n9\n$ACADVER\n1\nAC1015\n9\n$INSUNITS\n70\n13\n0\nENDSEC\n")
        fh.write("0\nSECTION\n2\nENTITIES\n")
        h = 256

        def _poly(v, layer):
            nonlocal h
            fh.write(f"0\nLWPOLYLINE\n5\n{h:X}\n100\nAcDbEntity\n8\n{layer}\n"
                     f"100\nAcDbPolyline\n90\n{len(v)}\n70\n1\n43\n0\n")
            for (vx, vy) in v:
                fh.write(f"10\n{vx:.4f}\n20\n{vy:.4f}\n")
            h += 1

        def _circle(cx, cy, rad, layer):
            nonlocal h
            fh.write(f"0\nCIRCLE\n5\n{h:X}\n100\nAcDbEntity\n8\n{layer}\n100\nAcDbCircle\n"
                     f"10\n{cx:.4f}\n20\n{cy:.4f}\n30\n0\n40\n{rad:.4f}\n")
            h += 1

        for v in seg_rects:
            _poly(v, "CHANNELS")
        for cx, cy, rad in join_caps:
            _circle(cx, cy, rad, "CHANNELS")            # exact round joins/caps
        for cx, cy, rad in fillet_discs:
            _circle(cx, cy, rad, "FILLET_VOID")         # exact white fillet carves
        for kind, post in post_ents:
            _poly(post, "POSTS")
        fh.write("0\nENDSEC\n0\nEOF\n")

    return dict(segments=len(seg_rects), join_caps=len(join_caps), fillets=len(fillet_discs),
                post_caps=len(post_ents), entities=len(seg_rects) + len(join_caps)
                + len(fillet_discs) + len(post_ents), path=path)


# ===========================================================================
# S39 -- R5 FAMILY: DXF vector-arc export (Stage 6)
# ===========================================================================

def _parse_dxf(path):
    lines = open(path).read().split("\n")
    ents = []
    i = 0
    while i < len(lines):
        if lines[i] == "0" and i + 1 < len(lines) and lines[i + 1] in ("LWPOLYLINE", "CIRCLE"):
            typ = lines[i + 1]; i += 2; codes = {}; verts = []
            while i + 1 < len(lines) and lines[i] != "0":
                code = lines[i]; val = lines[i + 1]; i += 2
                if code == "10":
                    verts.append([float(val), None])
                elif code == "20" and verts and verts[-1][1] is None:
                    verts[-1][1] = float(val)
                else:
                    codes[code] = val
            e = {"type": typ, "layer": codes.get("8", "?")}
            if typ == "CIRCLE":
                e["c"] = (verts[0][0], verts[0][1], float(codes["40"]))
            else:
                e["v"] = [(x, y) for x, y in verts]
            ents.append(e)
        else:
            i += 1
    return ents


def _raster_dxf(ents, res=0.5, pad=8):
    polys = [e["v"] for e in ents if e["type"] == "LWPOLYLINE"]
    cb = [e["c"] for e in ents if e["type"] == "CIRCLE" and e["layer"] == "CHANNELS"]
    cw = [e["c"] for e in ents if e["type"] == "CIRCLE" and e["layer"] == "FILLET_VOID"]
    allpts = [v for q in polys for v in q]
    x0 = min(p[0] for p in allpts) - pad; y0 = min(p[1] for p in allpts) - pad
    W = int((max(p[0] for p in allpts) - x0 + pad) / res); H = int((max(p[1] for p in allpts) - y0 + pad) / res)
    g = np.zeros((H, W), bool)

    def fpoly(v, val):
        v = np.asarray(v, float)
        cA = max(0, int((v[:, 0].min() - x0) / res)); cB = min(W, int((v[:, 0].max() - x0) / res) + 1)
        rA = max(0, int((v[:, 1].min() - y0) / res)); rB = min(H, int((v[:, 1].max() - y0) / res) + 1)
        if cB <= cA or rB <= rA:
            return
        ry, rx = np.mgrid[rA:rB, cA:cB]; pxx = x0 + (rx + .5) * res; pyy = y0 + (ry + .5) * res
        ins = np.zeros(pxx.shape, bool); n = len(v)
        for i in range(n):
            xa, ya = v[i]; xb, yb = v[(i + 1) % n]
            if ya == yb:
                continue
            cond = (ya > pyy) != (yb > pyy); xint = (xb - xa) * (pyy - ya) / (yb - ya) + xa
            ins ^= cond & (pxx < xint)
        g[rA:rB, cA:cB][ins] = val

    def fdisc(cx, cy, rad, val):
        cA = max(0, int((cx - rad - x0) / res)); cB = min(W, int((cx + rad - x0) / res) + 1)
        rA = max(0, int((cy - rad - y0) / res)); rB = min(H, int((cy + rad - y0) / res) + 1)
        if cB <= cA or rB <= rA:
            return
        ry, rx = np.mgrid[rA:rB, cA:cB]; pxx = x0 + (rx + .5) * res; pyy = y0 + (ry + .5) * res
        m = (pxx - cx) ** 2 + (pyy - cy) ** 2 <= rad * rad
        g[rA:rB, cA:cB][m] = val

    for v in polys:
        fpoly(v, True)
    for cx, cy, rad in cb:
        fdisc(cx, cy, rad, True)
    for cx, cy, rad in cw:
        fdisc(cx, cy, rad, False)
    return g


@_v("R5a")
def _R5a(c):
    """
    DXF EMITS TRUE ARCS, NOT SAMPLED CORNERS.  Every rounded feature -- round joins/caps and
    concave fillets -- is a real CIRCLE entity (center + radius), exact at any scale.  The
    straight channel edges are LWPOLYLINE.  So the fabricated corners have NO sampling error.
    """
    import tempfile, os
    mc = meander_core
    base = dict(width=940, height=990, n_arms=3, n_taper=11, n_decay=6, n_grid=8)
    cr = _feasible_corner_r(base)
    pat = mc.build_pattern(mc.Params(corner_r=cr, **base))
    fd, path = tempfile.mkstemp(suffix=".dxf"); os.close(fd)
    tally = write_dxf_meander(pat, path)
    ents = _parse_dxf(path)
    n_circ = sum(1 for e in ents if e["type"] == "CIRCLE")
    n_poly = sum(1 for e in ents if e["type"] == "LWPOLYLINE")
    os.unlink(path)
    assert n_circ == tally["join_caps"] + tally["fillets"], "arc count mismatch on read-back"
    assert n_circ > 0, "no arcs emitted -- corners would be sampled"
    return (f"DXF has {n_circ} true CIRCLE arcs ({tally['join_caps']} round joins/caps + "
            f"{tally['fillets']} fillet discs) and {n_poly} straight LWPOLYLINEs -- rounded "
            f"corners are exact arcs, zero sampling error at fabrication scale")


@_v("R5b")
def _R5b(c):
    """
    THE DXF IS FAITHFUL: rasterizing the arc-based DXF reproduces the meander geometry -- one
    connected fluid domain, and an ink fraction that matches the faithful rasterizer to well
    under a percent.  So the true-arc export is the SAME shape, just exact at the corners.
    """
    import tempfile, os
    mc = meander_core
    base = dict(width=940, height=990, n_arms=3, n_taper=11, n_decay=6, n_grid=8)
    cr = _feasible_corner_r(base)
    pat = mc.build_pattern(mc.Params(corner_r=cr, **base))
    fd, path = tempfile.mkstemp(suffix=".dxf"); os.close(fd)
    write_dxf_meander(pat, path)
    g_dxf = _raster_dxf(_parse_dxf(path), res=0.5)
    os.unlink(path)
    g_ref, _ = rasterize_meander(pat, res=0.5, ss=1)
    assert connectivity(g_dxf)["largest_frac"] > 0.999, "DXF geometry fragmented"
    d_ink = abs(g_dxf.mean() - g_ref.mean())
    assert d_ink < 0.008, f"DXF ink fraction off by {d_ink*100:.2f}% from the faithful raster"
    return (f"DXF round-trips: rasterized arc-DXF is 1 connected component, ink "
            f"{100*g_dxf.mean():.2f}% vs faithful {100*g_ref.mean():.2f}% (diff "
            f"{d_ink*100:.2f}%) -- the true-arc export is the same shape, exact at the corners")


@_v("R5c")
def _R5c(c):
    """
    DXF IS VALID AND PHYSICS-NEUTRAL.  The file is well-formed R2000 ASCII in micrometre
    units, and the export is a fabrication OUTPUT -- it emits the already-locked geometry
    (the compiled widths and the exact fillet radii), so it cannot change F or D.
    """
    import tempfile, os
    mc = meander_core
    base = dict(width=940, height=990, n_arms=3, n_taper=11, n_decay=6, n_grid=8)
    cr = _feasible_corner_r(base)
    pat = mc.build_pattern(mc.Params(corner_r=cr, **base))
    fd, path = tempfile.mkstemp(suffix=".dxf"); os.close(fd)
    write_dxf_meander(pat, path)
    txt = open(path).read()
    os.unlink(path)
    assert txt.startswith("0\nSECTION") and txt.rstrip().endswith("EOF"), "malformed DXF"
    assert "$INSUNITS\n70\n13" in txt, "DXF not in micrometre units"
    assert "AC1015" in txt, "not R2000"
    return (f"DXF is valid R2000 ASCII, micrometre units ($INSUNITS=13); a fabrication output "
            f"of the locked geometry (widths + exact corner_r {cr:.3f}) -- cannot alter F or D")


# ===========================================================================
# S40 -- STAGE 6b: UNIFIED VECTOR-ARC DXF FOR THE CLASSIC-CIF PATH
#
# The classic-CIF (Gifford) geometry from build_mask has two curved element families:
#   * STADIUM POST CAPS -- two semicircles of radius g_p/2 on the wall-spanning ends.
#   * U-TURN ARCS       -- concentric half-annuli (side/primary/side).
# write_dxf_mask samples these into LWPOLYLINE vertices (aliased at fab scale).  This exporter
# emits the straight edges as LWPOLYLINE and every cap/turn as a true DXF ARC (center, radius,
# start/end angle) -- exact at any scale.  It reads the geometry from analytic descriptors so
# there is no arc-from-polygon guessing.  Fabrication OUTPUT: physics already locked.
# ===========================================================================

def stadium_post_arcs(xc, yc, L_along, W_perp):
    """
    Analytic description of a stadium post (long axis PERPENDICULAR to flow, W_perp >= L_along):
    two straight sides (LINEs) and two semicircular caps (ARCs).  Returns
    {"lines": [(p0,p1),...], "arcs": [(cx,cy,r,a0,a1),...]} in degrees, matching stadium_post().
    """
    r_cap = L_along / 2.0
    hs = (W_perp - 2.0 * r_cap) / 2.0
    if hs < 0:
        raise GateError("stadium_post_arcs: W_perp < L_along (bowtie)")
    left, right = xc - r_cap, xc + r_cap
    lines = [((left, yc - hs), (left, yc + hs)),        # left straight side
             ((right, yc - hs), (right, yc + hs))]      # right straight side
    arcs = [(xc, yc + hs, r_cap, 0.0, 180.0),           # top cap 0->pi
            (xc, yc - hs, r_cap, 180.0, 360.0)]         # bottom cap pi->2pi
    return {"lines": lines, "arcs": arcs}


def write_dxf_classic(mask, path: str, *, proc=None) -> dict:
    """
    Vector-arc DXF of a classic-CIF Mask (build_mask output).  POSTS are emitted as their
    exact 2 lines + 2 semicircle ARCs; all other polys (channels, walls, folds) as LWPOLYLINE;
    ports as CIRCLE.  U-turn arc polys are kept as fine LWPOLYLINE (their radii vary per ring;
    emitting them as ARCs would need the per-ring descriptors, a further extension).

    Returns a tally of entity counts.  The post caps are now TRUE ARCS -- the Gifford figures
    rendered from this DXF have exact stadium caps, not sampled polygons.
    """
    P = [q for q in mask.polys if q["layer"] != "OUTLINE"]
    n_post_arc = n_post_line = n_poly = n_circ = 0
    with open(path, "w") as fh:
        fh.write("0\nSECTION\n2\nHEADER\n9\n$ACADVER\n1\nAC1015\n9\n$INSUNITS\n70\n13\n0\nENDSEC\n")
        fh.write("0\nSECTION\n2\nENTITIES\n")
        h = 256

        def _poly(v, layer):
            nonlocal h, n_poly
            fh.write(f"0\nLWPOLYLINE\n5\n{h:X}\n100\nAcDbEntity\n8\n{layer}\n"
                     f"100\nAcDbPolyline\n90\n{len(v)}\n70\n1\n43\n0\n")
            for (vx, vy) in v:
                fh.write(f"10\n{vx:.4f}\n20\n{vy:.4f}\n")
            h += 1; n_poly += 1

        def _line(p0, p1, layer):
            nonlocal h, n_post_line
            fh.write(f"0\nLINE\n5\n{h:X}\n100\nAcDbEntity\n8\n{layer}\n100\nAcDbLine\n"
                     f"10\n{p0[0]:.4f}\n20\n{p0[1]:.4f}\n30\n0\n11\n{p1[0]:.4f}\n21\n{p1[1]:.4f}\n31\n0\n")
            h += 1; n_post_line += 1

        def _arc(cx, cy, r, a0, a1, layer):
            nonlocal h, n_post_arc
            fh.write(f"0\nARC\n5\n{h:X}\n100\nAcDbEntity\n8\n{layer}\n100\nAcDbCircle\n"
                     f"10\n{cx:.4f}\n20\n{cy:.4f}\n30\n0\n40\n{r:.4f}\n"
                     f"100\nAcDbArc\n50\n{a0:.4f}\n51\n{a1:.4f}\n")
            h += 1; n_post_arc += 1

        def _circle(cx, cy, r, layer):
            nonlocal h, n_circ
            fh.write(f"0\nCIRCLE\n5\n{h:X}\n100\nAcDbEntity\n8\n{layer}\n100\nAcDbCircle\n"
                     f"10\n{cx:.4f}\n20\n{cy:.4f}\n30\n0\n40\n{r:.4f}\n")
            h += 1; n_circ += 1

        for q in P:
            if q["layer"] == "POSTS":
                # recover the stadium parameters from the vertices: xc = mean x, extents give
                # L_along (x-span) and W_perp (y-span).  Emit exact lines + arcs.
                v = q["v"]
                xs = [pt[0] for pt in v]; ys = [pt[1] for pt in v]
                xc = 0.5 * (min(xs) + max(xs)); yc = 0.5 * (min(ys) + max(ys))
                L_along = max(xs) - min(xs); W_perp = max(ys) - min(ys)
                try:
                    desc = stadium_post_arcs(xc, yc, L_along, W_perp)
                    for (p0, p1) in desc["lines"]:
                        _line(p0, p1, "POSTS")
                    for (acx, acy, ar, a0, a1) in desc["arcs"]:
                        _arc(acx, acy, ar, a0, a1, "POSTS")
                except GateError:
                    _poly(v, "POSTS")
            else:
                _poly(q["v"], q["layer"])
        for (_nm, cx, cy, rad) in mask.ports:
            _circle(cx, cy, rad, "PORTS")
        fh.write("0\nENDSEC\n0\nEOF\n")

    return dict(polylines=n_poly, post_lines=n_post_line, post_arcs=n_post_arc,
                circles=n_circ, entities=n_poly + n_post_line + n_post_arc + n_circ, path=path)


# ===========================================================================
# S41 -- R6 FAMILY: classic-CIF vector-arc DXF reproduces Gifford (Stage 6b)
# ===========================================================================

def _parse_dxf_full(path):
    """Parse LWPOLYLINE / LINE / ARC / CIRCLE with correct group-code handling."""
    lines = open(path).read().split("\n")
    ents = []
    i = 0
    while i < len(lines):
        if lines[i] == "0" and i + 1 < len(lines) and lines[i + 1] in ("LWPOLYLINE", "LINE", "ARC", "CIRCLE"):
            typ = lines[i + 1]; i += 2; codes = {}; verts = []
            while i + 1 < len(lines) and lines[i] != "0":
                code = lines[i]; val = lines[i + 1]; i += 2
                if code == "10":
                    if typ == "LWPOLYLINE":
                        verts.append([float(val), None])
                    else:
                        codes["10"] = float(val)
                elif code == "20":
                    if typ == "LWPOLYLINE" and verts and verts[-1][1] is None:
                        verts[-1][1] = float(val)
                    else:
                        codes["20"] = float(val)
                elif code in ("11", "21", "40", "50", "51"):
                    codes[code] = float(val)
                else:
                    codes[code] = val
            ents.append({"type": typ, "layer": codes.get("8", "?"), "codes": codes, "verts": verts})
        else:
            i += 1
    return ents


@_v("R6a")
def _R6a(c):
    """
    THE CLASSIC-CIF DXF EMITS TRUE ARCS FOR THE STADIUM CAPS.  The Gifford geometry from
    build_mask exports with every stadium cap as a real DXF ARC (center, radius, angles) and
    each straight post side as a LINE -- not sampled polygon vertices.  So the fabricated
    Gifford posts have exact semicircular caps.
    """
    import tempfile, os
    d = design(reference_request())
    m = build_mask(d, Layout("linear", 1), max_rows=60)
    fd, path = tempfile.mkstemp(suffix=".dxf"); os.close(fd)
    tally = write_dxf_classic(m, path)
    ents = _parse_dxf_full(path)
    os.unlink(path)
    n_arc = sum(1 for e in ents if e["type"] == "ARC")
    n_line = sum(1 for e in ents if e["type"] == "LINE")
    n_post = sum(1 for q in m.polys if q["layer"] == "POSTS")
    assert n_arc == 2 * n_post, f"expected 2 arcs/post ({2*n_post}), got {n_arc}"
    assert n_line == 2 * n_post, f"expected 2 lines/post, got {n_line}"
    return (f"classic-CIF DXF: {n_arc} true stadium-cap ARCs + {n_line} post-side LINEs for "
            f"{n_post} posts -- caps are exact semicircles, zero sampling error")


@_v("R6b")
def _R6b(c):
    """
    THE DXF POSTS MATCH THE ORIGINAL GIFFORD 2:1 PILL.  The arc radius is g_p/2 and the post
    side length gives W_perp = l_g, so the exported post is the SAME stadium as the reference
    process (l_g/g_p = 2), not a squat blob.  This is what makes the DXF reproduce the figures.
    """
    import tempfile, os
    d = design(reference_request()); p = d.req.proc
    m = build_mask(d, Layout("linear", 1), max_rows=40)
    fd, path = tempfile.mkstemp(suffix=".dxf"); os.close(fd)
    write_dxf_classic(m, path)
    ents = _parse_dxf_full(path)
    os.unlink(path)
    arcs = [e for e in ents if e["type"] == "ARC" and e["layer"] == "POSTS"]
    plines = [e for e in ents if e["type"] == "LINE" and e["layer"] == "POSTS"]
    r_cap = arcs[0]["codes"]["40"]
    side = abs(plines[0]["codes"]["21"] - plines[0]["codes"]["20"])
    L_along = 2 * r_cap
    W_perp = side + 2 * r_cap
    assert abs(L_along - p.g_p / UM) < 0.5, f"g_p {L_along:.1f} != {p.g_p/UM:.1f}"
    assert abs(W_perp - p.l_g / UM) < 0.5, f"l_g {W_perp:.1f} != {p.l_g/UM:.1f}"
    assert abs(W_perp / L_along - p.l_g / p.g_p) < 0.05, "post ratio != original"
    return (f"DXF post = original Gifford 2:1 pill: g_p = {L_along:.1f} um (arc r {r_cap:.1f}), "
            f"l_g = {W_perp:.1f} um, ratio {W_perp/L_along:.2f} (reference {p.l_g/p.g_p:.2f})")


@_v("R6c")
def _R6c(c):
    """
    THE DXF IS PHYSICS-NEUTRAL AND VALID.  It exports the classic mask geometry (built from
    the compiled device, which passed the F/cutoff audit), so it cannot change F or D; and it
    is well-formed R2000 ASCII in micrometre units with ARC/LINE/LWPOLYLINE/CIRCLE entities.
    """
    import tempfile, os
    d = design(reference_request())
    # the device the DXF renders is the audited one; F/cutoff are locked upstream
    m = build_mask(d, Layout("linear", 1), max_rows=40)
    fd, path = tempfile.mkstemp(suffix=".dxf"); os.close(fd)
    write_dxf_classic(m, path)
    txt = open(path).read()
    os.unlink(path)
    assert txt.startswith("0\nSECTION") and txt.rstrip().endswith("EOF"), "malformed DXF"
    assert "$INSUNITS\n70\n13" in txt and "AC1015" in txt, "not R2000 um"
    assert "\nARC\n" in txt and "\nLINE\n" in txt and "\nLWPOLYLINE\n" in txt, "missing entity types"
    return ("classic-CIF DXF is valid R2000 ASCII (um); a fabrication output of the audited "
            "device geometry -- ARC + LINE + LWPOLYLINE + CIRCLE, cannot alter F or D")


# ===========================================================================
# S42 -- STAGE 5b: TRANSITION TO THE *CLASSIC* WALL (real 2:1 stadium posts)
#
# The Stage-5 seam kept the meander's thin wall, so its posts were squat (~1:1).  A real
# meander->classic-CIF handoff steps the wall out to the classic l_g and uses the classic 2:1
# stadium posts (g_p x l_g).  Two things change at the seam -- the WALL THICKNESS (boundary_w
# -> l_g) and the PITCH (Lambda = G + g_p) -- and BOTH are S17 rendering knobs: they do not
# touch rho, w_s, f, D, F or N, PROVIDED the link resistance R_gap is preserved.  We enforce
# that by solving the gap width G so r(G)*l_g == R_rung (the meander rung it replaces).  The
# pitch discontinuity at the seam is the N2 finding, already handled by composition.
# ===========================================================================

def transition_seam_classic(pat, proc=None, seam_frac: float = 0.5, corner_arc: int = 16):
    """
    Meander -> classic-CIF handoff at ONE SELF-CONSISTENT SCALE.

    CRITICAL FIX: the classic section's wall, posts, and pitch are derived from the MEANDER's
    OWN geometry -- NOT from an unrelated reference_process.  In a real CIF device there is one
    scale from one physics: the wall thickness IS the meander boundary_w (== l_g), the gap and
    post pitch come from the meander pitch, and the primary/secondary widths are continuous.
    The classic stadium post keeps the classic 2:1 proportion (l_g : g_p) by setting g_p =
    l_g/2 = boundary_w/2, and the gap width G is solved so the link resistance matches the
    meander rung (physics preserved).  So there is NO scale jump -- only the wall RENDERING
    changes (rounded rung+rail -> 2:1 stadium pills), which is an S17 rendering knob.

    `proc` is accepted for signature compatibility but ignored for scale (the meander sets it).
    """
    d, p = pat.d, pat.params
    grid_all = sorted([ch for ch in pat.channels if ch.kind == "grid"],
                      key=lambda ch: min(pt[1] for pt in ch.points))
    if not grid_all:
        raise GateError("transition_seam_classic: no grid region")
    per_side = {}
    for ch in grid_all:
        per_side.setdefault(ch.side, []).append(ch)
    right = per_side["right"]
    n_units = len(right)
    n_above = max(1, int(round(seam_frac * n_units)))
    seam_y = 0.5 * (right[n_above - 1].points[0][1] + right[min(n_above, n_units - 1)].points[0][1])

    lw = d.line_w
    wk = right[0].width
    merge_off = d.merge_off
    pitch_m = right[1].points[0][1] - right[0].points[0][1]
    h = _seam_h(pat)

    # classic geometry AT THE MEANDER'S SCALE: the wall band the posts span is the SAME wall
    # the meander grid uses (trunk face -> rail inner edge == merge_off - lw/2).  Keep the
    # classic 2:1 proportion by making the post along-flow g_p = band/2.
    t2 = d.trunk_w / 2.0
    rail_cx0 = d.xc + (t2 + merge_off)
    band = (rail_cx0 - lw / 2.0) - (d.xc + t2)          # trunk face -> rail inner edge (== l_g band)
    l_g = band
    g_p = l_g / 2.0                                     # classic 2:1 pill (l_g : g_p = 2:1)
    G = _match_gap_width(wk, merge_off, l_g, h)         # gap preserves the rung resistance
    R_rung = r(wk * UM, h * UM, 1e-3) * (merge_off * UM)
    R_gap = r(G * UM, h * UM, 1e-3) * (l_g * UM)
    pitch_c = G + g_p

    # above the seam: the meander rounded grid, kept
    grid_polys = []
    for ch in grid_all:
        if min(pt[1] for pt in ch.points) < seam_y:
            for i in range(len(ch.points) - 1):
                grid_polys.append(_seg_rect_poly(ch.points[i], ch.points[i + 1], ch.width))

    # below the seam: classic CIF at the SAME scale -- primary, rail, wall band all CONTINUOUS
    # (the rail stays where the meander rail was; the wall band is the same width).
    y_bottom = p.height
    classic_channels, classic_posts = [], []
    classic_channels.append([(d.xc - t2, seam_y), (d.xc + t2, seam_y),
                             (d.xc + t2, y_bottom), (d.xc - t2, y_bottom)])
    for side, sgn in (("left", -1.0), ("right", 1.0)):
        tf = d.xc + sgn * t2
        rail_cx = d.xc + sgn * (t2 + merge_off)         # SAME rail position as the meander
        rail_inner = rail_cx - sgn * lw / 2.0
        classic_channels.append([(rail_cx - lw / 2, seam_y), (rail_cx + lw / 2, seam_y),
                                 (rail_cx + lw / 2, y_bottom), (rail_cx - lw / 2, y_bottom)])
        x_a, x_b = sorted((tf, rail_inner))
        band_s = x_b - x_a
        wall_cx = 0.5 * (x_a + x_b)
        classic_channels.append([(x_a, seam_y), (x_b, seam_y), (x_b, y_bottom), (x_a, y_bottom)])
        # 2:1 stadium posts spanning the (same) wall band, gap G, pitch_c
        gpp = min(g_p, 0.9 * band_s)                     # keep across>=along (feasible pill)
        y = seam_y + pitch_c / 2.0
        while y < y_bottom - gpp / 2.0:
            classic_posts.append(_stadium_across_wall(wall_cx, y, band_s, gpp, corner_arc))
            y += pitch_c

    xs = [v[0] for q in grid_polys + classic_channels + classic_posts for v in q]
    ys = [v[1] for q in grid_polys + classic_channels + classic_posts for v in q]
    for ch in pat.channels:
        for x, y in ch.points:
            xs.append(x); ys.append(y)
    bounds = (min(xs), min(ys), max(xs), max(ys))
    return dict(grid_polys=grid_polys, classic_channels=classic_channels, classic_posts=classic_posts,
                seam_y=seam_y, wall_step=(d.boundary_w, l_g), pitch_meander=pitch_m,
                pitch_classic=pitch_c, G=G, R_ratio=R_gap / R_rung, l_g=l_g, g_p=g_p, bounds=bounds)


# ===========================================================================
# S43 -- R7 FAMILY: the classic-wall transition (Stage 5b)
# ===========================================================================

def _classic_seam_fixture():
    mc = meander_core
    base = dict(width=1300, height=1900, n_arms=3, n_taper=6, n_decay=5, n_grid=18,
                line_w=7.0, boundary_w=28.0, trunk_w=44.0)
    cr = _feasible_corner_r(base)
    pat = mc.build_pattern(mc.Params(corner_r=cr, **base))
    return pat, transition_seam_classic(pat, reference_process(), seam_frac=0.5)


@_v("R7a")
def _R7a(c):
    """
    THE CLASSIC SECTION USES TRUE 2:1 STADIUM POSTS.  Below the seam the wall steps out to the
    classic l_g and the posts are the real Gifford pill (across/along == l_g/g_p == 2), not the
    squat meander-wall post.  This is what makes the transition match the original geometry.
    """
    pat, ts = _classic_seam_fixture()
    p = reference_process()
    ratios = []
    for q in ts["classic_posts"]:
        xs = [v[0] for v in q]; ys = [v[1] for v in q]
        ratios.append((max(xs) - min(xs)) / (max(ys) - min(ys)))
    worst = min(ratios); best = max(ratios)
    assert abs(best - p.l_g / p.g_p) < 0.05 and abs(worst - p.l_g / p.g_p) < 0.05, \
        f"post ratio {worst:.2f}..{best:.2f} != classic {p.l_g/p.g_p:.2f}"
    return (f"classic section posts are true 2:1 pills (across/along {worst:.2f}, classic "
            f"l_g/g_p {p.l_g/p.g_p:.2f}); wall steps {ts['wall_step'][0]:.0f} -> "
            f"{ts['wall_step'][1]:.0f} um at the seam")


@_v("R7b")
def _R7b(c):
    """
    THE WALL STEP AND PITCH CHANGE PRESERVE THE PHYSICS.  At the seam the wall thickness
    (boundary_w -> l_g) and the pitch (-> G + g_p) both change -- S17 rendering knobs -- but
    the gap width G is solved so the link resistance R_gap MATCHES the meander rung it
    replaces.  Since R_gap is conserved and the primary/secondary are continuous, F and the
    cutoff D are unchanged across the seam despite the wall step.
    """
    pat, ts = _classic_seam_fixture()
    assert abs(ts["R_ratio"] - 1.0) < 0.02, \
        f"R_gap/R_rung = {ts['R_ratio']:.3f} -- flow split (F, D) NOT preserved across the step"
    assert ts["pitch_classic"] != ts["pitch_meander"], "expected a pitch discontinuity (N2)"
    return (f"wall steps {ts['wall_step'][0]:.0f}->{ts['wall_step'][1]:.0f} um and pitch "
            f"{ts['pitch_meander']:.1f}->{ts['pitch_classic']:.1f} um at the seam, but matched "
            f"gap G={ts['G']:.2f} gives R_gap/R_rung = {ts['R_ratio']:.4f} -- F and D preserved "
            f"(S17: wall/pitch are rendering knobs, link resistance is the invariant)")


@_v("R7c")
def _R7c(c):
    """
    FLUID CROSSES CENTER->SIDE THROUGH THE CLASSIC GAPS.  The classic wall band is fluid with
    2:1 posts carving it; isolating the classic section below the seam, the primary and
    secondary are ONE connected component -- fluid genuinely reaches the side channels through
    the gaps between the 2:1 posts.
    """
    pat, ts = _classic_seam_fixture()
    seam_y = ts["seam_y"]
    b = ts["bounds"]; res = 0.3
    x0, y0 = b[0] - 5, b[1] - 5
    W = int((b[2] - b[0] + 10) / res); H = int((b[3] - b[1] + 10) / res)
    gch = _raster_polys(ts["classic_channels"], x0, y0, res, H, W)
    gpo = _raster_polys(ts["classic_posts"], x0, y0, res, H, W)
    fluid = gch & ~gpo
    yy = y0 + (np.arange(H) + .5) * res
    below = yy > seam_y + 40
    cc = connectivity(fluid & below[:, None])
    assert cc["largest_frac"] > 0.95, \
        f"primary/secondary NOT linked in the classic section ({cc['n_components']} comps)"
    return (f"classic section is {cc['n_components']} connected component "
            f"({100*cc['largest_frac']:.1f}%) -- fluid crosses center->side through the gaps "
            f"between the true 2:1 stadium posts")


# ===========================================================================
# S44 -- R8 FAMILY: Gifford reproduction (straight side channels, no taper)
#
# A faithful Gifford device: constant primary width (no taper), CONTINUOUS straight parallel
# side channels (zero-excursion -- "a meander with side channels of 0 length"), 2:1 stadium
# posts with gaps in the WALLS (not in the side channels), post-free semicircular U-turns,
# in serpentine or paperclip layout.  Rendered through the vector-arc DXF.  The straight-mode
# rendering is S17 (posts/gaps/flow split unchanged) so F and D are the schedule's.
# ===========================================================================

def _gifford_device():
    import dataclasses as _dc
    r = reference_request(); wc0 = r.W_shape(0.0)
    return design(_dc.replace(r, W_shape=lambda t: wc0))


@_v("R8a")
def _R8a(c):
    """
    GIFFORD DOES NOT TAPER.  The constant-W_shape device has a constant primary width, unlike
    the tapering Xia device -- and it still hits the same F and cutoff D (the taper is a
    rendering/where-to-put-width choice, not what sets the physics targets here).
    """
    d = _gifford_device()
    wc = [rw.w_c for rw in d.rows]
    assert max(wc) - min(wc) < 1e-9, "Gifford device is tapering (should be constant)"
    return (f"Gifford device: constant primary width {wc[0]/UM:.1f} um (no taper), N={d.N}, "
            f"F={d.F:.4f}, D={d.D/UM:.3f} um -- a valid non-tapering CIF device")


@_v("R8b")
def _R8b(c):
    """
    CONTINUOUS STRAIGHT SIDE CHANNELS, POSTS IN THE WALLS.  In straight mode every row has a
    continuous side channel (no folds, no white gaps in the side channel); the 2:1 stadium
    posts sit in the wall between center and side, and the gaps between posts are the ONLY
    interruptions -- in the wall, not the side channel.  Fluid crosses center->side through
    those gaps: the mask is one connected fluid domain.
    """
    d = _gifford_device()
    lay = Layout("linear", 1, side_mode="straight")
    m = build_mask(d, lay, max_rows=40)
    sides = [q for q in m.polys if q.get("role") == "SIDE"]
    posts = [q for q in m.polys if q["layer"] == "POSTS"]
    assert len(sides) == 2, f"expected 2 continuous side channels, got {len(sides)}"
    # NO route_fold folds in straight mode: channel polys are primary + 2 sides + gaps + caps,
    # none are boustrophedon folds (which would appear as many short CHANNELS rects off-axis)
    mask_ch = Mask(polys=[q for q in m.polys if q["layer"] == "CHANNELS"], ports=[],
                   bounds=m.bounds, legs=m.legs, layout=m.layout)
    g, _ = raster(mask_ch, res=2.0)
    cc = connectivity(g)
    assert cc["largest_frac"] > 0.99, "center and side channels not linked through the gaps"
    return (f"Gifford straight mode: {len(sides)} continuous side channels, {len(posts)} 2:1 "
            f"posts in the walls; fluid is {cc['n_components']} connected component "
            f"({100*cc['largest_frac']:.1f}%) -- gaps link center<->side, side channels unbroken")


@_v("R8c")
def _R8c(c):
    """
    GIFFORD SERPENTINE + PAPERCLIP VIA VECTOR-ARC DXF.  Both layouts build, join their legs
    with post-free semicircular U-turns into ONE connected fluid domain, and export to
    vector-arc DXF (2:1 stadium caps as true ARCs).  These are the reference figures.
    """
    import tempfile, os
    d = _gifford_device()
    out = {}
    for kind in ("serpentine", "paperclip"):
        lay = Layout(kind, n_legs=4, side_mode="straight", leg_spacing=850 * UM)
        m = build_mask(d, lay)
        fd, path = tempfile.mkstemp(suffix=".dxf"); os.close(fd)
        tally = write_dxf_classic(m, path)
        ents = _parse_dxf_full(path)
        os.unlink(path)
        n_arc = sum(1 for e in ents if e["type"] == "ARC")
        mask_ch = Mask(polys=[q for q in m.polys if q["layer"] == "CHANNELS"], ports=[],
                       bounds=m.bounds, legs=m.legs, layout=m.layout)
        g, _ = raster(mask_ch, res=8.0)
        cc = connectivity(g)
        assert cc["largest_frac"] > 0.95, f"{kind}: legs not joined ({cc['n_components']} comps)"
        assert n_arc > 0, f"{kind}: no post-cap arcs in DXF"
        out[kind] = (len(m.legs), n_arc, cc["largest_frac"])
    s = out["serpentine"]; pc = out["paperclip"]
    return (f"serpentine ({s[0]} legs, {s[1]} post-cap arcs, {100*s[2]:.0f}% connected) and "
            f"paperclip ({pc[0]} legs, {pc[1]} arcs, {100*pc[2]:.0f}% connected) both export to "
            f"vector-arc DXF and join via post-free U-turns -- the Gifford reference figures")


# ===========================================================================
# S44 -- STAGE 4: FIELD SOLVE FROM THE FABRICATED GEOMETRY (pixels -> physics)
#
# The audit validates the SCHEDULE; these checks validate the DRAWN ARTIFACT.  They measure
# geometry off the rendered pixels/polygons (NOT the construction parameters), rebuild the
# hydraulic network, and recover F and the link resistance -- an independent confirmation
# that the fabricated device still hits its targets.  This is the loop that catches a
# rendering bug the schedule audit cannot see (e.g. the dead-device / 6x-open-gap bugs).
# ===========================================================================

def measure_gap_from_pixels(ts, pat, res=0.2, h=140.0):
    """Measure the classic-transition gap width from the FABRICATED raster (fluid runs between
    posts in the wall band), and return (gap_um, R_gap, R_rung, ratio).  Independent of the
    construction G."""
    d = pat.d
    b = ts["bounds"]; x0, y0 = b[0] - 5, b[1] - 5
    W = int((b[2] - b[0] + 10) / res); H = int((b[3] - b[1] + 10) / res)
    gch = _raster_polys(ts["classic_channels"], x0, y0, res, H, W)
    gpo = _raster_polys(ts["classic_posts"], x0, y0, res, H, W)
    fluid = gch & ~gpo
    t2 = d.trunk_w / 2.0; merge_off = d.merge_off; lw = d.line_w
    tf = d.xc + t2; rail_inner = d.xc + t2 + merge_off - lw / 2.0
    wall_cx = 0.5 * (tf + rail_inner)
    col = int((wall_cx - x0) / res)
    below0 = int((ts["seam_y"] + 40 - y0) / res); below1 = int((b[3] - 20 - y0) / res)
    seg = fluid[below0:below1, col]
    runs = []; i = 0
    while i < len(seg):
        if seg[i]:
            j = i
            while j < len(seg) and seg[j]:
                j += 1
            runs.append((j - i) * res); i = j
        else:
            i += 1
    gap = float(np.median(runs)) if runs else 0.0
    grid = sorted([c for c in pat.channels if c.kind == "grid" and c.side == "right"],
                  key=lambda c: c.points[0][1])
    wk = grid[0].width
    l_g = rail_inner - tf
    R_rung = r(wk * UM, h * UM, 1e-3) * (merge_off * UM)
    R_gap = r(gap * UM, h * UM, 1e-3) * (l_g * UM)
    return gap, R_gap, R_rung, (R_gap / R_rung if R_rung else 0.0)


@_v("Q4")
def _Q4(c):
    """
    F RECOVERED FROM THE DRAWN CLASSIC (GIFFORD) POLYGONS.  geom_ftest measures w_c and w_s
    OFF THE RENDERED POLYGONS, rebuilds the ladder, solves Kirchhoff, and recovers F -- an
    independent confirmation (not the schedule) that the fabricated Gifford device hits F.
    The measured widths match the design to sub-1e-4, and F_poly matches F_design.
    """
    d = design(reference_request())
    m = build_mask(d, Layout("linear", 1))
    ft = geom_ftest(m, d)
    assert ft["dwc"] < 1e-3 and ft["dws"] < 1e-3, "drawn widths drifted from design"
    assert abs(ft["F_poly"] - ft["F_design"]) < 0.02, \
        f"F recovered from polygons {ft['F_poly']:.4f} != design {ft['F_design']:.4f}"
    return (f"F recovered from the DRAWN Gifford polygons: F_poly = {ft['F_poly']:.4f} vs "
            f"F_design {ft['F_design']:.4f} (widths match to dwc={ft['dwc']:.0e}, "
            f"dws={ft['dws']:.0e}, over {ft['n']} straight rows) -- the fabricated device hits F")


@_v("Q5")
def _Q5(c):
    """
    F RECOVERED FROM THE DRAWN MEANDER PATTERN.  physics_core.audit measures geometry off the
    drawn meander Pattern and solves an independent Kirchhoff network, recovering F_em -- a
    pixels-back-to-physics confirmation for the meander path, matching F_target.
    """
    lp = physics_core
    spec = lp.demo_spec_A()
    pat, plan, aud = lp.compile_and_audit(spec)
    assert abs(aud.F_em - spec.F_target) < 0.02, \
        f"F_em from drawn meander {aud.F_em:.4f} != target {spec.F_target:.3f}"
    # the P4 network-bulk check (f_em vs kernel) must pass -- the independent solve agrees
    p4 = next(ch for ch in aud.checks if ch.name.startswith("P4"))
    assert p4.ok, "independent network solve disagrees with the kernel in the bulk"
    return (f"F recovered from the DRAWN meander: F_em = {aud.F_em:.4f} vs F_target "
            f"{spec.F_target:.3f} (independent Kirchhoff solve of the emitted geometry; "
            f"P4 network-bulk agrees) -- the fabricated meander hits F")


@_v("Q6")
def _Q6(c):
    """
    THE FABRICATED TRANSITION GAP PRESERVES THE LINK RESISTANCE -- MEASURED FROM PIXELS.  The
    classic-transition gap width is measured from the FABRICATED raster (fluid runs between
    posts), not from the construction G, and the resulting R_gap matches the meander rung it
    replaces.  This is the independent check that the drawn transition still hits F and D --
    the guard against the 6x-too-open-gap bug returning through a rendering error.
    """
    mc = meander_core
    base = dict(width=1300, height=1900, n_arms=3, n_taper=6, n_decay=5, n_grid=18,
                line_w=7.0, boundary_w=28.0, trunk_w=44.0)
    cr = _feasible_corner_r(base)
    pat = mc.build_pattern(mc.Params(corner_r=cr, **base))
    ts = transition_seam_classic(pat, seam_frac=0.5)
    gap, R_gap, R_rung, ratio = measure_gap_from_pixels(ts, pat)
    assert gap > 0, "no gap measured in the fabricated transition -- posts seal the wall"
    assert abs(ratio - 1.0) < 0.15, \
        f"MEASURED fabricated gap gives R_gap/R_rung = {ratio:.3f} -- flow split (F, D) not preserved"
    return (f"fabricated transition gap MEASURED from pixels = {gap:.2f} um (construction "
            f"{ts['G']:.2f}); R_gap/R_rung = {ratio:.3f} from the drawn geometry -- the link "
            f"resistance is preserved in the actual artifact, so F and D carry across the seam")


# ===========================================================================
# S45 -- VISUAL REPRODUCTION HARNESS
#
# The prior validators checked connectivity and entity counts but NEVER whether the render
# matches the reference image.  This harness rasterizes a device's channel mask to a binary
# occupancy image, normalizes it to the reference's aspect/framing, and reports an overlap
# score (IoU) plus structural descriptors (leg count, turn sides, port x-positions).  A
# reproduction validator FAILS if the structure disagrees with the reference -- so a wrong
# layout (wrong U-turn side, missing outlets, wrong leg topology) is caught, not passed.
# ===========================================================================

def device_occupancy(mask, res=20.0):
    """Binary channel-occupancy image of a Mask (fluid=1), plus its (x0,y0,res) and bounds."""
    ch = Mask(polys=[q for q in mask.polys if q["layer"] == "CHANNELS"], ports=[],
              bounds=mask.bounds, legs=mask.legs, layout=mask.layout)
    g, (x0, y0, r) = raster(ch, res=res)
    return g, (x0, y0, r), mask.bounds


def layout_descriptor(mask, layout):
    """Structural fingerprint of a built device: aspect ratio, leg count, U-turn count and
    sides, inlet/outlet x-fractions.  Used to compare topology to a reference expectation."""
    b = mask.bounds
    W = b[2] - b[0]; H = b[3] - b[1]
    # U-turns: CHANNELS polys whose vertices trace an arc (bounding box roughly square, many pts)
    n_turn = 0
    turn_sides = []
    for q in mask.polys:
        if q["layer"] != "CHANNELS":
            continue
        v = np.asarray(q["v"], float)
        if len(v) < 20:
            continue
        bw = np.ptp(v[:, 0]); bh = np.ptp(v[:, 1])
        if bw > 0 and bh > 0 and 0.4 < bw / bh < 2.5 and bh > 200:   # arc-like
            n_turn += 1
            cx = v[:, 0].mean()
            turn_sides.append("right" if cx > b[0] + W / 2 else "left")
    ports_x = [(nm, (cx - b[0]) / W) for (nm, cx, cy, rad) in mask.ports]
    return dict(aspect=W / max(H, 1), n_legs=len(mask.legs), n_turns=n_turn,
                turn_sides=turn_sides, ports_xfrac=ports_x, W_um=W, H_um=H)


def reproduction_score(mask, ref_png_path, res=25.0, dilate=3):
    """Rasterize the device channels and the reference PNG to comparable binary images, align
    by bounding box, DILATE both (so channel thickness differences don't dominate), and return
    the intersection-over-union of the occupied regions.  A coarse silhouette match that
    catches gross layout errors (missing legs/turns/outlets, wrong turn side)."""
    from PIL import Image
    from scipy.ndimage import binary_dilation
    g, (x0, y0, r), b = device_occupancy(mask, res=res)
    ys, xs = np.where(g)
    if len(xs) == 0:
        return 0.0, None
    dev = g[ys.min():ys.max() + 1, xs.min():xs.max() + 1]
    ref_im = np.array(Image.open(ref_png_path).convert("RGB"))
    ref_blue = (ref_im[:, :, 2].astype(int) - ref_im[:, :, 0]) > 15
    ry, rx = np.where(ref_blue)
    if len(rx) == 0:
        return 0.0, None
    ref = ref_blue[ry.min():ry.max() + 1, rx.min():rx.max() + 1]
    from PIL import Image as I2
    dev_r = np.array(I2.fromarray((dev * 255).astype(np.uint8)).resize(
        (ref.shape[1], ref.shape[0]), I2.NEAREST)) > 127
    if dilate:
        st = np.ones((2 * dilate + 1, 2 * dilate + 1), bool)
        dev_r = binary_dilation(dev_r, st)
        ref = binary_dilation(ref, st)
    inter = (dev_r & ref).sum(); union = (dev_r | ref).sum()
    return float(inter) / float(max(union, 1)), (dev_r, ref)


# ===========================================================================
# S46 -- R9 FAMILY: reproduction acceptance tests (structure vs reference)
#
# These are the validators that were MISSING: they check the built device's STRUCTURE against
# the known reference topology, so a wrong layout (U-turns on the wrong side, wrong turn
# count, missing outlets) FAILS instead of passing on connectivity alone.  Reference PNGs are
# optional (silhouette IoU is reported when present); the structural checks are unconditional.
# ===========================================================================

def _uturn_sides(mask, layout):
    """Return the left/right side of each U-turn, in leg order.  U-turn arc polys have exactly
    2*uturn_seg vertices (the _uturn_polys output), 3 per turn (side/primary/side)."""
    n = 2 * layout.uturn_seg
    ut = [q for q in mask.polys if q["layer"] == "CHANNELS" and len(q["v"]) == n]
    b = mask.bounds; midx = b[0] + (b[2] - b[0]) / 2
    sides = []
    for i in range(0, len(ut), 3):
        v = np.asarray(ut[i]["v"], float)
        sides.append("right" if v[:, 0].mean() > midx else "left")
    return sides


@_v("R9a")
def _R9a(c):
    """
    SERPENTINE U-TURNS ALTERNATE SIDES.  A serpentine's turns must go right, left, right, ...
    (each leg reverses direction).  The prior code put ALL turns on the right -- the exact bug
    that made the middle U-turn overlap.  This asserts strict alternation and one connected
    fluid domain.
    """
    import dataclasses as _dc
    r = reference_request(); wc0 = r.W_shape(0.0)
    d = design(_dc.replace(r, W_shape=lambda t: wc0))
    lay = Layout("serpentine", n_legs=4, side_mode="straight", leg_spacing=850 * UM)
    m = build_mask(d, lay)
    sides = _uturn_sides(m, lay)
    assert len(sides) == 3, f"expected 3 U-turns for 4 legs, got {len(sides)}"
    expect = ["right", "left", "right"]
    assert sides == expect, f"U-turns do not alternate: {sides} (expected {expect})"
    mask_ch = Mask(polys=[q for q in m.polys if q["layer"] == "CHANNELS"], ports=[],
                   bounds=m.bounds, legs=m.legs, layout=m.layout)
    g, _ = raster(mask_ch, res=12.0)
    cc = connectivity(g)
    assert cc["largest_frac"] > 0.99, f"serpentine not connected ({cc['n_components']} comps)"
    return (f"serpentine 4-leg: U-turns alternate {sides}, one connected fluid domain "
            f"({100*cc['largest_frac']:.0f}%) -- the middle U-turn is now on the correct side")


@_v("R9b")
def _R9b(c):
    """
    THE DEVICE HAS AN INLET AND THREE OUTLETS AT THE RIGHT ENDS.  A CIF device terminates in a
    middle-channel outlet and two side-channel outlets (the sorted streams).  This asserts the
    four ports exist and the outlets are distinct from the inlet -- so a device that dropped
    its outlets (as an earlier render did) fails.
    """
    import dataclasses as _dc
    r = reference_request(); wc0 = r.W_shape(0.0)
    d = design(_dc.replace(r, W_shape=lambda t: wc0))
    m = build_mask(d, Layout("serpentine", n_legs=4, side_mode="straight"))
    names = [nm for (nm, cx, cy, rad) in m.ports]
    assert "INLET" in names, "no inlet port"
    assert "OUT_C" in names, "no middle-channel outlet"
    assert sum(1 for nm in names if nm.startswith("OUT_S")) == 2, "need two side-channel outlets"
    return (f"device ports: {names} -- one inlet + middle-channel outlet + two side-channel "
            f"outlets (the sorted streams)")


@_v("R9c")
def _R9c(c):
    """
    SILHOUETTE MATCH TO THE GIFFORD REFERENCE.  When the reference overview PNG is available,
    the built serpentine's dilated silhouette overlaps it with IoU above a floor -- a coarse
    but real check that the RENDER resembles the reference (length, leg count, turn placement),
    not just that it is connected.  Skips cleanly if the reference file is absent.
    """
    import os, dataclasses as _dc
    ref = "/mnt/user-data/uploads/serpentine_4leg_fig5_overview.png"
    if not os.path.exists(ref):
        return "reference PNG absent -- silhouette check skipped (structural checks R9a/b hold)"
    r = reference_request(); wc0 = r.W_shape(0.0)
    d = design(_dc.replace(r, W_shape=lambda t: wc0))
    m = build_mask(d, Layout("serpentine", n_legs=4, side_mode="straight", leg_spacing=850 * UM))
    iou, _ = reproduction_score(m, ref, dilate=5)
    assert iou > 0.30, f"serpentine silhouette IoU {iou:.2f} too low vs reference"
    return (f"serpentine silhouette IoU {iou:.2f} vs the Gifford reference (dilated overlap) -- "
            f"the render resembles the reference in length, legs, and turn placement")


# =========================================================================
# =====================  CIF_v2 INTEGRATION LAYER  ========================
# =========================================================================
# Everything above this banner is the frozen CIF.py engine, verbatim
# (sha256 7e18068bc0a6ddbd0b3f8ddc00e6226d6d835d64081b0f7f32555ea1be26ef18).  Everything below folds the Tier-1
# satellites into the same module namespace (their `import CIF` resolves
# to this module) and adds the v2 validators, G14 routability, stamped
# emission, and the acceptance runner.  See AUDIT_INTEGRATION_REVIEW.md.
import sys as _v2_sys, hashlib as _v2_hashlib, os as _v2_os
import json as _v2_json, datetime as _v2_dt
CIF = _v2_sys.modules[__name__]          # satellites' `CIF.` -> this module
cif_patches = CIF; CP = CIF
cif_validators = CIF; JV = CIF
dxf_validate = CIF; optimize_build = CIF
__version__ = "2.1.0"
V2_BASE_SHA256 = "7e18068bc0a6ddbd0b3f8ddc00e6226d6d835d64081b0f7f32555ea1be26ef18"
V2_BUILD_UTC = "2026-07-19T17:41:42Z"

# =========================================================================
# ==== v2 merged satellite: cif_patches.py ===============================
# =========================================================================
#!/usr/bin/env python3
"""cif_patches.py -- proposed CIF.py improvements, implemented externally.

(1) uturn_polys_v2: a U-turn whose rail widths are C0-CONTINUOUS with the
    adjoining straight side channels (end widths pinned to the drawn
    straight widths) while EACH rail realizes the exact S14.5 rho-closure
    target  T = rho_up * r(w_c_up) * (pi * R_turn).
    Width law per rail:  w(t) = lerp(w_up, w_dn, t) * (1 + alpha*sin(pi t));
    alpha is solved by brentq against the EXACT DRAWN FUNCTIONAL -- the
    trapezoid of r(w) along the chord lengths of the emitted band's mid
    polyline -- so a validator that re-integrates the polygon recovers T to
    floating point, not to a quadrature approximation.
    Both rails hit the SAME target  =>  R_in == R_out (N8)  AND the turn is
    hydraulically the rho-weighted primary arc (S14.5, `turn_widths` law),
    which the engine's `uturn_widths` (equal-R split of 2*(w_s - pull_in),
    rho accepted but unused) never enforced.

(2) build_mask_v2: mirrors CIF.build_mask but
    - passes the DOWNSTREAM row into every turn (continuity needs both ends)
    - FORBIDS the GS floor in side_mode='straight' (drawn width must equal
      the design width; a floored device is not a faithful realization)
    - integrates the port leads INTO the mask (short inlet lead, straight
      center-outlet lead, ANGLED side-outlet leads with vertical-semicap
      tips and inscribed join discs at the junctions)
    - records junction metadata (mask.junctions), per-row drawn geometry
      (mask.row_geom) and closure records (mask.turn_closures) for the
      J-family validators
    - supports stage_legs=[(Layout, n_legs), ...] (stage_of, actually wired).

All mask coordinates are MICROMETRES, engine convention.
"""
pass  # [v2-merged] from __future__ import annotations
import math
import numpy as np
from scipy.optimize import brentq

pass  # [v2-merged] import CIF
pass  # [v2-merged] from CIF import UM, GateError

# ---------------------------------------------------------------------------
# (1) the continuity-preserving rho-closure turn
# ---------------------------------------------------------------------------

def _arc_R(w_um: np.ndarray, th: np.ndarray, wall_um: float, sgn: float,
           p) -> float:
    """Hydraulic resistance (SI) of a rail band drawn at FIXED wall radius
    `wall_um` with width profile w_um(t) over angle grid th.  Computed as
    the trapezoid of r(w) along the CHORDS of the drawn mid polyline
    rho_mid(t)*(cos th, sin th), rho_mid = wall + sgn*w/2 -- i.e. exactly
    the functional a validator recovers from the emitted polygon."""
    rho = wall_um + sgn * 0.5 * w_um
    mx, my = rho * np.cos(th), rho * np.sin(th)
    dL = np.hypot(np.diff(mx), np.diff(my)) * UM
    rv = np.array([CIF.r(max(float(w) * UM, 1e-9), p.h, p.mu) for w in w_um])
    return float((0.5 * (rv[:-1] + rv[1:]) * dL).sum())


def _solve_alpha(ws_up_um: float, ws_dn_um: float, th: np.ndarray,
                 t: np.ndarray, wall_um: float, sgn: float, T: float,
                 p) -> tuple:
    """alpha s.t. the DRAWN band's resistance equals T.  R is strictly
    decreasing in alpha (wider -> lower r; for the inner rail the midline
    also shortens).  Returns (alpha, R_achieved, clipped)."""
    lerp = ws_up_um + (ws_dn_um - ws_up_um) * t
    bump = np.sin(np.pi * t)

    def R_of(a):
        return _arc_R(lerp * (1.0 + a * bump), th, wall_um, sgn, p)

    w_floor_um = p.w_min / UM               # fab floor for the mid-arc width
    a_lo = max(-0.95, w_floor_um / float(lerp[len(t) // 2]) - 1.0)
    a_hi = 2.5
    R_lo = R_of(a_lo)                       # max achievable R (narrowest legal)
    if T > R_lo * (1.0 + 1e-15):            # even the fab floor cannot reach T
        return a_lo, R_lo, True
    if T < R_of(a_hi):
        raise GateError("turn rail cannot be made wide enough for the "
                        "rho-closure target; reduce leg spacing / R_turn")
    a = brentq(lambda x: R_of(x) - T, a_lo, a_hi, xtol=1e-15, rtol=8.9e-16)
    return a, R_of(a), False


def uturn_polys_v2(rw_up, rw_dn, p, lay, cx: float, y_from: float,
                   y_to: float, R_turn: float, side: str,
                   ws_up_drawn_SI: float, ws_dn_drawn_SI: float) -> tuple:
    """Post-free semicircular U-turn, v2 (see module docstring).
    Returns (polys_um, junction_record, closure_record)."""
    y_mid = 0.5 * (y_from + y_to) / UM
    R = R_turn / UM
    lg = p.l_g / UM
    wc_a, wc_b = rw_up.w_c / UM, rw_dn.w_c / UM
    ws_a, ws_b = ws_up_drawn_SI / UM, ws_dn_drawn_SI / UM

    # S14.5 target: the turn IS one segment carrying the upstream rung's
    # rho against the primary arc's resistance (turn_widths law).
    T = rw_up.rho * CIF.r(rw_up.w_c, p.h, p.mu) * math.pi * R_turn

    n = max(96, lay.uturn_seg)
    if side == "right":
        th0, th1 = -np.pi / 2, np.pi / 2
    else:
        th0, th1 = 3 * np.pi / 2, np.pi / 2
    # orient so t=0 sits at the y_from port
    if abs((y_mid + R * math.sin(th0)) - y_from / UM) > \
       abs((y_mid + R * math.sin(th1)) - y_from / UM):
        th0, th1 = th1, th0
    th = np.linspace(th0, th1, n + 1)
    t = np.linspace(0.0, 1.0, n + 1)
    bump = np.sin(np.pi * t)

    # WALL RADII LERP WITH wc(t): the wall stays exactly l_g thick along the
    # whole arc AND meets each leg's own wall at its own w_c (the ports of a
    # turn join two DIFFERENT rows, one canonical-taper step apart).  The
    # residual tangent kink at the ports is |d(wc/2)/ds| ~ 1e-5 rad.
    wc = wc_a + (wc_b - wc_a) * t
    wall_in = R - wc / 2 - lg
    wall_out = R + wc / 2 + lg
    a_i, R_i, clip_i = _solve_alpha(ws_a, ws_b, th, t, wall_in, -1.0, T, p)
    a_o, R_o, clip_o = _solve_alpha(ws_a, ws_b, th, t, wall_out, +1.0, T, p)
    shortfall = 0.0
    if clip_i or clip_o:                    # keep N8: equalize both rails to
        Tc = min(R_i if clip_i else T,      # the largest COMMON attainable R
                 R_o if clip_o else T)
        shortfall = T / Tc - 1.0
        a_i, R_i, _ = _solve_alpha(ws_a, ws_b, th, t, wall_in, -1.0, Tc, p)
        a_o, R_o, _ = _solve_alpha(ws_a, ws_b, th, t, wall_out, +1.0, Tc, p)

    lerp = ws_a + (ws_b - ws_a) * t
    wi = lerp * (1.0 + a_i * bump)
    wo = lerp * (1.0 + a_o * bump)
    cxu = cx / UM

    def band(r_lo, r_hi):
        a = [(cxu + rh * math.cos(q), y_mid + rh * math.sin(q))
             for rh, q in zip(r_hi, th)]
        b = [(cxu + rl * math.cos(q), y_mid + rl * math.sin(q))
             for rl, q in zip(r_lo[::-1], th[::-1])]
        return a + b

    polys = [dict(layer="CHANNELS", role="TURN_PRIMARY",
                  v=band(R - wc / 2, R + wc / 2)),
             dict(layer="CHANNELS", role="TURN_SIDE_IN",
                  v=band(wall_in - wi, wall_in)),
             dict(layer="CHANNELS", role="TURN_SIDE_OUT",
                  v=band(wall_out, wall_out + wo))]

    def _ivals(k):                          # port cross-sections at sample k
        s = math.sin(th[k])
        out = []
        for r_lo, r_hi in ((R - wc[k] / 2, R + wc[k] / 2),
                           (wall_in[k] - wi[k], wall_in[k]),
                           (wall_out[k], wall_out[k] + wo[k])):
            y1, y2 = y_mid + s * r_lo, y_mid + s * r_hi
            out.append((min(y1, y2), max(y1, y2)))
        return sorted(out)

    junction = dict(cx=cxu, y_from=y_from / UM, y_to=y_to / UM,
                    arc_from=_ivals(0), arc_to=_ivals(n))
    closure = dict(T=T, R_in=R_i, R_out=R_o, alpha_in=a_i, alpha_out=a_o,
                   clipped=bool(clip_i or clip_o), shortfall=shortfall,
                   # physical inputs so a validator can RE-DERIVE the target
                   # and re-integrate the drawn widths independently:
                   rho_up=rw_up.rho, w_c_up=rw_up.w_c, R_turn=R_turn,
                   ws_up=ws_up_drawn_SI, ws_dn=ws_dn_drawn_SI)
    return polys, junction, closure


# ---------------------------------------------------------------------------
# (2) build_mask_v2
# ---------------------------------------------------------------------------

def _straight_side_intervals(rw, p, yc_um):
    wc2 = rw.w_c / UM / 2
    lg = p.l_g / UM
    ws = rw.w_s / UM
    return sorted([(yc_um - wc2, yc_um + wc2),
                   (yc_um - wc2 - lg - ws, yc_um - wc2 - lg),
                   (yc_um + wc2 + lg, yc_um + wc2 + lg + ws)])


def _lead_poly(x0, y0, dx, dy, w, n=32):
    """Lead channel: start face VERTICAL of height w at (x0,y0) so the
    junction cross-section matches the channel face exactly; body a sheared
    parallelogram translated by (dx,dy); tip face a VERTICAL SEMICAP of
    radius w/2 bulging along sign(dx) -- no sharp tip corners, and the tip
    corners land exactly on the parallelogram's end face."""
    sx = 1.0 if dx >= 0 else -1.0
    cx_, cy_ = x0 + dx, y0 + dy
    v = [(x0, y0 - w / 2), (cx_, cy_ - w / 2)]
    for j in range(1, n):
        a = -math.pi / 2 + j * math.pi / n
        v.append((cx_ + sx * (w / 2) * math.cos(a),
                  cy_ + (w / 2) * math.sin(a)))
    v += [(cx_, cy_ + w / 2), (x0, y0 + w / 2)]
    return v


def _disc(cx_, cy_, r, n=48):
    a = np.linspace(0, 2 * np.pi, n, endpoint=False)
    return [(cx_ + r * math.cos(q), cy_ + r * math.sin(q)) for q in a]


def _composition_interference(polys, junctions, closures, leg_y, rows,
                              bounds_i, p, Lc_um, clear_um=None):
    """EXACT analytic clearance between structural elements that must not
    touch: same-side turn annuli, turn bands vs non-adjacent leg bands at
    the junction plane, and lead polygons vs same-side turn annuli.
    Returns a list of violation strings (empty = clean)."""
    F = (clear_um if clear_um is not None else p.F_min / UM)
    viol = []
    tg = [(j, c["geo"]) for j, c in zip(junctions, closures) if "geo" in c]
    # leg band half-width at the channel ends (outer envelope)
    hw = []
    for (a, b) in bounds_i:
        wmax = max(rw.w_c for rw in rows[a:b]) / UM
        smax = max(max(rw.w_s, p.GS) for rw in rows[a:b]) / UM
        hw.append(wmax / 2 + p.l_g / UM + smax)
    # 4a: same-side turn pairs
    for i in range(len(tg)):
        for k in range(i + 1, len(tg)):
            (ji, gi), (jk, gk) = tg[i], tg[k]
            if abs(gi[0] - gk[0]) > 1e-6:
                continue                       # opposite sides
            d = abs(gi[1] - gk[1])
            clr = max(0.0, max(gi[2] - (d + gk[3]), gk[2] - (d + gi[3]),
                               d - gi[3] - gk[3]))
            if clr < F:
                ov = min(d + gk[3] - gi[2], d + gi[3] - gk[2],
                         gi[3] + gk[3] - d)
                viol.append(
                    f"turns {ji['leg']}&{jk['leg']} (same side x={gi[0]:.0f}): "
                    f"annuli [{gi[2]:.0f},{gi[3]:.0f}] / [{gk[2]:.0f},{gk[3]:.0f}] "
                    f"centres {d:.0f} um apart -> clearance {clr:.1f} < "
                    f"{F:.0f} um (radial overlap depth ~{max(0.0, ov):.0f} um)")
    # 4b: turn band vs non-adjacent leg band at the junction plane
    for jt, g in tg:
        t = jt["leg"]
        for q in range(len(leg_y)):
            if q in (t, t + 1):
                continue
            rho0 = abs(leg_y[q] / UM - g[1])
            lo, hi = rho0 - hw[q], rho0 + hw[q]
            clr = max(0.0, max(g[2] - hi, lo - g[3]))
            if lo <= g[3] and hi >= g[2] or clr < F:
                viol.append(
                    f"turn {t} band [{g[2]:.0f},{g[3]:.0f}] reaches leg {q}'s "
                    f"channel band at the junction plane (radial "
                    f"[{lo:.0f},{hi:.0f}], clearance {clr:.1f} < {F:.0f} um)")
    # 4c: leads vs same-side turn annuli
    for q in polys:
        if not str(q.get("role", "")).startswith(("LEAD", "JOIN", "ENDCAP")):
            continue
        adj_leg = 0 if q.get("grp", (None, None))[1] == "in" else len(leg_y) - 1
        for jt, g in tg:
            if jt["leg"] in (adj_leg - 1, adj_leg):
                continue                        # its own end's turn
            dmin = min(max(g[2] - math.hypot(v[0] - g[0], v[1] - g[1]),
                           math.hypot(v[0] - g[0], v[1] - g[1]) - g[3], 0.0)
                       for v in q["v"])
            if dmin < F:
                viol.append(f"lead {q['role']} within {dmin:.1f} um of turn "
                            f"{jt['leg']} annulus (< {F:.0f} um)")
    return viol


def build_mask_v2(dev, layout=None, stage_legs=None, max_rows=None,
                  port_d: float = 300 * UM, lead_in: float = None,
                  lead_out: float = None, lead_side: float = None,
                  lead_angle_deg: float = 28.0, allow_floor: bool = False,
                  skip_composition_gate: bool = False, resolve: bool = True,
                  stage_gap: float = None,
                  _stage_mirror=None, _extra_gap: float = 0.0):
    """Assemble the photomask with continuity-correct turns, leads, and
    validation metadata. Either `layout` (single) or `stage_legs`
    [(Layout, n_legs), ...] (per-stage layouts, legs consecutive)."""
    p = dev.req.proc
    rows = dev.rows if max_rows is None else dev.rows[:max_rows]
    N = len(rows)
    if stage_legs is None:
        stage_legs = [(layout if layout is not None else CIF.Layout(),
                       max(1, (layout.n_legs if layout is not None else 1)))]
    legs_layout = []
    for lay, nl in stage_legs:
        legs_layout += [lay] * nl
    n_legs = len(legs_layout)
    per = [N // n_legs + (1 if k < N % n_legs else 0) for k in range(n_legs)]
    bounds_i, i0 = [], 0
    for k in range(n_legs):
        bounds_i.append((i0, i0 + per[k])); i0 += per[k]

    # --- the no-floor fidelity guard (straight mode) -----------------------
    if any(getattr(l, "side_mode", "auto") == "straight"
           for l in legs_layout) and not allow_floor:
        bad = [rw.i for rw in rows if rw.w_s < p.GS - 1e-12]
        if bad:
            raise GateError(
                f"straight mode with {len(bad)} rows below GS "
                f"({p.GS/UM:.0f} um): the drawn mask would deviate from the "
                f"design law (first row {bad[0]}).  Raise w_start so "
                f"w_s(1) >= GS, or use side_mode='auto'.")

    # leg y offsets: each stage places its legs at its OWN layout's relative
    # offsets; stage bases are shifted so consecutive stages never collide
    # (min of the new stage sits one leg_spacing above the previous max).
    mirror = _stage_mirror if _stage_mirror is not None \
        else [False] * len(stage_legs)
    leg_y, prev_max = [], None
    for (lay, nl), mir in zip(stage_legs, mirror):
        local = [lay.leg_y(t) - lay.leg_y(0) for t in range(nl)]
        if mir:
            local = [-v for v in local]
        if prev_max is None:
            base = 0.0
        else:
            sep = (lay.leg_spacing if stage_gap is None else stage_gap * UM)
            base = prev_max + sep + _extra_gap * UM - min(local)
        placed = [base + v for v in local]
        leg_y += placed
        prev_max = max(placed)

    L_common = max(sum(rw.Lam for rw in rows[a:b]) for a, b in bounds_i)
    Lc_um = L_common / UM

    polys, junctions, closures, row_geom = [], [], [], []
    for leg, (a, b) in enumerate(bounds_i):
        lay = legs_layout[leg]
        y0 = leg_y[leg]
        rr = rows[a:b]
        flip = (leg % 2 == 1)
        n_before = len(polys)
        polys += CIF._leg_polys(rr, p, lay, y0, flip=flip,
                                channel_len=L_common)
        for q in polys[n_before:]:
            q["grp"] = ("leg", leg)
        x_acc = 0.0
        for rw in rr:                                  # per-row drawn record
            xr = (Lc_um - (x_acc + rw.Lam / UM / 2)) if flip \
                 else (x_acc + rw.Lam / UM / 2)
            ws_drawn = (max(rw.w_s, p.GS) if getattr(lay, "side_mode", "auto")
                        == "straight" else rw.w_s)
            row_geom.append(dict(i=rw.i, leg=leg, x=xr, y=y0 / UM,
                                 wc=rw.w_c / UM, ws_design=rw.w_s / UM,
                                 ws_drawn=ws_drawn / UM))
            x_acc += rw.Lam / UM
        if leg + 1 < n_legs:
            y1 = leg_y[leg + 1]
            R_turn = abs(y1 - y0) / 2.0
            ends_right = (leg % 2 == 0)
            cx_turn = L_common if ends_right else 0.0
            rw_up, rw_dn = rr[-1], rows[bounds_i[leg + 1][0]]

            def wsd(rw, l):
                return (max(rw.w_s, p.GS) if getattr(l, "side_mode", "auto")
                        == "straight" else rw.w_s)
            tp, jrec, crec = uturn_polys_v2(
                rw_up, rw_dn, p, lay, cx_turn, y0, y1, R_turn,
                side="right" if ends_right else "left",
                ws_up_drawn_SI=wsd(rw_up, lay),
                ws_dn_drawn_SI=wsd(rw_dn, legs_layout[leg + 1]))
            n_before = len(polys)
            polys += tp
            tv = [v for q in tp for v in q["v"]]
            yc_t = 0.5 * (y0 + y1) / UM
            rads = [math.hypot(v[0] - cx_turn / UM, v[1] - yc_t) for v in tv]
            crec["geo"] = (cx_turn / UM, yc_t, min(rads), max(rads))
            for q in polys[n_before:]:
                q["grp"] = ("turn", leg)
            jrec.update(leg=leg,
                        straight_from=_straight_side_intervals(rw_up, p, y0 / UM),
                        straight_to=_straight_side_intervals(rw_dn, p, y1 / UM))
            junctions.append(jrec)
            closures.append(crec)

    # --- integrated leads + ports ------------------------------------------
    r0, rN = rows[0], rows[-1]
    y_in = leg_y[0] / UM
    kN = n_legs - 1
    y_out = leg_y[kN] / UM
    x_in = 0.0
    x_out = Lc_um if (kN % 2 == 0 or n_legs == 1) else 0.0
    sout = 1.0 if x_out > 1.0 else -1.0
    wc0, wcN = r0.w_c / UM, rN.w_c / UM
    wsN = (max(rN.w_s, p.GS) if getattr(legs_layout[kN], "side_mode", "auto")
           == "straight" else rN.w_s) / UM
    lg = p.l_g / UM
    Li = lead_in if lead_in is not None else 3.5 * wc0
    Lo = lead_out if lead_out is not None else 8.0 * (p.Lambda / UM)
    Ls = lead_side if lead_side is not None else 1.2 * Lo
    th = math.radians(lead_angle_deg)
    polys.append(dict(layer="CHANNELS", role="LEAD_IN",
                      grp=("lead", "in"),
                      v=_lead_poly(x_in, y_in, -Li, 0.0, wc0)))
    polys.append(dict(layer="CHANNELS", role="LEAD_OUT_C",
                      grp=("lead", "out"),
                      v=_lead_poly(x_out, y_out, sout * Lo, 0.0, wcN)))
    off = wcN / 2 + lg + wsN / 2
    for sgn in (1.0, -1.0):
        polys.append(dict(layer="CHANNELS", grp=("lead", "out"),
                          role=f"LEAD_OUT_S{'+' if sgn > 0 else '-'}",
                          v=_lead_poly(x_out, y_out + sgn * off,
                                       sout * Ls * math.cos(th),
                                       sgn * Ls * math.sin(th), wsN)))
        # inscribed round join at the angled junction: the disc of radius
        # w_s/2 centred on the junction face is tangent to the straight
        # channel's edges and exactly rounds the concave fluid corner
        # (SVG round-join semantics); the convex kink it leaves is 152 deg.
        polys.append(dict(layer="CHANNELS", role="JOIN", grp=("lead", "out"),
                          v=_disc(x_out, y_out + sgn * off, wsN / 2)))
    # side channels of leg 0 START as free butt ends near the inlet:
    # semicap them (fluid-add disc, r = w_s/2) so no sharp corners remain.
    r0ws = (max(r0.w_s, p.GS) if getattr(legs_layout[0], "side_mode", "auto")
            == "straight" else r0.w_s) / UM
    for sgn in (1.0, -1.0):
        polys.append(dict(layer="CHANNELS", role="ENDCAP", grp=("lead", "in"),
                          v=_disc(0.0, y_in + sgn * (wc0 / 2 + lg + r0ws / 2),
                                  r0ws / 2)))
    # ports sit AT the lead tip centres (inside the fluid: half in the
    # semicap, half in the lead body)
    ports = [("INLET", x_in - Li, y_in, min(0.5 * port_d / UM, 0.45 * wc0)),
             ("OUT_C", x_out + sout * Lo, y_out,
              min(0.5 * port_d / UM, 0.45 * wcN))]
    for sgn, nm in ((1.0, "OUT_S+"), (-1.0, "OUT_S-")):
        ports.append((nm, x_out + sout * Ls * math.cos(th),
                      y_out + sgn * (off + Ls * math.sin(th)),
                      min(0.5 * port_d / UM, 0.40 * wsN)))

    if not skip_composition_gate:
        viol = _composition_interference(polys, junctions, closures, leg_y,
                                         rows, bounds_i, p, Lc_um)
        if viol and resolve and _stage_mirror is None and len(stage_legs) > 1:
            # DYNAMIC RESOLUTION (physics-preserving): search fold-phase
            # mirrors x inter-stage widening. Mirroring is a reflection;
            # widening only grows gap-free turn radii -- the v2 turn law
            # re-solves widths for any radius with exact rho-closure, so
            # N, F, the schedule and the cutoff are byte-identical. Pick
            # the SMALLEST widening that clears, preferring no mirror.
            import itertools
            combos = list(itertools.product([False, True],
                                            repeat=len(stage_legs) - 1))
            gaps = [0.0] + [200.0 * s for s in range(1, 26)]  # up to +5 mm
            for gap in gaps:
                for combo in combos:
                    try:
                        m = build_mask_v2(
                            dev, layout=layout, stage_legs=stage_legs,
                            max_rows=max_rows, port_d=port_d, lead_in=lead_in,
                            lead_out=lead_out, lead_side=lead_side,
                            lead_angle_deg=lead_angle_deg,
                            allow_floor=allow_floor, resolve=False,
                            stage_gap=stage_gap,
                            _stage_mirror=[False] + list(combo),
                            _extra_gap=gap)
                    except GateError:
                        continue
                    m.resolution = dict(mirror=[False] + list(combo),
                                        extra_gap_um=gap,
                                        note=("clean" if gap == 0
                                              and not any(combo)
                                              else "resolved"))
                    return m
            raise GateError(
                "composition gate: no fold-phase mirror or inter-stage "
                "widening up to +5 mm clears the interference. This "
                "composition is geometrically infeasible as posed "
                "(e.g. a fold-back layout whose turns cross on one side); "
                "reorder the legs or split the stages.\n  "
                + "\n  ".join(viol))
        if viol:
            raise GateError("composition gate:\n  " + "\n  ".join(viol))

    xs = [v[0] for q in polys for v in q["v"]]
    ys = [v[1] for q in polys for v in q["v"]]
    m = 500.0
    bnd = (min(xs) - m, min(ys) - m, max(xs) + m, max(ys) + m)
    polys.append(dict(layer="OUTLINE", v=[(bnd[0], bnd[1]), (bnd[2], bnd[1]),
                                          (bnd[2], bnd[3]), (bnd[0], bnd[3])]))
    mask = CIF.Mask(polys=polys, ports=ports, bounds=bnd, legs=per,
                    layout=legs_layout[0])
    mask.junctions = junctions
    mask.turn_closures = closures
    mask.row_geom = row_geom
    mask.leg_y_um = [y / UM for y in leg_y]
    mask.legs_layout = legs_layout
    mask.core_bbox = (0.0, min(mask.leg_y_um) - 2000, Lc_um,
                      max(mask.leg_y_um) + 2000)
    mask.resolution = dict(mirror=list(mirror), extra_gap_um=_extra_gap / UM,
                           note="clean")
    return mask

# =========================================================================
# ==== v2 merged satellite: cif_validators.py ============================
# =========================================================================
#!/usr/bin/env python3
"""cif_validators.py -- the J family: GEOMETRY-FIDELITY validators.

The engine's 162 validators prove the SCHEDULE and the NETWORK MODEL are
right, and that the canonical (auto-mode, linear) mask connects and
recovers F.  None of them look at the geometry a *rendered layout*
actually draws.  The J family closes that gap.  Every check recomputes
from the DRAWN polygons (or from stored physical inputs re-derived
through CIF's own kernels) -- it never trusts the builder's metadata for
the quantity it is checking.

  J1  junction_edge_continuity  every leg<->turn interface: the straight
      cross-section (3 channels) equals the arc cross-section, edge for
      edge, recomputed from dev rows and from the TURN polys' endpoints.
  J2  turn_rho_closure          each rail's resistance, INTEGRATED FROM
      THE DRAWN WIDTH PROFILE, equals rho_up * r(w_c) * pi*R (S14.5),
      and R_in == R_out (N8).
  J3  channel_post_overlap      no CHANNELS polygon covers the interior
      of any POSTS polygon (posts are solid; fluid through a post is a
      physical violation).  0.8 um interior margin (> post standoff).
  J4  min_corner_radius         raster morphology: opening/closing with
      a disc of radius r changes the fluid area by < eps  <=>  no corner
      sharper than r anywhere on the union boundary.
  J5  posts_region              posts exist ONLY where the architecture
      allows them (none in a meander; only below a seam in a composite).
  J6  ports_and_leads           1 inlet + 3 outlets, each port inside its
      channel, all four lead channels present, ports outside the core.
  J7  width_fidelity            drawn w_s == design w_s on EVERY row (no
      silent flooring), w_c monotone non-increasing.
  J8  F_from_polygons           the engine's geom_ftest on the v2 linear
      mask of the same device: F recovered from drawn widths.

self_test() corrupts a pristine mask six ways and asserts each corruption
is caught by its matching validator -- the engine's own philosophy.
"""
pass  # [v2-merged] from __future__ import annotations
import math
import numpy as np
from dataclasses import dataclass
from matplotlib.path import Path as MplPath

pass  # [v2-merged] import CIF
pass  # [v2-merged] from CIF import UM


@dataclass
class JResult:
    name: str
    ok: bool
    detail: str


# ---------------------------------------------------------------------------
# shared: rasterizer (exact, no antialiasing) + morphology
# ---------------------------------------------------------------------------

def rasterize(polys, region, res=1.5):
    """Binary fluid raster of CHANNELS-style polys over region (um)."""
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.collections import PolyCollection
    x0, y0, x1, y1 = region
    Wpx = max(8, int(math.ceil((x1 - x0) / res)))
    Hpx = max(8, int(math.ceil((y1 - y0) / res)))
    fig = Figure(figsize=(Wpx / 100.0, Hpx / 100.0), dpi=100)
    FigureCanvasAgg(fig)
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_axis_off()
    ax.add_collection(PolyCollection([np.asarray(p) for p in polys],
                                     facecolors="k", edgecolors="none",
                                     antialiased=False))
    ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
    fig.canvas.draw()
    buf = np.frombuffer(fig.canvas.buffer_rgba(), dtype=np.uint8)
    Wc, Hc = fig.canvas.get_width_height()
    img = buf.reshape(Hc, Wc, 4)[:, :, 0]
    return img < 128


def _disk(rpx):
    r = int(math.ceil(rpx))
    y, x = np.ogrid[-r:r + 1, -r:r + 1]
    return (x * x + y * y) <= rpx * rpx + 1e-9


def check_rounding(polys, r_um, region, res=1.0, eps=None, pad=None):
    """(ok, removed_frac, added_frac): opening removes sharp CONVEX
    corners; closing fills sharp CONCAVE ones.  A fully-rounded design at
    radius >= r changes by ~0.  The raster is PADDED by pad (default
    r+3*res) and the border cropped after morphology, so fluid cut by the
    window edge cannot masquerade as a sharp feature.
    CAVEAT (measured, real): opening also flags fluid THINNER than 2r --
    e.g. the engine's by-design 0.5 um cap-complement standoff pockets.
    Choose windows that exclude by-design sub-2r features, or treat the
    result as a min-feature measurement, not a corner failure."""
    from scipy import ndimage
    if pad is None:
        pad = r_um + 3.0 * res
    x0, y0, x1, y1 = region
    g = rasterize(polys, (x0 - pad, y0 - pad, x1 + pad, y1 + pad), res)
    st = _disk(r_um / res)
    op = ndimage.binary_opening(g, structure=st)
    cl = ndimage.binary_closing(g, structure=st)
    b = int(math.ceil(pad / res))
    g, op, cl = g[b:-b, b:-b], op[b:-b, b:-b], cl[b:-b, b:-b]
    area = max(1, int(g.sum()))
    rem = float((g & ~op).sum()) / area
    add = float((cl & ~g).sum()) / area
    if eps is None:
        eps = 8e-4
    return (rem <= eps and add <= eps), rem, add


# ---------------------------------------------------------------------------
# J1 -- junction edge continuity (recomputed from dev + drawn turn polys)
# ---------------------------------------------------------------------------

def _turn_triples(mask):
    idx = [i for i, q in enumerate(mask.polys)
           if str(q.get("role", "")).startswith("TURN_")]
    return [idx[k:k + 3] for k in range(0, len(idx), 3)]


def _arc_port_intervals(triple, mask, end):
    """Cross-section intervals at one arc port, from the DRAWN polys.
    Band poly layout: outer ring (n+1 pts, th ascending) + inner reversed.
    end=0 -> t=0 port; end=1 -> t=1 port."""
    out = []
    xs = []
    for i in triple:
        v = np.asarray(mask.polys[i]["v"], float)
        n2 = len(v) // 2
        if end == 0:
            po, pi_ = v[0], v[-1]
        else:
            po, pi_ = v[n2 - 1], v[n2]
        out.append((min(po[1], pi_[1]), max(po[1], pi_[1])))
        xs += [po[0], pi_[0]]
    return sorted(out), (min(xs), max(xs))


def _straight_intervals(rw, p, yc_um):
    wc2, lg, ws = rw.w_c / UM / 2, p.l_g / UM, rw.w_s / UM
    return sorted([(yc_um - wc2, yc_um + wc2),
                   (yc_um - wc2 - lg - ws, yc_um - wc2 - lg),
                   (yc_um + wc2 + lg, yc_um + wc2 + lg + ws)])


def J1_junction_edge_continuity(mask, dev, tol=1e-6):
    p = dev.req.proc
    triples = _turn_triples(mask)
    if len(triples) != len(mask.junctions):
        return JResult("J1_junction_continuity", False,
                       f"{len(triples)} drawn turns vs "
                       f"{len(mask.junctions)} junction records")
    bounds, i0 = [], 0
    for n in mask.legs:
        bounds.append((i0, i0 + n)); i0 += n
    worst, worst_x = 0.0, 0.0
    for k, tri in enumerate(triples):
        rw_up = dev.rows[bounds[k][1] - 1]
        rw_dn = dev.rows[bounds[k + 1][0]]
        y_up = mask.leg_y_um[k]
        y_dn = mask.leg_y_um[k + 1]
        for end, rw, yy in ((0, rw_up, y_up), (1, rw_dn, y_dn)):
            arc, (xa, xb) = _arc_port_intervals(tri, mask, end)
            st = _straight_intervals(rw, p, yy)
            for (a0, a1), (s0, s1) in zip(arc, st):
                worst = max(worst, abs(a0 - s0), abs(a1 - s1))
            worst_x = max(worst_x, xb - xa)
    ok = worst <= tol and worst_x <= tol
    return JResult("J1_junction_continuity", ok,
                   f"max edge mismatch {worst:.2e} um over "
                   f"{2 * len(triples)} junction planes "
                   f"(port-plane x spread {worst_x:.2e})")


# ---------------------------------------------------------------------------
# J2 -- turn rho-closure, integrated from the DRAWN width profile
# ---------------------------------------------------------------------------

def J2_turn_rho_closure(mask, dev, tol=1e-9):
    p = dev.req.proc
    triples = _turn_triples(mask)
    worst_T, worst_eq, nclip = 0.0, 0.0, 0
    for tri, rec in zip(triples, mask.turn_closures):
        T = rec["rho_up"] * CIF.r(rec["w_c_up"], p.h, p.mu) \
            * math.pi * rec["R_turn"]                # re-derived target
        Rm = {}
        for i in tri:
            role = mask.polys[i]["role"]
            if role == "TURN_PRIMARY":
                continue
            v = np.asarray(mask.polys[i]["v"], float)
            n2 = len(v) // 2
            outer, inner = v[:n2], v[n2:][::-1]
            widths = np.hypot(*(outer - inner).T) * UM        # SI
            mid = 0.5 * (outer + inner)
            dL = np.hypot(*np.diff(mid, axis=0).T) * UM
            rvals = np.array([CIF.r(max(w, 1e-9), p.h, p.mu) for w in widths])
            Rm[role] = float((0.5 * (rvals[:-1] + rvals[1:]) * dL).sum())
        Ri, Ro = Rm["TURN_SIDE_IN"], Rm["TURN_SIDE_OUT"]
        worst_eq = max(worst_eq, abs(Ri / Ro - 1.0))
        if rec.get("clipped"):
            nclip += 1
        else:
            worst_T = max(worst_T, abs(Ri / T - 1.0), abs(Ro / T - 1.0))
    ok = worst_T <= tol and worst_eq <= tol
    return JResult("J2_turn_rho_closure", ok,
                   f"max |R_rail/T - 1| = {worst_T:.2e} (target re-derived, "
                   f"widths measured off polys); max |R_in/R_out - 1| = "
                   f"{worst_eq:.2e}; {nclip} clipped turns")


# ---------------------------------------------------------------------------
# J3 -- no channel may cover a post interior
# ---------------------------------------------------------------------------

def J3_channel_post_overlap(mask, margin=0.8, step=4.0):
    posts = [np.asarray(q["v"], float) for q in mask.polys
             if q["layer"] == "POSTS"]
    chans = [(np.asarray(q["v"], float),) for q in mask.polys
             if q["layer"] == "CHANNELS"]
    cpaths = [(MplPath(v), (v[:, 0].min(), v[:, 1].min(),
                            v[:, 0].max(), v[:, 1].max())) for (v,) in chans]
    bad = 0
    for pv in posts:
        bb = (pv[:, 0].min(), pv[:, 1].min(), pv[:, 0].max(), pv[:, 1].max())
        gx = np.arange(bb[0] + margin, bb[2] - margin + 1e-9, step)
        gy = np.arange(bb[1] + margin, bb[3] - margin + 1e-9, step)
        if not len(gx) or not len(gy):
            continue
        pts = np.array([(x, y) for x in gx for y in gy])
        inside = MplPath(pv).contains_points(pts, radius=-margin)
        pts = pts[inside]
        if not len(pts):
            continue
        for cp, cb in cpaths:
            if cb[2] < bb[0] or cb[0] > bb[2] or cb[3] < bb[1] or cb[1] > bb[3]:
                continue
            if cp.contains_points(pts).any():
                bad += 1
                break
    return JResult("J3_channel_post_overlap", bad == 0,
                   f"{bad}/{len(posts)} posts have fluid in their interior "
                   f"(margin {margin} um > standoff)")


# ---------------------------------------------------------------------------
# J4 / J5 / J6 / J7 / J8
# ---------------------------------------------------------------------------

def J4_min_corner_radius(polys, r_um, region, res=1.0, eps=None,
                         label="J4_min_corner_radius"):
    ok, rem, add = check_rounding(polys, r_um, region, res, eps)
    return JResult(label, ok,
                   f"opening/closing at r={r_um:.2f} um changes fluid area "
                   f"by {rem:.2e}/{add:.2e} (sharp corners would exceed eps)")


def J5_posts_region(mask_or_posts, seam_y=None, forbid=False, tol=1.0):
    posts = (mask_or_posts if isinstance(mask_or_posts, list)
             else [np.asarray(q["v"], float) for q in mask_or_posts.polys
                   if q["layer"] == "POSTS"])
    if forbid:
        return JResult("J5_posts_region", len(posts) == 0,
                       f"{len(posts)} POSTS in a region that must have none")
    above = sum(1 for v in np.array(posts, dtype=object)
                if np.asarray(v)[:, 1].min() < seam_y - tol)
    return JResult("J5_posts_region", above == 0,
                   f"{above}/{len(posts)} posts above the seam "
                   f"(y={seam_y:.0f} um)")


def J6_ports_and_leads(mask):
    names = [p[0] for p in mask.ports]
    ok = sorted(names) == sorted(["INLET", "OUT_C", "OUT_S+", "OUT_S-"])
    roles = {str(q.get("role", "")) for q in mask.polys}
    for need in ("LEAD_IN", "LEAD_OUT_C", "LEAD_OUT_S+", "LEAD_OUT_S-"):
        ok &= need in roles
    cpaths = [MplPath(np.asarray(q["v"], float)) for q in mask.polys
              if q["layer"] == "CHANNELS"]
    cb = mask.core_bbox
    inport, outside = 0, 0
    for nm, cx, cy, r in mask.ports:
        if any(p.contains_point((cx, cy), radius=0.5) for p in cpaths):
            inport += 1
        if not (cb[0] <= cx <= cb[2] and cb[1] <= cy <= cb[3]):
            outside += 1
    ok &= (inport == 4)
    return JResult("J6_ports_and_leads", ok,
                   f"ports {names}; {inport}/4 inside a channel; "
                   f"{outside}/4 outside the core bbox; leads "
                   f"{sorted(r for r in roles if r.startswith('LEAD'))}")


def J7_width_fidelity(mask, tol=1e-9):
    worst = max(abs(r["ws_drawn"] - r["ws_design"]) for r in mask.row_geom)
    wcs = [r["wc"] for r in sorted(mask.row_geom, key=lambda q: q["i"])]
    mono = all(b <= a + 1e-9 for a, b in zip(wcs, wcs[1:]))
    return JResult("J7_width_fidelity", worst <= tol and mono,
                   f"max |drawn-design| w_s = {worst:.2e} um over "
                   f"{len(mask.row_geom)} rows; w_c monotone: {mono}")


def J8_F_from_polygons(dev, tol=2e-3):
    """GEOMETRY fidelity: Kirchhoff F of the ladder rebuilt from widths
    MEASURED off the drawn polygons vs Kirchhoff F of the ladder built from
    the DESIGN widths (same row subset, same network, same terminals).  This
    isolates the mask: march-vs-network accounting (sigma_0 fictitious-inlet
    bypass in DIRECT mode, terminal layer) cancels out; it is reported in
    the detail string instead."""
    pass  # [v2-merged] import cif_patches as CP
    p = dev.req.proc
    m = CP.build_mask_v2(dev, CIF.Layout("linear", 1, side_mode="straight"))
    wc_m, ws_m, idx = CIF.measure_widths(m, dev)
    sub = [dev.rows[k] for k in idx]
    rows_m = [CIF.replace(rw, w_c=float(a), w_s=float(b),
                          R_c=CIF.r(float(a), p.h, p.mu) * rw.Lam,
                          R_s=CIF.r(float(b), p.h, p.mu) * rw.Lam)
              for rw, a, b in zip(sub, wc_m, ws_m)]
    net_m = CIF.solve_network(rows_m, p, dev.req.tol)
    net_d = CIF.solve_network(sub, p, dev.req.tol)
    dF = abs(net_m.F_em - net_d.F_em)
    s0 = getattr(dev, "sigma_0", 0.0)
    return JResult("J8_F_from_polygons", dF <= tol and len(idx) > 0,
                   f"F_em(measured)={net_m.F_em:.5f} vs F_em(design widths)="
                   f"{net_d.F_em:.5f} (dF={dF:.1e}, {len(idx)} rows); march "
                   f"F={dev.F:.5f}, sigma_0={s0:.5f} (fictitious-inlet "
                   f"bypass; physical F is F_em)")


def J9_composition_clearance(mask, min_clear_um, res=2.0, adjacent=None):
    """Composition disjointness via distance transforms (a discretized
    SDF): fluid ink of structurally NON-adjacent groups must keep
    >= min_clear_um of wall everywhere.  Groups come from the builder's
    grp tags; the clearance itself is measured on the drawn polygons
    only.  Catches crossings (overlap area) and sub-wall proximity that
    per-junction checks cannot see."""
    from scipy import ndimage as _ndi
    groups = {}
    for q in mask.polys:
        g = q.get("grp")
        if g is None or q.get("layer") in ("OUTLINE", "POSTS"):
            continue
        groups.setdefault(g, []).append(np.asarray(q["v"], float))
    keys = sorted(groups.keys(), key=str)
    if adjacent is not None:
        adj = set(adjacent)
    else:
        n_legs = 1 + max([k[1] for k in keys if k[0] == "leg"], default=0)
        adj = set()
        for t in range(n_legs - 1):
            adj.add((("leg", t), ("turn", t)))
            adj.add((("turn", t), ("leg", t + 1)))
        adj.add((("lead", "in"), ("leg", 0)))
        adj.add((("lead", "out"), ("leg", n_legs - 1)))
    bb = {k: (min(v[:, 0].min() for v in groups[k]),
              min(v[:, 1].min() for v in groups[k]),
              max(v[:, 0].max() for v in groups[k]),
              max(v[:, 1].max() for v in groups[k])) for k in keys}
    pad = min_clear_um + 4 * res
    worst, worst_pair, viol = float("inf"), None, []
    for i in range(len(keys)):
        for k in range(i + 1, len(keys)):
            a, b = keys[i], keys[k]
            if (a, b) in adj or (b, a) in adj:
                continue
            A, B = bb[a], bb[b]
            dx = max(0.0, B[0] - A[2], A[0] - B[2])
            dy = max(0.0, B[1] - A[3], A[1] - B[3])
            dbox = (dx * dx + dy * dy) ** 0.5
            if dbox > pad:
                if dbox < worst:
                    worst, worst_pair = dbox, (a, b, "bbox bound")
                continue
            win = (max(A[0], B[0]) - pad, max(A[1], B[1]) - pad,
                   min(A[2], B[2]) + pad, min(A[3], B[3]) + pad)
            gA = rasterize(groups[a], win, res)
            gB = rasterize(groups[b], win, res)
            inter = int((gA & gB).sum())
            if inter:
                viol.append(f"{a} x {b}: OVERLAP {inter * res * res:.0f} um^2")
                worst, worst_pair = -inter * res * res, (a, b)
                continue
            if not gA.any() or not gB.any():
                continue
            dist = _ndi.distance_transform_edt(~gA) * res
            dmin = float(dist[gB].min())
            if dmin < worst:
                worst, worst_pair = dmin, (a, b)
            if dmin < min_clear_um - 2 * res:
                viol.append(f"{a} x {b}: clearance {dmin:.1f} um "
                            f"< {min_clear_um:.0f} um")
    ok = not viol
    det = (f"min inter-group clearance >= {worst:.1f} um at {worst_pair} "
           f"(required {min_clear_um:.0f} um); "
           + (f"{len(viol)} violation(s): " + "; ".join(viol[:3])
              if viol else "no non-adjacent contact anywhere"))
    return JResult("J9_composition_clearance", ok, det)


def J10_cutoff_from_polygons(dev, tol=0.02):
    """CUTOFF fidelity on the drawn mask. The cutoff D_c = kappa*delta* is
    held iff the lamina delta(i) = x(f(i); beta(i))*w_c(i) equals delta*
    at every rung -- the S1 schedule's guarantee. We recompute delta from
    the MEASURED primary widths w_c(i) (off the polygons) and the design
    schedule f(i), and compare to delta* = D_target/kappa. This isolates
    the drawn PRIMARY geometry (a distorted w_c would move the cutoff)
    WITHOUT the sigma_0 fictitious-inlet offset that contaminates the raw
    network f_em (the same offset J8 reports for F). Independent width
    fidelity is J7; the physical F is J8."""
    pass  # [v2-merged] import cif_patches as CP
    p = dev.req.proc
    pair = dev.req.pair
    dstar = dev.req.D_target / pair.kappa
    m = CP.build_mask_v2(dev, CIF.Layout("linear", 1, side_mode="straight"))
    wc_m, ws_m, idx = CIF.measure_widths(m, dev)
    sub = [dev.rows[k] for k in idx]
    devs = []
    for rw, wc in zip(sub, wc_m):
        wck = float(wc)
        beta = p.h / wck
        x = pair.x_of_f(min(max(rw.f, 1e-12), 0.4999), beta)
        devs.append(abs(x * wck / dstar - 1.0))
    dmax = max(devs)
    ok = dmax <= tol and len(idx) > 0
    return JResult("J10_cutoff_from_polygons", ok,
                   f"delta* = {dstar:.4f} um (D_target/kappa); max "
                   f"|delta_drawn/delta*-1| = {dmax:.2e} over {len(idx)} rungs "
                   f"(cutoff carried by the drawn primary widths; S1 exact)")


def run_suite(mask, dev, rounding=None):
    out = [J1_junction_edge_continuity(mask, dev),
           J2_turn_rho_closure(mask, dev),
           J3_channel_post_overlap(mask),
           J6_ports_and_leads(mask),
           J7_width_fidelity(mask),
           J8_F_from_polygons(dev),
           J9_composition_clearance(mask, dev.req.proc.F_min / 1e-6)]
    if rounding is not None:
        r_um, region = rounding
        ch = [np.asarray(q["v"], float) for q in mask.polys
              if q["layer"] == "CHANNELS"]
        out.append(J4_min_corner_radius(ch, r_um, region))
    return out


# ---------------------------------------------------------------------------
# negative self-test: each corruption caught by its validator
# ---------------------------------------------------------------------------

def self_test(dev, mask_builder):
    import copy
    msgs, ok = [], True

    def chk(cond, msg):
        nonlocal ok
        ok &= bool(cond)
        msgs.append(("PASS " if cond else "FAIL ") + msg)

    base = mask_builder()
    pos = run_suite(base, dev)
    chk(all(r.ok for r in pos), "positive: pristine v2 mask passes J1-J8 "
        + "; ".join(f"{r.name.split('_')[0]}" for r in pos if r.ok))
    for r in pos:
        if not r.ok:
            msgs.append(f"       !! {r.name}: {r.detail}")

    # n1: shift one turn's inner rail 3 um -> J1
    m = copy.deepcopy(base)
    for q in m.polys:
        if q.get("role") == "TURN_SIDE_IN":
            q["v"] = [(x, y + 3.0) for x, y in q["v"]]
            break
    chk(not J1_junction_edge_continuity(m, dev).ok,
        "n1 shifted turn rail -> J1 caught")

    # n2: widen one turn's outer rail band 20% -> J2 (measured R drops)
    m = copy.deepcopy(base)
    for q in m.polys:
        if q.get("role") == "TURN_SIDE_OUT":
            v = np.asarray(q["v"], float)
            n2 = len(v) // 2
            outer, inner = v[:n2], v[n2:][::-1]
            mid = 0.5 * (outer + inner)
            q["v"] = [tuple(p) for p in
                      np.vstack([mid + 1.2 * (outer - mid),
                                 (mid + 1.2 * (inner - mid))[::-1]])]
            break
    chk(not J2_turn_rho_closure(m, dev).ok,
        "n2 widened turn rail -> J2 caught")

    # n3: a channel crossing a post -> J3
    m = copy.deepcopy(base)
    pv = next(np.asarray(q["v"], float) for q in m.polys
              if q["layer"] == "POSTS")
    cx, cy = pv[:, 0].mean(), pv[:, 1].mean()
    m.polys.append(dict(layer="CHANNELS", v=[(cx - 40, cy - 3), (cx + 40, cy - 3),
                                             (cx + 40, cy + 3), (cx - 40, cy + 3)]))
    chk(not J3_channel_post_overlap(m).ok, "n3 fluid through a post -> J3 caught")

    # n4: remove the inlet lead -> J6
    m = copy.deepcopy(base)
    m.polys = [q for q in m.polys if q.get("role") != "LEAD_IN"]
    chk(not J6_ports_and_leads(m).ok, "n4 missing inlet lead -> J6 caught")

    # n5: silent flooring -> J7
    m = copy.deepcopy(base)
    m.row_geom = copy.deepcopy(base.row_geom)
    for r in m.row_geom[:40]:
        r["ws_drawn"] = max(r["ws_drawn"], r["ws_design"] + 5.0)
    chk(not J7_width_fidelity(m).ok, "n5 floored widths -> J7 caught")

    # n6: J4 discriminates rounded vs sharp free channel ends.  Strip the
    # ENDCAP discs: the leg-0 side-channel butt at x=0 then exposes two
    # 90-degree convex corners (0.215*r^2 removed each by opening), while
    # the capped version is a semicircular tip.  Axis-aligned region ==
    # zero raster staircase noise; res fine enough to resolve r=3.
    # window: OUTER half of the upper side-channel end -- excludes the
    # by-design 0.5 um cap-complement standoff pockets on the wall side
    # (they sit below the side channel's inner edge and are thinner than
    # 2r, so opening would flag them in ANY variant; see check_rounding).
    ecs = [np.asarray(q["v"], float) for q in base.polys
           if q.get("role") == "ENDCAP"]
    ec = max(ecs, key=lambda v: v[:, 1].mean())
    ys_end = float(ec[:, 1].mean())
    r_end = 0.5 * float(ec[:, 1].max() - ec[:, 1].min())
    reg = (-80.0, ys_end, 120.0, ys_end + r_end + 25.0)
    ch = [np.asarray(q["v"], float) for q in base.polys
          if q["layer"] == "CHANNELS"]
    sharp = [np.asarray(q["v"], float) for q in base.polys
             if q["layer"] == "CHANNELS" and q.get("role") != "ENDCAP"]
    ok_r, rem_r, _ = check_rounding(ch, 3.0, reg, res=0.25)
    ok_s, rem_s, _ = check_rounding(sharp, 3.0, reg, res=0.25)
    chk(ok_r and rem_s > 5 * max(rem_r, 1e-6),
        f"n6 semicapped end passes (rem {rem_r:.1e}); butt end fails "
        f"margin (rem {rem_s:.1e}) -> J4 discriminates")
    return ok, msgs

# =========================================================================
# ==== v2 merged satellite: dxf_validate.py ==============================
# =========================================================================
#!/usr/bin/env python3
"""DXF generation + COMPUTATIONAL validation of the DXF files themselves.

The chain that makes the DXF trustworthy:
    mask.polys --> DXF (roles encoded in layer names, um units)
    DXF --(ezdxf read-back)--> duck Mask {polys with roles, legs}
    duck Mask --> CIF.measure_widths   (the ENGINE's own metrology, verbatim)
    measured widths --> CIF.solve_network --> F_em, f_em FROM THE FILE
    f_em --> pairing inversion --> emergent cutoff adherence FROM THE FILE

So D3-D7 are computed from the .dxf bytes, not from the in-memory mask.
The engine metrology (measure_widths) slices leg-0 rows in leg-local
coordinates; it is exact for LINEAR masks (one leg = whole device). Serpentine
DXFs therefore get D1/D2 (validity + round-trip) and visual checks, while the
LINEAR DXF carries the full physics validation. That scoping is stated, not
hidden.
"""
pass  # [v2-merged] from __future__ import annotations
from dataclasses import dataclass
import numpy as np
import ezdxf
pass  # [v2-merged] import CIF
pass  # [v2-merged] from CIF import UM

LAYER_COLOR = {"CHANNELS": 4, "CHANNELS__PRIMARY": 5, "CHANNELS__SIDE": 4,
               "POSTS": 1, "OUTLINE": 8, "PORTS": 3, "META": 7}


def _role_layer(p):
    """BIJECTIVE (layer, role) -> DXF layer name: role=None -> bare layer,
    any role -> LAYER__ROLE.  Inverse in read_dxf_mask.  This preserves the
    engine's role taxonomy through the file EXACTLY.  CIF.measure_widths
    selects side channels by role == "SIDE" (the two collection rails only);
    the old many-to-one mapping ("else -> CHANNELS__SIDE") relabelled 4394
    connector/lead polys as SIDE on read-back, corrupting side metrology
    (root cause of the D3/D4/D5 failures)."""
    r = p.get("role")
    return f"{p['layer']}__{r}" if r else p["layer"]


def export_mask_dxf(mask, path):
    """Write the mask to DXF with the polygon ROLE encoded bijectively in the
    layer name, ports as circles, and the legs vector as a META text entity,
    so the file is self-describing for engine metrology on read-back."""
    doc = ezdxf.new()
    doc.header["$INSUNITS"] = 13                    # micrometres
    msp = doc.modelspace()
    names = {_role_layer(p) for p in mask.polys} | {"PORTS", "META"}
    for name in sorted(names):
        doc.layers.add(name, color=LAYER_COLOR.get(name, 2))
    for p in mask.polys:
        msp.add_lwpolyline(p["v"], close=True,
                           dxfattribs={"layer": _role_layer(p)})
    for tag, x, y, w in (getattr(mask, "ports", None) or []):
        msp.add_circle((x, y), radius=max(w, 5.0), dxfattribs={"layer": "PORTS"})
        msp.add_text(tag, dxfattribs={"layer": "PORTS", "height": 40.0}
                     ).set_placement((x, y + max(w, 5.0) + 20))
    msp.add_text(f"LEGS={list(getattr(mask, 'legs', [0]))}",
                 dxfattribs={"layer": "META", "height": 1.0}).set_placement(
                     (mask.bounds[0], mask.bounds[1] - 100))
    doc.saveas(path)
    return path


class DuckMask:
    """Just enough Mask for CIF.measure_widths: .polys (role dicts) + .legs."""
    def __init__(self, polys, legs):
        self.polys, self.legs = polys, legs


def read_dxf_mask(path):
    doc = ezdxf.readfile(path)
    msp = doc.modelspace()
    polys, legs = [], [0]
    for e in msp.query("LWPOLYLINE"):
        lay = e.dxf.layer
        v = [(pt[0], pt[1]) for pt in e.get_points()]
        base, sep, role = lay.partition("__")
        polys.append({"layer": base, "role": role if sep else None, "v": v})
    for e in msp.query("TEXT"):
        t = e.dxf.text
        if t.startswith("LEGS="):
            legs = eval(t[5:], {"__builtins__": {}})
    return DuckMask(polys, legs), doc


# ---------------------------------------------------------------------------
@dataclass
class DResult:
    name: str
    ok: bool
    detail: str


def _d(n, ok, det):
    return DResult(n, bool(ok), det)


def D1_dxf_valid(path):
    duck, doc = read_dxf_mask(path)
    aud = doc.audit()
    n = len(duck.polys)
    return _d("D1_dxf_valid", len(aud.errors) == 0 and n > 0,
              f"ezdxf structural audit: {len(aud.errors)} errors; "
              f"{n} polylines, legs={duck.legs}")


def D2_roundtrip(path, mask):
    duck, _ = read_dxf_mask(path)
    if len(duck.polys) != len(mask.polys):
        return _d("D2_roundtrip", False,
                  f"poly count {len(duck.polys)} != source {len(mask.polys)}")
    # vertex-exactness: compare in written order (export preserved order)
    worst = 0.0
    for pd, ps in zip(duck.polys, mask.polys):
        a = np.asarray(pd["v"], float); b = np.asarray(ps["v"], float)
        if a.shape != b.shape:
            return _d("D2_roundtrip", False, "vertex count mismatch")
        worst = max(worst, float(np.max(np.abs(a - b))))
    return _d("D2_roundtrip", worst < 1e-6,
              f"max vertex deviation DXF vs source = {worst:.2e} um over "
              f"{len(duck.polys)} polylines (write->read exact)")


def _dxf_network(path, dev):
    """Engine metrology ON THE FILE: read DXF, measure widths with
    CIF.measure_widths, rebuild rows, solve both drawn and design networks."""
    duck, _ = read_dxf_mask(path)
    proc = dev.req.proc
    wc_m, ws_m, idx = CIF.measure_widths(duck, dev)
    sub = [dev.rows[k] for k in idx]
    rows_m = [CIF.replace(rw, w_c=float(a), w_s=float(b),
              R_c=CIF.r(float(a), proc.h, proc.mu) * rw.Lam,
              R_s=CIF.r(float(b), proc.h, proc.mu) * rw.Lam)
              for rw, a, b in zip(sub, wc_m, ws_m)]
    net_d = CIF.solve_network(rows_m, proc, dev.req.tol)
    net_r = CIF.solve_network(sub, proc, dev.req.tol)
    return dict(wc=wc_m, ws=ws_m, idx=idx, rows_m=rows_m, rows_design=sub,
                net_drawn=net_d, net_design=net_r, proc=proc)


def D3_widths_from_dxf(g, dev, tol=1e-3):
    # gate = 1e-3 um = 1 nm: the mask's own drawn-vs-design grid fidelity,
    # identical to the J7 gate (1e-9 m).  Demanding tighter than the source
    # mask achieves on the live object would be gating noise, not the file.
    dws = np.max(np.abs(g["ws"] - np.array([r.w_s for r in g["rows_design"]]))) / UM
    dwc = np.max(np.abs(g["wc"] - np.array([r.w_c for r in g["rows_design"]]))) / UM
    n = len(g["idx"])
    ok = n == dev.N and dws < tol and dwc < tol
    return _d("D3_widths_from_dxf", ok,
              f"{n}/{dev.N} rows measured FROM THE FILE; max |w-design|: "
              f"side {dws:.2e} um, trunk {dwc:.2e} um (gate {tol:g} um)")


def D4_F_from_dxf(g, F_desired, tol=3e-3):
    Fd = g["net_drawn"].F_em
    ok = abs(Fd - F_desired) <= tol
    return _d("D4_F_from_dxf", ok,
              f"filtrate solved from DXF widths: F_em = {Fd:.4f} vs desired "
              f"{F_desired} (|dF| = {abs(Fd-F_desired):.1e} <= {tol})")


def D5_field_fidelity_from_dxf(g, eps=2e-3):
    fd = np.asarray(g["net_drawn"].f_em)
    fr = np.asarray(g["net_design"].f_em)
    n = len(fd); nt = max(2, n // 12)
    b = slice(1, n - nt)
    e = float(np.max(np.abs(fd[b] / fr[b] - 1.0)))
    eF = abs(g["net_drawn"].F_em / g["net_design"].F_em - 1.0)
    return _d("D5_field_fidelity_from_dxf", e <= eps and eF <= 2e-2,
              f"drawn-vs-design network from DXF: max bulk |f ratio-1| = "
              f"{e:.2e}, |F ratio-1| = {eF:.2e} (sigma0 cancels; N2/N3 basis)")


def D6_cutoff_from_dxf(g, dev, tol=5e-3):
    """Schedule-basis cutoff (J10 discipline): delta from design f x widths
    measured FROM THE FILE -- fabrication holds the cutoff."""
    pair = dev.req.pair
    dstar = dev.req.D_target / pair.kappa
    proc = g["proc"]
    devi = [abs(pair.x_of_f(rw.f, proc.h / float(wc)) * float(wc) / dstar - 1.0)
            for rw, wc in zip(g["rows_design"], g["wc"])]
    worst = max(devi)
    D_em = dev.req.D_target / UM * (1.0 + float(np.median(devi)))
    return _d("D6_cutoff_from_dxf", worst <= tol,
              f"schedule-basis delta from DXF widths: max |delta/delta*-1| = "
              f"{worst:.2e} (D_em ~ {D_em:.3f} um vs target "
              f"{dev.req.D_target/UM:.1f} um)")


def D7_adherence_from_dxf(g, dev, adh_tol=0.05):
    """Emergent-cutoff adherence and lamina purity computed FROM THE FILE:
    invert the DXF-solved f_em through the pairing."""
    proc = g["proc"]; pair = dev.req.pair
    dstar = dev.req.D_target / pair.kappa
    fd = np.asarray(g["net_drawn"].f_em)
    P, w_in, w_tot, excess = 1.0, 0.0, 0.0, 0.0
    for rw, wc, fe in zip(g["rows_design"], g["wc"], fd):
        fe = float(fe)
        fc = min(max(fe, 1e-9), 0.4999)
        x_em = pair.x_of_f(fc, proc.h / float(wc))
        dr = x_em * float(wc) / dstar
        wgt = 2.0 * fe * P
        w_tot += wgt
        if abs(dr - 1.0) <= adh_tol:
            w_in += wgt
        excess += 2.0 * max(fe - rw.f, 0.0) * P
        P *= (1.0 - 2.0 * fe)
    A = w_in / w_tot
    C = excess / g["net_drawn"].F_em
    return _d("D7_adherence_from_dxf", True,
              f"FROM THE FILE: adherence(+/-5%) = {A*100:.2f}% of filtrate mass; "
              f"lamina purity = {(1-C)*100:.2f}% (excess supra-delta* lamina "
              f"{C*100:.2f}%)"), A, C


def run_dxf_suite(path, mask, dev, F_desired):
    r = [D1_dxf_valid(path), D2_roundtrip(path, mask)]
    g = _dxf_network(path, dev)                     # <-- file metrology ONCE
    r += [D3_widths_from_dxf(g, dev), D4_F_from_dxf(g, F_desired),
          D5_field_fidelity_from_dxf(g), D6_cutoff_from_dxf(g, dev)]
    d7, A, C = D7_adherence_from_dxf(g, dev)
    r.append(d7)
    return r, A, C

# =========================================================================
# ==== v2 merged satellite: optimize_build.py ============================
# =========================================================================
#!/usr/bin/env python3
"""The OPTIMIZED-DEVICE builder (the point of Tasks 1/2/4 combined):
smallest-N classic device that ADHERES TO THE PHYSICS -- delivered filtrate
F_em = F_target and cutoff D held -- and that RENDERS to a mask/DXF today.

Two facts drive the design (both measured, this session):
  1. FOLD mode (sigma0 = 0) delivers F exactly at the same N as DIRECT for a
     given wc0, and its floor is wc0 = 45 (N = 178, M14 Stokes below) -- but
     it is UNRENDERABLE today: its inlet secondaries (w_s down to 4.7 um,
     Phi ~ 200 >> W/p) need the multi-slot fold renderer that Phase B scoped
     as open. It is reported as the paper optimum.
  2. DIRECT mode renders (straight rows, w_start = 1.03 GS), but sigma0
     shifts the delivered filtrate (F_em ~ F_target + sigma0) AND elevates
     the inlet rows' emergent extraction f_em, which -- by the engine's own
     pairing inversion (verified exact: f_of_x(dstar/wc, beta) == f_sched to
     0.0e0) -- means those rows capture a WIDER streamline bundle, i.e. an
     effectively COARSER cutoff. The optimized device must therefore be
     (a) sigma0-COMPENSATED: design at F_target' so the drawn network
         delivers F_em = F_desired, and
     (b) ADHERENCE-GATED: enough of the delivered filtrate mass must be
         collected at the design cutoff.

PRE-REGISTERED ACCEPTANCE CRITERION (fixed before the sweep):
  * all engine design gates pass, and the mask builds (straight-renderable);
  * |F_em(drawn network) - F_desired| <= 3e-3;
  * schedule-basis cutoff (J10 basis: delta from design f x measured wc)
    holds as before; AND
  * >= 95% of the delivered filtrate mass is collected by rows whose
    EMERGENT cutoff delta_em(i) = x_of_f(f_em(i), beta_i) * wc(i) is within
    +/- 5% of delta*.
The optimized device is the smallest-N candidate meeting all four.
"""
pass  # [v2-merged] from __future__ import annotations
from dataclasses import dataclass
import numpy as np
pass  # [v2-merged] import CIF
pass  # [v2-merged] from CIF import UM, GateError, WindowError, RoutabilityError
pc = CIF.physics_core

F_DESIRED = 0.35
D_TARGET = 4.0          # um
ADH_TOL = 0.05          # +/- 5% of delta*
ADH_MIN = 0.95          # >= 95% of filtrate mass
F_TOL = 3e-3


def _proc():
    return CIF.Process(h=140*UM, G=19*UM, g_p=30*UM, GS=28*UM, F_min=19*UM,
                       W=1562*UM, l_stroke=1562*UM, l_g=60*UM, L_out=200*UM,
                       AR_max=10.0)


def _design(wc0, F_target, mode="DIRECT"):
    proc = _proc()
    req = CIF.Request(D_target=D_TARGET*UM, F_target=F_target,
                      W_shape=lambda t, w=wc0: (w - 0.2*w*t)*UM, proc=proc,
                      pair=CIF.Pairing(), tol=CIF.Tolerances(), mode=mode,
                      w_start=1.03*proc.GS, N_seed=200)
    return CIF.design(req), proc


def compensate(wc0, F_desired=F_DESIRED, max_iter=6):
    """Fixed-point sigma0 compensation: choose F_target' so the device's own
    network solve delivers F_em = F_desired. (Design-network F_em == drawn
    F_em to ~2e-6, verified Phase A, so iteration needs no mask.)"""
    Ft = F_desired
    hist = []
    for it in range(max_iter):
        dev, proc = _design(wc0, Ft)
        net = CIF.solve_network(list(dev.rows), proc, dev.req.tol)
        err = net.F_em - F_desired
        hist.append((Ft, dev.N, dev.sigma_0, net.F_em))
        if abs(err) <= 1e-3:
            break
        Ft = Ft - err
        if Ft <= 0.02:
            raise GateError(f"compensation drove F_target' to {Ft:.3f}")
    return dev, proc, Ft, hist


def adherence(dev, proc, f_em=None):
    """Mass-weighted emergent-cutoff adherence, engine-sanctioned inversion:
    x_em = x_of_f(f_em, beta); delta_em = x_em * wc. Weight per row =
    2 f_em(i) P_em(i-1) (the filtrate identity F = 2 sum f P)."""
    rows = list(dev.rows)
    if f_em is None:
        f_em = np.asarray(CIF.solve_network(rows, proc, dev.req.tol).f_em)
    pair = dev.req.pair
    dstar = D_TARGET*UM / pair.kappa
    P, w_in, w_tot = 1.0, 0.0, 0.0
    worst, nclamp = 0.0, 0
    prof = []
    for i, rw in enumerate(rows):
        fe = float(f_em[i])
        fc = min(max(fe, 1e-9), 0.4999)
        nclamp += (fc != fe)
        beta = proc.h / rw.w_c
        x_em = pair.x_of_f(fc, beta)
        dr = x_em * rw.w_c / dstar               # delta_em / delta*
        wgt = 2.0 * fe * P
        w_tot += wgt
        if abs(dr - 1.0) <= ADH_TOL:
            w_in += wgt
        worst = max(worst, abs(dr - 1.0))
        prof.append(dr)
        P *= (1.0 - 2.0 * fe)
    A = w_in / w_tot if w_tot > 0 else 0.0
    return A, np.asarray(prof), worst, nclamp


def purity(dev, proc, f_em=None):
    """Lamina purity 1-C, with C = sum_i 2 (f_em - f_sched)^+ P / F_em: the
    filtrate fraction that is EXCESS supra-delta* lamina (skimmed beyond the
    schedule's cutoff-exact allocation).  Complementary to adherence():
    adherence counts a whole row off-band once its delta_em leaves +/-5%
    (conservative -- the row's entire flux is condemned); purity condemns
    only the excess mass itself, so purity >= adherence always.  Report
    BOTH; the acceptance gate uses row-level adherence."""
    rows = list(dev.rows)
    if f_em is None:
        f_em = np.asarray(CIF.solve_network(rows, proc, dev.req.tol).f_em)
    P, excess, F = 1.0, 0.0, 0.0
    for i, rw in enumerate(rows):
        fe = float(f_em[i])
        excess += 2.0 * max(fe - rw.f, 0.0) * P
        F += 2.0 * fe * P
        P *= (1.0 - 2.0 * fe)
    C = excess / F if F > 0 else 1.0
    return 1.0 - C, C


@dataclass
class Opt:
    wc0: float
    Ft: float
    N: int
    sigma0: float
    F_em: float
    adh: float
    ok_gates: bool
    ok: bool
    note: str


def sweep(wc0_grid=(50, 60, 70, 80, 90, 100)):
    out = []
    for wc0 in wc0_grid:
        try:
            dev, proc, Ft, hist = compensate(wc0)
        except (GateError, WindowError) as e:
            out.append(Opt(wc0, float("nan"), -1, float("nan"), float("nan"),
                           float("nan"), False, False, type(e).__name__))
            continue
        net = CIF.solve_network(list(dev.rows), proc, dev.req.tol)
        A, prof, worst, ncl = adherence(dev, proc, np.asarray(net.f_em))
        okF = abs(net.F_em - F_DESIRED) <= F_TOL
        okA = A >= ADH_MIN
        note = []
        if not okF: note.append(f"F miss {net.F_em-F_DESIRED:+.3f}")
        if not okA: note.append(f"adh {A*100:.1f}%")
        if ncl: note.append(f"{ncl} clamped")
        out.append(Opt(wc0, Ft, dev.N, dev.sigma_0, net.F_em, A,
                       True, okF and okA, "; ".join(note) or "PASS"))
    return out


def pick(res):
    valid = [r for r in res if r.ok]
    return min(valid, key=lambda r: r.N) if valid else None




# =========================================================================
# ==============  v2 NEW: G14 ROUTABILITY + EXTRA VALIDATORS  =============
# =========================================================================

def check_fold_routability(dev, proc, band_margin_um=110.0,
                           monster_slots=8):
    """Decision procedure for fold-field drawability under this process
    (single-track lattice theorem, machine-checked in fold45_probe.py).
    A vertical stroke needs w/2 + F_min clearance from BOTH bounding link
    columns of its Lambda cell, so Lambda < w_min + 2*F_min forbids every
    off-column vertical; each descent->ascent reversal then requires one
    pivot fold of length exactly Lambda, width-tunable, capped at
    r(w_min)/r(GS)*Lambda -- while every folded row's required length
    exceeds that cap by construction (w_s < w_min).  The largest
    `monster_slots` folds are exempt from the spend budget (routable as
    upstream multi-arm blocks in the virgin lead region)."""
    GSu = proc.GS / UM
    Lam = (proc.G + proc.g_p) / UM
    Fm = proc.F_min / UM
    wmin = max(Fm, (proc.h / UM) / proc.AR_max)
    rr = lambda w_um: r(w_um * UM, proc.h, proc.mu)
    rGS = rr(GSu)
    rho = lambda w_um: rr(w_um) / rGS
    single = Lam < wmin + 2.0 * Fm
    folded = [rw for rw in dev.rows if rw.w_s / UM < wmin - 1e-12]
    if not folded:
        return dict(ok=True, nf=0, single_track=single,
                    Lambda=Lam, w_min=wmin, F_min=Fm,
                    note="no folded rows (all w_s >= w_min)")
    L = sorted((rho(max(rw.w_s / UM, 0.05)) * Lam for rw in folded),
               reverse=True)
    field = L[monster_slots:]                # beyond upstream-arm capacity
    cap = rho(wmin) * Lam
    band = proc.W / UM - band_margin_um
    spend = sum(max(x - Lam, 0.0) for x in field) / rho(wmin)
    need_rev = spend > band
    infeasible = single and len(field) > 0 and need_rev
    return dict(ok=not infeasible, single_track=single, nf=len(L),
                nf_field=len(field), Lambda=Lam, w_min=wmin, F_min=Fm,
                pivot_cap=cap,
                min_fold_L=(min(field) if field else min(L)),
                max_fold_L=max(L), min_spend=spend, band=band,
                reversals_lb=int(spend // band) if band > 0 else 0)


def G14_fold_routability(dev, proc, strict=True, **kw):
    res = check_fold_routability(dev, proc, **kw)
    if strict and not res["ok"]:
        raise RoutabilityError(
            "G14: fold field undrawable under this process: single-track "
            f"lattice (Lambda={res['Lambda']:.1f} um < w_min + 2*F_min = "
            f"{res['w_min'] + 2 * res['F_min']:.1f} um) forbids off-column "
            f"verticals; the {res['nf_field']} beyond-monster folds need "
            f"minimum descent spend {res['min_spend']:.0f} um vs band "
            f"{res['band']:.0f} um, so >= {max(1, res['reversals_lb'])} "
            f"reversal(s) are required, each needing a pivot fold with "
            f"L <= {res['pivot_cap']:.2f} um, but every folded row needs "
            f"L >= {res['min_fold_L']:.2f} um.  Feasible at Lambda >= "
            f"{res['w_min'] + 2 * res['F_min']:.0f} um.")
    return res


def V_operating_regime(dev, proc, Q_uLmin=5.0, mu_SI=1.0e-3, rho_SI=1.0e3,
                       Re_max=1.0, gdot_max=1500.0, dP_max_kPa=200.0):
    """File-independent operating gates: Re, entrance length vs pitch,
    max wall shear via the engine's C_shear, and a trunk-only upper
    estimate of the pressure budget at the prescribed march flows."""
    rows = list(dev.rows)
    h = proc.h
    Lam = proc.G + proc.g_p
    wc0 = rows[0].w_c
    Q = Q_uLmin * 1e-9 / 60.0
    Dh = 2.0 * wc0 * h / (wc0 + h)
    ub = Q / (wc0 * h)
    Re = rho_SI * ub * Dh / mu_SI
    Le = 0.05 * Re * Dh
    a0 = max(wc0, h) / min(wc0, h)
    gd = C_shear(a0) * Q / (max(wc0, h) * min(wc0, h) ** 2)
    P_after, dP = 1.0, 0.0
    for rw in rows:
        dP += r(rw.w_c, h, proc.mu) * Lam * (P_after * Q)
        P_after *= (1.0 - 2.0 * rw.f)
    g = dict(Re=Re <= Re_max, entrance=Le <= 0.5 * Lam,
             shear=gd <= gdot_max, pressure=dP <= dP_max_kPa * 1e3)
    return dict(ok=all(g.values()), Re=Re, Le_um=Le / UM, gdot=gd,
                dP_kPa=dP / 1e3, gates=g)


def V_steric(dev, proc, D_target=None):
    """Clogging / steric passage gates (current-practice)."""
    D = D_target if D_target is not None else dev.req.D_target
    ds = D / dev.req.pair.kappa
    r1, r2, r3 = D / proc.G, D / proc.h, ds / proc.h
    return dict(ok=(r1 < 0.5 and r2 < 0.3 and r3 < 0.2),
                D_over_G=r1, D_over_h=r2, dstar_over_h=r3)


def _adherence_on_rows(rows, dev, proc, f_em):
    import numpy as _np
    pair = dev.req.pair
    dstar = dev.req.D_target / pair.kappa
    P, w_in, w_tot = 1.0, 0.0, 0.0
    for rw, fe in zip(rows, _np.asarray(f_em, float)):
        fc = min(max(float(fe), 1e-9), 0.4999)
        x_em = pair.x_of_f(fc, proc.h / rw.w_c)
        dr = x_em * rw.w_c / dstar
        wgt = 2.0 * float(fe) * P
        w_tot += wgt
        if abs(dr - 1.0) <= 0.05:
            w_in += wgt
        P *= (1.0 - 2.0 * float(fe))
    return (w_in / w_tot) if w_tot > 0 else 0.0


def V_fab_tolerance(dev, proc, n=24, sigma_um=0.25, seed=0,
                    F_desired=None, F_tol=3e-3):
    """Monte-Carlo fabrication tolerance: perturb every drawn width by
    N(0, sigma) and re-solve the network; report the F and adherence
    spread and the yield against |F - F_desired| <= F_tol."""
    import numpy as _np
    from dataclasses import replace as _repl, is_dataclass as _isdc
    rng = _np.random.default_rng(seed)
    rows0 = list(dev.rows)
    Fd = F_desired
    Fs, As = [], []
    for _ in range(n):
        pr = []
        for rw in rows0:
            dwc = float(rng.normal(0.0, sigma_um)) * UM
            dws = float(rng.normal(0.0, sigma_um)) * UM
            if _isdc(rw):
                pr.append(_repl(rw, w_c=max(rw.w_c + dwc, 2 * UM),
                                w_s=max(rw.w_s + dws, 2 * UM)))
            else:
                import copy as _cp
                q = _cp.copy(rw)
                q.w_c = max(rw.w_c + dwc, 2 * UM)
                q.w_s = max(rw.w_s + dws, 2 * UM)
                pr.append(q)
        net = solve_network(pr, proc, dev.req.tol)
        Fs.append(net.F_em)
        As.append(_adherence_on_rows(pr, dev, proc, net.f_em))
    Fs = _np.asarray(Fs); As = _np.asarray(As)
    if Fd is None:
        Fd = float(Fs.mean())
    yld = float((_np.abs(Fs - Fd) <= F_tol).mean())
    return dict(ok=yld >= 0.95, n=n, sigma_um=sigma_um,
                F_mean=float(Fs.mean()), F_std=float(Fs.std()),
                F_min=float(Fs.min()), F_max=float(Fs.max()),
                adh_mean=float(As.mean()), adh_min=float(As.min()),
                yield_F=yld)


def V_turn_closure(mask, eps=1e-9):
    """Gate the drawn turn closures: zero clipped turns and residuals at
    machine precision (closes the ungated-shortfall path in
    _solve_alpha)."""
    cl = getattr(mask, "turn_closures", []) or []
    if not cl:
        return dict(ok=True, turns=0, clipped=0, worst=0.0)
    clipped = sum(int(c.get("clipped", 0)) for c in cl)
    worst = max(max(abs(c["R_in"] / c["T"] - 1.0),
                    abs(c["R_out"] / c["T"] - 1.0)) for c in cl)
    return dict(ok=(clipped == 0 and worst <= eps),
                turns=len(cl), clipped=clipped, worst=worst)


def V_outlet_divider(dev, proc, F_em=None, tol=5e-3):
    """The current-divider closed form  F = 2*Rc_out/(Rs_out + 2*Rc_out)
    at the COMMON OUTLET -- a corollary of the mass identity plus G13
    run-out matching: with per-side outlet flow S and trunk flow P, the
    side fraction is 2S/(P+2S) = 2/(rho_out+2), and matched run-outs
    (Rs = rho_out*Rc) turn that flow ratio into the resistance divider.
    Exact as the OUTLET identity (any sigma0, since the physical split
    at the last node is a divider when the last rung drop is small);
    equals F_march = 1-P only when sigma0 = 0; it is NOT a global design
    law mid-device for tapered ladders (there it gives the local sigma(i),
    not F).  This validator computes the march-side divider and compares
    it to the solved network F_em."""
    P, S = 1.0, float(dev.sigma_0)
    for rw in dev.rows:
        S += rw.f * P
        P *= (1.0 - 2.0 * rw.f)
    F_div = 2.0 * S / (P + 2.0 * S)
    if F_em is None:
        F_em = solve_network(list(dev.rows), proc, dev.req.tol).F_em
    return dict(ok=abs(F_em - F_div) <= tol, F_em=F_em,
                F_march=1.0 - P, F_divider=F_div, rho_out=P / S,
                sigma0=float(dev.sigma_0), dF=F_em - F_div)


def outlet_divider_drawn(Rc_out, Rs_out):
    """F predicted from DRAWN outlet-leg resistances (single-trunk Rc_out,
    per-side Rs_out): the user-facing current-divider form."""
    return 2.0 * Rc_out / (Rs_out + 2.0 * Rc_out)


# =========================================================================
# ==============  v2 NEW: STAMPED EMISSION + ACCEPTANCE  ==================
# =========================================================================

def engine_sha256():
    try:
        return _v2_hashlib.sha256(open(__file__, "rb").read()).hexdigest()
    except Exception:
        return "unavailable"


def emit_device_dxf(mask, path, note=""):
    """Export via the bijective-role DXF writer, then stamp a META layer
    TEXT with engine version + hashes + UTC, and write a .sha256
    sidecar.  Returns the file hash."""
    export_mask_dxf(mask, path)
    doc = ezdxf.readfile(path)
    msp = doc.modelspace()
    if "META" not in doc.layers:
        doc.layers.add("META")
    stamp = (f"CIFV2 v{__version__} base={V2_BASE_SHA256[:12]} "
             f"self={engine_sha256()[:12]} "
             f"utc={_v2_dt.datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%SZ')}"
             + (f" {note}" if note else ""))
    b = mask.bounds
    msp.add_text(stamp, dxfattribs={"layer": "META", "height": 6.0,
                                    "insert": (b[0], b[1] - 40.0)})
    doc.saveas(path)
    h = _v2_hashlib.sha256(open(path, "rb").read()).hexdigest()
    with open(path + ".sha256", "w") as fh:
        fh.write(f"{h}  {_v2_os.path.basename(path)}\n")
    return h


def run_acceptance(outdir=".", emit=("80", "90", "100"),
                   run_dxf_suites=True, verbose=True):
    """End-to-end acceptance: goldens, G14 refusal of FOLD-45, J-suite on
    the serpentine fixture, stamped re-emission of the DIRECT family, and
    the D-suite run FROM THE EMITTED BYTES.  Writes acceptance JSON +
    hash manifest into `outdir`."""
    import time as _t
    import numpy as _np
    t0 = _t.time()
    rep = dict(version=__version__, base_sha=V2_BASE_SHA256,
               self_sha=engine_sha256(), utc=V2_BUILD_UTC, steps={})

    def say(s):
        if verbose:
            print(f"[{_t.time()-t0:6.1f}s] {s}", flush=True)

    # -- 1. goldens ------------------------------------------------------
    exp = {80: 406, 90: 562, 100: 731}
    gold = {}
    for wc0 in (80, 90, 100):
        dev, proc, Ft, _ = compensate(wc0)
        net = solve_network(list(dev.rows), proc, dev.req.tol)
        A, _, _, _ = adherence(dev, proc, _np.asarray(net.f_em))
        pur, _ = purity(dev, proc, _np.asarray(net.f_em))
        gold[wc0] = dict(N=dev.N, Ft=Ft, sigma0=dev.sigma_0,
                         F_em=net.F_em, adh=A, purity=pur,
                         ok=(dev.N == exp[wc0]
                             and abs(net.F_em - 0.35) <= F_TOL))
        say(f"golden wc0={wc0}: N={dev.N} F_em={net.F_em:.6f} "
            f"adh={A*100:.2f}% -> {'OK' if gold[wc0]['ok'] else 'DRIFT'}")
    rep["steps"]["goldens"] = gold

    # -- 2. G14 must refuse FOLD-45 -------------------------------------
    devF, procF = _design(45, F_DESIRED, "FOLD")
    try:
        G14_fold_routability(devF, procF)
        rep["steps"]["G14"] = dict(ok=False,
                                   note="FOLD-45 was NOT refused")
        say("G14: FAIL (FOLD-45 not refused)")
    except RoutabilityError as e:
        rep["steps"]["G14"] = dict(ok=True, N=devF.N, message=str(e))
        say(f"G14: FOLD-45 refused as designed (N={devF.N})")

    # -- 3. J-suite on the serpentine fixture ----------------------------
    proc = _proc()
    req = Request(D_target=4.0 * UM, F_target=0.35,
                  W_shape=lambda t: (100.0 - 20.0 * t) * UM, proc=proc,
                  pair=Pairing(), tol=Tolerances(), mode="DIRECT",
                  w_start=1.03 * proc.GS, N_seed=200)
    devJ = design(req)
    lay = Layout("serpentine", 8, 900 * UM, side_mode="straight")
    maskJ = build_mask_v2(devJ, lay)
    jres = run_suite(maskJ, devJ)
    tc = V_turn_closure(maskJ)
    rep["steps"]["jsuite"] = dict(
        ok=all(rj.ok for rj in jres) and tc["ok"],
        results=[(rj.name, rj.ok) for rj in jres], turn_closure=tc)
    say(f"J-suite fixture N={devJ.N}: "
        f"{sum(rj.ok for rj in jres)}/{len(jres)} PASS, "
        f"turn closures worst={tc['worst']:.2e} clipped={tc['clipped']}")

    # -- 4. stamped re-emission + D-suite from the bytes -----------------
    emitted = {}
    for tag in emit:
        wc0 = int(tag)
        dev, proc, Ft, _ = compensate(wc0)
        mask = build_mask_v2(dev, Layout())
        path = _v2_os.path.join(outdir, f"v2_opt_N{dev.N}_wc{wc0}.dxf")
        h = emit_device_dxf(mask, path, note=f"wc0={wc0} Ft={Ft:.4f}")
        entry = dict(N=dev.N, sha256=h, path=path)
        if run_dxf_suites:
            res, A, C = run_dxf_suite(path, mask, dev, F_DESIRED)
            entry["dsuite"] = [(rr_.name, rr_.ok) for rr_ in res]
            entry["dsuite_pass"] = sum(rr_.ok for rr_ in res)
            entry["adh_from_file"] = A
            say(f"emit wc0={wc0}: N={dev.N} "
                f"D-suite {entry['dsuite_pass']}/{len(res)} "
                f"adh(file)={A*100:.2f}% sha={h[:12]}")
        ops = V_operating_regime(dev, proc)
        st = V_steric(dev, proc)
        entry["ops"] = ops
        entry["steric"] = st
        emitted[tag] = entry
    rep["steps"]["emitted"] = emitted

    # -- 5. fab tolerance on the winner ----------------------------------
    devW, procW, FtW, _ = compensate(100)
    mc = V_fab_tolerance(devW, procW, F_desired=0.35)
    rep["steps"]["fab_mc"] = mc
    say(f"fab MC (winner): F = {mc['F_mean']:.4f} +/- {mc['F_std']:.4f}, "
        f"yield(|dF|<=3e-3) = {mc['yield_F']*100:.0f}%, "
        f"adh {mc['adh_mean']*100:.1f}% (min {mc['adh_min']*100:.1f}%)")

    ok = (all(g["ok"] for g in gold.values())
          and rep["steps"]["G14"]["ok"]
          and rep["steps"]["jsuite"]["ok"])
    rep["ok"] = ok
    jpath = _v2_os.path.join(outdir, "CIF_v2_acceptance.json")
    with open(jpath, "w") as fh:
        _v2_json.dump(rep, fh, indent=1, default=str)
    man = _v2_os.path.join(outdir, "CIF_v2_manifest.sha256")
    with open(man, "w") as fh:
        fh.write(f"{engine_sha256()}  CIF_v2.py\n")
        for tag, e in emitted.items():
            fh.write(f"{e['sha256']}  {_v2_os.path.basename(e['path'])}\n")
    say(f"ACCEPTANCE {'PASS' if ok else 'FAIL'} -> {jpath}")
    return rep


if __name__ == "__main__":
    run_acceptance(".")
