"""opencif.physics.kernels: the few formulas on which every other part of opencif is built.

Each function below is one result of the dissertation's Appendix A ("Step-by-Step Derivations Behind Controlled
Incremental Filtration") or one definition of the Ladder-Array Design Specification v5 (spec v5), written out so that
it can be read next to the derivation.

    duct_bracket, r              the resistance of a rectangular duct (Appendix A, Result 14, eq. 88; Gifford eqn 5)
    r_inverse                    the width of a duct with a given resistance per length (Result 14, "Monotonicity")
    g_cubic                      the fraction of the flow carried by a wall lamina (Appendix A, Result 1, eq. 10)
    x_of_f                       the exact inverse of the cubic (Appendix A, Result 2, eq. 20)
    Pairing                      the rule that kappa goes with the kernel it was fitted with (Appendix A, Result 9
                                 corollary, eq. 59), and kappa_g = 4.9138 [4.75, 5.08] on the cubic kernel
                                 (Appendix B, claim C3b)
    c2, Pairing.transfer_bias    the laminar kernel's leading coefficient c2(beta), used only for the design code's
                                 transfer-bias report (spec v5 S4 E-d and S5); opencif's kernel is the cubic, c2 = 3
    C_shear                      the peak wall shear of a rectangular duct (spec v5 S16-b)
    entrance_length              the laminar entrance length, reported and never gated (spec v5 S16-c)
    stadium_post                 the outline of a stadium-shaped post (spec v5 S14)
    Tolerances, Process          the tolerances (spec v5 S2.4) and the process inputs (fabrication floors, fluid)

These are the design code's own formulas. The design code, ``reference/CIF_v2.py`` (version 2.1.0), is kept in the
repository byte for byte, because Appendix A names it, and opencif does not import it. The test
``tests/test_kernels.py`` proves that every function here returns exactly the same floating-point numbers as the
design code (bit for bit, over grids of inputs), so the photomasks that opencif draws are those that the design code's
kernels would give. Each formula below therefore keeps the design code's order of operations, because another order
can change the last digit of a width, and with it the bytes of a photomask.

The functions take SI units (meters, Pa s) unless a docstring says otherwise.
"""
import math
from dataclasses import dataclass

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq

UM = 1e-6                      # meters per micrometer

#: The number of odd terms of the duct series (n = 1, 3, ..., 99), as in Gifford et al. 2014, p. 4499.
N_ODD_TERMS = 50

_ODD_CACHE = {}


def _odd(n_terms):
    """Return the odd integers 1, 3, ..., 2 n_terms - 1 as a float array, cached because r() is called many thousand
    times."""
    a = _ODD_CACHE.get(n_terms)
    if a is None:
        a = np.arange(1, 2 * n_terms, 2, dtype=float)
        _ODD_CACHE[n_terms] = a
    return a


# ---------------------------------------------------------------------------------------------------------------------
# Appendix A, Result 14: the resistance of a rectangular duct
# ---------------------------------------------------------------------------------------------------------------------

def duct_bracket(alpha, n_terms=N_ODD_TERMS):
    """Return the bracket of eq. 88 as a function of the aspect ratio alpha = (long side)/(short side) >= 1:

        B(alpha) = 1 - (192 / (pi^5 alpha)) * sum over odd n of tanh(n pi alpha / 2) / n^5

    Sorting the sides so that alpha >= 1 keeps B between 0.4217310 (a square duct) and 1 (a thin slot), so the
    subtraction never loses digits. (The design code calls this B.)"""
    n = _odd(n_terms)
    s = np.sum(np.tanh(np.minimum(n * np.pi * alpha / 2.0, 350.0)) / n ** 5)
    return 1.0 - (192.0 / (np.pi ** 5 * alpha)) * s


def r(w, h, mu, n_terms=N_ODD_TERMS):
    """Return the resistance per unit length of a duct of width w and depth h (SI, Pa s m^-4):

        r = 12 mu / (a b^3 B(a/b)),   a = max(w, h), b = min(w, h)            (Appendix A, eq. 88)

    A segment of length L has resistance R = r L. The resistance r decreases strictly as the width grows (Result 14,
    "Monotonicity"), so each resistance belongs to exactly one width (see r_inverse). For a square duct,
    r h^4 / mu = 28.4542 (Result 14, "Limits")."""
    if not (w > 0 and h > 0 and mu > 0):
        raise ValueError("r: w, h, mu must be strictly positive")
    a, b = (w, h) if w >= h else (h, w)
    br = duct_bracket(a / b, n_terms)
    if br <= 0:
        raise ValueError("r: non-positive bracket -- increase n_terms or check inputs")
    return 12.0 * mu / (a * b ** 3 * br)


def r_inverse(y, h, mu, n_terms=N_ODD_TERMS):
    """Return the width w of a duct of depth h whose resistance per length is y, which is the unique root of r(w) = y.

    The search starts from the thin-slot estimate w0 = (12 mu / (h y))^(1/3), widens the bracket until it holds the
    root, and finishes with Brent's method. Gifford increased w step by step, and the result is the same width."""
    seed = (12.0 * mu / (h * y)) ** (1.0 / 3.0)
    lo, hi = seed * 0.25, seed * 4.0
    fn = lambda w: r(w, h, mu, n_terms) - y
    k = 0
    while fn(lo) < 0:                # r(lo) is too small, so lo is too wide and is halved
        lo *= 0.5
        k += 1
        if k > 200 or lo < 1e-12:
            raise ValueError("r_inverse: cannot bracket below")
    k = 0
    while fn(hi) > 0:                # r(hi) is too big, so hi is too narrow and is doubled
        hi *= 2.0
        k += 1
        if k > 200:
            raise ValueError("r_inverse: cannot bracket above")
    return brentq(fn, lo, hi, xtol=1e-16, rtol=1e-13)


# ---------------------------------------------------------------------------------------------------------------------
# Appendix A, Results 1 and 2: the skim cubic and its inverse
# ---------------------------------------------------------------------------------------------------------------------

def g_cubic(x):
    """Return the fraction of the center channel's flow carried by a lamina of width delta next to one wall, with
    x = delta / w_c. The velocity profile between two walls is a parabola, and its integral across the lamina gives

        f = g(x) = 3 x^2 - 2 x^3                                              (Appendix A, Result 1, eq. 10)

    Thus g(1/2) = 1/2, and half the channel carries half the flow."""
    return 3.0 * x ** 2 - 2.0 * x ** 3


def x_of_f(f, beta=1.0, kernel="cubic"):
    """Return the lamina x = delta / w_c that carries the fraction f, which is the root of g(x) = f on [0, 1/2].

    The root is unique because g increases strictly (Appendix A, Result 2), and eq. 20 gives it in closed form
    (Viete's formula). For small f the root is close to sqrt(f/3) (eq. 18), and at the calibration cap f = 5.76e-4 the
    two differ by 0.47 %. opencif uses the exact root. Only the cubic kernel is provided, because the cutoff prefactor
    kappa_g was fitted with it (see Pairing); ``beta`` is accepted for the design code's signature and does not enter
    the cubic."""
    if not (0.0 < f < 0.5):
        raise ValueError(f"x_of_f: need 0 < f < 1/2 (hard bound, every beta), got {f}")
    if kernel != "cubic":
        raise ValueError("x_of_f: opencif uses the cubic kernel only (the pairing rule; see Pairing)")
    return brentq(lambda x: g_cubic(x) - f, 0.0, 0.5, xtol=1e-16, rtol=1e-15)


# ---------------------------------------------------------------------------------------------------------------------
# Appendix A, Result 9 corollary and Appendix B, claim C3b: the cutoff prefactor and its kernel
# ---------------------------------------------------------------------------------------------------------------------

@dataclass(frozen=True)
class Pairing:
    """The cutoff law D = kappa * delta and the kernel with which kappa was fitted.

    A prefactor is meaningful only with the kernel it was fitted with, which is the "pairing rule" of Appendix A,
    Result 9 corollary (eq. 59). Appendix B (claim C3b) fits kappa = 4.91 [4.75, 5.08] on Gifford's calibration array
    with the cubic kernel, and the glossary's value is kappa_g = 4.9138. opencif uses the cubic kernel with kappa_g
    everywhere.

    The prefactor kappa is constant only on a plateau of the reduced size d~ = D / w_c >= u_lo = 0.035, and below a
    foot at 0.012 it is a different constant (about 1.34). The value eta = 0.903 is the measured ratio of Dinh's
    filtrate fraction to the 1-D prediction (a 3-D outlet effect); opencif reports F from the 1-D network and does not
    apply eta. The fields f_soft and gamma_GS are spec v5's soft bound on the per-gap fraction and the band for the
    meander channel width GS/G. Finally, beta_cal records the calibration array's assumed aspect ratios (unpublished),
    which enter only the transfer-bias report (transfer_bias) and never the design."""
    kernel: str = "cubic"
    kappa: float = 4.9138
    kappa_ci: tuple = (4.750, 5.078)
    beta_cal: tuple = (0.80, 0.64, 0.53)
    u_lo: float = 0.035
    u_ff: float = 0.012
    eta: float = 0.903
    f_soft: float = 5e-4
    gamma_GS: tuple = (1.4, 1.5)

    def D(self, delta):
        """Return the cutoff D = kappa * delta of a lamina delta (Appendix A, Result 5)."""
        return self.kappa * delta

    def transfer_bias(self, betas):
        """Return the design code's mandatory report (spec v5 S5, P3), the ratio
        D_true / D_predicted = sqrt(c2(beta_cal) / c2(beta)) at each aspect ratio beta = h / w_c, as (min, max) over
        the betas.

        Here c2(beta_cal) is the mean over the calibration array's assumed aspect ratios. A kappa fitted with the cubic
        carries c2(beta_cal), and used at another aspect ratio it biases D by this factor. The factor is reported and
        never applied, as in the design code."""
        c_cal = float(np.mean([c2(b) for b in self.beta_cal]))
        bias = [math.sqrt(c_cal / c2(b)) for b in betas]
        return min(bias), max(bias)


# ---------------------------------------------------------------------------------------------------------------------
# spec v5 S4: the laminar kernel's leading coefficient (only for the design code's transfer-bias report)
# ---------------------------------------------------------------------------------------------------------------------

def _t_minus_tanh(t):
    """Return t - tanh t without cancellation, by its Taylor series below t = 1/2 (t^3/3 - 2t^5/15 + ...) and by the
    direct form above (the design code's E-e.1 guard)."""
    t = np.asarray(t, dtype=float)
    out = np.empty_like(t)
    lo = t < 0.5
    ts = t[lo]
    out[lo] = (ts ** 3 / 3.0 - 2 * ts ** 5 / 15.0 + 17 * ts ** 7 / 315.0
               - 62 * ts ** 9 / 2835.0 + 1382 * ts ** 11 / 155925.0)
    th = t[~lo]
    out[~lo] = th - np.tanh(np.minimum(th, 350.0))
    return out


def _A_n(beta, n_terms):
    """Return the odd n and the duct series' coefficients A_n(beta) = (2 / (pi n^5)) (t_n - tanh t_n), with
    t_n = n pi beta / 2 (all > 0)."""
    n = _odd(n_terms)
    return n, (2.0 / (np.pi * n ** 5)) * _t_minus_tanh(n * np.pi * beta / 2.0)


def c2_exact(beta, n_terms=200_000):
    """Return the leading coefficient c2(beta) of the exact laminar skim fraction F_lam ~ c2 x^2 of a w x h duct, with
    beta = h / w (spec v5 S4 E-d):

        c2(beta) = lim_{x->0} F_lam / x^2 = (pi^2 / 4) sum n^2 A_n / sum A_n

    The sum of n^2 A_n converges only like 1/n^2, so its analytic tail beta / (2 n_last) is added (E-e.2). As
    beta -> infinity, c2 -> 3 (the cubic)."""
    n, A = _A_n(beta, n_terms)
    S0 = float(A.sum()) + beta / (6.0 * n[-1] ** 3)
    S2 = float((n ** 2 * A).sum()) + beta / (2.0 * n[-1])
    return (np.pi ** 2 / 4.0) * S2 / S0


_C2 = {}


def c2(beta):
    """Return c2(beta) from a monotone spline (PCHIP in log-log) through c2_exact at 240 aspect ratios from 1e-3 to 1e4,
    as the design code caches it. The spline agrees with c2_exact to about 1e-9 and is built on the first call."""
    if 'spline' not in _C2:
        grid = np.exp(np.linspace(math.log(1e-3), math.log(1e4), 240))
        vals = np.array([c2_exact(b) for b in grid])
        _C2['spline'] = PchipInterpolator(np.log(grid), np.log(vals), extrapolate=True)
    return float(np.exp(_C2['spline'](math.log(beta))))


# ---------------------------------------------------------------------------------------------------------------------
# spec v5: shear, the entrance length, the post outline, the tolerances, and the process inputs
# ---------------------------------------------------------------------------------------------------------------------

def C_shear(alpha, n_terms=2000):
    """Return the peak wall-shear coefficient C(alpha) of a rectangular duct, for which gamma_max = C(alpha) Q / (a b^2)
    (spec v5 S16-b), with

        C(alpha) = (48 / (pi^2 B(alpha))) * sum over odd n of (1/n^2) [1 - sech(n pi alpha / 2)]

    C tends to 6 only for a thin slot (alpha -> infinity), and at alpha = 2 it is 8.13. The sum's tail beyond the last
    term is added analytically, because the sum over odd n > n_max of 1/n^2 is about 1/(2 n_max)."""
    n = _odd(n_terms)
    s = np.sum((1.0 / n ** 2) * (1.0 - 1.0 / np.cosh(np.minimum(n * np.pi * alpha / 2.0, 350.0))))
    s += 1.0 / (2.0 * n[-1])
    return float(48.0 / (np.pi ** 2 * duct_bracket(alpha, n_terms)) * s)


def entrance_length(Re, D_h):
    """Return the laminar entrance length (spec v5 S16-c), which keeps a Stokes floor as Re -> 0:

        L_e / D_h = [0.619^1.6 + (0.0567 Re)^1.6]^(1/1.6)

    The textbook L_e = 0.05 Re D_h vanishes as Re -> 0 and is wrong in this regime. opencif reports L_e against the
    row pitch and never gates on it, and the perturbation it stands for is gated instead by the gap-dominance gate
    M13, as the glossary says. Spec v5's "O(f^2)" order for it is not used as a bound, because that order needs a
    first-order cancellation that no record shows (docs/ALIGNMENT.md AG9)."""
    return D_h * (0.619 ** 1.6 + (0.0567 * Re) ** 1.6) ** (1.0 / 1.6)


class GateError(RuntimeError):
    """The error raised for an input outside what the formulas allow, for example a post whose caps would cross."""


def stadium_post(xc, yc, L_along, W_perp, n_arc=16):
    """Return the outline of a stadium (pill-shaped) post centered at (xc, yc), with its long axis across the flow
    (spec v5 S14).

    L_along is the post's extent along the flow (g_p), and W_perp is its extent across the flow (l_g, the length of the
    gap that the skimmed fluid passes through). The two ends are semicircles of radius g_p / 2, so W_perp must be at
    least L_along, or the two caps would cross. The outline is returned as its vertices, with each cap drawn with n_arc
    segments."""
    r_cap = L_along / 2.0
    hs = (W_perp - 2.0 * r_cap) / 2.0
    if hs < 0:
        raise GateError(f"stadium_post: W_perp={W_perp:.1f} < L_along={L_along:.1f} -- the caps "
                        f"would CROSS (hs={hs:.1f}).  That is a bowtie, not a post.")
    v = []
    for j in range(n_arc + 1):                       # top cap, 0 -> pi
        a = j * math.pi / n_arc
        v.append((xc + r_cap * math.cos(a), yc + hs + r_cap * math.sin(a)))
    for j in range(n_arc + 1):                       # bottom cap, pi -> 2 pi
        a = math.pi + j * math.pi / n_arc
        v.append((xc + r_cap * math.cos(a), yc - hs + r_cap * math.sin(a)))
    return np.array(v)


@dataclass(frozen=True)
class Tolerances:
    """The design code's numerical tolerances and regime limits (spec v5 S2.4), of which opencif uses the regime limits:

        eps_lam   = 0.10   thin lamina: delta / h <= 0.10                      (spec v5 E11)
        eps_taper = 1e-2   quasi-developed taper: |w_c(i) - w_c(i-1)| / Lambda  (spec v5 M16)
        Re_max    = 10     channel Reynolds number                            (spec v5 M14; Appendix B)
        Re_stokes = 1e-2   gap and lamina Reynolds numbers                    (spec v5 M14)
        gap_dom   = 0.05   gap dominance R_gap f / R_c                        (spec v5 M13)"""
    eps_r: float = 1e-9
    eps_inv: float = 1e-12
    eps_cub: float = 1e-3
    eps_id: float = 1e-10
    eps_geo: float = 1e-9
    eps_net: float = 2e-3
    eps_lam: float = 0.10
    eps_taper: float = 1e-2
    Re_max: float = 10.0
    Re_stokes: float = 1e-2
    gap_dom: float = 0.05


@dataclass(frozen=True)
class Process:
    """The design code's default process (SI). opencif takes every process value from the specification and uses this
    class only for the default aspect-ratio limit AR_max = 10 (a channel must be at least h / AR_max wide)."""
    h: float = 140 * UM
    mu: float = 1.0e-3
    rho_f: float = 1000.0
    F_min: float = 19 * UM
    R_min: float = 0.0
    AR_max: float = 10.0
    G: float = 20.5 * UM
    g_p: float = 25.6 * UM
    l_g: float = 50.6 * UM
    c_ent: float = 1.0

    @property
    def w_min(self):
        """The narrowest buildable channel, max(F_min, h / AR_max) (spec v5 S12)."""
        return max(self.F_min, self.h / self.AR_max)

    @property
    def R_gap(self):
        """The resistance of one gap, r(G) (l_g + c_ent D_h) (spec v5 S14.4), which is the duct law plus an entrance
        allowance of one hydraulic diameter."""
        D_h = 2 * self.G * self.h / (self.G + self.h)
        return r(self.G, self.h, self.mu) * (self.l_g + self.c_ent * D_h)
