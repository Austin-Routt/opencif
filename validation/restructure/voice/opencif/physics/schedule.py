"""opencif.physics.schedule: the row schedule of a CIF device, that is, the march at a constant lamina, the head, and
the exact closure.

The schedule uses the frozen physics of the v114 record. Its kernels are those of opencif.physics.kernels, which return
the design code's numbers bit for bit (tests/test_kernels.py). The cubic skim kernel is f = g(x) = 3x^2 - 2x^3 with
x = delta/w_c, and the cutoff is D = kappa_g * delta* with kappa_g = 4.9138 (the cubic pair, which is never mixed with
kappa_F). A constant cutoff requires f(i) = g(delta*/w_c(i)), which is the corrected Xia 2016 eqn 1 in its exact cubic
form. The march is P(i) = (1 - 2 f(i)) P(i-1), with F = 1 - P(N). The duct law r(w, h) is Gifford eqn 5 with 50 odd
terms, sorted (kernels.r, which equals CIF_v2.r). The gap resistance is R_gap = r(w_gap) [l_g + c_ent D_h(w_gap, h)]
(spec v5 section 14.4, the official Process.R_gap form).

The closure is the exact flow-prescribed pressure field with gap drops (spec v5 section 13), which generalizes App. A
Result 10 (eq. 69) and Gifford eqn 4 (App. A eq. 64) (internal review patch RP-B5):

      P_c(i+1) = P_c(i) - R_c(i) P(i)        center node pressures, where P(i) is the center flow after gap i
      P_s(i)   = P_c(i) - R_gap(i) Q_gap(i)   side node of gap i
      R_s(i)   = (P_s(i) - P_s(i+1)) / S(i)   side segment i -> i+1, which carries S(i) on each side

The closure reduces to Gifford eqn 4 (pressure equilibration) when R_gap Q_gap is equal at consecutive rows, and it
stays exact where the gap geometry or the row pitch changes (the arm rows, the seam between the meander and the linear
section, and the U-turn loops). Flows are per unit inlet flow, resistances are in SI units (Pa s m^-3), and lengths are
in um unless marked otherwise.

An 'A' row is an arm row, a meander unit drawn as an arm, with channel width w_arm, mouth w_arm, and pitch 2 s_a. A
'T' row is a taper row, a meander unit drawn as a perpendicular hairpin, with width GS, mouth GS, and pitch 2 s. An
'L' row is a linear row, a G-wide gap between stadium posts, with a straight side channel of width w_s(i) and pitch
G + g_p. The schedule also marks an 'N' row, a taper row of a narrow head drawn at w_narrow, and an 'R' row, one of the
n_rect rectilinear rows after the seam, with a GS-wide opening beside a wall block and pitch 2 s (the last one in the
papers' transition excepted). The last meander row
is the transition segment, whose unit's drop turns into the straight side channel.
"""
import math
from dataclasses import dataclass, field
from typing import Optional
import numpy as np

import os
PKG = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))     # the package directory, .../opencif
SP = os.path.join(PKG, 'data')                                         # the measured paper data (r4/*.json)
from opencif.physics import kernels as C    # the kernels, bit-identical to reference/CIF_v2.py (tests/test_kernels.py)

UM = 1e-6
MU = 1.0e-3
D_TR = 24.0          # extra drop length (um) of the transition unit, which makes room for the GS -> w_s(in) width taper
KAPPA_G = C.Pairing().kappa          # 4.9138, the cubic pair
KAPPA_CI = C.Pairing().kappa_ci      # (4.750, 5.078)
U_LO = C.Pairing().u_lo              # 0.035
U_FF = C.Pairing().u_ff              # 0.012, the foot of the kappa curve, below which kappa is another constant
# The kappa(u) collapse curve below the knee is the Hill fit (kmax, u_half, n) frozen in the glossary (\hillKmax,
# \hillUhalf, \hillN) and used by App. B P-C7.1 for Dinh's band. The test tests/test_appendix_trace.py checks it against
# the macros.
HILL = (5.568, 0.01651, 1.791)
# The calibration range of kappa_g in f_gap ends at 5.76e-4, the cap of Gifford's 33-device array (glossary
# \fgapDesignCap). Above it the cutoff law is extrapolated, which the record reports rather than gates.
FGAP_CAP = 5.76e-4


def rr(w, h):
    """Return the resistance per meter (SI) of a w x h duct, with w and h in um."""
    return C.r(w * UM, h * UM, MU)


_BEND_CACHE = {}


def clear_caches():
    """Empty the module's memo caches (bend_ratio, depth_lambda).

    Their keys are rounded inputs, so a cached value is the one computed at the first exact inputs that reached the
    key. A later build in the same process would therefore start from an earlier build's values and could differ from
    a build on its own in the last printed digit (BUILD-1). pipeline.run calls this first, so that every photomask is a
    function of its input alone."""
    _BEND_CACHE.clear()
    _LAM_CACHE.clear()


def bend_ratio(r1, r2, h, n_modes=150):
    """Return the exact Stokes conductance of a curved rectangular duct (walls at radii r1 < r2, depth h, all in um)
    relative to the 1-D rule used everywhere else, which is the straight-duct conductance divided by the centerline
    length rm = (r1 + r2)/2.

    At Re -> 0, a duct curved in its own plane with a purely azimuthal pressure gradient carries a purely azimuthal flow
    u(r, z), with u = 0 on the four walls, that solves

        d2u/dr2 + (1/r) du/dr - u/r^2 + d2u/dz2 = -1/r      (mu = 1, unit pressure drop per radian)

    Expanding in sin(n pi z / h) (odd n) gives, for each mode, the modified Bessel equation of order 1 with the
    particular solution c_n / (k_n^2 r), and the flux follows in closed form from I_0 and K_0. The straight duct of the
    same section is the same series with u'' - k^2 u = -c_n. Two exact limits are checked in tests/test_units.py, in
    which an infinitely deep duct gives the planar annulus and a shallow one gives the Hele-Shaw value rm ln(r2/r1)/w.

    Dean (inertial) corrections to the flux (Dean 1928, 1 - 0.03058 (K/576)^2 with K = 2 Re^2 a/R) are not included.
    At Re <= 10 they stay below about 1.5e-4 for a/R <= 0.2 (about 1e-3 at a/R 0.5, and 3.7e-3 at a/R 1, the most any
    drawable arc reaches), where a is the half-width. build.gates reports the bound for every bend from a Reynolds
    bound that holds in every channel, where that bound is <= 10 (a circular-pipe result, used as an order of
    magnitude)."""
    key = (round(r1, 4), round(r2, 4), round(h, 4), n_modes)
    v = _BEND_CACHE.get(key)
    if v is not None:
        return v
    from scipy.special import ive, kve
    r1, r2, h = float(r1), float(r2), float(h)
    if r1 <= 0.01 * (r2 - r1):                      # only reached by root-finder probes far outside any drawable loop
        r1 = 0.01 * (r2 - r1)
    w, rm = r2 - r1, 0.5 * (r1 + r2)
    n = np.arange(1, 2 * n_modes, 2, dtype=float)
    k = n * math.pi / h
    c = 4.0 / (n * math.pi)
    e = np.exp(-k * w)
    # The homogeneous part is A ive(1,kr) e^{k(r-r2)} + B kve(1,kr) e^{-k(r-r1)}, scaled so that both terms stay O(1).
    i1a, i1b = ive(1, k * r1), ive(1, k * r2)
    k1a, k1b = kve(1, k * r1), kve(1, k * r2)
    ra, rb = -c / (k * k * r1), -c / (k * k * r2)
    det = (i1a * e) * (k1b * e) - k1a * i1b
    A = (ra * (k1b * e) - k1a * rb) / det
    B = ((i1a * e) * rb - i1b * ra) / det
    I0a, I0b = ive(0, k * r1), ive(0, k * r2)
    K0a, K0b = kve(0, k * r1), kve(0, k * r2)
    integ = (c / (k * k)) * math.log(r2 / r1) + A * (I0b - I0a * e) / k + B * (K0a - K0b * e) / k
    Q = np.sum(2.0 * h / (n * math.pi) * integ)
    qs = np.sum(2.0 * h / (n * math.pi) * (c / (k * k)) * (w - (2.0 / k) * np.tanh(k * w / 2.0)))
    v = float(Q / (qs / rm))
    if len(_BEND_CACHE) < 500000:
        _BEND_CACHE[key] = v
    return v


# Junction terms (reviews R3-1, R4-1). Where a 1-D path meets a geometric change, the flow is still adjusting, and the
# centerline convention miscounts a length of order w. 2-D Stokes solves in the two limits of a rectangular duct, the
# planar limit (deep, h >> w) and the Hele-Shaw limit (shallow, h << w), give the two terms below
# (rev6/r3_junction2d.py, rev7/r4_tjunc2d*.py). First, a bend's two straight-arc junctions remove kappa_j * w of the arc
# from the developed correction, with
#       kappa_j = 0.50 (deep; 0.46 to 0.52 over 90 and 180 degrees, r_m/w 0.98 to 4.8)  ...  0.33 (shallow)
# Second, a meander head's T-junction (the drop of unit k turning into the stem of unit k+1, with the gap's mouth below)
# is shorter by c_T * w than the convention of a centerline through the junction's center, with
#       c_T = -0.28 (deep)  ...  -0.47 (shallow)            (plain 90-degree corner: -0.24 / -0.44)
# A real duct lies between the limits. The interpolation weight lam(h/w) in [0, 1] is taken from the exact curved-duct
# solution of a tight reference bend (r_m = w), whose conductance moves from its deep to its shallow value as h/w falls,
#       lam = (ratio(h) - ratio_deep) / (ratio_shallow - ratio_deep)      (0.05 at h/w = 5, 0.33 at 1, 0.58 at 0.5)
# That weight is a proxy rather than a 3-D result. The record's model band therefore re-solves the drawn device with
# lam moved by +-max(lam, 0.1) (both junction terms together), and the 3-D calibration of the terms is future work
# (docs/LIMITATIONS.md, item 8).
KAPPA_J_DEEP, KAPPA_J_SHALLOW = 0.50, 0.33
C_T_DEEP, C_T_SHALLOW = -0.28, -0.47
# The drawn T-junction is filleted (Gm.smooth, radius r_fillet) at both re-entrant corners. 2-D solves with fillets
# (review R4b, rev7/r4b_fillet*.py) move c_T by about S_T * (r_fillet / w). At r_f/w = 0.18 the shift is -0.055
# (shallow, converged) and -0.15 +- 0.03 (deep, staircase walls), and at r_f/w = 0.26 it is -0.09 (shallow) and -0.21
# (deep). The term is used up to r_f/w = 0.3, and the band adds a quarter of the fillet term plus 0.02 w on either side.
S_T_DEEP, S_T_SHALLOW = -0.83, -0.33
KAPPA_J = KAPPA_J_DEEP                  # a legacy name for the deep value
KAPPA_J_RANGE = (KAPPA_J_SHALLOW, KAPPA_J_DEEP)
_LAM_CACHE = {}


def depth_lambda(h, w):
    """Return the interpolation weight of the junction terms at depth h and width w, from 0 (deep) to 1 (shallow).

    The weight comes from the exact curved-duct solution of the reference bend described above."""
    key = (round(h, 3), round(w, 3))
    v = _LAM_CACHE.get(key)
    if v is None:
        r1, r2 = 0.5 * w, 1.5 * w
        c = -0.5
        Mx = np.array([[r1, 1 / r1], [r2, 1 / r2]])
        A, B = np.linalg.solve(Mx, -np.array([c * r1 * math.log(r1), c * r2 * math.log(r2)]))
        F_ = lambda r: c * (r * r / 2 * math.log(r) - r * r / 4) + A * r * r / 2 + B * math.log(r)
        deep = (F_(r2) - F_(r1)) / ((1.0 / w) * w ** 3 / 12.0)
        shallow = math.log(3.0)
        v = min(max((bend_ratio(r1, r2, h) - deep) / (shallow - deep), 0.0), 1.0)
        _LAM_CACHE[key] = v
    return v


def junction_lambda(spec, w):
    """Return the junction weight used for this spec. It is the depth proxy, the proxy moved to an end of its band in
    the model scenarios (spec.junction_band = 'lo' or 'hi'), or the deep or shallow limit itself ('deep' or
    'shallow')."""
    lam = depth_lambda(spec.h, w)
    band = getattr(spec, 'junction_band', None)
    if band == 'lo':
        return max(lam - max(lam, 0.1), 0.0)
    if band == 'hi':
        return min(lam + max(lam, 0.1), 1.0)
    if band == 'deep':
        return 0.0
    if band == 'shallow':
        return 1.0
    return lam


def kappa_j(spec, w):
    """Return the bend junction term kappa_j at width w, placed between its deep and shallow values by the junction
    weight. The term is 0 when spec.junction_band is 'none'."""
    if getattr(spec, 'junction_band', None) == 'none':
        return 0.0
    lam = junction_lambda(spec, w)
    return KAPPA_J_DEEP + lam * (KAPPA_J_SHALLOW - KAPPA_J_DEEP)


# The width change at the end of the outlet stretch (reviews of rc4). Where the side channel (width w_s(N)) meets a lead
# of another width, the 1-D rule (each piece at its own r(w)) misses the flow's adjustment, as it does at a bend or a
# T-junction. 2-D solves in the two limits of a rectangular duct of uniform depth give the excess as a length of the
# narrower channel, c * w_narrow, at the width ratio R = w_wide / w_narrow. The solves
# (provenance/rc5_scripts/expfem.py) use Taylor-Hood finite elements in scikit-fem, with the planar limit by the
# dissipation integral and the Hele-Shaw limit by Laplace's equation, checked against the conformal
# (Schwarz-Christoffel) closed form to 4 digits at every ratio. They cover three cases. First, a lead wider than the
# side channel, with its inner (axis-facing) wall aligned with the side channel's, as every outlet draws a wider lead,
# is a one-sided step; its deep limit is tabulated below, and its shallow limit is the closed form
# step_shallow(R, True). Second, a lead narrower than the side channel and centered on it (Xia's ring) is a symmetric
# step; its deep limit is tabulated, and its shallow limit, step_shallow(R, False), is half the one-sided value at the
# same R. Finally, a linear flare of the outer wall at flare_deg per side (Dinh 2024 Fig. 2B), with the inner wall
# straight, is tabulated deep and shallow over the angle and R. The flare's own widening is integrated by the exact
# wall rule, so its table holds only the excess over that integral.
# Each value is placed between the limits by the depth weight at the narrow width, like kappa_j and c_T. It is banded by
# moving that weight (junction_band 'lo' / 'hi') and by +-OUTLET_EXP_WIDEN of the value, which covers the
# planar-versus-duct proxy, the lead turning outward right after the step, and the 5 um fillets, which are sharp corners
# in these 2-D solves. A 3-D calibration is future work.
STEP_ONE_DEEP = ((1.0, 0.0), (1.25, 0.1320), (1.50, 0.2652), (2.00, 0.4225), (2.40, 0.4883), (2.71, 0.5213), (3.00, 0.5432), (4.00, 0.5850), (6.00, 0.6146), (8.00, 0.6249))
STEP_SYM_DEEP = ((1.0, 0.0), (1.25, 0.0763), (1.50, 0.1576), (2.00, 0.2591), (2.40, 0.3040), (3.00, 0.3429), (4.00, 0.3739), (6.00, 0.3972), (8.00, 0.4056))
FLARE_RATIOS = (1.0, 1.25, 1.50, 2.00, 2.40, 3.00, 4.00)
FLARE_DEGS = (0.0, 5.0, 10.0, 14.5, 20.0, 30.0, 45.0)
FLARE_DEEP = ((0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
              (0.0, 0.0123, 0.0242, 0.0344, 0.0459, 0.0641, 0.0860),
              (0.0, 0.0192, 0.0380, 0.0546, 0.0745, 0.1090, 0.1553),
              (0.0, 0.0260, 0.0517, 0.0746, 0.1023, 0.1521, 0.2246),
              (0.0, 0.0287, 0.0571, 0.0825, 0.1132, 0.1687, 0.2511),
              (0.0, 0.0309, 0.0615, 0.0889, 0.1222, 0.1823, 0.2723),
              (0.0, 0.0326, 0.0649, 0.0939, 0.1291, 0.1929, 0.2887))
FLARE_SHALLOW = ((0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
              (0.0, 0.0059, 0.0107, 0.0142, 0.0177, 0.0228, 0.0288),
              (0.0, 0.0112, 0.0213, 0.0294, 0.0382, 0.0519, 0.0688),
              (0.0, 0.0196, 0.0380, 0.0537, 0.0718, 0.1018, 0.1415),
              (0.0, 0.0249, 0.0487, 0.0692, 0.0932, 0.1340, 0.1899),
              (0.0, 0.0314, 0.0617, 0.0881, 0.1193, 0.1736, 0.2501),
              (0.0, 0.0398, 0.0784, 0.1125, 0.1531, 0.2247, 0.3284))
OUTLET_EXP_WIDEN = 0.3


def step_shallow(R, one_sided=True):
    """Return the Hele-Shaw (potential-flow) excess of an abrupt width change of ratio R >= 1, in units of the narrow
    width.

    The value is the conformal-map closed form for a symmetric step of half-widths ratio 1/R, which is also the
    one-sided step of the same R, because the symmetry line is a no-flux line in potential flow. A symmetric step in
    units of its narrow full width is half of it."""
    if R <= 1.0 + 1e-12:
        return 0.0
    r = 1.0 / R
    c = ((r + 1.0 / r) * math.log((1.0 + r) / (1.0 - r)) - 2.0 * math.log(4.0 * r / (1.0 - r * r))) / math.pi
    return c if one_sided else 0.5 * c


def _interp_logR(tab, R):
    """Interpolate a table of (ratio, value) pairs linearly in ln R, clamped at its ends (R below 1 is read as 1)."""
    xs = [math.log(a) for a, _ in tab]
    ys = [b for _, b in tab]
    return float(np.interp(math.log(max(R, 1.0)), xs, ys))


def _flare_tab(tab, deg, R):
    """Interpolate a flare table bilinearly in (angle, ln R) over the solved grid, clamped at its edges."""
    lr = [math.log(a) for a in FLARE_RATIOS]
    col = [float(np.interp(deg, FLARE_DEGS, row)) for row in tab]
    return float(np.interp(math.log(max(R, FLARE_RATIOS[0])), lr, col))


def outlet_exp_c(spec, w_a, w_b, parts=False):
    """Return the excess length of the outlet stretch's end, in units of the narrower width, where the side channel
    w_a meets the lead w_b (an abrupt step, or the flare when spec.flare_deg is set and the lead is wider).

    The excess is 0 for equal widths or junction_band 'none'. With parts=True the call returns the tuple (c, c_deep,
    c_shallow, depth weight, kind)."""
    band = getattr(spec, 'junction_band', None)
    wn, ww = min(w_a, w_b), max(w_a, w_b)
    if band == 'none' or ww <= wn * (1.0 + 1e-6):
        return (0.0, 0.0, 0.0, 0.0, 'none') if parts else 0.0
    R = ww / wn
    flare = getattr(spec, 'flare_deg', None) is not None and w_b > w_a
    if flare:
        deg = float(spec.flare_deg)
        deep, shal, kind = _flare_tab(FLARE_DEEP, deg, R), _flare_tab(FLARE_SHALLOW, deg, R), 'flare'
    elif w_b > w_a:
        deep, shal, kind = _interp_logR(STEP_ONE_DEEP, R), step_shallow(R, True), 'step, one-sided'
    else:
        deep, shal, kind = _interp_logR(STEP_SYM_DEEP, R), step_shallow(R, False), 'step, centred'
    lam = junction_lambda(spec, wn)
    c = deep + lam * (shal - deep)
    if band == 'lo':
        c *= 1.0 + OUTLET_EXP_WIDEN
    elif band == 'hi':
        c *= 1.0 - OUTLET_EXP_WIDEN
    return (c, deep, shal, lam, kind) if parts else c


def outlet_exp_len(spec, w_a, w_b):
    """Return the same excess as a length (um) of the narrower channel."""
    return outlet_exp_c(spec, w_a, w_b) * min(w_a, w_b)


def outlet_exp_domain(spec, w_a, w_b):
    """Return None where the stretch-end term lies inside the solved 2-D tables, and otherwise a note of what was held
    at a table's end, because np.interp clamps silently (a review of rc5).

    The deep step tables run to R 8, and the shallow closed form has no end. The flare tables run to R 4 and 45 deg,
    and the spec limits flare_deg to [1, 45]."""
    wn, ww = min(w_a, w_b), max(w_a, w_b)
    if getattr(spec, 'junction_band', None) == 'none' or ww <= wn * (1.0 + 1e-6):
        return None
    R = ww / wn
    if getattr(spec, 'flare_deg', None) is not None and w_b > w_a:
        return ('flare width ratio %.2f beyond the tables\' %.1f: deep and shallow held there' % (R, FLARE_RATIOS[-1])
                if R > FLARE_RATIOS[-1] + 1e-9 else None)
    tab = STEP_ONE_DEEP if w_b > w_a else STEP_SYM_DEEP
    return ('step width ratio %.2f beyond the deep table\'s %.1f: deep held there (shallow: closed form)' % (R, tab[-1][0])
            if R > tab[-1][0] + 1e-9 else None)


def c_t(spec, w):
    """Return the length correction of the head's T-junction in units of w, which is negative because the junction is
    shorter than the centerline convention. The correction is the sharp-corner value plus the fillet term, each placed
    between its deep and shallow 2-D values by the depth weight."""
    band = getattr(spec, 'junction_band', None)
    if band == 'none' or not getattr(spec, 't_junction', True):
        return 0.0
    lam = junction_lambda(spec, w)
    rf = min(max(float(getattr(spec, 'r_fillet', 0.0) or 0.0), 0.0) / w, 0.3)
    fillet = rf * (S_T_DEEP + lam * (S_T_SHALLOW - S_T_DEEP))
    c = C_T_DEEP + lam * (C_T_SHALLOW - C_T_DEEP) + fillet
    widen = 0.25 * abs(fillet) + 0.02          # spread of the fillet solves (+-0.03 deep at r_f/w 0.18), corner details
    if band == 'lo':
        c += widen
    elif band == 'hi':
        c -= widen
    return c


# The head's square quarter-turn (meander.turns: square) is a right-angle corner whose outer (convex) corner is rounded
# by r_out and whose inner (re-entrant) corner, the end of the wall between the lanes, is rounded by r_in (DINH-1). The
# tables below give its excess over the 1-D centerline through the corner point, in units of w, from 2-D solves in the
# two limits of a rectangular duct (provenance/rc7_scripts/corner/corner2d.py, corner_table.json). The solves use
# Hele-Shaw finite volumes and planar Stokes psi-omega by the dissipation integral, with n = 40 and 80 cells per w and
# first-order extrapolation, and the solver was checked against rc6's sharp-corner values and its exact curved-duct
# solution of the arc turn. The outer rounding moves c by under 2e-4 over r_out/w from 0 to 0.31, because the corner's
# outer region is nearly stagnant, so the table runs over r_in/w at r_out/w 0.25. Beyond its last row, c is held there
# and the record says so.
CORNER_RIN = (0.0, 0.18, 0.32, 0.44)
# The corner is solved up to this outer rounding r_out / w (corner_table.json; code review R14).
CORNER_ROUT_SOLVED = 0.307
CORNER_DEEP = (-0.2352, -0.3882, -0.5287, -0.6370)
CORNER_SHALLOW = (-0.4417, -0.4892, -0.5567, -0.6105)


def c_corner(spec, w, parts=False):
    """Return the length correction of the square quarter-turn in units of w, which is negative because the flow cuts
    the corner. The correction is placed between its deep and shallow 2-D values by the depth weight and banded like
    c_t (junction_band 'lo' / 'hi')."""
    band = getattr(spec, 'junction_band', None)
    ri = max(float(getattr(spec, 'r_corner_in', 0.0) or 0.0), 0.0) / w
    held = ri > CORNER_RIN[-1] + 1e-9
    deep = float(np.interp(ri, CORNER_RIN, CORNER_DEEP))
    shal = float(np.interp(ri, CORNER_RIN, CORNER_SHALLOW))
    if band == 'none' or getattr(spec, 'bend_model', 'stokes') == 'none':
        # With no junction terms, or with the plain 1-D rule at bends (the corner's centerline through its point, as an
        # arc's is its centerline arc under bend_model: none), there is no term. The parts still report the table's
        # values and the depth weight, which is not applied (codereview2).
        return (0.0, deep, shal, depth_lambda(spec.h, w), held) if parts else 0.0
    lam = junction_lambda(spec, w)
    c = deep + lam * (shal - deep)
    widen = 0.25 * abs(c - (CORNER_DEEP[0] + lam * (CORNER_SHALLOW[0] - CORNER_DEEP[0]))) + 0.02
    if band == 'lo':
        c += widen
    elif band == 'hi':
        c -= widen
    return (c, deep, shal, lam, held) if parts else c


def arc_excess_len(spec, rm, theta, w):
    """Return the extra 1-D length (um) of one arc (centerline radius rm, |angle| theta in rad, width w) under the bend
    model. The extra length is (rm theta - kappa_j w) (factor - 1), the developed correction less its two junctions,
    with the length rm theta - kappa_j w never taken below zero."""
    if getattr(spec, 'bend_model', 'stokes') == 'none':
        return 0.0
    return max(rm * theta - kappa_j(spec, w) * w, 0.0) * (bend_factor(spec, rm - w / 2.0, rm + w / 2.0) - 1.0)


def bend_factor(spec, r1, r2):
    """Return the multiplier on the 1-D resistance of a curved piece, which is 1/bend_ratio under the Stokes bend model
    and 1 under 'none'."""
    if getattr(spec, 'bend_model', 'stokes') == 'none':
        return 1.0
    return 1.0 / bend_ratio(r1, r2, spec.h)


def rinv(rpm, h):
    """Return the width (um) whose resistance per meter is rpm, by the exact official inverse."""
    return C.r_inverse(rpm, h * UM, MU) / UM


def Dh(w, h):
    """Return the hydraulic diameter 2 w h / (w + h) of a w x h duct."""
    return 2.0 * w * h / (w + h)


def r_gap(w, l_g, h, c_ent=1.0):
    """Return the resistance (SI) of a gap of width w, link length l_g, and depth h (um), which is
    r(w) [l_g + c_ent D_h(w, h)] (spec v5 section 14.4)."""
    return rr(w, h) * (l_g + c_ent * Dh(w, h)) * UM


def f_of(delta, wc):
    """Return the cubic skim kernel f = g(x) with x = delta/w_c. The kernel is defined only for a lamina inside the half
    channel, x < 1/2."""
    x = np.asarray(delta / wc, float)
    if np.any(x >= 0.5) or np.any(x < 0):
        raise ValueError('lamina delta/w_c = %.4g outside [0, 1/2): the cubic kernel is not defined there' % float(np.max(x)))
    return C.g_cubic(delta / wc)


def x_of(f):
    """Return the lamina x = delta/w_c that carries the fraction f, the exact inverse of the cubic kernel."""
    return C.x_of_f(f, 1.0, 'cubic')


SHEAR_FORMULAS = ('exact', 'parallel_plates')


def shear(Q, wc, h, formula='exact'):
    """Return the maximum wall shear gamma = C(alpha) Q / (a b^2) (official S16-b), with Q relative.

    The default formula, 'exact', uses the exact Stokes coefficient C(alpha) of the rectangular duct (spec v5 S16-b).
    The formula 'parallel_plates' holds C at its wide-channel limit 6, so that gamma = 6 Q / (a b^2) (XIA-1). That form
    is a departure from spec v5 S16-b (S:456, constant shear with C; S:606, "max gamma_dot via C(alpha) [NEVER 6]"),
    allowed in reproduce mode only, and it is the form that Xia 2016's Fig. 3 panels d and e read (w_c proportional to
    the center flow). Xia's text does not state the shear coefficient. The formulas it cites are Yang et al.'s
    rectangular-duct formulas, whose Eqn 3 the glossary identifies as the resistance series (G:2036-2038). The record
    reports the true peak shear via C(alpha) along the taper (validation.gates.shear_peak)."""
    a, b = (wc, h) if wc >= h else (h, wc)
    if formula == 'parallel_plates':
        return 6.0 * Q / (a * b * b)
    return C.C_shear(a / b) * Q / (a * b * b)


class ShearTable:
    """A monotone table of shear(1, wc, h) for the constant-shear (Xia RS) taper. A lookup in the table ends with an
    exact refinement by bisection on `shear` with the table's formula (C_shear for 'exact')."""
    def __init__(self, h, wlo, whi, n=600, formula='exact'):
        self.h = h
        self.formula = formula
        self.w = np.linspace(wlo, whi, n)
        self.g = np.array([shear(1.0, w, h, formula) for w in self.w])     # decreasing in w

    def wc_for(self, target):
        """Return the width wc at which shear(1, wc) = target, clamped to the table's range."""
        g = self.g
        if target >= g[0]:
            return self.w[0]
        if target <= g[-1]:
            return self.w[-1]
        j = int(np.searchsorted(-g, -target))
        lo, hi = self.w[j - 1], self.w[j]
        for _ in range(40):
            mid = 0.5 * (lo + hi)
            if shear(1.0, mid, self.h, self.formula) > target:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)


@dataclass
class Spec:
    """The engine's specification of one device, with lengths in um. spec.to_engine builds it from a user spec, and
    opencif.published builds the papers' devices with it directly."""
    name: str
    h: float                      # depth (um)
    G: float                      # linear-section gap (um)
    g_p: float                    # post width along the flow (um), so that Lambda_L = G + g_p
    l_g: float                    # post or wall thickness, which is the link length (um), also of the meander wall to
                                  # the trunk
    GS: float                     # meander (taper-unit) channel width (um)
    wc0: float                    # center width at the inlet and through the meander (um)
    delta_star: float             # lamina width (um), which sets f(i) = g(delta*/w_c(i))
    F_target: float               # terminal filtrate fraction
    ws_in: float                  # side width (um) required at the first linear row (the seam criterion)
    W: float                      # device width envelope (um), both sides
    F_min: float = 19.0           # minimum feature, open and solid (um)
    wall_M: Optional[float] = None    # wall between meander strokes, F_min by default
    w_arm: Optional[float] = None     # arm channel width, by default w_narrow when that is set and GS otherwise
                                      # (spec 14.2 "width == GS")
    w_narrow: Optional[float] = None  # narrow-head option, in which the rows that would need arms at GS are drawn at
                                      # this width (arms and taper units) and the rest of the head at GS, because the
                                      # handover needs GS
    taper: str = 'shear'          # 'const' | 'shear' (Xia RS, a constant maximum wall shear after the seam)
    wc_min: float = 0.0           # floor of the shear taper (150 um for Xia; Dinh's published end is 52 um)
    shear_formula: str = 'exact'  # shear that the constant-shear taper holds, 'exact' C(alpha) | 'parallel_plates'
                                  # (6 Q/(a b^2), a departure from spec v5 S16-b allowed in reproduce mode only; XIA-1)
    turns: str = 'arc'            # the head's quarter turns, 'arc' (radius s/2) | 'square' (corners) (DINH-1)
    r_corner_out: float = 0.0     # outer-corner rounding of a square turn (um)
    r_corner_in: float = 0.0      # inner-corner rounding of a square turn (um), at the end of the wall between lanes
    r_junction: float = 0.0       # rounding of the wall corners where each head unit's stem meets the center channel
                                  # (the head gaps' entrance), for the drawing only because c_ent is external (DINH-7)
    leg_order: str = 'nested'     # paperclip legs 'nested' (0, +d, -d, +2d, ...) or side by side in 'path' order
                                  # (0, d, 2d, ..., with loops that then alternate their turning direction) (GIF4-1, B)
    rho_f: float = 1000.0         # the fluid's density (kg/m^3), used by the gates only (Reynolds numbers; ADJ-C5)
    mu_f: float = 1.0e-3          # its viscosity (Pa s), used by the gates only, because the schedule's resistance
                                  # ratios do not depend on mu for a Newtonian fluid (the resistances keep MU, a common
                                  # factor; ADJ-C5)
    c_ent: float = 1.0
    margin: float = 60.0          # lateral margin inside W (um)
    Q_in_uL_min: float = 30.0     # operating flow (reports only)
    head: str = 'pack'            # 'pack' = maximum packing | 'dinh3' = three goal-post arms (Dinh Fig. 2A) |
                                  # 'none' = no meander head (classic Gifford)
    n_arms_fixed: int = 3
    cap_frac: float = 1.0          # head planning uses cap_frac x the taper-unit capacity, which leaves headroom for
                                  # the fail-safe margin that lengthens the units after it (build.py lowers it if
                                  # needed)
    n_rect: int = 2                # rectilinear rows after the seam (Dinh 2024 Fig. 2B bottom, Xia 2016 Fig. 1c(ii)),
                                  # with a straight side channel beside wall blocks, GS-wide openings, and pitch 2s;
                                  # the stadium posts start after them
    stem_half: float = 300.0      # half-width (um) of the inlet stem or aggregate filter beside which the arms rise
    L_out: float = 200.0          # straight outlet stretch after the last gap (um)
    flare_deg: Optional[float] = None  # flared outlet lead, in which the side channel's outer wall diverges at this
                                     # angle (deg, per side) from w_s(N) to the lead width after `flare_start`;
                                     # None = the right-angle step at the end of L_out
    flare_start: float = 65.4        # straight run (um) from the last gap's center to the flare's onset (Dinh 2024
                                     # Fig. 2B, fourth panel), 65.4 +- 2 um at the bar read with its anti-aliased ends
                                     # (provenance/rc6_scripts/measure_rc6.py), whereas rc5 read 67 at the short bar
    filter: object = False        # aggregate filter: False | True ('grid' stand-in) | 'dinh' (Dinh 2024 Fig. 2A, two
                                  # stages)
    # ---- composition options (meander3) ----------------------------------------------------------------------
    f_head: Optional[tuple] = None   # prescribed gap fractions f(1..n) for the first n rows, for a head drawn as
                                     # depicted, that is, the f that makes each unit as long as the paper's figure
                                     # shows (see geometry/depicted.py)
    NM_force: Optional[int] = None   # fixed number of meander rows (head units, including the transition unit)
    R0_arms: Optional[float] = None  # lateral centerline of the innermost arm lane (um); None = beside the stem.
                                     # Dinh 2024 Fig. 2A puts the arms at the chip's outer edge (about 1.28 mm)
    transition: str = 'exact'        # 'exact' = rectilinear rows at their solved widths | 'paper' = the side channel
                                     # beside the wall blocks is ws_rect wide (GS in Dinh Fig. 2B and Xia Fig. 1c(ii))
                                     # and narrows to w_s(in) at the first post (the papers' deliberate pull-in)
    ws_rect: Optional[float] = None  # side width beside the blocks for transition 'paper' (GS by default)
    outlet: str = 'separate'         # 'separate' = each side to its own port (Gifford Fig. 1A/4) | 'combined' = both
                                     # sides wrap around the center port into one filtrate port (Dinh, Xia)
    port_r: float = 750.0            # access-port radius (um) of the side (filtrate) ports
    port_r_c: Optional[float] = None # center (retentate) port radius (port_r by default)
    port_r_in: Optional[float] = None  # inlet port radius (port_r by default)
    r_fillet: float = 5.0            # fillet radius (um) of the head and outlet smoothing
    rmin_ok: float = 4.5             # V2 curvature criterion (um)
    loop_pull_in: float = 0.0        # amount (um) subtracted from each loop's side widths (Xia 2016 p. 4)
    loop_pull_in_leg: bool = False   # also pull in the next leg's side channels, 'as well as the downstream side
                                     # channels of the next device leg' (p. 4)
    loop_pull_mode: str = 'width'    # 'width' = subtract loop_pull_in from the solved widths of both loop sides (Xia's
                                     # wording; the narrower inner side then gains more resistance than the outer, so
                                     # flow crosses between the sides at every loop, val/network_studies.py S2) |
                                     # 'resistance' = close the loop rule on the narrowed width w_s - loop_pull_in, so
                                     # that both sides gain the same fraction of resistance as the narrowed leg
    loop_ramp: Optional[float] = None  # center-path length (um) at each arc end over which the loop's side widths ramp
                                     # from the leg's to the loop's; None = a quarter of the arc. Xia's printed loop-1
                                     # widths are those of a constant-width loop (val/xia_loop_check.py)
    pc_last_short: float = 0.0       # paperclip distance (um) by which the last leg's post array ends before the
                                     # others' end, with its outlet stretch running on past the loop (in the Gifford
                                     # 2014 Fig. 4 inset, the array of leg 4 ends beside the still-straight leg 2)
    R_ext_c: float = 0.0             # off-chip resistance (SI) on the center (retentate) outlet, such as tubing
    R_ext_s: float = 0.0             # off-chip resistance (SI) on each filtrate port (a combined port carries both
                                     # sides)
    seam_at: str = 'segment'         # where ws_in is required, 'segment' = the first straight segment after the head
                                     # (meander2) | 'post' = the first stadium post, after the n_rect rectilinear rows
                                     # (Xia 2016 Fig. 1c(ii) has w_s(i=1) = G 'at the first set of pill-shaped posts')
    # ---- measured layouts (replicas drawn as the papers draw them) ---------------------------------------------
    leg_Y: Optional[tuple] = None    # paperclip leg axes (um, leg 1 = 0, +Y = leg 2's side); None = uniform pitch d.
                                     # Loop j still joins legs j and j+1 with radius |dY|/2 (unequal in Xia 2016
                                     # Fig. 3a)
    inlet_gap: Optional[float] = None  # gap (um) from the inlet port's edge to the head's top; None = 150 um. In
                                     # Xia Fig. 3a the head's top lies 1.45 mm below the port's center
    pc_inlet_offset: float = 0.0     # paperclip distance (um) by which the inlet port lies downstream of the loop-2
                                     # centers (2.55 mm in Xia Fig. 3a at rc6's reading, which the preset sets)
    pc_out_ext: Optional[float] = None  # paperclip ring outlet, the distance (um) by which the ring's top lies beyond
                                     # the centers of the loops at the last leg's end (3.44 mm in Xia Fig. 3a at rc6's
                                     # reading, which the preset sets); None = the clearance minimum
    outlet_style: str = 'teardrop'   # combined outlet, 'teardrop' (Dinh 2024 Fig. 2A) | 'ring' (Xia 2016 Fig. 3a), in
                                     # which both side leads follow one circle around the center port into a filtrate
                                     # port on the axis beyond it and enter it from opposite sides
    ring_R: float = 4179.0           # 'ring' centerline radius of the circle (um). These three defaults are rc5's
                                     # Fig. 3a reading, and no shipped spec uses them (spec.to_engine sets every ring
                                     # value, and rc6's Xia values are in data/r6/rc6_values.json)
    ring_dc: float = 1187.0          # 'ring' offset (um) of the circle's center beyond the center port's center
    ring_dfp: float = 4762.0         # 'ring' offset (um) of the filtrate port's center beyond the center port's center
    ring_rq: float = 250.0           # 'ring' radius (um) of the turn from the side channel onto the circle
    note: str = ''
    trim_F: bool = False             # trim the last gap so that F = F_target exactly (opencif F10)
    bend_model: str = 'stokes'       # 'stokes' = every arc carries the exact Stokes curved-duct correction
                                     # (bend_ratio) less its junctions; 'none' = the plain 1-D rule r(w) ds
                                     # (opencif R2-4)
    t_junction: bool = True          # the head's T-junction term (c_t; review R4-1)
    S0: float = 0.0                  # side flow assumed at gap 1 by the design closure (Gifford 2014 p. 4498's start
                                     # rule, which chooses w_s(1) at the fabrication floor instead of the sealed
                                     # value). The drawn device stays sealed, so its first gaps over-extract, and the
                                     # measured network reports by how much. 0 = the sealed-start design.

    def __post_init__(self):
        if self.wall_M is None:
            self.wall_M = self.F_min
        if self.port_r_c is None:
            self.port_r_c = self.port_r
        if self.port_r_in is None:
            self.port_r_in = self.port_r
        if self.w_arm is None:
            self.w_arm = self.w_narrow if self.w_narrow else self.GS

    def L_cap_w(self, w):
        """Return the capacity of a taper unit of channel width w at its full reach, which includes the corner term
        and is scaled by cap_frac."""
        s = w + self.wall_M
        R = self.y_max - s / 2.0 - w / 2.0
        yJ = self.yJ(w)
        rq = s / 2.0
        return self.cap_frac * ((R - yJ) + math.pi * s / 2.0 + (R - yJ - rq) + self.qturn(rq) + (s - rq) +
                                self.corner_term(w))

    def corner_term(self, w):
        """Return the resistance length c_corner w that the one corner of a square-turn taper unit adds (code review
        R7).

        The term is negative, because the corner's centerline through its point overstates its resistance, so the
        unit's capacity in effective length is that much smaller than its centerline. The term is 0 for arc turns,
        whose capacity leaves out the arcs' and the T-junction's terms (docs/LIMITATIONS.md, item 20)."""
        return c_corner(self, w) * w if self.turns == 'square' else 0.0

    @property
    def s(self):        # taper-unit stroke pitch
        return self.GS + self.wall_M

    @property
    def s_a(self):      # arm lane pitch
        return self.w_arm + self.wall_M

    @property
    def Lam_L(self):    # linear row pitch G + g_p
        return self.G + self.g_p

    def yJ(self, w):    # lateral centerline of the T-junction drop for a unit of width w (right side)
        return self.wc0 / 2.0 + self.l_g + w / 2.0

    @property
    def y_max(self):    # outermost allowed channel edge (right side)
        return self.W / 2.0 - self.margin

    @property
    def R_max(self):    # largest taper reach (U-turn center) whose outer edge stays inside y_max
        return self.y_max - self.s / 2.0 - self.GS / 2.0

    def L_taper(self, R):
        """Return the centerline length J_k -> J_k+1 of a taper unit whose outer U-turn is centered at reach R
        (geometry/shapes.py)."""
        s, yJ, rq = self.s, self.yJ(self.GS), self.s / 2.0
        return (R - yJ) + math.pi * s / 2.0 + (R - yJ - rq) + self.qturn(rq) + (s - rq)

    def qturn(self, r):
        """Return the centerline length of the unit's quarter turn, which is pi r / 2 for an arc and 2 r through a
        square corner's point (DINH-1; shapes.qturn_len)."""
        return 2.0 * r if self.turns == 'square' else math.pi * r / 2.0

    def R_of_L(self, L):
        """Return the reach R at which a taper unit of width GS has the centerline length L, the inverse of
        L_taper."""
        s, yJ, rq = self.s, self.yJ(self.GS), self.s / 2.0
        const = math.pi * s / 2.0 - rq + self.qturn(rq) + (s - rq)
        return (L - const) / 2.0 + yJ

    @property
    def L_cap(self):
        """The capacity of a GS-wide taper unit at full reach, which includes the corner term and is scaled by
        cap_frac."""
        return self.cap_frac * (self.L_taper(self.R_max) + self.corner_term(self.GS))

    @property
    def R0(self):       # centerline of the innermost arm lane (rise), clear of the stem by one wall (or as prescribed)
        r_min = max(self.stem_half + self.wall_M + self.w_arm / 2.0, self.yJ(self.w_arm) + self.s_a)
        return max(self.R0_arms, r_min) if self.R0_arms is not None else r_min

    @property
    def A_max(self):    # maximum number of arms (>= 2 lanes each) between R0 and y_max
        n_lanes = int(math.floor((self.y_max - self.w_arm / 2.0 - self.R0) / self.s_a)) + 1
        return n_lanes // 2


@dataclass
class Sched:
    """The planned schedule of one device, which holds the march, the row kinds and positions, and the closure's
    resistances, pressures, and side targets, row by row."""
    spec: Spec
    f_seed: float
    A: int                        # arm rows (1..A)
    N: int
    NM: int                       # meander rows (gaps 1..NM are meander mouths); the NM-th segment is the transition
    kind: np.ndarray              # 'A' / 'N' / 'T' / 'R' / 'L'
    wM: np.ndarray                # meander channel width per row (w_arm, w_narrow, or GS; nan for R and L)
    wL: np.ndarray                # link width per row, the mouth (w_arm / w_narrow / GS), the opening GS ('R'), or
                                  # the gap G ('L')
    x: np.ndarray                 # gap position along the leg path (um), gap 1 at 0
    ell: np.ndarray               # center distance gap i -> i+1 (um); the last is the outlet stretch
    wc: np.ndarray                # center width at gap i (um)
    f: np.ndarray
    P: np.ndarray
    Qg: np.ndarray
    S: np.ndarray
    Rc: np.ndarray
    Rg: np.ndarray
    Pc: np.ndarray
    Ps: np.ndarray
    Rs: np.ndarray                # side target i -> i+1 (len N-1)
    Rc_out: float
    Rs_out: float
    L_req: np.ndarray             # required centerline length (um) of meander rows 1..NM-1 at their own width
    ws: np.ndarray                # linear side width (um) for segments NM+1..N-1 (nan elsewhere)
    delta: np.ndarray             # lamina per row (um)
    F: float
    info: dict = field(default_factory=dict)


def build(spec: Spec, f1: Optional[float] = None, A: int = 0, st: Optional[ShearTable] = None,
          f_mult: Optional[dict] = None, NM_force: Optional[int] = None, K_narrow: int = 0) -> Sched:
    """Run the march, classify the rows, and compute the exact closure targets.

    The optional arguments are

        f1        the gap-1 extraction (None = f*)
        A         the number of arm rows
        st        the ShearTable of a constant-shear taper (built here when None)
        f_mult    {row index (0-based): multiplier} applied to f (the fail-safe seam margin)
        NM_force  the fixed number of meander rows, so that a margin does not move the seam
        K_narrow  the last narrow row, so that rows A+1..K_narrow are narrow taper rows ('N')"""
    h = spec.h
    fstar = f_of(spec.delta_star, spec.wc0)
    f1 = fstar if f1 is None else f1
    rcw0 = rr(spec.wc0, h)
    RgA = r_gap(spec.w_arm, spec.l_g, h, spec.c_ent)
    RgT = r_gap(spec.GS, spec.l_g, h, spec.c_ent)
    RgN = r_gap(spec.w_narrow, spec.l_g, h, spec.c_ent) if spec.w_narrow else RgT
    RgL = r_gap(spec.G, spec.l_g, h, spec.c_ent)
    r_wsin = rr(spec.ws_in, h)
    if spec.taper == 'shear' and st is None:
        st = ShearTable(h, max(spec.wc_min, 0.3 * spec.wc0), spec.wc0, formula=spec.shear_formula)
    kinds, wcs, fs, Ps, Qgs, Ss = [], [], [], [], [], []
    S0_ = float(getattr(spec, 'S0', 0.0) or 0.0)
    P_prev, S_prev = 1.0, S0_
    lin = False
    gamma_ref = None
    if spec.head == 'none':                     # classic Gifford, with no meander and straight side channels from gap 1
        lin = True
        gamma_ref = shear(1.0, spec.wc0, h, spec.shear_formula)
    if NM_force is None and spec.NM_force is not None:
        NM_force = spec.NM_force
    fh = list(spec.f_head) if spec.f_head is not None else []
    i = 1
    while True:
        if not lin:
            kind = 'A' if i <= A else ('N' if i <= K_narrow else 'T')
            wc = spec.wc0
        else:
            kind = 'L'
            if spec.taper == 'shear':
                wc = st.wc_for(gamma_ref / P_prev)
            else:
                wc = spec.wc0
        f = f1 if i == 1 else f_of(spec.delta_star, wc)
        if i <= len(fh):                        # prescribed fractions of a head drawn as depicted
            f = fh[i - 1]
        if f_mult and (i - 1) in f_mult:
            f *= f_mult[i - 1]
        Qg = f * P_prev
        P = P_prev * (1.0 - 2.0 * f)
        S = S_prev + Qg
        kinds.append(kind); wcs.append(wc); fs.append(f); Ps.append(P); Qgs.append(Qg); Ss.append(S)
        # The march stops on the delivered filtrate fraction. With an assumed start side flow S0 (a floor start), the
        # sealed device carries the design flows scaled by 1/(1 + 2 S0), so F_delivered = 1 - P/(1 + 2 S0), which is
        # 1 - P when sealed.
        if 1.0 - P / (1.0 + 2.0 * S0_) >= spec.F_target:
            break
        if not lin and NM_force is not None:
            if i == NM_force:
                lin = True
                gamma_ref = shear(P, spec.wc0, h, spec.shear_formula)
        elif not lin and i > max(A, K_narrow, len(fh)):
            # Test whether this row can be the first linear row, that is, whether a straight segment at the linear
            # pitch from here would be at least ws_in wide. The relabel leaves f unchanged, because w_c = wc0 at the
            # first linear row.
            Qg_next = f_of(spec.delta_star, spec.wc0) * P
            Rs_lin = (rcw0 * spec.Lam_L * UM * P - RgL * (Qg - Qg_next)) / S
            if spec.seam_at == 'post' and spec.n_rect > 0:
                # Project the march over the rectilinear rows (constant cutoff at w_c0) to the first stadium post.
                fq = f_of(spec.delta_star, spec.wc0)
                Pp, Sp, Qgp = P, S, Qg
                for _ in range(spec.n_rect):
                    Qgp = fq * Pp
                    Sp += Qgp
                    Pp *= (1.0 - 2.0 * fq)
                Rs_lin = (rcw0 * spec.Lam_L * UM * Pp - RgL * (Qgp - fq * Pp)) / Sp
            if i >= 2 and Rs_lin > 0 and Rs_lin / (spec.Lam_L * UM) <= r_wsin:
                kinds[-1] = 'L'
                lin = True
                gamma_ref = shear(P_prev, spec.wc0, h, spec.shear_formula)
        P_prev, S_prev = P, S
        i += 1
        if i > 400000:
            raise RuntimeError('march did not reach F_target')
    if getattr(spec, 'trim_F', False) and len(fs) >= 2:
        # The march stops at the first row whose delivered fraction reaches F_target, overshooting by up to one row's
        # extraction (opencif F10). The last gap's fraction is therefore trimmed so that the delivered fraction
        # 1 - P(N)/(1 + 2 S0), which is 1 - P(N) when sealed, equals F_target exactly. The trimmed f is smaller than
        # the row's constant-cutoff f, so that gap's cutoff is smaller than D (retentive), which is safe for the
        # cutoff.
        f_last = 0.5 * (1.0 - (1.0 - spec.F_target) * (1.0 + 2.0 * S0_) / P_prev)
        fs[-1] = f_last
        Qgs[-1] = f_last * P_prev
        Ps[-1] = P_prev * (1.0 - 2.0 * f_last)
        Ss[-1] = S_prev + Qgs[-1]
    N = len(fs)
    kinds = np.array(kinds)
    NM = int(np.sum(kinds != 'L'))
    if spec.head != 'none':
        for j in range(NM, min(NM + spec.n_rect, N - 1)):      # the rectilinear end of the head
            kinds[j] = 'R'
        if NM < 1 or NM < max(A, K_narrow) + 1:
            raise RuntimeError('head shorter than the requested arms')
    wc = np.array(wcs, float)
    f = np.array(fs); P = np.array(Ps); Qg = np.array(Qgs); S = np.array(Ss)
    wM = np.where(kinds == 'A', spec.w_arm, np.where(kinds == 'N', spec.w_narrow or spec.GS,
                                                      np.where(kinds == 'T', spec.GS, np.nan)))
    ell = np.zeros(N)
    for k in range(N - 1):
        if kinds[k] == 'L':
            ell[k] = spec.Lam_L
        elif kinds[k] == 'R':                               # wall block (2s - GS) + opening GS
            ell[k] = 2.0 * spec.s
        elif k == NM - 1:                                   # the transition unit, whose drop is extended by D_TR
            ell[k] = 2.0 * spec.s + D_TR
        else:                                               # from the return line of unit k to the out line of unit k+1
            wk, wn = wM[k], wM[k + 1]
            sk = wk + spec.wall_M
            ell[k] = sk + 0.5 * (wk + wn) + spec.wall_M
    if getattr(spec, 'transition', 'exact') == 'paper' and spec.head != 'none' and spec.n_rect >= 1:
        # The papers' transition ends 'opening (GS) - pill - gap (G) - pill' (Dinh 2024 Fig. 2B second panel; Xia 2016
        # Fig. 1c(ii)), so the last rectilinear opening is followed directly by a stadium post, and that row is
        # GS/2 + g_p + G/2 long. With Dinh's GS 28, g_p 26.1, and G 19 the row is 49.6 um (rc6), whereas rc5's g_p 27
        # gave 50.5. Measured on Dinh's panel, it is 49 um at rc5's short-bar reading and 48 um at the bar read with
        # its anti-aliased ends.
        k_last = NM + spec.n_rect - 1
        if k_last < N - 1 and kinds[k_last] == 'R' and kinds[k_last + 1] == 'L':
            ell[k_last] = spec.GS / 2.0 + spec.g_p + spec.G / 2.0
    ell[N - 1] = spec.L_out
    x = np.concatenate([[0.0], np.cumsum(ell[:-1])])
    Rc = np.zeros(N)
    for k in range(N):
        w0 = wc[k]
        w1 = wc[k + 1] if k + 1 < N else wc[k]
        if abs(w1 - w0) < 1e-12:
            Rc[k] = rr(w0, h) * ell[k] * UM
        else:
            wm = 0.5 * (w0 + w1)
            Rc[k] = (rr(w0, h) + 4 * rr(wm, h) + rr(w1, h)) / 6.0 * ell[k] * UM
    Rg = np.where(kinds == 'A', RgA, np.where(kinds == 'N', RgN, np.where((kinds == 'T') | (kinds == 'R'), RgT, RgL)))
    wL = np.where(kinds == 'L', spec.G, np.where(kinds == 'R', spec.GS, wM))
    Pc = np.zeros(N)
    for k in range(N - 1):
        Pc[k + 1] = Pc[k] - Rc[k] * P[k]
    Psd = Pc - Rg * Qg
    Rs = (Psd[:-1] - Psd[1:]) / S[:-1]
    Rc_out = Rc[N - 1]
    P_out = Pc[N - 1] - Rc_out * P[N - 1]
    Rs_out = (Psd[N - 1] - P_out) / S[N - 1]
    L_req = np.array([Rs[k] / (rr(wM[k], h) * UM) for k in range(NM - 1)])
    ws = np.full(N, np.nan)
    for k in range(NM, N - 1):
        ws[k] = rinv(Rs[k] / (ell[k] * UM), h)
    delta = np.array([x_of(fi) * wi if 0 < fi < 0.5 else float('nan') for fi, wi in zip(f, wc)])
    return Sched(spec=spec, f_seed=f1, A=A, N=N, NM=NM, kind=kinds, wM=wM, wL=wL, x=x, ell=ell, wc=wc, f=f, P=P, Qg=Qg,
                 S=S, Rc=Rc, Rg=Rg, Pc=Pc, Ps=Psd, Rs=Rs, Rc_out=Rc_out, Rs_out=Rs_out, L_req=L_req, ws=ws,
                 delta=delta, F=1.0 - P[-1] / (1.0 + 2.0 * S0_),
                 info=dict(fstar=fstar, RgA=RgA, RgT=RgT, RgN=RgN, RgL=RgL, gamma_ref=gamma_ref, st=st,
                           K_narrow=K_narrow))


def over_cap(sc: Sched) -> int:
    """Return the number of taper rows (narrow or GS) whose required length exceeds one taper unit at full reach."""
    n = 0
    for k in range(sc.NM - 1):
        if sc.kind[k] in ('T', 'N') and sc.L_req[k] > sc.spec.L_cap_w(sc.wM[k]):
            n += 1
    return n


def arms_needed(sc: Sched) -> int:
    """Return the length of the prefix of rows whose required length exceeds the taper cap, for a schedule planned
    with no arm rows. The over-cap rows must form a prefix of the head."""
    over = sc.L_req > sc.spec.L_cap
    if not over.any():
        return 0
    last = int(np.max(np.nonzero(over)[0]))
    if not over[:last + 1].all():
        raise RuntimeError('over-cap rows are not a prefix of the head')
    return last + 1


def plan(spec: Spec, tol: float = 2e-3) -> Sched:
    """Choose the arm count and the gap-1 extraction.

    For the 'pack' head, A is the number of rows that need arms with a sealed inlet when A <= A_max; otherwise
    A = A_max, and f1 is the smallest gap-1 extraction for which every later row fits a taper unit. For the 'dinh3'
    head, A = n_arms_fixed (3, as in Dinh Fig. 2A), and f1 is the smallest for which rows A+1 onward fit a taper unit.
    A design without a head ('none') and a head drawn as depicted (f_head) are built directly, and a narrow head
    (w_narrow) is planned by its own iteration."""
    st = ShearTable(spec.h, max(spec.wc_min, 0.3 * spec.wc0), spec.wc0, formula=spec.shear_formula) if spec.taper == 'shear' else None
    if spec.head == 'none':
        sc = build(spec, None, 0, st)
        sc.info.update(A_need_sealed=0, A_cap=0)
        return sc
    if spec.f_head is not None:                 # a head as depicted, with the first n_arms_fixed rows as arms (if any)
        A = spec.n_arms_fixed if spec.head == 'dinh3' else 0
        sc = build(spec, spec.f_head[0], A, st)
        sc.info.update(A_need_sealed=None, A_cap=A, depicted_rows=len(spec.f_head))
        return sc
    sc0 = build(spec, None, 0, st)
    A_need = arms_needed(sc0)
    if spec.w_narrow:
        # In a narrow head, the rows that would need arms at GS are drawn narrow, and of those, the ones over the
        # narrow taper cap become (narrow) arms. The plan iterates, because narrow rows change the pitch and the
        # mouths.
        K = max(A_need, 1)
        for _ in range(20):
            scn = build(spec, None, 0, st, K_narrow=K)
            over_gs = [k for k in range(K, scn.NM - 1) if scn.L_req[k] > spec.L_cap]
            if over_gs:
                K = max(over_gs) + 1
                continue
            capn = spec.L_cap_w(spec.w_narrow)
            An = sum(1 for k in range(K) if scn.L_req[k] > capn)
            if An > spec.A_max:
                raise RuntimeError('narrow head still needs %d arms > %d' % (An, spec.A_max))
            for _ in range(10):
                sc = build(spec, None, An, st, K_narrow=K)
                if over_cap(sc) == 0:
                    sc.info['A_need_sealed'] = A_need
                    sc.info['A_cap'] = spec.A_max
                    sc.info['K_narrow'] = K
                    return sc
                An += 1
            K += 1
        raise RuntimeError('narrow head planning did not converge')
    if spec.head == 'pack':
        A_cap = spec.A_max
    else:
        A_cap = spec.n_arms_fixed
    if spec.head == 'pack':
        A = min(A_need, A_cap)
        while True:                       # arm rows change the pitch and mouth of rows 1..A, so re-check and grow A
            sc = build(spec, None, A, st)
            if over_cap(sc) == 0:
                sc.info['A_need_sealed'] = A_need
                sc.info['A_cap'] = A_cap
                return sc
            if A >= A_cap:
                break
            A += 1
    else:
        A = A_cap
        sc = build(spec, None, A, st)
        if over_cap(sc) == 0:
            sc.info['A_need_sealed'] = A_need
            sc.info['A_cap'] = A_cap
            return sc
    lo, hi = sc0.info['fstar'], 80.0 * sc0.info['fstar']
    if over_cap(build(spec, hi, A, st)) > 0:
        raise RuntimeError('even f1 = 80 f* leaves rows over the taper cap')
    while hi / lo - 1 > tol:
        mid = math.sqrt(lo * hi)
        if over_cap(build(spec, mid, A, st)) == 0:
            hi = mid
        else:
            lo = mid
    sc = build(spec, hi, A, st)
    sc.info['A_need_sealed'] = A_need
    sc.info['A_cap'] = A_cap
    return sc


def reclose(sc: Sched, turns=(), R_turn=0.0, Rc_lead=0.0) -> Sched:
    """Redo the closure (positions, center resistances, pressures, and side targets) for a layout.

    The layout has a U-turn of axis radius R_turn (a number, or a dict of radii by row) after each row in `turns`, and
    an extra center outlet lead resistance Rc_lead (SI). Each U-turn is a gap-free semicircle, and its segment takes the
    exact target of its row (spec v5 section 14.5). The march (f, P, S, w_c) does not depend on the layout and is
    kept."""
    spec, h = sc.spec, sc.spec.h
    N, NM = sc.N, sc.NM
    ell = sc.ell.copy()
    turns = sorted(int(t) for t in turns)
    Rt = R_turn if isinstance(R_turn, dict) else {k: R_turn for k in turns}
    for k in turns:
        if not (NM <= k < N - 1):
            raise ValueError('turns must be after linear rows')
        ell[k] = spec.Lam_L + math.pi * Rt[k]
    x = np.concatenate([[0.0], np.cumsum(ell[:-1])])
    wc = sc.wc
    ell[N - 1] = spec.L_out
    Rc = np.zeros(N)
    tset = set(turns)
    for k in range(N):
        w0 = wc[k]
        w1 = wc[k + 1] if k + 1 < N else wc[k]
        L_eff = ell[k]
        if k in tset:                                     # the arc's resistance with the curved-duct correction
            wm_ = 0.5 * (w0 + w1)
            L_eff = spec.Lam_L + math.pi * Rt[k] + arc_excess_len(spec, Rt[k], math.pi, wm_)
        if abs(w1 - w0) < 1e-12:
            Rc[k] = rr(w0, h) * L_eff * UM
        else:
            wm = 0.5 * (w0 + w1)
            Rc[k] = (rr(w0, h) + 4 * rr(wm, h) + rr(w1, h)) / 6.0 * L_eff * UM
    Pc = np.zeros(N)
    for k in range(N - 1):
        Pc[k + 1] = Pc[k] - Rc[k] * sc.P[k]
    Psd = Pc - sc.Rg * sc.Qg
    Rs = (Psd[:-1] - Psd[1:]) / sc.S[:-1]
    Rc_out = Rc[N - 1] + Rc_lead
    P_out = Pc[N - 1] - Rc_out * sc.P[N - 1]
    Rs_out = (Psd[N - 1] - P_out) / sc.S[N - 1]
    ws = np.full(N, np.nan)
    for k in range(NM, N - 1):
        if k not in turns:
            ws[k] = rinv(Rs[k] / (ell[k] * UM), h)
    import copy
    sc2 = copy.copy(sc)
    sc2.ell, sc2.x, sc2.Rc, sc2.Pc, sc2.Ps, sc2.Rs, sc2.Rc_out, sc2.Rs_out, sc2.ws = ell, x, Rc, Pc, Psd, Rs, Rc_out, Rs_out, ws
    sc2.info = dict(sc.info, turns=turns, R_turn=(max(Rt.values()) if Rt else 0.0), R_turns=dict(Rt), Rc_lead=Rc_lead)
    return sc2


def R_of(sc, k):
    """Return the axis radius of the U-turn after row k, which is set per turn for a paperclip and is one value for a
    serpentine."""
    Rt = sc.info.get('R_turns')
    if Rt and k in Rt:
        return Rt[k]
    return sc.info['R_turn']


def turn_side_width(sc: Sched, k, inner: bool):
    """Return the width of one side channel through the U-turn after row k, chosen so that
    r(w) [Lam + pi (R -/+ y_c)] = R_s(k)."""
    spec, h = sc.spec, sc.spec.h
    R = sc.info['R_turn']
    wc_avg = 0.5 * (sc.wc[k] + sc.wc[k + 1])
    target = sc.Rs[k]

    def res(w):
        yc = wc_avg / 2.0 + spec.l_g + w / 2.0
        L = spec.Lam_L + math.pi * (R - yc if inner else R + yc)
        return rr(w, h) * L * UM - target
    from scipy.optimize import brentq
    if res(2.0) * res(5000.0) > 0.0:
        raise ValueError('U-turn after row %d: no %s side width between 2 and 5000 um meets its resistance target'
                         % (k + 1, 'inner' if inner else 'outer'))
    return brentq(res, 2.0, 5000.0, xtol=1e-12)
