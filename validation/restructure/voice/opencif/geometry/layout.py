"""opencif.geometry.layout: the assembly of a whole CIF device.

A device consists of the head (the arms meander), the linear legs, the U-turn loops, the inlet (with an optional
aggregate filter), and the outlets. The outlets are the center port and either a side port for each side channel (the
separate outlet) or one filtrate port for both side channels (the combined outlet). There are three layouts:

    'linear'      one straight leg
    'serpentine'  n legs side by side, joined by gap-free 180-degree loops of one radius at alternating ends. This is
                  a new layout, because the papers draw only linear and paperclip devices. Legs 2 to n run beside
                  the linear part of leg 1 only, so they stop short of the wide head.
    'paperclip'   n nested legs (Gifford 2014 Fig. 4; Xia 2016 Figs. 1a, 3a), with leg j at Y = 0, +d, -d, +2d, -2d,
                  and so on. Every loop turns the same way with radius |Y(j+1) - Y(j)| / 2, so the loops at one end
                  are concentric. The inlet (and the head) sit inside loop 2, and the outlet leaves the outermost leg
                  beyond the loops. With `layout.order: path`, the legs lie side by side in path order instead
                  (Y = 0, d, 2d, ...), the loops alternate their turning direction, and the inlet does not sit inside
                  loop 2 (Gifford Fig. 4's reading B, register GIF4-1).

Every loop is one segment of the ladder (spec v5 section 14.5). It takes the exact target R_s(k) of its row, and the
inner and outer side channels get their own widths so that each matches that target on its own path length.

The geometry follows these conventions. The path coordinate is the schedule's x (gap 1 at 0). In each leg, local x
runs along the flow and local +y is the left-hand side of the flow ('side A'), whereas side B is local -y. In the
linear and serpentine layouts, leg j's axis lies at global y = 2 R_turn (j-1). Odd legs flow +x and even legs flow -x,
which is a 180-degree rotation of the local frame."""
import math
import numpy as np
from shapely.geometry import Polygon, Point, box, LineString
from shapely.ops import unary_union
from shapely import affinity
import shapely
from opencif.physics import schedule as S
from opencif.geometry import shapes as Gm

# The default access pad radius (um). The per-device radii are spec.port_r (side), port_r_c, and port_r_in.
PORT_R = 750.0
INLET_CLEAR = 150.0     # clearance between the top of the head and the inlet port (um)
from opencif.geometry import dinh_inlet as DI


class OutletOutsideTables(Exception):
    """Raised in strict mode when the end of the outlet stretch lies outside the solved 2-D tables. It is a plain
    Exception, so the leg fitter treats it neither as an OutletInfeasible ('try the next leg count') nor as a
    ValueError ('n legs not drawable'), and strict mode refuses at the first leg count whose settled lead lies beyond a
    table. That is the conservative choice, because another leg count could move the lead width through the stretch
    length."""


def check_exp_domain(spec, w_last, w_lead):
    """Return the domain note of the outlet stretch's end (None inside the solved tables), and raise
    OutletOutsideTables when the spec refuses an end outside them."""
    dom = S.outlet_exp_domain(spec, w_last, w_lead)
    if dom and getattr(spec, 'exp_domain_refuse', False):
        raise OutletOutsideTables('the outlet stretch\'s end lies outside the solved 2-D tables (%s)' % dom)
    return dom


def check_plan_domain(sc2):
    """Raise OutletOutsideTables when the plan's stretch end (sc2.info['exp_domain']) lies outside the solved tables
    and the spec refuses such an end."""
    dom = sc2.info.get('exp_domain')
    if dom and getattr(sc2.spec, 'exp_domain_refuse', False):
        raise OutletOutsideTables('the outlet stretch\'s end lies outside the solved 2-D tables (%s)' % dom)


class OutletInfeasible(ValueError):
    """Raised when the requested outlet cannot close this design, because its lead would need a resistance, or a
    clearance, that the outlet's geometry cannot give. This is a property of the request rather than a drawing bug. The
    pipeline turns it into a named refusal (exit 5) before drawing, in every mode, and the leg fitter tries the next leg
    count instead of stopping."""


def _root(fn, lo, hi, what, **kw):
    """Run brentq on [lo, hi] after checking that the root is bracketed, so that a failure names what was being solved
    and no library text reaches the user. The result is identical to brentq's whenever the root is bracketed."""
    from scipy.optimize import brentq
    a, b = fn(lo), fn(hi)
    if a * b > 0.0:
        raise ValueError('%s: no solution between %g and %g um' % (what, lo, hi))
    return brentq(fn, lo, hi, **kw)


def head_top_x(sc):
    """Return the upstream extent of the drawn head in leg-1 local x (gap 1 at 0). Without a head (classic Gifford),
    this is the upstream edge of the first post; otherwise it is the most upstream point of any unit's ribbon (the arm
    tops, or row 1's out line)."""
    spec = sc.spec
    if sc.NM == 0:
        return -spec.Lam_L / 2.0 - spec.g_p / 2.0
    units = Gm.plan_head(sc)
    Gm.solve_transition(sc, units[-1])
    xmin = 0.0
    for u in units:
        Gm.exact_unit(spec, u)
        p = Gm.unit_path(spec, u)
        xmin = min(xmin, min(x for x, y in p.p) - u['w'] / 2.0)
    return xmin


def inlet_port(sc):
    """Return the center x and the radius of the inlet port in leg-1 local coordinates."""
    spec = sc.spec
    if spec.filter == 'dinh':
        return DI.PORT_X, DI.PORT_R
    _, inl = inlet_local(spec, head_top_x(sc), filt=bool(spec.filter))
    return inl['x_port_c'], spec.port_r_in


# ------------------------------------------------------------------------------------------------ helpers
def wc_at(sc, xp):
    """Return the center channel's width at path coordinate xp, interpolated linearly between the gap nodes."""
    return np.interp(xp, sc.x, sc.wc)


def band_poly(fx_lo, fx_hi, xa, xb, n):
    """Return the polygon between fx_lo(x) and fx_hi(x) for x in [xa, xb], sampled at n points (at least 2)."""
    xs = np.linspace(xa, xb, max(2, n))
    lo = [(x, fx_lo(x)) for x in xs]
    hi = [(x, fx_hi(x)) for x in xs]
    return Polygon(lo + hi[::-1])


def annulus_piece(cx, cy, th0, dth, r_in_fn, r_out_fn, ext=0.0, n=None, L=None):
    """Return the band between two radius functions of u in [0, 1] along the angle th0 + u dth (radians).

    At both ends the band is continued straight along the end tangents by `ext` um (u is extrapolated as -d/L and
    1 + d/L, with d = ext and L the center-path length of the arc). In other words, it overlaps the adjoining straight
    leg exactly on that leg's own width profile, so the union carries no step. Earlier versions extended the arc
    itself with u clamped, which left a step of slope x overlap (0.02 to 0.03 um) where a loop's width ramp meets its
    leg, flagged by the curvature validator."""
    # The chords are at most about 2 um long on the outer edge, because a loop's short width step (loop_ramp) must be
    # resolved. Chords of 0.25 degree are 16 um long on Xia's 3.6 mm loop and turned a 40 um step into a 2 um corner
    # (V2).
    n = n or max(16, int(abs(math.degrees(dth)) / 0.25), int(abs(dth) * abs(r_out_fn(0.5)) / 2.0))
    sgn = 1.0 if dth > 0 else -1.0
    th1 = th0 + dth
    ths = np.linspace(th0, th1, n + 1)
    us = (ths - th0) / dth

    def ring(rfn):
        pts = [(cx + rfn(u) * math.cos(t), cy + rfn(u) * math.sin(t)) for u, t in zip(us, ths)]
        if ext > 0 and L:
            t0 = (-sgn * math.sin(th0), sgn * math.cos(th0))
            t1 = (-sgn * math.sin(th1), sgn * math.cos(th1))
            r0, r1 = rfn(-ext / L), rfn(1.0 + ext / L)
            a = (cx + r0 * math.cos(th0) - ext * t0[0], cy + r0 * math.sin(th0) - ext * t0[1])
            b = (cx + r1 * math.cos(th1) + ext * t1[0], cy + r1 * math.sin(th1) + ext * t1[1])
            pts = [a] + pts + [b]
        return pts
    inner, outer = ring(r_in_fn), ring(r_out_fn)
    return Polygon(inner + outer[::-1]).buffer(0)


def loop_R(sc, k, inner, wfun, n_str=400, n_arc=4000):
    """Return the integral of r(w) ds along one side's centerline through the loop after row k. The three pieces (a
    straight half segment, the arc, and a straight half segment) are integrated separately with Simpson's rule, so that
    the jump of ds/dp at the arc ends is resolved exactly. wfun(p) gives the width at center-path coordinate p from
    node k."""
    from scipy.integrate import simpson
    spec, h = sc.spec, sc.spec.h
    R = S.R_of(sc, k)
    Lam, lg = spec.Lam_L, spec.l_g
    a, b = Lam / 2.0, Lam / 2.0 + math.pi * R
    PL = Lam + math.pi * R
    tot = 0.0
    for p0, p1, arc, n in ((0.0, a, False, n_str), (a, b, True, n_arc), (b, PL, False, n_str)):
        ps = np.linspace(p0, p1, n + 1)
        w = np.asarray(wfun(ps), float)
        rr_ = np.array([S.rr(ww, h) for ww in w])
        if arc:
            hw = wc_at(sc, sc.x[k] + ps) / 2.0
            yc = hw + lg + w / 2.0
            rho = (R - yc if inner else R + yc)
            ds = rho / R
            if getattr(spec, 'bend_model', 'stokes') != 'none':   # exact Stokes curved-duct correction (R2-4)
                isub = np.unique(np.linspace(0, len(ps) - 1, 41).astype(int))   # 41 exact values, smooth along the arc
                bf = np.array([S.bend_factor(spec, rho[i] - w[i] / 2.0, rho[i] + w[i] / 2.0) for i in isub])
                rr_ = rr_ * np.interp(ps, ps[isub], bf)
            if getattr(spec, 'bend_model', 'stokes') != 'none':   # the arc's two junctions (S.KAPPA_J)
                im = len(ps) // 2
                bfm = S.bend_factor(spec, rho[im] - w[im] / 2.0, rho[im] + w[im] / 2.0)
                tot -= min(S.kappa_j(spec, w[im]) * w[im], rho[im] * math.pi) * (bfm - 1.0) * S.rr(w[im], h)
        else:
            ds = 1.0
        tot += simpson(rr_ * ds, x=ps)
    return tot * S.UM


def loop_eq_width(sc, k):
    """Return the straight width whose r(w) over the loop row's center-path length gives the row's target R_s(k).
    This is the width that the loop rule stands in for, because the loop row itself has no straight side width."""
    return _root(lambda w: S.rr(w, sc.spec.h) * sc.ell[k] * S.UM - sc.Rs[k], 1.0, 5000.0,
                 'U-turn after row %d: the equivalent straight side width' % (k + 1), xtol=1e-10)


def loop_profile(sc, k, inner):
    """Return the side-channel width w_t through the gap-free loop after row k (spec v5 14.5) and the width profile
    as a function of the center path p in [0, PL] (p = 0 at node k, and PL = Lam + pi R at node k+1).

    The profile is smooth. It goes from the width drawn on segment k-1 to w_t over the first quarter of the arc (over
    loop_ramp when that is set, lengthened where a short step would bend the wall too tightly), stays at w_t, and
    returns to the width drawn on segment k+1 over the same length at the arc's end. w_t is solved so that the
    integral of r(w) along this side's own centerline equals the row's exact target R_s(k), and under Xia's pull-in in
    'width' mode, loop_pull_in is then subtracted from it."""
    spec = sc.spec
    R = S.R_of(sc, k)
    Lam = spec.Lam_L
    PL = Lam + math.pi * R
    # The loop's ends meet the widths drawn on the legs, which Xia's leg pull-in narrows after a loop.
    wa, wb = drawn_side_width(sc, k - 1), drawn_side_width(sc, k + 1)
    target = sc.Rs[k]
    pull = getattr(spec, 'loop_pull_in', 0.0) or 0.0
    if pull > 0.0 and getattr(spec, 'loop_pull_mode', 'width') == 'resistance':
        # The loop rule closes on the narrowed width, so both sides' R scale by the narrowed leg's factor and the
        # solved widths are drawn as they are, with no further subtraction.
        weq = loop_eq_width(sc, k)
        target = target * S.rr(weq - pull, spec.h) / S.rr(weq, spec.h)
        pull = 0.0

    # The width ramps lie entirely on the arc (its first and last quarters), with zero slope at both arc ends,
    # because a ramp that straddled the arc end left a kink where the arc met the straight leg. On the arc the
    # boundary slope is dw/dp x R/r, whereas on the straight it is dw/dp, and the kink read as a 3.5 um radius on the
    # tight Dinh loops.
    a0, qa = Lam / 2.0, math.pi * R / 4.0
    ramp = getattr(spec, 'loop_ramp', None)
    # Xia 2016 Fig. 1e(ii) draws constant loop widths with short smooth steps at the ends.
    if ramp:
        qa = min(qa, float(ramp))

    def wfun(p, wt):
        p = np.asarray(p, float)
        t1 = np.clip((p - a0) / qa, 0, 1)
        t2 = np.clip((p - (PL - a0 - qa)) / qa, 0, 1)
        w = wa + (wt - wa) * Gm.smoothstep(t1)
        return np.where(p > PL - a0 - qa, wt + (wb - wt) * Gm.smoothstep(t2), w)
    def solve_wt():
        wt_ = _root(lambda wt: loop_R(sc, k, inner, lambda p: wfun(p, wt), 200, 1500) - target, 3.0, 3000.0,
                    'U-turn after row %d: the %s side channel\'s width' % (k + 1, 'inner' if inner else 'outer'), xtol=1e-12)
        # One secant step with the fine integral refines the root, which is enough because the function is smooth.
        f0 = loop_R(sc, k, inner, lambda p: wfun(p, wt_)) - target
        d = 1e-4
        f1 = loop_R(sc, k, inner, lambda p: wfun(p, wt_ + d)) - target
        return wt_ - f0 * d / (f1 - f0)
    wt = solve_wt()
    if ramp:
        # A short step must stay smooth. A smoothstep of height dw over a center-path length q bends the side
        # channel's wall with a radius of about q^2 / (6 dw) (less on the inner wall), so we keep q >= sqrt(90 dw), that
        # is, a radius of at least about 15 um. A 40 um ramp under a 30 um step at Xia's loop 6 left a 1.7 um radius
        # (V2).
        for _ in range(3):
            dw = max(abs(wt - pull - wa), abs(wt - pull - wb))
            # On the inner side the step runs along walls closer to the loop center, so its length there is
            # shorter by the factor r_wall / R (380 / 663 on Xia's loop 1), and the center-path ramp is lengthened
            # by R / r_wall.
            r_wall = R - (sc.wc[k] / 2.0 + spec.l_g + max(wa, wb, wt)) if inner else R
            q_need = min(math.sqrt(90.0 * dw) * R / max(r_wall, 1.0), math.pi * R / 4.0)
            if qa >= q_need - 1e-9:
                break
            qa = q_need
            wt = solve_wt()
    if pull > 0.0:                         # Xia 2016 p. 4: "subtract 1-5 um from the predicted side channel widths"
        wt = wt - pull
    return wt, (lambda p, wt=wt: wfun(p, wt))


# ------------------------------------------------------------------------------------------------ plan
def head_clear_x(sc, lim):
    """Return the last path x where the head's half-width exceeds lim, with the arms counted as infinitely wide."""
    spec = sc.spec
    x_clear = 0.0
    for u in Gm.plan_head(sc):
        if u['kind'] == 'A':
            x_clear = max(x_clear, u['X'] + u['ell'])
            continue
        R = u['R'] if u['R'] is not None else spec.yJ(spec.GS) + spec.s
        reach = R + spec.s / 2.0 + spec.GS / 2.0
        if reach > lim:
            x_clear = max(x_clear, u['X'] + u['ell'])
    return x_clear


def plan_layout(sc, kind='linear', n_legs=1, R_turn=None, near_clear=200.0, Rc_lead_len=600.0):
    """Plan the layout, growing the loop radius until every inner loop channel is at least F_min + 1 um wide and the
    solid disc inside each loop's innermost channel is at least F_min + 1 um across."""
    R = R_turn
    for _ in range(40):
        sc2, legs = _plan_layout(sc, kind, n_legs, R, near_clear, Rc_lead_len)
        ok = True
        for k in sc2.info['turns']:
            wi, _ = loop_profile(sc2, k, inner=True)
            r_disc = S.R_of(sc2, k) - (max(sc2.wc[k], sc2.wc[k + 1]) / 2 + sc2.spec.l_g + wi)
            if wi < sc2.spec.F_min + 1.0 or 2 * r_disc < sc2.spec.F_min + 1.0:
                ok = False
        if ok:
            check_plan_domain(sc2)                            # only the plan that the radius search accepts
            return sc2, legs
        R = (sc2.info['leg_d'] if kind == 'paperclip' else sc2.info['R_turn']) * 1.1
    raise RuntimeError('could not find a loop radius meeting the floor')


def _plan_layout(sc, kind='linear', n_legs=1, R_turn=None, near_clear=200.0, Rc_lead_len=600.0):
    """Choose the rows after which loops are placed and return (sc_layout, legs), where each leg is a dict with these
    entries:

        rows        (ka, kb), the leg's first and last rows
        p0          the path start (None for leg 1)
        p1          the path end
        Y           the axis y
        dirn        +1 or -1
        turn_after  the row after which the leg's loop is placed, or None

    A serpentine has an odd n_legs, so the outlet leaves at +x, away from the head. Its legs 2 to n lie between x_near
    and x_far beside leg 1's linear part, with x_near clearing the head."""
    spec = sc.spec
    N, NM = sc.N, sc.NM
    Lam = spec.Lam_L
    hw_max = float(np.nanmax(sc.wc[NM:N - 1] / 2.0 + spec.l_g + sc.ws[NM:N - 1]))
    if kind == 'paperclip' and n_legs > 1:
        return _plan_paperclip(sc, n_legs, R_turn)        # R_turn here is the leg spacing d (None for automatic)
    if R_turn is None:
        R_turn = hw_max + spec.F_min / 2.0 + 20.0
    turns = []
    if kind == 'serpentine' and n_legs > 1:
        if n_legs % 2 == 0:
            raise ValueError('serpentine: use an odd number of legs (the outlet must leave away from the head)')
        lim = 2 * R_turn - hw_max - spec.F_min
        x_near = max(head_clear_x(sc, lim), sc.x[NM]) + R_turn + hw_max + near_clear
        x_end = sc.x[N - 1]
        x_far = (x_end + (n_legs - 1) * x_near) / n_legs
        if x_far <= x_near + 4 * Lam:
            raise ValueError('device too short for %d legs' % n_legs)
        Lj = x_far - x_near
        kk = int(np.searchsorted(sc.x, x_far - Lam / 2.0, side='right')) - 1
        turns.append(kk)
        for j in range(2, n_legs):
            start = sc.x[kk + 1] - Lam / 2.0
            kk = int(np.searchsorted(sc.x, start + Lj - Lam / 2.0, side='right')) - 1
            turns.append(kk)
        if turns[-1] >= N - 3:
            raise ValueError('loop placement ran past the end')
    if turns:
        # The last leg ends beside the far-end loop of legs 1 and 2, so its outlet stretch is lengthened for the side
        # leads to start beyond that loop's extent (R_turn + half-width), with room for their quarter arcs.
        import copy
        sc = copy.copy(sc)
        sc.spec = copy.copy(spec)
        # The stretch also carries the extension that build.plan_and_draw adds when the outlet would meet the inlet
        # end (the clearance guarantee).
        sc.spec.L_out = R_turn + hw_max + 400.0 + float(getattr(spec, 'L_out_extra', 0.0) or 0.0)
        spec = sc.spec
    elif getattr(spec, 'L_out_extra', 0.0):
        # A layout without turns lengthens its own stretch when the outlet meets the inlet (code review 4, F7).
        import copy
        sc = copy.copy(sc)
        sc.spec = copy.copy(spec)
        sc.spec.L_out = spec.L_out + float(spec.L_out_extra)
        spec = sc.spec
    if flare_on(spec):
        sc = _flare_fixed_point(sc, turns, R_turn, kind, spec.L_out)
        spec = sc.spec
    sc_c = S.reclose(sc, turns, R_turn, spec.R_ext_c)
    Rc_lead, Lc, fr = outlet_Rc_lead(sc_c)
    dom = S.outlet_exp_domain(sc_c.spec, sc_c.ws[sc_c.N - 2], float(fr['w'] if 'w' in fr else fr['wl']))
    sc2 = S.reclose(sc, turns, R_turn, Rc_lead + spec.R_ext_c)
    sc2.info['exp_domain'] = dom                              # checked in plan_layout on the accepted plan
    sc2.info['R_turn'] = R_turn
    sc2.info['Rc_lead_len'] = Lc
    sc2.info['Rc_lead'] = Rc_lead
    sc2.info['layout'] = kind
    if flare_on(spec):
        sc2.info.update(L_out_base=sc.info['L_out_base'], w_lead=sc.info['w_lead'], flare_len=sc.info['flare_len'],
                        flare_kind=sc.info['flare_kind'])
    if getattr(spec, 'outlet', 'separate') == 'combined':
        sc2.info['ring'] = fr
    return sc2, make_legs(sc2)


def head_halfwidth(sc):
    """Return the largest lateral extent of the head on one side (um), counting the units, the arms, the inlet stem
    or filter, and the port."""
    spec = sc.spec
    hw = max(spec.wc0 / 2.0 + spec.l_g, spec.stem_half if spec.filter else spec.wc0 / 2.0, inlet_port(sc)[1])
    for u in Gm.plan_head(sc):
        if u['kind'] == 'A':
            hw = max(hw, u['R'] + (2 * u['n'] - 1) * spec.s_a + u['w'] / 2.0)
        elif u['R'] is not None:
            hw = max(hw, u['R'] + spec.s / 2.0 + u['w'] / 2.0)
    if spec.filter == 'dinh':
        hw = max(hw, DI.BULB_HALF, DI.PORT_R)
    return hw


def paperclip_Y(j, d):
    """Return the axis of leg j (0-based) in a paperclip, in the sequence 0, +d, -d, +2d, -2d, ..."""
    if j == 0:
        return 0.0
    m = (j + 1) // 2
    return m * d if j % 2 == 1 else -m * d


def _plan_paperclip(sc, n_legs, d=None, clear=150.0):
    """Plan the legs of a paperclip. The leg spacing d clears the widest leg and the head and inlet, and loop j (after
    leg j) has radius |Y(j+1) - Y(j)| / 2.

    With spec.leg_order 'nested' (the default), the legs nest (Y = 0, +d, -d, +2d, ...), every loop turns the same
    way, legs 2 to n run between x_lo (the inlet port's x, the center of loop 2) and x_hi, and leg 1 runs from gap 1 to
    x_hi. With 'path' (GIF4-1, reading B), the legs lie side by side in path order (Y = j d, or a reproduction.leg_Y
    stepping one way), the loops alternate their turning direction, and the inlet does not sit inside loop 2. Either
    way, the last leg ends at the far end of its own run, and its outlet stretch leaves beyond the outermost loop at
    that end."""
    import copy
    spec = sc.spec
    N, NM, Lam = sc.N, sc.NM, spec.Lam_L
    hw_max = float(np.nanmax(sc.wc[NM:N - 1] / 2.0 + spec.l_g + sc.ws[NM:N - 1]))
    hw_head = head_halfwidth(sc)
    if d is None:
        d = max(2.0 * hw_max, hw_head + hw_max) + clear
    if getattr(spec, 'leg_Y', None) is not None:          # measured leg axes (Xia 2016 Fig. 3a)
        Y = [float(v) for v in spec.leg_Y]
        if len(Y) != n_legs or Y[0] != 0.0:
            raise ValueError('leg_Y must list %d axes starting at 0' % n_legs)
        if getattr(spec, 'leg_order', 'nested') == 'path':
            # The legs lie side by side in path order, so the axes step one way (0, +, ++, ... or 0, -, --, ...), and
            # the loops then alternate their turning direction, each loop's geometry and closure taking its own
            # (GIF4-1, B).
            st_ = [Y[j] - Y[j - 1] for j in range(1, n_legs)]
            if not (all(s_ > 0 for s_ in st_) or all(s_ < 0 for s_ in st_)):
                raise ValueError('leg_Y with layout.order: path must step one way (0, +, ++, ...)')
        else:
            for j in range(1, n_legs):                     # the alternation of paperclip_Y, so the loops all turn alike
                if (Y[j] - Y[j - 1]) * (1 if j % 2 == 1 else -1) <= 0:
                    raise ValueError('leg_Y must alternate as 0, +, -, +, ... with growing spread')
        d = None
    elif getattr(spec, 'leg_order', 'nested') == 'path':
        Y = [j * d for j in range(n_legs)]                 # legs side by side in path order
    else:
        Y = [paperclip_Y(j, d) for j in range(n_legs)]
    # The inlet port sits at the center of loop 2 (x_lo), or pc_inlet_offset downstream of it (Xia Fig. 3a).
    x_lo = inlet_port(sc)[0] - (getattr(spec, 'pc_inlet_offset', 0.0) or 0.0)
    x_end = sc.x[N - 1]                                  # pre-layout path of the last gap (no loops yet)
    dshort = getattr(spec, 'pc_last_short', 0.0) or 0.0
    x_hi = (x_end + dshort + (n_legs - 1) * x_lo) / n_legs
    Lj = x_hi - x_lo
    if Lj < 10 * Lam or x_hi <= sc.x[NM] + Lam:
        raise ValueError('device too short for %d paperclip legs' % n_legs)
    turns, Rt = [], {}
    kk = int(np.searchsorted(sc.x, x_hi - Lam / 2.0, side='right')) - 1
    turns.append(kk)
    for j in range(1, n_legs - 1):
        start = sc.x[kk + 1] - Lam / 2.0
        kk = int(np.searchsorted(sc.x, start + Lj - Lam / 2.0, side='right')) - 1
        turns.append(kk)
    if turns[-1] >= N - 3:
        raise ValueError('loop placement ran past the end')
    for j, k in enumerate(turns):
        Rt[k] = abs(Y[j + 1] - Y[j]) / 2.0
    # The outlet stretch reaches beyond the outermost loop at the last leg's end, with room for the outlet leads.
    end_hi = (n_legs - 1) % 2 == 0                       # last leg flows +x (ends at x_hi)
    R_end = max([Rt[k] for j, k in enumerate(turns) if (j % 2 == 0) == end_hi] or [0.0])
    sc = copy.copy(sc)
    sc.spec = copy.copy(spec)
    # The stretch also carries the extension that build.plan_and_draw adds when the outlet would meet the inlet end,
    # as in a paperclip whose last leg ends at the inlet end with no loop there to clear, for example two nested legs
    # (the clearance guarantee).
    sc.spec.L_out = R_end + hw_max + clear + 200.0 + dshort + float(getattr(spec, 'L_out_extra', 0.0) or 0.0)
    if getattr(spec, 'outlet', 'separate') == 'combined' and getattr(spec, 'outlet_style', 'teardrop') == 'ring':
        # The ring is centered on the last leg's axis beyond the stretch end. It must clear every loop at that end by
        # `clear` (the loops are semicircles beyond their turn rows, which sit at x_hi, close to the last gap), and its
        # top lies pc_out_ext beyond the loop centers when that is given (3.44 mm in Xia Fig. 3a). The lead width from
        # a preliminary outlet plan sets the frame when it is wider than the side channel.
        try:
            _, _, fr0 = outlet_Rc_lead(S.reclose(sc, turns, Rt, sc.spec.R_ext_c))
            w_est = float(fr0.get('w', 0.0))
        except ValueError:
            w_est = 0.0
        g0 = ring2_geom(sc.spec, sc.wc[-1], sc.ws[N - 2], w_lead=w_est if w_est > sc.ws[N - 2] else None)
        w_g = max(sc.ws[N - 2], w_est)
        need = 0.0
        for j, k in enumerate(turns):
            if (j % 2 == 0) == end_hi:
                Yc = 0.5 * (Y[j] + Y[j + 1]) - Y[-1]
                hw_k = sc.wc[k] / 2.0 + spec.l_g + sc.ws[k]
                reach = spec.ring_R + w_g / 2.0 + Rt[k] + hw_k + clear
                if reach > abs(Yc):
                    need = max(need, math.sqrt(reach ** 2 - Yc ** 2) - g0['dx_rc'])
        L_ring = need + dshort
        if getattr(spec, 'pc_out_ext', None) is not None:
            L_want = spec.pc_out_ext - (g0['dx_rc'] - spec.ring_R) + dshort
            if L_want < L_ring - 1e-6:
                raise ValueError('ring outlet: pc_out_ext %.0f um collides with the loops (needs %.0f)' % (
                    spec.pc_out_ext, L_ring - dshort + g0['dx_rc'] - spec.ring_R))
            L_ring = L_want
        sc.spec.L_out = max(L_ring, hw_max + clear)
    if flare_on(spec):
        sc = _flare_fixed_point(sc, turns, Rt, 'paperclip', sc.spec.L_out)
    sc_c = S.reclose(sc, turns, Rt, sc.spec.R_ext_c)
    Rc_lead, Lc, fr = outlet_Rc_lead(sc_c)
    dom = S.outlet_exp_domain(sc_c.spec, sc_c.ws[sc_c.N - 2], float(fr['w'] if 'w' in fr else fr['wl']))
    sc2 = S.reclose(sc, turns, Rt, Rc_lead + sc.spec.R_ext_c)
    sc2.info['exp_domain'] = dom                              # checked in plan_layout on the accepted plan
    sc2.info['Rc_lead_len'] = Lc
    sc2.info['Rc_lead'] = Rc_lead
    if flare_on(spec):
        sc2.info.update(L_out_base=sc.info['L_out_base'], w_lead=sc.info['w_lead'], flare_len=sc.info['flare_len'],
                        flare_kind=sc.info['flare_kind'])
    if getattr(spec, 'outlet', 'separate') == 'combined':
        sc2.info['ring'] = fr
    sc2.info['leg_Y'] = Y
    sc2.info['leg_d'] = d
    sc2.info['layout'] = 'paperclip'
    return sc2, make_legs(sc2)


def make_legs(sc2):
    """Return the legs of a closed schedule whose info carries the turns and R_turn (or the per-turn radii and the leg
    axes)."""
    spec = sc2.spec
    N, Lam = sc2.N, spec.Lam_L
    turns = list(sc2.info.get('turns', []))
    R_turn = sc2.info['R_turn']
    Ylist = sc2.info.get('leg_Y')
    legs = []
    starts = [0] + [k + 1 for k in turns]
    ends = turns + [N - 1]
    for j in range(len(ends)):
        ka, kb = starts[j], ends[j]
        p0 = None if j == 0 else sc2.x[ka] - Lam / 2.0
        p1 = sc2.x[kb] + Lam / 2.0 if j < len(ends) - 1 else sc2.x[N - 1] + spec.L_out
        Yj = Ylist[j] if Ylist is not None else 2 * R_turn * j
        legs.append(dict(j=j, rows=(ka, kb), p0=p0, p1=p1, Y=Yj, dirn=(1 if j % 2 == 0 else -1),
                         turn_after=(kb if j < len(ends) - 1 else None)))
    return legs


# ------------------------------------------------------------------------------------------------ leg pieces
def transition_params(sc):
    """Return (paper, ws_rect, kL), where paper says whether the papers' transition is drawn. That transition draws the
    side channel of the rectilinear rows ws_rect wide (GS by default) and narrows it to w_s(kL) at the first stadium
    post, where kL is the first linear ('L') row after the head."""
    spec = sc.spec
    kL = next((k for k in range(sc.NM, sc.N) if sc.kind[k] == 'L'), sc.NM)
    paper = (getattr(spec, 'transition', 'exact') == 'paper') and sc.NM > 0
    wR = spec.ws_rect if getattr(spec, 'ws_rect', None) is not None else spec.GS
    return paper, wR, kL


def pulled_rows(sc):
    """Return the rows whose straight side channels Xia's pull-in narrows (loop_pull_in_leg), which are every linear
    row of every leg that follows a loop. The outlet stretch is not narrowed, so the outlet closure stays exact."""
    spec = sc.spec
    if not getattr(spec, 'loop_pull_in_leg', False) or not (getattr(spec, 'loop_pull_in', 0.0) or 0.0):
        return set()
    turns = sorted(sc.info.get('turns', []))
    if not turns:
        return set()
    return set(range(turns[0] + 1, sc.N - 1)) - set(turns)


def drawn_side_width(sc, k, paper=None, wR=None, pulled=None):
    """Return the side width drawn on straight segment k. This is the closure width, except in the papers' transition
    (ws_rect on the rectilinear rows) and under Xia's leg pull-in (the closure width minus loop_pull_in on the legs
    after a loop)."""
    if paper is None:
        paper, wR, _ = transition_params(sc)
    if pulled is None:
        pulled = pulled_rows(sc)
    if paper and sc.kind[k] == 'R':
        return wR
    if k in pulled:
        return sc.ws[k] - sc.spec.loop_pull_in
    return sc.ws[k]


# ---- straight side-channel walls: resistance-exact profile (opencif fix F1) ------------------------------------
# The legacy rule drew each straight wall piecewise linear through (node, mean of the two neighboring widths) and
# (segment midpoint, w_s(k)), assuming that the segment's resistance is then r(w_s(k)) ell(k) to second order. That
# assumption fails wherever the width or the pitch changes fast (the sealed start, the seam, and the last rows), and
# Gifford row 2 drew -14.8% in f. Now the node widths stay the means (for continuity), and each segment's midpoint
# width is solved so that the drawn piecewise-linear profile has exactly the closure's resistance r(w_s(k)) ell(k)
# (6-point Gauss-Legendre per half, exact to about 1e-12 for these smooth profiles). WALL_RULE = 'legacy' restores
# the old drawing for the golden comparisons.
WALL_RULE = 'exact'
SIDE_PROFILES = {}            # (leg index, side 0/1) -> [(path x, intended side width)], filled by leg_linear_local
_GLX, _GLW = np.polynomial.legendre.leggauss(6)


def lin_profile_R(a, b, L, h):
    """Return the integral of r(w) ds (SI per m x um) over a straight run of length L (um) whose width goes linearly
    from a to b."""
    t = 0.5 * (_GLX + 1.0)
    return 0.5 * L * sum(wi * S.rr(a + (b - a) * ti, h) for ti, wi in zip(t, _GLW))


def seg_profile_R(n0, m, n1, ell, h):
    """Return the integral of r(w) ds over a segment of length ell whose width goes linearly from n0 to m over its
    first half and from m to n1 over its second half."""
    return lin_profile_R(n0, m, 0.5 * ell, h) + lin_profile_R(m, n1, 0.5 * ell, h)


def solve_mid_width(n0, n1, w_target, ell, h):
    """Return the midpoint width m for which the profile n0 -> m -> n1 over ell has the resistance of a constant width
    w_target. Under the legacy wall rule, or when both ends are already at w_target, m is w_target itself."""
    if WALL_RULE == 'legacy' or (abs(n0 - w_target) < 1e-12 and abs(n1 - w_target) < 1e-12):
        return w_target
    from scipy.optimize import brentq
    target = S.rr(w_target, h) * ell
    fn = lambda m: seg_profile_R(n0, m, n1, ell, h) - target
    lo, hi = 0.5 * w_target, 2.0 * w_target
    while fn(lo) < 0:
        lo *= 0.5
    while fn(hi) > 0:
        hi *= 2.0
    return brentq(fn, lo, hi, xtol=1e-12, rtol=1e-14)


def leg_linear_local(sc, leg, side_w_start=(None, None), side_w_end=(None, None)):
    """Draw the trunk, both side channels, and both post rows of one leg in local coordinates (x = path - p0, or the
    path itself on leg 1). side_w_start and side_w_end give (w_A, w_B) for the half loop segments at the leg's ends
    (None where there is no loop)."""
    spec = sc.spec
    Lam, lg = spec.Lam_L, spec.l_g
    N, NM = sc.N, sc.NM
    ka, kb = leg['rows']
    k_first = max(ka, NM)
    p0 = leg['p0'] if leg['p0'] is not None else 0.0
    p1 = leg['p1']
    off = p0 if leg['p0'] is not None else 0.0

    def yin(xp):
        return wc_at(sc, xp) / 2.0 + lg

    pieces = []
    paper, wR, kL = transition_params(sc)

    wR_open = spec.GS                                   # width of a rectilinear opening

    pulled = pulled_rows(sc)

    def wdraw(k):
        return drawn_side_width(sc, k, paper, wR, pulled)
    # ---- the trunk over the leg (on leg 1 from gap 1 on, because the stem is drawn with the head)
    pa = p0 if leg['p0'] is not None else 0.0
    xs = np.concatenate([[pa], sc.x[(sc.x > pa) & (sc.x < p1)], [p1]])
    hw = wc_at(sc, xs) / 2.0
    pieces.append(Polygon([(x - off, -h) for x, h in zip(xs, hw)] + [(x - off, h) for x, h in zip(xs[::-1], hw[::-1])]))
    # ---- the side channels, one side at a time
    for sgn, iw in ((+1, 0), (-1, 1)):
        pts_in, pts_out = [], []
        ramp = 1.0

        prof = []

        def add(xp, w):
            y = yin(xp)
            pts_in.append((xp - off, sgn * y))
            pts_out.append((xp - off, sgn * (y + w)))
            prof.append((float(xp), float(w)))
        if leg['p0'] is not None:                         # half loop segment at the leg start (width profile)
            fw = side_w_start[iw]
            kprev = ka - 1
            for xp in np.linspace(p0, sc.x[k_first], 60):
                add(xp, float(fw(xp - sc.x[kprev])))
        elif NM == 0:                                     # no head (classic Gifford), so the side channel starts closed
            add(sc.x[0] - Lam / 2.0, sc.ws[0])            # behind the first post (its end is capped round below)
        # On the straight segments the outer wall is piecewise linear through (node, mean of the two neighboring
        # widths) and (segment midpoint, m_k), with m_k solved so that the segment's resistance is exactly
        # r(w_s(k)) ell(k) (fix F1, see solve_mid_width). Constant widths joined by 2 um ramps left 1 um jogs per row
        # on Xia's fast-widening channels, which drew 1.9 um radii (V2).
        k_last_seg = kb if leg['turn_after'] is not None else N - 1
        segs = list(range(k_first, k_last_seg))
        # The first and last nodes of the leg carry their segment's width, and the inner nodes carry the mean of their
        # two neighbors. A sealed start's first segment is drawn flat at w_s(1), its end node included, so that the
        # solved midpoints never dip below the narrowest width that the design asks for (the process floor at a floor
        # start).
        nodes = [wdraw(segs[0])] if segs else []
        for j_ in range(1, len(segs)):
            nodes.append(0.5 * (wdraw(segs[j_ - 1]) + wdraw(segs[j_])))
        if segs:
            nodes.append(wdraw(segs[-1]))
        if segs and k_first == 0 and NM == 0 and len(segs) > 1:
            nodes[1] = wdraw(segs[0])
        for j, k in enumerate(segs):
            w = wdraw(k)
            if paper and k == kL - 1 and j > 0 and sc.kind[k] == 'R':
                # In the papers' transition, the side channel beside the wall blocks is ws_rect wide and narrows to
                # w_s at the first stadium post (Xia 2016 Fig. 1c(ii); Dinh 2024 Fig. 2B, second panel, where it
                # narrows from the end of the last block to the center of the first pill). A smooth ramp runs from
                # the start of the last opening (x - GS/2) to the first post's center, and w_s(kL) follows to the next
                # node.
                xa_ = sc.x[k] - wR_open / 2.0
                xb_ = sc.x[k + 1] - spec.G / 2.0 - spec.g_p / 2.0
                wn = sc.ws[k + 1]
                for t in np.linspace(0.0, 1.0, 21):
                    add(xa_ + t * (xb_ - xa_), wR + (wn - wR) * Gm.smoothstep(t))
                continue
            if paper and k == kL and j > 0:
                add(sc.x[k], w)
                add(0.5 * (sc.x[k] + sc.x[k + 1]), w)
                continue
            n0, n1 = (nodes[j], nodes[j + 1]) if WALL_RULE != 'legacy' else \
                (w if j == 0 else 0.5 * (wdraw(k - 1) + w), None)
            # The midpoint is resistance-exact (F1), whereas the papers' transition rows (reproductions only) keep the
            # legacy drawing.
            m_ = w if (paper and NM <= k <= kL) or WALL_RULE == 'legacy' else solve_mid_width(n0, n1, w, sc.ell[k], spec.h)
            add(sc.x[k], n0)
            add(0.5 * (sc.x[k] + sc.x[k + 1]), m_)
        if segs:
            add(sc.x[segs[-1] + 1], wdraw(segs[-1]))
        if leg['turn_after'] is not None:                 # half loop segment at the leg end (width profile)
            fw = side_w_end[iw]
            for xp in np.linspace(sc.x[kb], p1, 60):
                add(xp, float(fw(xp - sc.x[kb])))
        else:                                             # outlet stretch
            w_end = wdraw(segs[-1]) if segs else sc.ws[N - 2]
            if abs(w_end - sc.ws[N - 2]) > 0.05:
                # Xia's leg pull-in stops at the last gap, so that the outlet stretch stays exact for the outlet
                # closure, and a smoothstep over at most 30 um returns the pulled width to w_s. A 1 um jog drew a
                # 0.95 um radius (V2).
                Lr_ = min(30.0, 0.3 * (p1 - sc.x[N - 1]))
                for t in np.linspace(0.0, 1.0, 16)[1:]:
                    add(sc.x[N - 1] + t * Lr_, w_end + (sc.ws[N - 2] - w_end) * Gm.smoothstep(t))
            else:
                add(sc.x[N - 1] + ramp, sc.ws[N - 2])
            if flare_on(spec) and sc.info.get('flare_len', 0.0) > 0.0:
                # After the straight run, the outer wall diverges linearly to the lead width, whereas the inner wall
                # (pts_in at yin) stays straight (Dinh 2024 Fig. 2B, fourth panel). Tangent arcs round both kinks.
                for s_, w_ in flare_profile(spec, sc.ws[N - 2], sc.info['w_lead'], spec.L_out)[1:]:
                    if s_ > ramp + 1e-9:
                        add(sc.x[N - 1] + s_, w_)
            else:
                add(p1, sc.ws[N - 2])
        # The intended outer-wall profile of this side, as (path x, width) drawn here before any union or smoothing,
        # is the reference with which the drawing check compares the measured polygon (opencif, review R2-1).
        SIDE_PROFILES[(leg['j'], iw)] = prof
        pieces.append(Polygon(pts_in + pts_out[::-1]).buffer(0))
        if leg['p0'] is None and NM == 0:                 # round cap on the closed upstream end
            w0 = sc.ws[0]
            pieces.append(Point(sc.x[0] - Lam / 2.0 - off, sgn * (yin(sc.x[0] - Lam / 2.0) + w0 / 2.0)).buffer(
                w0 / 2.0, quad_segs=32))
    # ---- the post rows (the band minus the solid islands) on both sides. Between two 'L' gaps there is a stadium
    # post (g_p along the flow, l_g across). Any interval touching an 'R' opening holds a wall block, a rounded
    # rectangle filling the interval between the two openings (Dinh 2024 Fig. 2B bottom; Xia 2016 Fig. 1c(ii)). The
    # leg's upstream end has a stadium before an 'L' gap and a block before an 'R' opening.
    kind = sc.kind
    wL = sc.wL
    if kind[k_first] == 'R':
        Lb_up = 2.0 * spec.s - spec.GS
        xa = sc.x[k_first] - wL[k_first] / 2.0 - Lb_up / 2.0
        if k_first == NM and NM > 0 and kind[NM - 1] != 'A':
            # The band overlaps the side channel by EPS above its blocks, but on leg 1 nothing lies above the first
            # block upstream of the transition unit's drop (x0 = X + s + s/2, where the hairpin's quarter arc ends,
            # tangent to y_in). A band starting upstream of x0 left a 1 um dead-end sliver under the arc, which the
            # fillet closing turned into a pocket of about 10 um2 against the arc's outer edge. In the wbc example,
            # where ell < 2.5 s, the transition unit read +7e-4 long and the seam row moved 1%. The band therefore
            # starts under the drop.
            x0_drop = sc.x[NM - 1] + 1.5 * spec.s
            xa = min(max(xa, x0_drop + 0.5), sc.x[k_first] - wL[k_first] / 2.0 - 1.0)
    else:
        xa = sc.x[k_first] - Lam / 2.0
    xb = sc.x[kb] + Lam / 2.0
    xx = np.linspace(xa, xb, max(2, int((xb - xa) / 25.0)))
    yc = wc_at(sc, xx) / 2.0

    def block(x0, x1, sgn):
        """Return a wall block whose two long sides lie exactly on the band edges (the trunk edge w_c(x)/2 and the side
        channel's inner edge w_c(x)/2 + l_g, both linear within one row interval) and whose four corners are rounded
        by an opening with radius rb, so that it meets the trunk and the side channel tangentially. A block
        overlapping the band edges by EPS left its rounded corners crossing the trunk edge at an angle, with a 2 um
        chord radius."""
        rb = min(10.0, 0.5 * (x1 - x0) - 0.5, lg / 2.0 - 0.5)
        if x1 > xb and rb > xb - x0:
            # The leg's last block runs on past the leg's end xb, where the band stops and the U-turn's solid wall
            # takes over. Its upstream corners' arcs must end on the band edges before xb, because an arc cut at xb
            # left a 0.2 um step against the turn's wall (a corner radius of about 0). Only a leg whose last opening
            # ends closer to xb than rb changes (composition sweep, case 12).
            rb = max(xb - x0 - 0.5, 0.0)
        y0a, y0b = wc_at(sc, x0) / 2.0, wc_at(sc, x1) / 2.0
        if sgn > 0:
            pts = [(x0, y0a), (x1, y0b), (x1, y0b + lg), (x0, y0a + lg)]
        else:
            pts = [(x0, -y0a - lg), (x1, -y0b - lg), (x1, -y0b), (x0, -y0a)]
        b = Polygon([(x - off, y) for x, y in pts])
        return b.buffer(-rb, quad_segs=16).buffer(rb, quad_segs=16)
    for sgn in (+1, -1):
        lower = [(x - off, sgn * (y - Gm.EPS)) for x, y in zip(xx, yc)]
        upper = [(x - off, sgn * (y + lg + Gm.EPS)) for x, y in zip(xx, yc)]
        band = Polygon(lower + upper[::-1]).buffer(0)
        posts = []
        for k in range(k_first, kb + 1):
            if kind[k] == 'L' and (k + 1 >= N or kind[k + 1] == 'L' or k == kb):
                xc = sc.x[k] + Lam / 2.0
                ycen = sgn * (wc_at(sc, xc) / 2.0 + lg / 2.0)
                posts.append(Polygon(S.C.stadium_post(xc - off, ycen, spec.g_p, lg, 24)))
            elif paper and kind[k] == 'R' and k + 1 < N and kind[k + 1] == 'L':
                xc = 0.5 * (sc.x[k] + wL[k] / 2.0 + sc.x[k + 1] - wL[k + 1] / 2.0)  # first pill after the last opening
                ycen = sgn * (wc_at(sc, xc) / 2.0 + lg / 2.0)
                posts.append(Polygon(S.C.stadium_post(xc - off, ycen, spec.g_p, lg, 24)))
            else:
                posts.append(block(sc.x[k] + wL[k] / 2.0, sc.x[k + 1] - wL[k + 1] / 2.0, sgn))
        if kind[k_first] == 'R':
            x1 = sc.x[k_first] - wL[k_first] / 2.0
            posts.append(block(x1 - (2.0 * spec.s - spec.GS), x1, sgn))
        else:
            xc = sc.x[k_first] - Lam / 2.0
            ycen = sgn * (wc_at(sc, xc) / 2.0 + lg / 2.0)
            posts.append(Polygon(S.C.stadium_post(xc - off, ycen, spec.g_p, lg, 24)))
        pieces.append(band.difference(unary_union(posts)))
    return unary_union(pieces)


def to_global(g, leg, gx0):
    """Map a leg's local geometry to global coordinates, as (gx0 + x, Y + y) on odd legs (dirn +1) and as
    (gx0 - x, Y - y) on even legs."""
    if leg['dirn'] > 0:
        return affinity.translate(g, gx0, leg['Y'])
    g2 = affinity.rotate(g, 180.0, origin=(0, 0))
    return affinity.translate(g2, gx0, leg['Y'])


def loop_pieces(sc, k, cx, cy, ccw, fA, fB, th0=-math.pi / 2.0):
    """Draw the gap-free 180-degree loop after row k, centered at (cx, cy), as the trunk band, side A, and side B.

    The path fraction u in [0, 1] maps to the center path from the leg end to the next leg's start. th0 is the angle of
    the leg end seen from the center, which is -90 deg when the next leg lies at larger y (the center above the leg
    end) and +90 deg when it lies below."""
    spec = sc.spec
    R = S.R_of(sc, k)
    lg = spec.l_g
    p_a = sc.x[k] + spec.Lam_L / 2.0
    L = math.pi * R

    def hw(u):
        return wc_at(sc, p_a + u * L) / 2.0
    # The band is extended by 1 um at both ends, straight along the end tangents, to overlap the straight legs.
    ext = 1.0
    dth = math.pi if ccw else -math.pi
    # Local +y (side A) is toward the center for a ccw turn and away from it for a cw turn.
    sA = -1.0 if ccw else +1.0
    trunk = annulus_piece(cx, cy, th0, dth, lambda u: R - hw(u), lambda u: R + hw(u), ext, L=L)
    pA = lambda u: spec.Lam_L / 2.0 + u * L            # center path from node k (u may leave [0, 1] by ext)
    A_in = lambda u: R + sA * (hw(u) + lg)
    A_out = lambda u: R + sA * (hw(u) + lg + float(fA(pA(u))))
    B_in = lambda u: R - sA * (hw(u) + lg)
    B_out = lambda u: R - sA * (hw(u) + lg + float(fB(pA(u))))
    sideA = annulus_piece(cx, cy, th0, dth, lambda u: min(A_in(u), A_out(u)), lambda u: max(A_in(u), A_out(u)), ext, L=L)
    sideB = annulus_piece(cx, cy, th0, dth, lambda u: min(B_in(u), B_out(u)), lambda u: max(B_in(u), B_out(u)), ext, L=L)
    return unary_union([trunk, sideA, sideB])


# ------------------------------------------------------------------------------------------------ inlet / outlet
def inlet_local(spec, x_top, filt=True, L_f=1500.0, clear=150.0):
    """Return the stem from the inlet port down to gap 1 (in leg-1 local coordinates) and its info, with an optional
    aggregate filter between the arms or above the head. The filter is a bulb of half-width up to 150 um that holds a
    staggered array of round posts (d 30, pitch 60 um)."""
    wc0 = spec.wc0
    PR = spec.port_r_in
    pieces = []
    if getattr(spec, 'inlet_gap', None) is not None:     # measured gap between the port's edge and the head's top
        clear = spec.inlet_gap
    x_port_c = x_top - clear - PR
    pieces.append(Point(x_port_c, 0).buffer(PR, quad_segs=64))
    if filt:
        # The bulb narrows to the trunk over 400 um, and that narrowing ends 50 um before anything else is drawn,
        # that is, before gap 1's stem when the filter sits between the arms and before the head's top (or the first
        # post) when it sits above the head. The narrowing used to end up to 100 um past the head's top, so that with a
        # hairpin or narrow head the channel just upstream of the first head mouth was up to 3 um a side wider than
        # the trunk. Also, the bulb was stem_half - 2 = 57 um wide, which is too narrow for one 30 um post, so every
        # matrix "filter" inlet was an empty bulb.
        x_f0 = x_top + 200.0                          # the filter sits between the arms, below their tops
        x_f1 = min(x_f0 + L_f, -450.0)
        between = x_f1 - x_f0 >= 400.0
        if not between:                               # arms too short to host a filter, so it goes above the head
            x_f1 = x_top - 450.0
            x_f0 = x_f1 - L_f
            x_port_c = x_f0 - 400.0 - PR
            pieces = [Point(x_port_c, 0).buffer(PR, quad_segs=64)]
        # Between the arms the bulb keeps the floor (at least 20 um) from the innermost lane, and above the head its
        # room is 150 um.
        room = (spec.R0 - spec.w_arm / 2.0 - max(spec.F_min, 20.0)) if between else 150.0
        Hf = max(spec.stem_half - 2.0, min(150.0, room))
        xs = np.linspace(x_port_c, 1.0, 2000)         # 1 um into the trunk, since with no head nothing else joins them

        def half(x):
            if x <= x_f0 - 300:                         # the port neck, at Hf
                return Hf
            if x <= x_f1:
                return Hf
            if x <= x_f1 + 400:                         # smooth narrowing to the stem
                t = (x - x_f1) / 400.0
                return Hf + (wc0 / 2.0 - Hf) * Gm.smoothstep(t)
            return wc0 / 2.0
        up = [(x, half(x)) for x in xs]
        stem = Polygon([(x, -y) for x, y in up] + [(x, y) for x, y in up[::-1]]).buffer(0)
        posts = []
        pitch, d = 60.0, 30.0
        j = 0
        x = x_f0 + 60.0
        while x < x_f1 - 30.0:
            yoff = 0.0 if j % 2 == 0 else pitch / 2.0
            y = -Hf + 2 * d + yoff
            while y < Hf - 2 * d:
                posts.append(Point(x, y).buffer(d / 2.0, quad_segs=32))
                y += pitch
            x += pitch * 0.866
            j += 1
        if not posts:
            raise ValueError('inlet filter: its bulb (half-width %.0f um) has no room for one 30 um post' % Hf)
        pieces.append(stem.difference(unary_union(posts)))
        info = dict(x_f0=x_f0, x_f1=x_f1, Hf=Hf, n_posts=len(posts))
    else:
        pieces.append(box(x_port_c, -wc0 / 2.0, 1.0, wc0 / 2.0))
        info = dict()
    info['x_port_c'] = x_port_c
    info['port_r'] = PR
    return unary_union(pieces), info


OUT_GAP = 150.0         # clearance between outlet features (um)


def outlet_frame(spec, wce, w_last):
    """Return the fixed geometry of the separate outlet in the last leg's local frame (flow +x, axis y = 0).

    The product (center) lead runs on the axis to the center port. Each filtrate (side) lead turns outward with a
    quarter arc, runs laterally to y_run, turns downstream with a second quarter arc, and runs to its own side port.
    y_run >= r_side + r_centre + gap puts the side ports clear of the center port wherever it lands, so the outlets can
    never merge."""
    wl = max(w_last, 40.0)
    rq = max(2.0 * wl, 60.0)
    y_side_c = wce / 2.0 + spec.l_g + wl / 2.0           # inner edges aligned when the lead is the wider
    y_run = max(y_side_c + 2.0 * rq, spec.port_r + spec.port_r_c + OUT_GAP + max(wl, wce) / 2.0)
    arcs = [(rq, math.pi / 2.0), (rq, math.pi / 2.0)]
    return dict(wl=wl, rq=rq, y_side_c=y_side_c, y_run=y_run, Ly=y_run - y_side_c - 2.0 * rq,
                Lc_min=rq + wl / 2.0 + OUT_GAP, arcs=arcs, excess=arcs_excess(spec, arcs, wl))


def arcs_excess(spec, arcs, w):
    """Return the extra 1-D length of a lead's arcs [(centerline radius, |angle| rad)] under the bend model, which is
    the exact Stokes curved-duct correction less the junctions (schedule.arc_excess_len). The lead's resistance is
    r(w) (geometric length + this excess)."""
    return float(sum(S.arc_excess_len(spec, r, a, w) for r, a in arcs))


def side_lead_path(fr, x_e, Lr, sg=1.0):
    """Return the side lead's centerline from (x_e - 1, sg y_side_c), which consists of a 1 um overlap into the outlet
    stretch, a quarter arc outward, the lateral run, a quarter arc downstream, and the downstream run Lr. The port's
    center lies at the end + (r_side - 1, 0)."""
    p = Gm.Path(x_e - 1.0, sg * fr['y_side_c'], 0.0)
    p.line(1.0)
    p.arc(fr['rq'], sg * 90.0)
    p.line(fr['Ly'])
    p.arc(fr['rq'], -sg * 90.0)
    p.line(Lr)
    return p


# ---- combined outlet (Dinh 2024 Fig. 2A; Xia 2016 Figs. 1a, 3a, 3g-i) ------------------------------------------
RING_SEG = 3.0          # max chord (um) of the ring's arcs
LR_GOAL = 1000.0        # planned run before the filtrate port (um); the final re-solve absorbs small changes


def ring_geom(spec, wce, w_last, w, Rr=None):
    """Return the geometry of the combined outlet, a teardrop.

    Both side channels leave the outlet stretch with an outward quarter arc (radius rq), wrap over the center port
    along a half circle of radius Rr centered on the port's x (a 'D' around it), come back toward the axis, turn
    downstream with a quarter arc (radius rho), and run side by side (with a wall 2 yb - w between them) into one
    filtrate port beyond the center port. The two paths meet only inside the filtrate port, so each side keeps its own
    outlet resistance and the closure is that of separate ports at one pressure. The center lead runs on the axis into
    the center port inside the ring. The run Lr before the filtrate port is the tuning length, solved for the exact
    closure target."""
    rq = max(2.0 * w, 60.0)
    rho = max(1.5 * w, 60.0)
    # A lead wider than the side channel starts with its inner edge on the side channel's inner edge, so that it does
    # not eat into the wall between the center lead and the side channel, which left an 18 um wall tip on Dinh.
    y_sc = wce / 2.0 + spec.l_g + max(w, w_last) / 2.0
    rc, rs = spec.port_r_c, spec.port_r
    yb = w / 2.0 + max(spec.F_min, 20.0) / 2.0 + 5.0           # branch axis offset at the filtrate port
    Rr_min = rc + OUT_GAP + w / 2.0
    Rr = Rr_min if Rr is None else max(Rr, Rr_min)
    Ly = 0.0
    y0 = y_sc + rq + Ly
    if y0 < yb + rho:                                          # room for the descent and the quarter arc
        Ly = yb + rho - y0
        y0 = yb + rho
    Lc = rq + Rr - rc                                          # center lead, stretch end to the center port's edge
    sx = math.sqrt(max(rs ** 2 - yb ** 2, 0.0))
    # Lr_min keeps Lr >= 0 and the filtrate port clear of the center port.
    Lr_min = max(0.0, rs + rc + OUT_GAP - Rr - rho - sx)
    const = math.pi * rq / 2.0 + Ly + math.pi * Rr + (y0 - yb - rho) + math.pi * rho / 2.0
    arcs = [(rq, math.pi / 2.0), (Rr, math.pi), (rho, math.pi / 2.0)]
    return dict(w=w, rq=rq, rho=rho, y_sc=y_sc, Rr=Rr, Ly=Ly, y0=y0, yb=yb, Lc=Lc, const=const, rc=rc, rs=rs,
                sx=sx, Lr_min=Lr_min, arcs=arcs, excess=arcs_excess(spec, arcs, w))


def ring_path(g, x_e, Lr, sg=1.0):
    """Return the side lead of the combined outlet, which ends 1 um inside the filtrate port (whose center lies at
    x_fp = ring_fp(g, x_e, Lr), y = 0)."""
    p = Gm.Path(x_e - 1.0, sg * g['y_sc'], 0.0)
    p.line(1.0)
    p.arc(g['rq'], sg * 90.0, seg_max=RING_SEG)
    p.line(g['Ly'])
    p.arc(g['Rr'], -sg * 180.0, seg_max=RING_SEG)
    p.line(g['y0'] - g['yb'] - g['rho'])
    p.arc(g['rho'], sg * 90.0, seg_max=RING_SEG)
    p.line(Lr + 1.0)
    return p


def ring_fp(g, x_e, Lr):
    """Return the x of the filtrate port's center for the teardrop geometry g, the stretch end x_e, and the run Lr."""
    return x_e + g['rq'] + 2.0 * g['Rr'] + g['rho'] + Lr + g['sx']


# ---- Xia's ring (outlet_style 'ring', Xia 2016 Fig. 3a; r4/meas_xia_fig3a.py) -----------------------------------
def ring2_geom(spec, wce, w_last, w_lead=None):
    """Return the ring outlet's frame, which is independent of the lead width w except when the lead (w_lead) is wider
    than the side channel.

    The frame is the last leg's local frame, with the end of the outlet stretch at x = 0 and flow +x. Each side lead
    leaves its side channel (center y_sc) with a turn of radius rq onto one circle (centerline radius R, center on the
    axis at dx_rc), follows it around the center port, and enters the filtrate port (radius rs, on the axis at x_fp)
    from its own side, where the circle crosses the port's edge. The center port (radius rc) lies at
    x_cp = dx_rc - ring_dc, and the filtrate port lies at x_cp + ring_dfp (Xia's offsets). The frame fixes the lead's
    path length L, so the closure fixes its width alone."""
    R, rq, dc, dfp = spec.ring_R, spec.ring_rq, spec.ring_dc, spec.ring_dfp
    rc, rs = spec.port_r_c, spec.port_r
    # A lead wider than the side channel starts with its axis-facing edge on the side channel's inner edge, like the
    # teardrop's, because a lead centered on the side channel would overlap the stretch wall and the center lead, and
    # the center and filtrate flows would meet (Gifford Fig. 6 has a 313 um lead from a 140 um side channel).
    y_sc = wce / 2.0 + spec.l_g + max(w_last, w_lead or 0.0) / 2.0
    ca = (y_sc + rq) / (R + rq)
    if not 0.0 < ca < 1.0:
        raise OutletInfeasible('ring outlet: side channel at %.0f um outside the ring (R %.0f)' % (y_sc, R))
    al = math.acos(ca)                                    # the turn onto the circle
    dx_rc = (R + rq) * math.sin(al)
    x_cp = dx_rc - dc
    Lc = x_cp - rc
    x_fp = x_cp + dfp
    D = dfp - dc
    a = (R ** 2 - rs ** 2 + D ** 2) / (2.0 * D)
    if abs(a) >= R:
        raise OutletInfeasible('ring outlet: the circle misses the filtrate port')
    psi = math.acos(a / R)
    ph_j, ph_e = math.pi / 2.0 - al, math.pi - psi       # angles on the circle from its upstream-most point
    L = rq * al + R * (ph_e - ph_j)
    return dict(style='ring', R=R, rq=rq, y_sc=y_sc, al=al, dx_rc=dx_rc, x_cp=x_cp, Lc=Lc, x_fp=x_fp, ph_j=ph_j,
                ph_e=ph_e, L=L, rc=rc, rs=rs, dc=dc, dfp=dfp, y_inner=wce / 2.0 + spec.l_g, w_frame=w_lead)


def ring2_check(g, w):
    """Check the clearances of the ring outlet for lead width w, and raise OutletInfeasible on a violation."""
    errs = []
    if g['Lc'] < 50.0:
        errs.append('centre lead %.0f um < 50' % g['Lc'])
    if g['R'] - w / 2.0 - (g['dc'] + g['rc']) < OUT_GAP:
        errs.append('ring inner edge %.0f um from the centre port' % (g['R'] - w / 2.0 - g['dc'] - g['rc']))
    if g['dfp'] - g['rc'] - g['rs'] < OUT_GAP:
        errs.append('ports %.0f um apart' % (g['dfp'] - g['rc'] - g['rs']))
    if g['rq'] - w / 2.0 < 10.0:
        errs.append('turn radius %.0f um leaves an inner edge radius below 10 um at w %.0f' % (g['rq'], w))
    if g['y_sc'] - w / 2.0 < g.get('y_inner', 0.0) - 0.5:
        errs.append('the lead (%.0f um) overlaps the stretch wall: its axis-facing edge is %.1f um inside the side '
                    "channel's inner edge" % (w, g.get('y_inner', 0.0) - (g['y_sc'] - w / 2.0)))
    if errs:
        raise OutletInfeasible('ring outlet: ' + '; '.join(errs))


def ring2_path(g, x_e, sg=1.0, ext=0.0):
    """Return the side lead of the ring outlet, which consists of a 1 um overlap into the stretch, the turn onto the
    circle, and the circle to the filtrate port's edge, followed by 1 um more (the resistance path) or `ext` rad more
    (the drawn ribbon) inside the port."""
    p = Gm.Path(x_e - 1.0, sg * g['y_sc'], 0.0)
    p.line(1.0)
    p.arc(g['rq'], sg * math.degrees(g['al']), seg_max=RING_SEG)
    dph = g['ph_e'] - g['ph_j'] + (ext if ext > 0.0 else 1.0 / g['R'])
    p.arc(g['R'], -sg * math.degrees(dph), seg_max=RING_SEG)
    return p


def ring2_ext(g, w, x_e):
    """Return the extra angle (rad) along the circle that puts both corners of the lead's flat end 1 um inside the
    filtrate port, so that the port disc and the ribbon leave no solid notch between them."""
    xr, xf, R, rs = x_e + g['dx_rc'], x_e + g['x_fp'], g['R'], g['rs']
    t = 1.0 / R
    while t < 0.6:
        ph = g['ph_e'] + t
        ok = all(math.hypot(xr - rr_ * math.cos(ph) - xf, rr_ * math.sin(ph)) <= rs - 1.0 for rr_ in (R - w / 2.0, R + w / 2.0))
        if ok:
            return t
        t += 0.0005
    raise OutletInfeasible('ring outlet: lead end not inside the filtrate port')


def ring2_arcs(g):
    """Return the ring lead's arcs as [(centerline radius, |angle| rad)], which are the turn onto the circle and
    the circle itself."""
    return [(g['rq'], g['al']), (g['R'], g['ph_e'] - g['ph_j'])]


def ring2_width(g, target, h, spec=None):
    """Return the lead width whose resistance r(w) (L + bend excess(w)) is the target, and raise OutletInfeasible when
    the target is not positive or no width from 1 um to 100 mm gives it."""
    from scipy.optimize import brentq
    arcs = ring2_arcs(g)
    ex = (lambda ww: arcs_excess(spec, arcs, ww)) if spec is not None else (lambda ww: 0.0)
    fn = lambda ww: S.rr(ww, h) * (g['L'] + ex(min(ww, 1.9 * g['rq']))) * S.UM - target
    if target <= 0.0:
        raise OutletInfeasible('ring outlet: the closure leaves the side lead a resistance of %.3g Pa s/m^3 (not positive): '
                               "the side channel's straight run to the outlet already exceeds the side outlet resistance "
                               "that the ring's fixed centre lead allows" % target)
    if fn(1.0e5) > 0.0:
        raise OutletInfeasible("ring outlet: the side lead would need %.3g Pa s/m^3, less than even a 100 mm wide lead "
                               "along the ring's %.0f um path gives" % (target, g['L']))
    if fn(1.0) < 0.0:
        raise OutletInfeasible("ring outlet: the side lead would need %.3g Pa s/m^3, more than a 1 um wide lead along the "
                               "ring's %.0f um path gives" % (target, g['L']))
    return brentq(fn, 1.0, 1.0e5, xtol=1e-9)


def combined_plan(sc):
    """Plan the combined outlet for a schedule closed with Rc_lead = 0, and return (Rc_lead, Lc, geometry).

    The teardrop takes the smallest ring, and the ring sets the center lead's length (Lc = rq + Rr - r_c), hence
    R_c,lead. R_s,out is linear in it:

        R_s,out = R_s,out(0) + (P_N / S_N) R_c,lead

    Lr follows from the exact target. The lead width w (at least 40 um) is chosen so that Lr is LR_GOAL longer than its
    clearance minimum Lr_min, or is 40 um when a 40 um lead already leaves a longer run. For Xia's ring
    (outlet_style 'ring'), the frame and the lead width are iterated to a fixed point instead."""
    from scipy.optimize import brentq
    spec, h, N = sc.spec, sc.spec.h, sc.N
    wce, w_last = sc.wc[-1], sc.ws[N - 2]
    slope = sc.P[N - 1] / sc.S[N - 1]
    if getattr(spec, 'outlet_style', 'teardrop') == 'ring':
        # The frame depends on the lead width once the lead is wider than the side channel, so we iterate to a fixed
        # point.
        w_frame, w_end = None, w_last
        for _ in range(400):                              # damped, so the error halves per step from any start
            g = ring2_geom(spec, wce, w_last, w_lead=w_frame)
            Rc_lead = S.rr(wce, h) * g['Lc'] * S.UM
            # The stretch end's excess is taken at the current width estimate w_end (its ratio to w_s(N)).
            target = sc.Rs_out + slope * Rc_lead - stretch_R(spec, w_last, w_end, h) - R_ext_side(spec)
            w = ring2_width(g, target, h, spec)
            end_ok = abs(w - w_end) < 1e-7
            if w <= w_last:                               # the frame does not depend on a narrower lead
                if end_ok:
                    break
                w_end = 0.5 * (w + w_end)
                continue
            if w_frame is not None and abs(w - w_frame) < 1e-7 and end_ok:
                break
            w_frame = w if w_frame is None else 0.5 * (w + w_frame)
            w_end = w_frame
        else:
            raise OutletInfeasible('ring outlet: the lead width does not settle (frame %s, width %.3f um)' % (w_frame, w))
        g['w'] = w
        ring2_check(g, g['w'])
        return Rc_lead, g['Lc'], g

    def Lr_for(w):
        g = ring_geom(spec, wce, w_last, w)
        Rc_lead = S.rr(wce, h) * g['Lc'] * S.UM
        target = sc.Rs_out + slope * Rc_lead - stretch_R(spec, w_last, w, h) - R_ext_side(spec)
        L_res = target / (S.rr(w, h) * S.UM)
        return L_res - g['const'] - g['excess'] - g['Lr_min'] - LR_GOAL, g, Rc_lead
    # The lead width is chosen so that the run before the filtrate port is LR_GOAL longer than its clearance minimum
    # Lr_min (a compact teardrop), with wider leads for a small target and narrower ones (down to 40 um) for a large
    # one, for example Xia's, which the off-chip tubing sets.
    w_lo, w_hi = 40.0, 20.0 * max(w_last, 40.0) + 2000.0
    if Lr_for(w_lo)[0] >= 0.0:
        w = w_lo
    elif Lr_for(w_hi)[0] < 0.0:
        raise OutletInfeasible('teardrop outlet: the side lead\'s resistance target leaves a run shorter than the '
                               'teardrop\'s clearance minimum at every lead width from %.0f to %.0f um' % (w_lo, w_hi))
    else:
        w = brentq(lambda ww: Lr_for(ww)[0], w_lo, w_hi, xtol=1e-6)
    d, g, Rc_lead = Lr_for(w)
    g['Lr_plan'] = d + g['Lr_min'] + LR_GOAL
    return Rc_lead, g['Lc'], g


def R_ext_side(spec):
    """Return the per-side off-chip filtrate resistance, which is doubled for a combined port because its tube carries
    both sides (2 S), whereas a separate port's tube carries one."""
    return spec.R_ext_s * (2.0 if getattr(spec, 'outlet', 'separate') == 'combined' else 1.0)


def flare_on(spec):
    """Return True when the spec asks for a flared outlet stretch (flare_deg is set)."""
    return getattr(spec, 'flare_deg', None) is not None


def flare_len(spec, w_last, w_lead):
    """Return the length (um) of the flared widening of the outlet stretch from w_last to w_lead at spec.flare_deg per
    side, or 0 when the lead is not wider than the side channel or the flare is off. The inner wall stays straight and
    only the outer wall diverges, as in Dinh 2024 Fig. 2B."""
    if not flare_on(spec) or w_lead <= w_last + 1e-9:
        return 0.0
    return (w_lead - w_last) / math.tan(math.radians(float(spec.flare_deg)))


FLARE_KINK_R = 25.0     # radius (um) of the tangent arc at each kink of the flared outer wall, drawn and integrated


def flare_profile(spec, w_last, w_lead, L_out):
    """Return the outer wall of a flared outlet stretch as (s, w) points from the last gap node (s = 0) to the lead
    (s = L_out).

    The wall runs straight at w_last, turns along a tangent arc of radius FLARE_KINK_R into the linear flare at
    flare_deg per side, follows the flare, and turns back to parallel along a tangent arc that ends at w_lead at L_out.
    The wall is piecewise linear between the points. The drawing uses them as its outer-wall polyline, and stretch_R
    integrates the same polyline, so the two agree exactly."""
    th = math.radians(float(spec.flare_deg))
    Lf = flare_len(spec, w_last, w_lead)
    if Lf <= 0.0:
        return [(0.0, w_last), (L_out, w_last)]
    R, tt = flare_kink(spec, Lf)
    s2 = L_out - tt                                   # kink 2 (the end of the linear flare), whose arc ends at L_out
    s1 = s2 - Lf                                      # kink 1 (the end of the straight run)
    if s1 - tt < 0.0:
        raise ValueError('flared outlet: the straight run %.1f um is too short for the kink arc' % (s1))
    pts = [(0.0, w_last)]
    # Kink 1 at (s1, w_last) is a convex solid corner with its center at (s1 - tt, w_last + R).
    for ph in np.linspace(0.0, th, 13):
        pts.append((s1 - tt + R * math.sin(ph), w_last + R - R * math.cos(ph)))
    # Kink 2 at (s2, w_lead) is a concave solid corner with its center at (s2 + tt, w_lead - R).
    for ph in np.linspace(th, 0.0, 13):
        pts.append((s2 + tt - R * math.sin(ph), w_lead - R + R * math.cos(ph)))
    return pts


def flare_kink(spec, Lf):
    """Return the radius and the tangent length of the kink arcs. The radius is FLARE_KINK_R unless the flare is too
    short for two arcs."""
    th = math.radians(float(spec.flare_deg))
    R = min(FLARE_KINK_R, 0.45 * Lf / math.tan(th / 2.0))
    return R, R * math.tan(th / 2.0)


def stretch_R(spec, w_last, w_lead, h, n=401):
    """Return the 1-D resistance (SI) of one side channel's outlet stretch, from the last gap node to the start of the
    lead path.

    Without a flare, the stretch is straight at w_last over L_out. With a flare, it is straight at w_last over
    L_out - L_f and then widens linearly to w_lead over L_f = flare_len, and the exact wall rule integrates r(w(s))
    along the polyline of flare_profile, its tangent arcs included, with the trapezoid rule on stations about 0.1 um
    apart on the sloped pieces. The excess of the stretch's end (schedule.outlet_exp_len) is added, which accounts
    for the abrupt width change into a lead of another width, or for the flare's excess over its integral."""
    R_end = S.rr(min(w_last, w_lead), h) * S.outlet_exp_len(spec, w_last, w_lead)
    if not flare_on(spec):
        return (S.rr(w_last, h) * spec.L_out + R_end) * S.UM
    Lf = flare_len(spec, w_last, w_lead)
    if Lf <= 0.0 or spec.L_out - Lf <= 0.0:
        return (S.rr(w_last, h) * spec.L_out + R_end) * S.UM
    pts = flare_profile(spec, w_last, w_lead, spec.L_out)
    R = R_end
    for (s0, w0), (s1, w1) in zip(pts[:-1], pts[1:]):        # each straight piece of the polyline, integrated
        if s1 <= s0:
            continue
        # Stations 0.1 um apart keep the trapezoid error below 1e-6.
        ss = np.linspace(s0, s1, 2 if abs(w1 - w0) < 1e-12 else max(41, int((s1 - s0) / 0.1) + 1))
        ww = w0 + (w1 - w0) * (ss - s0) / (s1 - s0)
        R += float(np.trapezoid([S.rr(x, h) for x in ww], ss))
    return R * S.UM


def exp_info(spec, w_last, w_lead, h):
    """Return the record's account of the excess term at the stretch's end."""
    c, deep, shal, lam, kind = S.outlet_exp_c(spec, w_last, w_lead, parts=True)
    wn = min(w_last, w_lead)
    d = dict(kind=kind, w_side_um=float(w_last), w_lead_um=float(w_lead), ratio=float(max(w_last, w_lead) / wn),
             c=float(c), c_deep=float(deep), c_shallow=float(shal), depth_weight=float(lam),
             excess_len_um=float(c * wn), R=float(S.rr(wn, h) * c * wn * S.UM),
             table_domain=S.outlet_exp_domain(spec, w_last, w_lead) or 'inside the solved tables')
    if flare_on(spec) and w_lead > w_last:
        Lf = flare_len(spec, w_last, w_lead)
        d.update(flare_deg=float(spec.flare_deg), flare_len_um=float(Lf),
                 straight_run_um=float(spec.L_out - Lf - flare_kink(spec, Lf)[1]) if Lf > 0 else float(spec.L_out))
    return d


def stretch_total(spec, kind, base_L_out, w_last, w_lead):
    """Return the outlet stretch's total length with the flare.

    The straight run is spec.flare_start where no U-turn precedes the outlet (kind 'linear', which is a linear layout
    or one leg of any layout), or the layout's clearance-based stretch where that is longer (kind 'clearance', where
    the last leg ends beside a loop), and the flare and its end arc follow. With no widening (the lead not wider than
    the side channel), the stretch is base_L_out, as without a flare, because a flare on a lead that does not widen
    must not shorten the stretch (rc5 review)."""
    Lf = flare_len(spec, w_last, w_lead)
    if Lf <= 0.0:
        return base_L_out
    # The outlet-stretch extension (spec.L_out_extra, included in base_L_out) lengthens the straight run, and the
    # flare's start must not absorb it. The extension is zero for a design that is clear at once, which is then
    # unchanged (code review 4, F7).
    ex = float(getattr(spec, 'L_out_extra', 0.0) or 0.0)
    Ls = (spec.flare_start if kind == 'linear' else max(base_L_out - ex, spec.flare_start)) + ex
    return Ls + Lf + flare_kink(spec, Lf)[1]          # the straight run, the flare, the end arc's tangent length


def _flare_fixed_point(sc, turns, Rt, kind, base_L_out):
    """Iterate (damped) to the fixed point of a flared outlet. With a flared lead, the outlet stretch's length depends
    on the lead width (L_out = straight run + L_f(w)), and the lead width depends on the closure, which depends on
    L_out through the center stretch. The function returns the schedule copy with spec.L_out set, and its sc.info
    carries the base length and the lead width."""
    import copy
    spec = sc.spec
    N = sc.N
    w_last = sc.ws[N - 2]
    kind = 'clearance' if turns else 'linear'         # one leg of a serpentine or paperclip is linear here (rc5 review)
    sc = copy.copy(sc)
    sc.spec = copy.copy(spec)
    L_prev = None
    for it in range(60):
        Rc_lead, Lc, fr = outlet_Rc_lead(S.reclose(sc, turns, Rt, sc.spec.R_ext_c))
        w_lead = float(fr['w'] if 'w' in fr else fr['wl'])
        L_new = stretch_total(sc.spec, kind, base_L_out, w_last, w_lead)
        if L_prev is not None and abs(L_new - L_prev) < 1e-7:
            break
        L_prev = L_new if L_prev is None else 0.5 * (L_new + L_prev)
        sc.spec.L_out = L_prev
    else:
        raise OutletInfeasible('flared outlet: the stretch length does not settle (%.4f vs %.4f um)' % (L_new, L_prev))
    sc.spec.L_out = L_new
    sc.info = dict(sc.info, L_out_base=base_L_out, w_lead=w_lead, flare_len=flare_len(sc.spec, w_last, w_lead),
                   flare_kind=kind)
    return sc


def outlet_Rc_lead(sc, Lr_plan=300.0):
    """Choose the center lead's resistance at the plan stage, and return (Rc_lead, Lc, frame).

    The closure fixes only the pair (R_c,out, R_s,out) through R_s,out = (R_c,out P_N - R_g,N Q_g,N) / S_N. For the
    separate outlet, the side lead is fixed at width max(w_s,last, 40) with a downstream run of Lr_plan (slack for the
    fail-safe re-plan), and the center lead's length follows. If that length is shorter than the clearance minimum, the
    center lead is set to the minimum, and the side lead's run absorbs the difference. The combined outlet is planned
    by combined_plan."""
    if getattr(sc.spec, 'outlet', 'separate') == 'combined':
        return combined_plan(sc)
    spec, h, N = sc.spec, sc.spec.h, sc.N
    wce, w_last = sc.wc[-1], sc.ws[N - 2]
    fr = outlet_frame(spec, wce, w_last)
    L_side = side_lead_path(fr, 0.0, Lr_plan).length() - 2.0      # to the port edge (both 1 um overlaps removed)
    Rs_out = stretch_R(spec, w_last, fr['wl'], h) + S.rr(fr['wl'], h) * (L_side + fr['excess']) * S.UM + R_ext_side(spec)
    Rc_out = (Rs_out * sc.S[N - 1] + sc.Rg[N - 1] * sc.Qg[N - 1]) / sc.P[N - 1]
    Lc = (Rc_out - sc.Rc[N - 1] - spec.R_ext_c) / (S.rr(wce, h) * S.UM)
    if Lc < fr['Lc_min'] + 50.0:
        Lc = fr['Lc_min'] + 50.0
    return S.rr(wce, h) * Lc * S.UM, Lc, fr


def _l2g_pts(pts, leg, gx0):
    """Map points in a leg's local frame to global coordinates, as to_global does for geometry."""
    if leg['dirn'] > 0:
        return [(gx0 + x, leg['Y'] + y) for x, y in pts]
    return [(gx0 - x, leg['Y'] - y) for x, y in pts]


def ribbon_into_port(p, w, r_port=PORT_R):
    """Return the ribbon of a lead whose path ends 1 um inside a port. The drawn ribbon runs on by the port's sagitta
    at the lead's half-width (+1 um), so the lead's flat end lies wholly inside the disc. A flat end at the port edge
    leaves two solid notches between the lead's corners and the circle (12 um wide for a 300 um lead into a 750 um
    port, flagged by the DRC). The resistance path (and its length) still ends at the port edge."""
    sag = r_port - math.sqrt(max(r_port ** 2 - (w / 2.0) ** 2, 0.0)) + 1.0
    (x0, y0), (x1, y1) = p.p[-2], p.p[-1]
    d = math.hypot(x1 - x0, y1 - y0)
    ext = (x1 + (x1 - x0) / d * sag, y1 + (y1 - y0) / d * sag)
    return LineString(p.p + [ext]).buffer(w / 2.0, cap_style='flat', join_style='round', quad_segs=8)


def outlet_global(sc, leg, gx0):
    """Draw the outlet at the end of the last leg in global coordinates, combined or separate as the spec asks, and
    return it with its info."""
    if getattr(sc.spec, 'outlet', 'separate') == 'combined':
        return outlet_global_combined(sc, leg, gx0)
    return outlet_global_separate(sc, leg, gx0)


def outlet_global_separate(sc, leg, gx0):
    """Draw the separate outlet, with the center lead and the center port and two side leads, each to its own side
    port (no merge and no enclosed solid). The center lead's length is the plan-stage value
    (sc.info['Rc_lead']), and the side leads' downstream run is solved so that each side's outlet resistance equals
    the closure target R_s,out exactly (minus the drawn outlet stretch)."""
    spec = sc.spec
    h, N = spec.h, sc.N
    rs, rc = spec.port_r, spec.port_r_c
    x_e = sc.x[N - 1] - (leg['p0'] or 0.0) + spec.L_out          # local end of the outlet stretch
    wce, w_last = sc.wc[-1], sc.ws[N - 2]
    fr = outlet_frame(spec, wce, w_last)
    Lc = sc.info['Rc_lead'] / (S.rr(wce, h) * S.UM)
    target = sc.Rs_out - stretch_R(spec, w_last, fr['wl'], h) - R_ext_side(spec)   # per side, stretch end to port
    L0 = side_lead_path(fr, x_e, 0.0).length() - 2.0         # minus the 1 um overlaps at both ends (to the port edge)
    Lr = target / (S.rr(fr['wl'], h) * S.UM) - L0 - fr['excess']
    if Lr < 0.0:
        raise OutletInfeasible('separate outlet: the side lead\'s resistance target is below its fixed path (run %.1f um)' % Lr)
    pc = Gm.Path(x_e - 1.0, 0.0, 0.0)
    pc.line(Lc + 2.0)                                         # Lc to the port edge, the path ending 1 um inside
    xcp = x_e + Lc + rc
    centre = unary_union([ribbon_into_port(pc, wce, rc), Point(xcp, 0).buffer(rc, quad_segs=64)])
    pieces = [centre]
    paths, ports = {}, {}
    for side, sg in (('A', 1.0), ('B', -1.0)):
        p = side_lead_path(fr, x_e, Lr, sg)
        xe_, ye_ = p.xy
        pieces += [ribbon_into_port(p, fr['wl'], rs), Point(xe_ + rs - 1.0, ye_).buffer(rs, quad_segs=64)]
        paths[side] = _l2g_pts(p.p, leg, gx0)
        ports[side] = _l2g_pts([(xe_ + rs - 1.0, ye_)], leg, gx0)[0]
    paths['C'] = _l2g_pts(pc.p, leg, gx0)
    ports['C'] = _l2g_pts([(xcp, 0.0)], leg, gx0)[0]
    g = unary_union(pieces)
    L_lead = side_lead_path(fr, x_e, Lr).length() - 2.0          # to the port edge
    check_exp_domain(spec, w_last, fr['wl'])                     # the domain of the stretch-end tables
    info = dict(kind='separate', w_lead=fr['wl'], L_lead=L_lead, R_side_target=target,
                expansion_term=exp_info(spec, w_last, fr['wl'], h), R_stretch_design=stretch_R(spec, w_last, fr['wl'], h),
                R_side_drawn=S.rr(fr['wl'], h) * (L_lead + fr['excess']) * S.UM, Lc_lead=Lc, w_centre_lead=wce,
                bend_arcs=dict(A=fr['arcs'], B=fr['arcs'], C=[]), bend_excess=fr['excess'],
                Rc_lead=sc.info['Rc_lead'], Lr_run=Lr, y_run=fr['y_run'], x_e_local=x_e, holes=0,
                R_ext_c=spec.R_ext_c, R_ext_side=R_ext_side(spec),
                paths=paths, ports=ports, port_r=dict(A=rs, B=rs, C=rc))
    return to_global(g, leg, gx0), info


def outlet_global_combined(sc, leg, gx0):
    """Draw the combined outlet (ring_geom). The lead width w and the ring are the plan-stage values
    (sc.info['ring']), and Lr is re-solved here so that each side path takes the final closure target exactly. Xia's
    ring style is drawn by outlet_global_ring."""
    spec = sc.spec
    h, N = spec.h, sc.N
    x_e = sc.x[N - 1] - (leg['p0'] or 0.0) + spec.L_out
    wce, w_last = sc.wc[-1], sc.ws[N - 2]
    if getattr(spec, 'outlet_style', 'teardrop') == 'ring':
        return outlet_global_ring(sc, leg, gx0, x_e, wce, w_last)
    gp = sc.info['ring']
    g = ring_geom(spec, wce, w_last, gp['w'], gp['Rr'])
    Lc = sc.info['Rc_lead'] / (S.rr(wce, h) * S.UM)
    if abs(Lc - g['Lc']) > 1e-6:
        raise ValueError('combined outlet: centre lead %.4f != ring %.4f' % (Lc, g['Lc']))
    if flare_on(spec) and abs(spec.L_out - stretch_total(spec, sc.info.get('flare_kind', 'linear'), sc.info.get('L_out_base', spec.L_out),
                                                       w_last, g['w'])) > 1e-6:
        raise ValueError('flared outlet: the stretch length %.4f um is not the planned %.4f um' % (
            spec.L_out, stretch_total(spec, sc.info.get('flare_kind', 'linear'), sc.info.get('L_out_base', spec.L_out), w_last, g['w'])))
    target = sc.Rs_out - stretch_R(spec, w_last, g['w'], h) - R_ext_side(spec)
    L_res = target / (S.rr(g['w'], h) * S.UM)
    Lr = L_res - g['const'] - g['excess']
    if Lr < g['Lr_min'] - 1e-6:
        raise OutletInfeasible('teardrop outlet: run %.2f um below its clearance minimum %.2f um' % (Lr, g['Lr_min']))
    rc, rs = g['rc'], g['rs']
    pc = Gm.Path(x_e - 1.0, 0.0, 0.0)
    pc.line(Lc + 2.0)                                         # Lc to the port edge, the path ending 1 um inside
    xcp = x_e + Lc + rc                                       # the ring's center, as ring_geom defines it
    pieces = [ribbon_into_port(pc, wce, rc), Point(xcp, 0).buffer(rc, quad_segs=64)]
    paths = {}
    x_fp = ring_fp(g, x_e, Lr)
    for side, sg in (('A', 1.0), ('B', -1.0)):
        p = ring_path(g, x_e, Lr, sg)
        # The drawn ribbon is extended so that its flat end (corners at yb +- w/2) lies inside the port disc.
        yo = min(g['yb'] + g['w'] / 2.0, 0.999 * rs)
        ext = g['sx'] - math.sqrt(max(rs ** 2 - yo ** 2, 0.0)) + 1.0
        xb, yb_ = p.p[-1]
        pieces.append(LineString(p.p + [(xb + ext, yb_)]).buffer(g['w'] / 2.0, cap_style='flat', join_style='round',
                                                                 quad_segs=8))
        paths[side] = _l2g_pts(p.p, leg, gx0)
    pieces.append(Point(x_fp, 0.0).buffer(rs, quad_segs=64))
    fp = _l2g_pts([(x_fp, 0.0)], leg, gx0)[0]
    ports = dict(A=fp, B=fp, C=_l2g_pts([(xcp, 0.0)], leg, gx0)[0])
    paths['C'] = _l2g_pts(pc.p, leg, gx0)
    gg = unary_union(pieces)
    L_lead = ring_path(g, x_e, Lr).length() - 2.0          # minus the 1 um overlaps at both ends
    check_exp_domain(spec, w_last, g['w'])                       # the domain of the stretch-end tables
    exp_term = exp_info(spec, w_last, g['w'], h)
    info = dict(kind='combined', w_lead=g['w'], L_lead=L_lead, R_side_target=target, expansion_term=exp_term,
                R_stretch_design=stretch_R(spec, w_last, g['w'], h),
                R_side_drawn=S.rr(g['w'], h) * (L_lead + g['excess']) * S.UM, Lc_lead=Lc, w_centre_lead=wce,
                bend_arcs=dict(A=g['arcs'], B=g['arcs'], C=[]), bend_excess=g['excess'],
                Rc_lead=sc.info['Rc_lead'], Lr_run=Lr, Rr=g['Rr'], Ly=g['Ly'], yb=g['yb'], x_e_local=x_e, holes=1,
                R_ext_c=spec.R_ext_c, R_ext_side=R_ext_side(spec),
                paths=paths, ports=ports, port_r=dict(A=rs, B=rs, C=rc))
    return to_global(gg, leg, gx0), info


def outlet_global_ring(sc, leg, gx0, x_e, wce, w_last):
    """Draw Xia's ring outlet (ring2_geom). The frame is fixed, and the lead width is re-solved here so that each side
    path takes the final closure target exactly. The center lead's length does not depend on the lead width."""
    spec = sc.spec
    h = spec.h
    g = ring2_geom(spec, wce, w_last, w_lead=(sc.info.get('ring') or {}).get('w_frame'))   # the plan's frame
    Lc = sc.info['Rc_lead'] / (S.rr(wce, h) * S.UM)
    if abs(Lc - g['Lc']) > 1e-6:
        raise ValueError('ring outlet: centre lead %.4f != frame %.4f' % (Lc, g['Lc']))
    if flare_on(spec):
        raise ValueError('flared outlet lead: not available with the ring outlet (teardrop and separate only)')
    # The stretch end's excess depends on the lead width (its ratio to w_s(N)), so the width is found as a fixed point
    # that starts from the plan's width.
    w = float((sc.info.get('ring') or {}).get('w', w_last))
    for _ in range(200):
        target = sc.Rs_out - stretch_R(spec, w_last, w, h) - R_ext_side(spec)
        w_new = ring2_width(g, target, h, spec)
        if abs(w_new - w) < 1e-9:
            w = w_new
            break
        w = w_new
    else:
        raise OutletInfeasible('ring outlet: the lead width does not settle with the stretch end term (%.6f um)' % w)
    target = sc.Rs_out - stretch_R(spec, w_last, w, h) - R_ext_side(spec)
    ring2_check(g, w)
    ex_ring = arcs_excess(spec, ring2_arcs(g), w)
    rc, rs = g['rc'], g['rs']
    pc = Gm.Path(x_e - 1.0, 0.0, 0.0)
    pc.line(Lc + 2.0)                                         # Lc to the port edge, the path ending 1 um inside
    xcp, xfp = x_e + g['x_cp'], x_e + g['x_fp']
    pieces = [ribbon_into_port(pc, wce, rc), Point(xcp, 0).buffer(rc, quad_segs=64),
              Point(xfp, 0.0).buffer(rs, quad_segs=96)]
    ext = ring2_ext(g, w, x_e)
    paths = {}
    for side, sg in (('A', 1.0), ('B', -1.0)):
        pieces.append(LineString(ring2_path(g, x_e, sg, ext=ext).p).buffer(w / 2.0, cap_style='flat',
                                                                         join_style='round', quad_segs=8))
        paths[side] = _l2g_pts(ring2_path(g, x_e, sg).p, leg, gx0)
    fp = _l2g_pts([(xfp, 0.0)], leg, gx0)[0]
    ports = dict(A=fp, B=fp, C=_l2g_pts([(xcp, 0.0)], leg, gx0)[0])
    paths['C'] = _l2g_pts(pc.p, leg, gx0)
    L_lead = ring2_path(g, x_e).length() - 2.0          # minus the 1 um overlaps at both ends
    check_exp_domain(spec, w_last, w)                            # the domain of the stretch-end tables
    info = dict(kind='combined', style='ring', w_lead=w, L_lead=L_lead, R_side_target=target,
                expansion_term=exp_info(spec, w_last, w, h), R_stretch_design=stretch_R(spec, w_last, w, h),
                R_side_drawn=S.rr(w, h) * (L_lead + ex_ring) * S.UM, Lc_lead=Lc, w_centre_lead=wce, Rc_lead=sc.info['Rc_lead'],
                bend_arcs=dict(A=ring2_arcs(g), B=ring2_arcs(g), C=[]), bend_excess=ex_ring,
                Rr=g['R'], ring_centre_local=x_e + g['dx_rc'], x_cp_local=xcp, x_fp_local=xfp, turn_deg=math.degrees(g['al']),
                arc_deg=math.degrees(g['ph_e'] - g['ph_j']), x_e_local=x_e, holes=1,
                R_ext_c=spec.R_ext_c, R_ext_side=R_ext_side(spec), paths=paths, ports=ports, port_r=dict(A=rs, B=rs, C=rc))
    return to_global(unary_union(pieces), leg, gx0), info


# ------------------------------------------------------------------------------------------------ assemble
def build(sc0, kind='linear', n_legs=1, R_turn=None, filt=None, r_f=None, tile=1500.0):
    """Plan the layout of a closed schedule and draw it, returning (fluid, sc, info) as assemble does."""
    sc, legs = plan_layout(sc0, kind, n_legs, R_turn)
    return assemble(sc, legs, filt=filt, r_f=r_f, tile=tile)


def junction_roundings(spec, units):
    """Round the wall corners where each head unit's stem leaves the center channel (the entrance of the head gaps) by
    spec.r_junction, which is about 11 um in Dinh Fig. 2B (DINH-7). Each corner on the right side (the head is mirrored
    afterward) gains a square of side r in the wall less the disc that rounds it. The rounding is drawing only, because
    the gap entrance term c_ent is external. The function returns the polygons to add to the fluid ([] when r_junction
    is 0)."""
    rj = float(getattr(spec, 'r_junction', 0.0) or 0.0)
    if rj <= 0.0:
        return []
    from shapely.geometry import Point as _Pt
    yc = spec.wc0 / 2.0
    out = []
    for u in units:
        for sgn in (-1.0, 1.0):
            xc = u['X'] + sgn * u['w'] / 2.0
            sq = box(min(xc, xc + sgn * rj), yc, max(xc, xc + sgn * rj), yc + rj)
            out.append(sq.difference(_Pt(xc + sgn * rj, yc + rj).buffer(rj, quad_segs=32)))
    return out


def assemble(sc, legs, filt=None, r_f=None, tile=1500.0):
    """Draw a planned (and possibly fail-safe re-planned) schedule with its legs, and return (fluid, sc, info)."""
    spec = sc.spec
    filt = spec.filter if filt is None else filt
    r_f = spec.r_fillet if r_f is None else r_f
    units = Gm.plan_head(sc)
    if units:
        Gm.solve_transition(sc, units[-1])
    for u in units:
        Gm.exact_unit(spec, u)
    # ---- the head (right side), mirrored, and the stem or inlet
    parts = []
    paths = []
    for u in units:
        p = Gm.unit_path(spec, u)
        paths.append(p)
        if u['kind'] == 'X':
            p.line(Gm.EPS)                    # overlap straight into the drop polygon (see shapes.transition_drop)
        parts.append(Gm.ribbon(p, u['w'], spec))
        if u['kind'] == 'X':
            parts.append(u['drop_poly'])
    if parts:
        parts.extend(junction_roundings(spec, units))          # none unless meander.junction_radius_um is set (DINH-7)
    right = unary_union(parts) if parts else Polygon()
    head = unary_union([right, affinity.scale(right, 1.0, -1.0, origin=(0, 0))]) if parts else Polygon()
    x_top = head_top_x(sc)
    if filt == 'dinh':
        inlet, inl = DI.inlet(spec)
    else:
        inlet, inl = inlet_local(spec, x_top, filt=bool(filt))
    head = unary_union([head, inlet] + ([] if filt else [box(inl['x_port_c'], -spec.wc0 / 2, 1.0, spec.wc0 / 2)]))
    # ---- legs
    turns = sc.info['turns']
    tw = {}
    for j, k in enumerate(turns):
        dY = legs[j + 1]['Y'] - legs[j]['Y']
        ccw = legs[j]['dirn'] * dY > 0                  # a left turn, whose inner side is side A (local +y)
        wi, fi = loop_profile(sc, k, inner=True)
        wo, fo = loop_profile(sc, k, inner=False)
        tw[k] = dict(A=(wi, fi) if ccw else (wo, fo), B=(wo, fo) if ccw else (wi, fi), ccw=ccw,
                     cyoff=dY / 2.0, th0=(-math.pi / 2.0 if dY > 0 else math.pi / 2.0))
    pieces = [head]
    part_names = ['inlet and head']
    gx = [0.0]
    loops = []
    for j, leg in enumerate(legs):
        def fpair(k):
            return (tw[k]['A'][1], tw[k]['B'][1])
        sw_start = fpair(legs[j - 1]['turn_after']) if j > 0 else (None, None)
        sw_end = fpair(leg['turn_after']) if leg['turn_after'] is not None else (None, None)
        loc = leg_linear_local(sc, leg, sw_start, sw_end)
        g = to_global(loc, leg, gx[j])
        pieces.append(g)
        part_names.append('leg %d' % (j + 1))
        if leg['turn_after'] is not None:
            Lloc = leg['p1'] - (leg['p0'] if leg['p0'] is not None else 0.0)
            gx_end = gx[j] + Lloc if leg['dirn'] > 0 else gx[j] - Lloc
            t_ = tw[leg['turn_after']]
            ccw, cy = t_['ccw'], leg['Y'] + t_['cyoff']
            lp = loop_pieces(sc, leg['turn_after'], gx_end, cy, ccw, t_['A'][1], t_['B'][1], th0=t_['th0'])
            loops.append(dict(k=leg['turn_after'], cx=gx_end, cy=cy, ccw=ccw, th0=t_['th0'], R=S.R_of(sc, leg['turn_after']),
                              wA=t_['A'][0], wB=t_['B'][0], fA=t_['A'][1], fB=t_['B'][1]))
            pieces.append(lp)
            part_names.append('U-turn %d' % (j + 1))
            gx.append(gx_end)
    outlet, outl = outlet_global(sc, legs[-1], gx[-1])
    pieces.append(outlet)
    part_names.append('outlet')
    fluid = unary_union(pieces)
    parts = list(zip(part_names, pieces))               # in flow order, so each part joins only its neighbors
    # The fillets go everywhere except over the posts of the linear legs, because the stadium posts are round already
    # and smoothing their polygonal caps would bias them by about 0.03 um, which the ladder amplifies.
    x_cut = sc.x[sc.NM] - spec.Lam_L / 2.0 - spec.g_p / 2.0 - 2.0
    # The head zone covers the head's own lateral extent and the inlet only, because in a paperclip legs 2 and 3 run
    # beside it.
    hz = head_halfwidth(sc) + 60.0
    x0z = min(inl['x_port_c'] - inl.get('port_r', spec.port_r_in), x_top) - 10.0
    big = box(x0z - 8 * r_f, -hz - 8 * r_f, x_cut + 8 * r_f, hz + 8 * r_f)
    headzone = box(x0z, -hz, x_cut, hz)
    head_s = Gm.smooth(fluid.intersection(big), r_f, tile).intersection(headzone)
    rest = fluid.difference(headzone)
    # The outlet zone and the loop junction edges are smooth by construction, except for the lead and collector T and
    # the port necks, so the outlet region is filleted separately.
    ob = outlet.bounds
    ozone = box(ob[0] - 10, ob[1] - 10, ob[2] + 10, ob[3] + 10)
    olarge = box(ob[0] - 10 - 8 * r_f, ob[1] - 10 - 8 * r_f, ob[2] + 10 + 8 * r_f, ob[3] + 10 + 8 * r_f)
    # The outlet is smoothed from the whole fluid rather than from what lies outside the head zone. When the outlet
    # reaches into the head zone (a two-leg paperclip ends beside its inlet, and m12's ring lead crossed the zone's
    # edge), smoothing the cut pieces rounded the cut itself, leaving 0.5 um corners and a step in the lead.
    out_s = Gm.smooth(fluid.intersection(olarge), r_f, 1e9).intersection(ozone).difference(headzone)
    rest = rest.difference(ozone)
    fluid = unary_union([head_s, rest, out_s])
    fluid = shapely.set_precision(fluid, 1e-5)
    info = dict(units=units, legs=legs, loops=loops, gx=gx, inlet=inl, outlet=outl, x_top=x_top,
                parts=[(n, outline(p)) for n, p in parts])
    return fluid, sc, info


def outline(g):
    """Return the filled outline of a drawn part, ignoring its holes (the posts)."""
    gs = list(g.geoms) if hasattr(g, 'geoms') else [g]
    return unary_union([Polygon(q.exterior) for q in gs if not q.is_empty and q.geom_type == 'Polygon'])


PORT_GAP_IN_OUT = 2000.0   # um, edge to edge, between the inlet port and every outlet port


def port_spacing(info, spec):
    """Check that the device's two ends are kept apart, with the inlet port at least PORT_GAP_IN_OUT (edge to edge)
    from every outlet port. A folded layout with an even number of legs ends beside its inlet, and with this rule it
    may still do so, although its outlet ports may not crowd the inlet port (0.3 mm for m10 drawn on two legs). Tubing
    or connectors sit on every port."""
    import math as _m
    ip = (float(info['inlet']['x_port_c']), 0.0)
    ir = float(info['inlet'].get('port_r', spec.port_r_in))
    o = info['outlet']
    pr = o.get('port_r') or {}
    gaps = {k: _m.hypot(float(x) - ip[0], float(y) - ip[1]) - ir - float(pr.get(k, spec.port_r))
            for k, (x, y) in o['ports'].items()}
    k, g = min(gaps.items(), key=lambda kv: kv[1])
    return dict(ok=bool(g >= PORT_GAP_IN_OUT - 1e-6), inlet_to_nearest_outlet_port_um=round(float(g), 1), port=k,
                required_um=PORT_GAP_IN_OUT)


def part_clearance(parts, floor):
    """Check that no parts intersect, that is, that every drawn part (the inlet with the head, each leg, each U-turn,
    and the outlet) keeps at least the process floor from every part it is not joined to.

    The parts are listed in flow order, and each joins only its neighbors (the inlet and head join leg 1, a U-turn
    joins the legs on either side, and the outlet joins the last leg). The check measures the drawn outlines before
    the fillets, so a crossing, an overlap, or a wall thinner than the floor between unrelated parts is named
    rather than only inferred from the hole count or the DRC."""
    found = []
    for i in range(len(parts)):
        for j in range(i + 2, len(parts)):
            (na, a), (nb, b) = parts[i], parts[j]
            ba, bb = a.bounds, b.bounds
            if (ba[0] - floor > bb[2] or bb[0] - floor > ba[2] or ba[1] - floor > bb[3] or bb[1] - floor > ba[3]):
                continue
            found.append((float(a.distance(b)), na, nb))
    found.sort()
    ok = all(d >= floor - 0.1 for d, _, _ in found)
    return dict(ok=bool(ok), floor_um=float(floor), parts=[n for n, _ in parts], pairs_within_reach=len(found),
                closest=[dict(gap_um=round(d, 3), parts=[x, y]) for d, x, y in found[:3]])
