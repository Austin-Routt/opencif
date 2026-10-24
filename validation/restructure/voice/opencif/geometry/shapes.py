"""opencif.geometry.shapes: the meander head, with its unit planning, exact unit lengths, transition unit, unit
centerlines, and fillet operator, in leg-local coordinates (opencif.geometry.layout assembles and draws the device).

The coordinate x runs along the flow (gap 1 at x = 0) and y is lateral, so the right side is y > 0 and the left side
is its mirror. Every channel is a centerline of straight runs and circular arcs, buffered to its width with flat
caps, and every remaining corner (the T-junctions, stub mouths, and channel ends) is filleted by a closing and an
opening of radius r_f. On a square path (turns = 'square', DINH-1), each quarter turn is a corner instead of an arc
(Path.turn).

Meander unit k (row k+1, with its gap at X_k) on the right side has the topology of Dinh 2024 Fig. 2A and Xia 2016
Fig. 1b. Its out line starts on the trunk face (the gap mouth, of width w) and runs outward (+y) through the
T-junction J_k at (X_k, yJ). The two kinds of unit then continue as follows:

    taper unit   on to the reach R, a 180-degree U-turn of centerline radius s/2 (downstream), back (-y) to the
                 stop-short line, a quarter arc of radius s/2, and the drop (+x, at y = yJ) onto the next unit's out
                 line at J_k+1
    arm unit     on to its lane, a quarter arc turning upstream (-x), then n lane pairs (rise -x, descend +x) joined
                 by 180-degree arcs at the top (and at the bottom between pairs), the last lane coming down past the
                 out row to the return row, a quarter arc, back (-y), a quarter arc, and the drop

Arms nest outward, so the lanes of arm k+1 lie outside those of arm k, and all turns are arcs of radius
(lane pitch)/2. The transition unit (the last meander row) is a taper unit whose drop narrows smoothly from GS to the
next row's side width and continues as the straight side channel over the wall blocks of the rectilinear rows (at
least one). The stadium posts start at the first linear gap.
"""
import math
import numpy as np
from shapely.geometry import LineString, Polygon, box
from shapely.ops import unary_union
from opencif.physics import schedule as S

ARC_DEG = 1.5          # arc discretization (degrees per segment)
EPS = 1.0              # overlap of a stub into the trunk (um)


class Path:
    """A centerline of straight runs, arcs, and (on a square path) corners, starting at (x, y) with a heading in
    degrees."""

    def __init__(self, x, y, heading_deg, square=False):
        self.p = [(float(x), float(y))]
        self.th = math.radians(heading_deg)
        self.arcs = []                         # (centerline radius, |angle| rad) of every arc, for the bend correction
        self.square = bool(square)             # quarter turns drawn as corners by turn() (DINH-1)
        self.corners = []                      # the index into p of every square corner's centerline point

    @property
    def xy(self):
        """The current end point of the path."""
        return self.p[-1]

    def line(self, L):
        """Extend the path straight ahead by L (um), refusing a negative length."""
        if L < -1e-9:
            raise ValueError('negative line length %.6f' % L)
        if L > 1e-12:
            x, y = self.xy
            self.p.append((x + L * math.cos(self.th), y + L * math.sin(self.th)))
        return self

    def arc(self, r, ang_deg, seg_max=None):
        """Turn by ang_deg degrees along an arc of centerline radius r, to the left (CCW) for a positive angle and to
        the right (CW) for a negative one. The optional seg_max (um) caps the chord length of large-radius arcs,
        because the 1.5-degree chords of a 2 mm arc are 51 um long, and the opening test of the DRC reads the
        chord-and-round-join outline as a thin crescent."""
        x, y = self.xy
        th = self.th
        sgn = 1.0 if ang_deg > 0 else -1.0
        cx = x - sgn * r * math.sin(th)
        cy = y + sgn * r * math.cos(th)
        a0 = math.atan2(y - cy, x - cx)
        ang = math.radians(ang_deg)
        self.arcs.append((float(r), abs(ang)))
        n = max(4, int(math.ceil(abs(ang_deg) / ARC_DEG)))
        if seg_max:
            n = max(n, int(math.ceil(r * abs(ang) / seg_max)))
        for j in range(1, n + 1):
            a = a0 + ang * j / n
            self.p.append((cx + r * math.cos(a), cy + r * math.sin(a)))
        self.th = th + ang
        return self

    def turn(self, r, ang_deg):
        """Make a quarter turn of centerline radius r, which is an arc, or, on a square path, a corner whose centerline
        runs r to the corner point and r on (DINH-1). The 1-D length of a corner is then the centerline through the
        corner point, and its 2-D excess is schedule.c_corner."""
        if self.square and abs(abs(ang_deg) - 90.0) < 1e-9:
            self.line(r)
            self.corners.append(len(self.p) - 1)
            self.th += math.radians(ang_deg)
            self.line(r)
            return self
        return self.arc(r, ang_deg)

    def length(self):
        """Return the polygonal length of the path (um)."""
        P = np.asarray(self.p)
        return float(np.sum(np.hypot(*np.diff(P, axis=0).T)))


def qturn_len(spec, r):
    """Return the centerline length of a quarter turn of radius r, which is pi r / 2 for an arc and 2 r for a square
    corner (through its point)."""
    return 2.0 * r if getattr(spec, 'turns', 'arc') == 'square' else math.pi * r / 2.0


def corner_cut(C, u, v, w, r_out, r_in):
    """Return the region to remove at the outer corner (rounding r_out) and the region to add at the inner corner
    (fillet r_in) of a square corner with centerline point C and unit headings u in and v out (DINH-1). With
    d = v - u, the outer corner is C - (w/2) d and the inner one is C + (w/2) d, and each rounding is the square
    spanned by -u and v at that corner less the disc of radius r centered at corner + r d."""
    from shapely.geometry import Point
    C, u, v = np.asarray(C, float), np.asarray(u, float), np.asarray(v, float)
    d = v - u
    out = []
    for j_, (P0, r) in enumerate(((C - 0.5 * w * d, r_out), (C + 0.5 * w * d, r_in))):
        if r <= 0:
            out.append(None)
            continue
        if j_ == 0:
            # The cut reaches 1e-3 um beyond both outer edges (outside the ribbon, where it removes nothing), so that
            # the miter vertex, which the buffer places within the last bits of P0, never survives as a zero-width
            # spike. The closing in smooth() grew such a spike and partly refilled the rounding (rc7 code review R13).
            e_ = 1e-3
            Q0 = P0 + e_ * u - e_ * v
            sq = Polygon([tuple(Q0), tuple(Q0 - (r + e_) * u), tuple(Q0 - (r + e_) * u + (r + e_) * v),
                          tuple(Q0 + (r + e_) * v)])
        else:
            sq = Polygon([tuple(P0), tuple(P0 - r * u), tuple(P0 - r * u + r * v), tuple(P0 + r * v)])
        out.append(sq.difference(Point(tuple(P0 + r * d)).buffer(r, quad_segs=32)))
    return out[0], out[1]


def corner_round_area(spec):
    """Return the net area (um^2) that the roundings give one square corner as ribbon() draws them, that is, the inner
    fillet added less the outer rounding removed (DINH-1). schedule.c_corner is the excess of the 2-D solve over the
    centerline through the corner point on the rounded corner, so it already counts this area. Therefore a measurement
    that converts a unit's drawn area into length (measure) subtracts it, or the corner would be counted twice."""
    ro = float(getattr(spec, 'r_corner_out', 0.0) or 0.0)
    ri = float(getattr(spec, 'r_corner_in', 0.0) or 0.0)
    if ro <= 0.0 and ri <= 0.0:
        return 0.0
    w_ = 4.0 * max(ro, ri)
    o_, i_ = corner_cut((0.0, 0.0), (0.0, 1.0), (1.0, 0.0), w_, ro, ri)
    # The outer cut reaches 1e-3 um beyond the ribbon's edges, so only its part inside the ribbon's corner square
    # counts.
    inside = Polygon([(-w_ / 2.0, -w_ / 2.0), (w_ / 2.0, -w_ / 2.0), (w_ / 2.0, w_ / 2.0), (-w_ / 2.0, w_ / 2.0)])
    return (i_.area if i_ is not None else 0.0) - (o_.intersection(inside).area if o_ is not None else 0.0)


def ribbon(path, w, spec=None):
    """Buffer the path to a channel of width w with flat caps, rounding the outer and filleting the inner side of each
    square corner with the spec's radii."""
    if not getattr(path, 'corners', None):
        return LineString(path.p).buffer(w / 2.0, cap_style='flat', join_style='round', quad_segs=8)
    # A path with square corners is buffered with mitered joins, which give square outer corners and leave the arcs'
    # 1.5-degree joins unaffected. Each corner's outer rounding is then removed and its inner fillet added.
    g = LineString(path.p).buffer(w / 2.0, cap_style='flat', join_style='mitre', mitre_limit=10.0)
    r_out = float(getattr(spec, 'r_corner_out', 0.0) or 0.0) if spec is not None else 0.0
    r_in = float(getattr(spec, 'r_corner_in', 0.0) or 0.0) if spec is not None else 0.0
    P = np.asarray(path.p, float)
    for ci in path.corners:
        u = P[ci] - P[ci - 1]; u /= np.hypot(*u)
        v = P[ci + 1] - P[ci]; v /= np.hypot(*v)
        cut, add = corner_cut(P[ci], u, v, w, r_out, r_in)
        if cut is not None:
            g = g.difference(cut)
        if add is not None:
            g = g.union(add)
    return g


# ----------------------------------------------------------------------------------------------------------
# Head planning: the reaches, lanes, and tops
# ----------------------------------------------------------------------------------------------------------
def taper_fixed(spec, ell, w=None):
    """Return the fixed part fixed(ell) of a taper unit's length, which is 2 (R - yJ) + fixed(ell)."""
    w = spec.GS if w is None else w
    s = w + spec.wall_M
    rq = s / 2.0
    return math.pi * s / 2.0 - rq + qturn_len(spec, rq) + (ell - s - rq)


def arm_length(spec, X, ell, R, n, t, w):
    """Return the centerline length J_k -> J_k+1 of an arm with n lane pairs, its first lane at y = R and all its tops
    at x = t."""
    sa = w + spec.wall_M
    ra = rr_ = rq = sa / 2.0
    yJ = spec.yJ(w)
    xb = X - sa / 2.0
    y_last = R + (2 * n - 1) * sa
    L = (R - ra - yJ) + qturn_len(spec, ra)
    for m in range(1, n + 1):
        bot_rise = (X - ra) if m == 1 else xb
        bot_desc = (X + sa - rr_) if m == n else xb
        L += (bot_rise - t) + (bot_desc - t) + math.pi * sa / 2.0          # two lanes + top U-turn
        if m < n:
            L += math.pi * sa / 2.0                                         # bottom U-turn
    L += qturn_len(spec, rr_) + (y_last - rr_ - yJ - rq) + qturn_len(spec, rq) + (ell - sa - rq)
    return L


def solve_arm_top(spec, X, ell, R, n, w, L_req):
    """Return the top t at which arm_length equals L_req; the length is linear in t, with slope -2n."""
    L0 = arm_length(spec, X, ell, R, n, 0.0, w)
    return (L0 - L_req) / (2.0 * n)


def plan_head(sc, greedy=True):
    """Return a geometry descriptor for each meander row. For an arm, it holds the lane positions and the top, which
    the 'pack' head sets by greedy skyline flattening."""
    spec = sc.spec
    units = []
    A = sc.A
    w_arm = spec.w_arm
    sa = w_arm + spec.wall_M
    # ---- arms
    n = [1] * A

    def layout_arms(n):
        R = []
        r = spec.R0
        for k in range(A):
            R.append(r)
            r += 2 * n[k] * sa
        tops = [solve_arm_top(spec, sc.x[k], sc.ell[k], R[k], n[k], w_arm, sc.L_req[k]) for k in range(A)]
        return R, tops

    R, tops = layout_arms(n)
    lanes_total = int(math.floor((spec.y_max - w_arm / 2.0 - spec.R0) / sa)) + 1
    if A and 2 * sum(n) > lanes_total:
        raise RuntimeError('arms do not fit laterally: %d lanes > %d' % (2 * sum(n), lanes_total))
    if spec.head == 'pack' and greedy and A:
        # Skyline flattening gives spare lane pairs to the tallest arm while that lowers the head top.
        while 2 * sum(n) + 2 <= lanes_total:
            k = int(np.argmin(tops))
            n2 = list(n); n2[k] += 1
            R2, t2 = layout_arms(n2)
            ok = all(t2[j] < (sc.x[j] - sa / 2.0) - sa for j in range(A))
            if ok and min(t2) > min(tops) + 1e-9:
                n, R, tops = n2, R2, t2
            else:
                break
    for k in range(sc.NM):
        kind = sc.kind[k]
        X, ell = float(sc.x[k]), float(sc.ell[k])
        if kind == 'A':
            units.append(dict(k=k, kind='A', X=X, ell=ell, w=w_arm, R=R[k], n=n[k], t=tops[k], L=sc.L_req[k]))
        elif k < sc.NM - 1:
            L = sc.L_req[k]
            w = float(sc.wM[k])
            Rk = spec.yJ(w) + (L - taper_fixed(spec, ell, w)) / 2.0
            units.append(dict(k=k, kind='T', X=X, ell=ell, w=w, R=Rk, L=L))
        else:
            units.append(dict(k=k, kind='X', X=X, ell=ell, w=spec.GS, R=None, L=None))   # the transition, solved below
    return units


def bend_excess(spec, path, w):
    """Return the extra 1-D length that the path's arcs and square corners add under the bend model. The arcs' term is
    schedule.arc_excess_len, the exact Stokes curved-duct correction less the junctions, and the corners' term is
    schedule.c_corner per corner, against the centerline through the corner point (DINH-1). A unit's resistance is
    r(w) (geometric length + bend_excess)."""
    arcs = float(sum(S.arc_excess_len(spec, r, a, w) for r, a in path.arcs))
    if getattr(path, 'corners', None):
        arcs += len(path.corners) * S.c_corner(spec, w) * w
    return arcs


def junction_excess(spec, u, w=None):
    """Return the length c_t w (c_t < 0, schedule.c_t) that the T-junction at the unit's end adds to the unit, where
    its drop turns into the next unit's stem with the gap's mouth below. Taper and arm units end in such a junction,
    whereas the transition unit ends in the straight side channel."""
    if u['kind'] not in ('T', 'A'):
        return 0.0
    w = u['w'] if w is None else w
    return float(S.c_t(spec, w) * w)


def unit_eff_length(spec, u):
    """Return the unit's resistance length J_k -> J_k+1, which is its drawn (polygonal) length plus the bend correction
    of its arcs plus its end T-junction."""
    p = unit_path(spec, u, from_J=True)
    return p.length() + bend_excess(spec, p, u['w']) + junction_excess(spec, u)


def exact_unit(spec, u, iters=4):
    """Make the unit's resistance length J_k -> J_k+1 (the drawn polygonal length plus the bend correction) equal the
    required length exactly, by adjusting the taper reach R (slope 2) or the arm top t (slope -2n) while the arcs stay
    unchanged. A taper unit cannot be shorter than its smallest loop (reach R = yJ + rq), so such a unit is refused
    with the numbers."""
    if u['kind'] in ('T', 'X') and u.get('R') is not None:
        w = u['w']
        R_min = spec.yJ(w) + (w + spec.wall_M) / 2.0
        if u['R'] < R_min:
            u['R'] = R_min                    # start the iteration from a drawable loop
    for _ in range(iters):
        dL = unit_eff_length(spec, u) - u['L']
        if u['kind'] in ('T', 'X'):
            u['R'] -= dL / 2.0
            R_min = spec.yJ(u['w']) + (u['w'] + spec.wall_M) / 2.0
            if u['R'] < R_min:                # it cannot be drawn shorter, so stop at the smallest loop (refused below)
                u['R'] = R_min - 1e-3
                break
        else:
            u['t'] += dL / (2.0 * u['n'])
    if u['kind'] in ('T', 'X') and u.get('R') is not None:
        w = u['w']
        R_min = spec.yJ(w) + (w + spec.wall_M) / 2.0
        if u['R'] < R_min - 1e-6:
            raise ValueError('head unit %d needs a resistance length of %.1f um of %.1f um channel, shorter than its smallest '
                             'drawable loop: widen meander.channel_um, narrow meander.wall_um, or end the head at a narrower '
                             'seam (meander.seam_side_width_um)' % (u['k'] + 1, u['L'], w))
    p = unit_path(spec, u, from_J=True)
    u['L_drawn'] = p.length()
    u['L_bend_excess'] = bend_excess(spec, p, u['w'])
    u['L_junction_excess'] = junction_excess(spec, u)
    return u


def smoothstep(t):
    """Return the cubic smoothstep t^2 (3 - 2 t), which rises from 0 to 1 with zero slope at both ends."""
    return t * t * (3 - 2 * t)


def transition_drop(spec, X, ell, w_next, n=200):
    """Draw the transition unit's drop, which runs in x from X+s+rq to X+ell with its inner edge fixed at wc0/2+l_g
    while its width goes from GS to w_next. Return (polygon, resistance in SI units, x0, x1)."""
    s = spec.s
    rq = s / 2.0
    x0, x1 = X + s + rq, X + ell
    y_in = spec.wc0 / 2.0 + spec.l_g
    xs = np.linspace(x0, x1, n + 1)
    ts = (xs - x0) / (x1 - x0)
    ws = spec.GS + (w_next - spec.GS) * smoothstep(ts)
    up = [(x, y_in + w) for x, w in zip(xs, ws)]
    # The drop starts exactly where the unit's quarter arc ends, and the arc's ribbon is extended straight by EPS into
    # it (layout.assemble). A drop polygon that reached back EPS under the arc left a step against the arc's outer
    # edge, which the fillet closing turned into a pocket of about 10 um2 beside the bend. In the wbc example, the
    # transition unit read +7e-4 long and the seam row moved 1%.
    poly = Polygon([(x0, y_in), (x1, y_in)] + up[::-1] + [(x0, y_in + spec.GS)])
    rs = np.array([S.rr(w, spec.h) for w in ws])
    R = float(np.trapezoid(rs, xs) * S.UM)
    return poly, R, x0, x1


def solve_transition(sc, u):
    """Solve for the reach R of the transition unit at which the resistance of its loop plus its drop equals R_s(NM)."""
    spec = sc.spec
    k = u['k']
    target = sc.Rs[k]
    w_next = sc.ws[k + 1]
    if getattr(spec, 'transition', 'exact') == 'paper' and sc.kind[k + 1] == 'R':
        # In the papers' transition, the drop continues at the rectilinear rows' drawn width (GS by default).
        w_next = spec.ws_rect if getattr(spec, 'ws_rect', None) is not None else spec.GS
    poly, R_drop, x0, x1 = transition_drop(spec, u['X'], u['ell'], w_next)
    s, rq = spec.s, spec.s / 2.0
    yJ = spec.yJ(spec.GS)
    # The loop's length is loop_fixed + 2 (R - yJ), and both arcs (radius s/2 = rq) carry the bend correction
    # (S.arc_excess_len).
    if getattr(spec, 'turns', 'arc') == 'square':          # the quarter turn is a corner (DINH-1)
        q_ex = S.c_corner(spec, spec.GS) * spec.GS
    else:
        q_ex = S.arc_excess_len(spec, rq, math.pi / 2.0, spec.GS)
    loop_fixed = (math.pi * s / 2.0 - rq + qturn_len(spec, rq) + S.arc_excess_len(spec, s / 2.0, math.pi, spec.GS)
                  + q_ex)
    L_loop = (target - R_drop) / (S.rr(spec.GS, spec.h) * S.UM)
    R = yJ + (L_loop - loop_fixed) / 2.0
    u.update(R=R, L=L_loop, R_drop=R_drop, w_next=w_next, drop_poly=poly)
    return u


# ----------------------------------------------------------------------------------------------------------
# Drawing
# ----------------------------------------------------------------------------------------------------------
def unit_path(spec, u, from_J=False):
    """Return the centerline of unit u, starting inside the trunk (for drawing) or at its T-junction J_k (for
    measurement)."""
    X, w = u['X'], u['w']
    s = w + spec.wall_M
    rq = s / 2.0
    yJ = spec.yJ(w)
    y0 = yJ if from_J else spec.wc0 / 2.0 - EPS
    p = Path(X, y0, 90.0, square=getattr(spec, 'turns', 'arc') == 'square')
    if u['kind'] in ('T', 'X'):
        R = u['R']
        p.line(R - y0)
        p.arc(s / 2.0, -180.0)                 # outer U-turn, downstream
        p.line(R - yJ - rq)
        p.turn(rq, 90.0)                       # turn to +x (a corner on a square path)
        if u['kind'] == 'T':
            p.line(u['ell'] - s - rq)          # drop onto the next out line
        return p
    # an arm unit
    sa = s
    ra = rr_ = sa / 2.0
    R, n, t = u['R'], u['n'], u['t']
    xb = X - sa / 2.0
    p.line(R - ra - y0)
    p.turn(ra, 90.0)                           # turn upstream (-x)
    for m in range(1, n + 1):
        bot_rise = (X - ra) if m == 1 else xb
        p.line(bot_rise - t)                   # rise to the top
        p.arc(sa / 2.0, -180.0)                # top U-turn (bulges upstream), now heading +x
        bot_desc = (X + sa - rr_) if m == n else xb
        p.line(bot_desc - t)                   # descend
        if m < n:
            p.arc(sa / 2.0, 180.0)             # bottom U-turn, heading -x again
    y_last = R + (2 * n - 1) * sa
    p.turn(rr_, -90.0)                         # +x -> -y
    p.line(y_last - rr_ - yJ - rq)             # return toward the trunk
    p.turn(rq, 90.0)                           # -y -> +x
    p.line(u['ell'] - sa - rq)                 # drop
    return p


def smooth(g, r, tile, margin=None):
    """Apply a closing and then an opening of radius r, tile by tile (in x-slabs), and reassemble the tiles."""
    if r <= 0:
        return g
    margin = margin if margin is not None else 6 * r
    x0, y0, x1, y1 = g.bounds
    out = []
    xa = x0
    while xa < x1:
        xb = min(xa + tile, x1)
        win = box(xa - margin, y0 - 10, xb + margin, y1 + 10)
        piece = g.intersection(win)
        piece = piece.buffer(r, quad_segs=16).buffer(-2 * r, quad_segs=16).buffer(r, quad_segs=16)
        out.append(piece.intersection(box(xa, y0 - 20, xb, y1 + 20)))
        xa = xb
    res = unary_union(out)
    return res
