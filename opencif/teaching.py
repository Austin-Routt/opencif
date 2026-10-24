"""opencif.teaching: the whole path from the two targets to an audited photomask, for the simplest CIF device.

The device is a straight ladder. It has a center channel of constant width w_c, a row of posts on each side of it with
a gap of width G between neighboring posts, and a side channel beyond each row of posts that collects what the gaps
skim off. There is no meander head, and the inlet is sealed, so all the sample enters the center channel. This
module does for that device what the full generator (opencif.generate) does, in about 300 lines meant to be read
beside Appendix A, and tests/test_teaching.py checks that its schedule is the generator's, number for number.

    from opencif import teaching as T
    lad = T.design(D_um=6.0, F=0.8, depth_um=80, gap_um=20.5, post_um=25.3, wall_um=59.5, centre_um=100)
    fluid = T.draw(lad)                       # a shapely polygon, in micrometers
    T.write_photomask(fluid, 'ladder.dxf')    # the DXF that a mask shop reads
    a = T.audit_photomask('ladder.dxf', depth_um=80)   # reads only the file (and the depth, which it cannot see)

Several parts are left to the full generator. These are a meander head (which lets a coarse process seal gap 1), a
tapered center channel, folded layouts, inlet filters, outlet leads and ports, rounded corners, the design-rule check,
the declared departures of the published devices, and the physics auditor's full set of checks (opencif.checks.auditor).
The drawing here has sharp corners, and its channels end at the chip's edge, so it is a model to learn from rather than
a mask to fabricate.

The symbols follow Appendix A:

    D               the cutoff diameter
    F               the filtrate fraction
    delta           the lamina width
    f               the fraction of the center flow that each gap removes
    P(i)            the center flow after gap i, in units of the inlet flow
    S(i)            one side channel's flow after gap i, in units of the inlet flow
    R_c, R_s, R_g   the resistances of a center segment, a side segment, and a gap
"""
import math
from dataclasses import dataclass

import numpy as np

from opencif.physics import kernels as K

UM = 1e-6        # meters per micrometer
# Every resistance is in Pa s m^-3 at water's viscosity MU, and the schedule does not depend on the viscosity.
MU = 1e-3        # water's viscosity, Pa s


def r_um(w_um, h_um):
    """Return the resistance per meter of a w x h duct, with w and h in micrometers (App. A Result 14, eq. 88).

    The result is in Pa s m^-4, so a segment's resistance is this value times the segment's length in meters."""
    return K.r(w_um * UM, h_um * UM, MU)


def gap_resistance(G_um, wall_um, h_um, c_ent=1.0):
    """Return the resistance of one gap, R_g = r(G) (l_g + c_ent D_h), in Pa s m^-3 (spec v5 sec. 14.4).

    A gap is a short duct of width G through the wall of thickness l_g, and an entrance allowance of c_ent hydraulic
    diameters D_h is added to its length."""
    D_h = 2.0 * G_um * h_um / (G_um + h_um)
    return r_um(G_um, h_um) * (wall_um + c_ent * D_h) * UM


@dataclass
class Ladder:
    """A designed straight ladder, which holds its inputs, the march, the closure, and the widths.

    Lengths are in micrometers, flows are in units of the inlet flow, and resistances are in Pa s m^-3."""
    D_um: float
    F_target: float
    h: float            # the channel depth
    G: float            # the gap width, along the flow
    g_p: float          # the post length along the flow
    l_g: float          # the post width across the flow, that is, the wall the skimmed fluid crosses
    w_c: float          # the center channel width
    L_out: float        # the length of the three outlet channels after the last gap
    kappa: float        # the cutoff prefactor kappa_g
    delta: float        # the lamina width delta* = D / kappa_g
    f: float            # the fraction per gap and side, f = g(delta*/w_c)
    N: int              # the number of gaps on each side
    P: np.ndarray       # the center flow after gap i, i = 1..N (P[0] is after gap 1)
    S: np.ndarray       # one side channel's flow after gap i
    Qg: np.ndarray      # the flow through gap i (one side)
    Rc: float           # a center segment, gap to gap
    Rg: float           # one gap
    Rs: np.ndarray      # side segment i -> i+1, i = 1..N-1 (the exact closure)
    Rs_pe: np.ndarray   # the same from the PE closure, for comparison
    Rc_out: float       # the center outlet channel
    Rs_out: float       # a side outlet channel
    ws: np.ndarray      # the side channel width of segment i -> i+1
    ws_out: float       # the side outlet channel width
    F: float            # the march's filtrate fraction, 1 - P(N)

    @property
    def pitch(self):
        """Return the pitch Lambda = G + g_p, from gap center to gap center, which spans one gap and one post."""
        return self.G + self.g_p


def design(D_um, F, depth_um, gap_um, post_um, wall_um, centre_um, outlet_um=200.0, c_ent=1.0):
    """Design a straight ladder from the two targets (the cutoff D and the filtrate fraction F) and the process.

    The design runs in six numbered steps, from the cutoff law to the side-channel widths, and returns a Ladder.

    D_um, F      the targets, the cutoff D in micrometers and the filtrate fraction F
    depth_um     the channel depth h (every channel has the same depth)
    gap_um       the gap width G between two posts
    post_um      the post's length along the flow, g_p (the pitch is G + g_p)
    wall_um      the post's width across the flow, l_g (the wall the skimmed fluid crosses)
    centre_um    the center channel width w_c (constant)
    outlet_um    the length of the three outlet channels after the last gap
    c_ent        the gaps' entrance allowance in hydraulic diameters (see gap_resistance)"""
    h, G, g_p, l_g, w_c = float(depth_um), float(gap_um), float(post_um), float(wall_um), float(centre_um)
    pitch = G + g_p
    # Step 1 is the cutoff law (App. A Result 5; App. B C3b). A particle is skimmed off when its center lies inside
    # the wall lamina, and the lamina that sets a cutoff D is delta* = D / kappa_g, where kappa_g = 4.9138 with the
    # cubic kernel.
    kappa = K.Pairing().kappa
    delta = D_um / kappa
    # Step 2 is the skim kernel (App. A Result 1). A lamina of width delta next to one wall of a channel of width w_c
    # carries the fraction f = g(x) = 3x^2 - 2x^3 of the channel's flow, where x = delta/w_c. Here w_c is constant,
    # so every gap removes the same fraction.
    x = delta / w_c
    if not 0.0 < x < 0.5:
        raise ValueError('the lamina %.3g um is not inside the half channel %.3g um' % (delta, w_c / 2))
    f = K.g_cubic(x)
    # Step 3 is the march (App. A eq. march). After gap i the center keeps P(i) = (1 - 2f) P(i-1), because there is a
    # gap on each side, and each side channel carries S(i) = S(i-1) + f P(i-1). The inlet is sealed, so P(0) = 1 and
    # S(0) = 0. The march stops at the first gap where the filtrate fraction F = 1 - P(N) (eq. Frec) reaches the
    # target, which gives N = ceil(ln(1-F)/ln(1-2f)) (eq. Nconst).
    P, S, Qg = [], [], []
    P_prev, S_prev = 1.0, 0.0
    while True:
        q = f * P_prev
        P_prev, S_prev = P_prev * (1.0 - 2.0 * f), S_prev + q
        P.append(P_prev)
        S.append(S_prev)
        Qg.append(q)
        if 1.0 - P_prev >= F:
            break
    P, S, Qg = np.array(P), np.array(S), np.array(Qg)
    N = len(P)
    # Step 4 computes the resistances (App. A Result 14) of a center segment one pitch long and of one gap.
    Rc = r_um(w_c, h) * pitch * UM
    Rg = gap_resistance(G, l_g, h, c_ent)
    # Step 5 is the closure. Between gaps i and i+1 the side channel and the center channel must lose the same
    # pressure, less the change of the pressure drop across the gaps (R_g Q_g(i) at gap i and R_g Q_g(i+1) at gap i+1).
    # Walking along the center channel, the pressures are
    #     p_c(i+1) = p_c(i) - R_c P(i)       in the center channel at gap i+1
    #     p_s(i)   = p_c(i) - R_g Q_g(i)     in the side channel beside gap i
    # so a side segment that carries S(i) needs R_s(i) = [p_s(i) - p_s(i+1)] / S(i). Without the gap terms this is
    # the PE closure of App. A (eq. closure), R_s = R_c P/S.
    p_c = np.zeros(N)
    for i in range(N - 1):
        p_c[i + 1] = p_c[i] - Rc * P[i]
    p_s = p_c - Rg * Qg
    Rs = (p_s[:-1] - p_s[1:]) / S[:-1]
    Rs_pe = Rc * P[:-1] / S[:-1]
    # At the outlet, three straight channels of length L_out end at the same pressure (0), and the side one must carry
    # S(N).
    Rc_out = r_um(w_c, h) * outlet_um * UM
    p_out = p_c[N - 1] - Rc_out * P[N - 1]
    Rs_out = (p_s[N - 1] - p_out) / S[N - 1]
    # Step 6 finds the widths. Each resistance belongs to exactly one width (App. A Result 14, "Monotonicity"), so
    # inverting the duct law gives the width of every side segment and of the side outlet.
    ws = np.array([K.r_inverse(R / (pitch * UM), h * UM, MU) / UM for R in Rs])
    ws_out = K.r_inverse(Rs_out / (outlet_um * UM), h * UM, MU) / UM
    return Ladder(D_um=float(D_um), F_target=float(F), h=h, G=G, g_p=g_p, l_g=l_g, w_c=w_c, L_out=float(outlet_um),
                  kappa=kappa, delta=delta, f=f, N=N, P=P, S=S, Qg=Qg, Rc=Rc, Rg=Rg, Rs=Rs, Rs_pe=Rs_pe, Rc_out=Rc_out,
                  Rs_out=Rs_out, ws=ws, ws_out=ws_out, F=float(1.0 - P[-1]))


# ---------------------------------------------------------------------------------------------------------------------
# drawing
# ---------------------------------------------------------------------------------------------------------------------

def draw(lad, inlet_um=300.0):
    """Return the fluid region of the ladder as one polygon in micrometers. The flow runs along +x, gap i is centered
    at x = (i-1) Lambda, the center channel's axis is y = 0, side A lies at y > 0, and side B is its mirror image.

    The polygon is the union of rectangles, which are the center channel (from inlet_um before gap 1's upstream edge to
    L_out after gap N), each gap's opening through the wall (G wide and l_g long), each side segment (from gap i to gap
    i+1 at width w_s(i)), and the side outlet channels. The first side segment starts at gap 1's upstream edge, so the
    side channel is sealed upstream. What is left between two gaps is a post, which is a hole in the polygon and is
    written to the DXF's ISLANDS layer."""
    from shapely.geometry import box
    from shapely.ops import unary_union
    lam, hw = lad.pitch, lad.w_c / 2.0
    xg = np.arange(lad.N) * lam                                      # the gap centers
    x_end = xg[-1] + lad.L_out
    y_in = hw + lad.l_g                                              # the side channel's inner edge
    parts = [box(-inlet_um - lad.G / 2.0, -hw, x_end, hw)]
    for sg in (1.0, -1.0):
        rects = []
        for i in range(lad.N):
            rects.append((xg[i] - lad.G / 2.0, hw, xg[i] + lad.G / 2.0, y_in))       # gap i
            if i < lad.N - 1:
                x0 = xg[i] - (lad.G / 2.0 if i == 0 else 0.0)
                rects.append((x0, y_in, xg[i + 1], y_in + lad.ws[i]))                # side segment i -> i+1
        rects.append((xg[-1], y_in, x_end, y_in + lad.ws_out))                       # the side outlet channel
        for x0, y0, x1, y1 in rects:
            parts.append(box(x0, min(sg * y0, sg * y1), x1, max(sg * y0, sg * y1)))
    return unary_union(parts)


def write_photomask(fluid, path):
    """Write the fluid polygon to a DXF file in micrometers with the generator's own writer (opencif.geometry.dxf).

    The layer CHANNELS holds the fluid's outline and the layer ISLANDS holds the posts. The function returns the
    number of outlines written."""
    from opencif.geometry import dxf
    return dxf.write_dxf(fluid, path)


# ---------------------------------------------------------------------------------------------------------------------
# the blind audit: everything from the photomask file
# ---------------------------------------------------------------------------------------------------------------------

def _cut(fluid, p0, p1):
    """Return the pieces of the straight cut from p0 to p1 through the fluid, as a sorted list [(start, end), ...] of
    distances from p0."""
    from shapely.geometry import LineString
    seg = fluid.intersection(LineString([p0, p1]))
    out = []
    for g in getattr(seg, 'geoms', [seg]):
        if g.is_empty or g.length == 0:
            continue
        cs = list(g.coords)
        d = [math.hypot(c[0] - p0[0], c[1] - p0[1]) for c in cs]
        out.append((min(d), max(d)))
    return sorted(out)


def audit_photomask(path, depth_um, c_ent=1.0, D_um=None):
    """Read a straight-ladder photomask back and work out what it does, knowing nothing of its design.

    From the file alone, the function measures the posts (the ISLANDS outlines, whose extent across the flow is the wall
    l_g), the gaps (the openings along a line through the middle of each wall), the center width (a cut across the
    channel between two gaps), each side segment's width (a cut beyond the wall), the side outlets' widths, and the
    outlets' common length; the center outlet takes the center width at the last gap. The depth cannot be seen on a mask
    and is an input. From these measurements it computes the resistances (App. A Result 14), the Kirchhoff network of
    the drawn device, each gap's emergent fraction f_i, its realized cutoff D_i = kappa_g x(f_i) w_c (App. A Results 2
    and 5), and the filtrate fraction F. It returns a dict, which also holds the largest relative deviation of the
    realized cutoff from D_um when D_um is given."""
    from opencif.checks import auditor as AU
    fluid, facts = AU.read_photomask(path)
    h = float(depth_um)
    holes = [np.array(r.coords) for r in fluid.interiors]
    A = sorted((p for p in holes if p[:, 1].min() > 0), key=lambda p: p[:, 0].min())
    B = sorted((p for p in holes if p[:, 1].max() < 0), key=lambda p: p[:, 0].min())
    if not A or len(A) != len(B):
        raise ValueError('not a straight ladder: %d posts on side A, %d on side B' % (len(A), len(B)))
    y_lo, y_hi = float(np.median([p[:, 1].min() for p in A])), float(np.median([p[:, 1].max() for p in A]))
    l_g = y_hi - y_lo                                                # the wall the skimmed fluid crosses
    xmin, ymin, xmax, ymax = fluid.bounds
    rows = {}
    for side, sg in (('A', 1.0), ('B', -1.0)):
        # The gaps are the openings along the middle of the wall.
        ym = sg * 0.5 * (y_lo + y_hi)
        gaps = [(xmin + a, xmin + b) for a, b in _cut(fluid, (xmin, ym), (xmax, ym))]
        xg = np.array([0.5 * (a + b) for a, b in gaps])
        G = np.array([b - a for a, b in gaps])
        # The center width between two gaps and each side segment's width are measured at the segment's middle.
        xm = 0.5 * (xg[:-1] + xg[1:])
        # Each cut starts or ends inside the post, where there is no fluid, so the fluid's own edges bound it.
        y_post = sg * 0.5 * (y_lo + y_hi)
        wc = np.array([sum(b - a for a, b in _cut(fluid, (x, 0.0), (x, y_post))) * 2.0 for x in xm])
        ys0, ys1 = y_post, sg * (ymax if sg > 0 else -ymin)
        ws = np.array([sum(b - a for a, b in _cut(fluid, (x, ys0), (x, ys1))) for x in xm])
        # Each outlet runs from the last gap to the end of its channel.
        x_last = xg[-1]
        L_c = xmax - x_last
        ws_out = sum(b - a for a, b in _cut(fluid, (x_last + 0.5 * L_c, ys0), (x_last + 0.5 * L_c, ys1)))
        rows[side] = dict(xg=xg, G=G, wc=wc, ws=ws, ws_out=ws_out, L_out=L_c)
    a_, b_ = rows['A'], rows['B']
    N = len(a_['xg'])
    lam = np.diff(a_['xg'])                                          # the pitches, gap to gap
    wc = 0.5 * (a_['wc'] + b_['wc'])                                  # each half measured on its own side
    # The resistances (App. A Result 14) feed the network, which opencif.checks.auditor.solve_ladder solves with
    # Kirchhoff's laws. Each list of segment resistances ends with a 0.0 that pads it to N entries, and solve_ladder
    # does not use that entry.
    Rc = np.array([r_um(w, h) * L * UM for w, L in zip(wc, lam)] + [0.0])
    RsA = np.array([r_um(w, h) * L * UM for w, L in zip(a_['ws'], lam)] + [0.0])
    RsB = np.array([r_um(w, h) * L * UM for w, L in zip(b_['ws'], lam)] + [0.0])
    RgA = np.array([gap_resistance(g, l_g, h, c_ent) for g in a_['G']])
    RgB = np.array([gap_resistance(g, l_g, h, c_ent) for g in b_['G']])
    wc_gap = np.concatenate([[wc[0]], 0.5 * (wc[:-1] + wc[1:]), [wc[-1]]])   # the center width at each gap
    Rc_out = r_um(wc_gap[-1], h) * a_['L_out'] * UM
    net = AU.solve_ladder(Rc, RsA, RsB, RgA, RgB, Rc_out, r_um(a_['ws_out'], h) * a_['L_out'] * UM,
                          r_um(b_['ws_out'], h) * b_['L_out'] * UM)
    kappa = K.Pairing().kappa
    cutoff = lambda f, w: kappa * K.x_of_f(f) * w if 0 < f < 0.5 else np.nan     # D = kappa_g x(f) w_c
    D = np.array([[cutoff(fa, w), cutoff(fb, w)] for fa, fb, w in zip(net['fA'], net['fB'], wc_gap)])
    out = dict(file=facts['file'], sha256=facts['sha256'], gaps_per_side=N, posts_per_side=len(A), wall_um=l_g,
               gap_um=[float(np.min(a_['G'])), float(np.max(a_['G']))],
               centre_um=[float(np.min(wc)), float(np.max(wc))],
               side_um=[float(np.min(a_['ws'])), float(np.max(a_['ws']))], F=float(net['F']),
               f=[float(np.min(net['fA'])), float(np.max(net['fA']))],
               D_um=[float(np.nanmin(D)), float(np.nanmax(D))], network=net, D_rows=D)
    if D_um is not None:
        out['D_rel_dev_max'] = float(np.nanmax(np.abs(D / D_um - 1.0)))
    return out
