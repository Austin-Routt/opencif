"""opencif.checks.ruler: the physics auditor's own ruler for the straight rows.

The auditor needs a second ruler because the generator measures its own drawing (opencif.checks.measure), reading the
widths of every part at the places its plan puts them and taking the lengths of the segments from the plan. An
auditor that used that measurement would share any error in it, so the auditor measures the straight rows again, from
the photomask file alone and with its own code. For every post of every ordinary straight row it finds the post
among the file's islands. From the posts' outlines and from cross-sections of the fluid it then reads (i) each gap's
position and width (the opening between two neighboring posts) and the wall it crosses (the posts' extent across the
flow), (ii) each row's pitch (from gap center to gap center), (iii) the center channel's width at every post tip (the
two posts' inner tips) and, at every gap, at the corner of the trunk's edge that the tips fix (wall_corners), and (iv) each side channel's width along its segment,
from the posts' outer tips to the outer wall.

These widths become resistances by Appendix A's duct law (Result 14, opencif.physics.kernels via schedule.rr) and
spec v5's gap resistance (sec. 14.4, schedule.r_gap). A center segment's resistance comes from Simpson's rule over
its three widths (gap, post tip, and gap). A side segment's resistance comes from integrating r(w) piecewise between
the polygon's own vertices, because on each piece both of the channel's edges are straight. A gap is a short duct
through the wall at the width of its mean resistance over five cross-sections.

The plan says which rows are ordinary straight rows (excluding the head, the head-to-ladder transition row, the rows
next to a U-turn, and the last row) and where on the mask to start looking for each post. No planned width,
position, or length enters a measurement, and nothing here imports opencif.checks.measure. The head units, the
U-turns, and the outlet are not measured here (docs/PHYSICS_AUDIT.md)."""
import numpy as np
from shapely.geometry import LineString, Point, Polygon, box
from shapely.strtree import STRtree

from opencif.physics import schedule as S

UM = S.UM


def _leg_of(legs, k):
    """Return the leg of the layout whose rows include row k."""
    for lg in legs:
        if lg['rows'][0] <= k <= lg['rows'][1]:
            return lg
    raise ValueError('row %d is on no leg' % k)


def _to_mask(leg, gx0, x, y):
    """Return the mask coordinates of a point given in a leg's own frame. In that frame x runs along the flow from
    the leg's start and y runs across it, positive on side A."""
    return (gx0 + x, leg['Y'] + y) if leg['dirn'] > 0 else (gx0 - x, leg['Y'] - y)


def _to_leg(leg, gx0, bounds):
    """Return a mask bounding box in the leg's frame, as (x_min, x_max, y_min, y_max)."""
    b = bounds
    if leg['dirn'] > 0:
        return b[0] - gx0, b[2] - gx0, b[1] - leg['Y'], b[3] - leg['Y']
    return gx0 - b[2], gx0 - b[0], leg['Y'] - b[3], leg['Y'] - b[1]


def straight_posts(sc):
    """Return the posts that the ruler reads. Post k lies between gaps k and k+1 (0-based) of a straight row, and
    neither it nor its two neighbors touches a U-turn."""
    kind = np.asarray(sc.kind)
    turns = set(int(t) for t in sc.info.get('turns', []))
    return [k for k in range(int(sc.N) - 1)
            if kind[k] == 'L' and kind[k + 1] == 'L' and not (turns & {k - 1, k, k + 1})]


def find_posts(pm, sc, info):
    """Find every straight-row post on the photomask. Each post is the island of the polygon read from the file that
    contains the point half a wall thickness outside the planned center channel, midway between the post's two gaps.
    That point is only a search point, so nothing of the plan is measured.

    The function returns {(k, side): (x_min, x_max, y_min, y_max) in the leg's frame} with k 0-based, the posts not
    found as a list [(k + 1, side)] numbered from 1, and the legs' frames by post."""
    spec = sc.spec
    lg = float(spec.l_g)
    holes = [Polygon(r) for p in (list(pm.geoms) if hasattr(pm, 'geoms') else [pm]) for r in p.interiors]
    tree = STRtree(holes)
    x, wc = np.asarray(sc.x, float), np.asarray(sc.wc, float)
    found, missing, frame = {}, [], {}
    for k in straight_posts(sc):
        leg = _leg_of(info['legs'], k)
        gx0 = info['gx'][leg['j']]
        off = leg['p0'] if leg['p0'] is not None else 0.0
        xm = 0.5 * (x[k] + x[k + 1])
        hw = 0.5 * float(np.interp(xm, x, wc))
        frame[k] = (leg, gx0)
        for sg, side in ((1.0, 'A'), (-1.0, 'B')):
            c = Point(_to_mask(leg, gx0, xm - off, sg * (hw + lg / 2.0)))
            hit = [holes[i] for i in tree.query(c) if holes[i].contains(c)]
            if not hit:
                missing.append((k + 1, side))
                continue
            found[(k, side)] = _to_leg(leg, gx0, hit[0].bounds)
    return found, missing, frame


def _first_run(poly, p0, p1):
    """Return the length of fluid from p0 along the segment p0 -> p1 before the first wall, or 0 when p0 is not in
    the fluid."""
    seg = LineString([p0, p1])
    inter = poly.intersection(seg)
    parts = [inter] if inter.geom_type == 'LineString' else [g for g in getattr(inter, 'geoms', [])
                                                            if g.geom_type == 'LineString']
    q0 = Point(p0)
    for g in parts:
        # GEOS returns the piece's ends to within a few 1e-6 um of the segment's own points, so the run is measured
        # from p0 to the piece's far end.
        if not g.is_empty and q0.distance(g) < 1e-3:
            return max(q0.distance(Point(c)) for c in g.coords)
    return 0.0


# the 3-point Gauss-Legendre rule on [0, 1], as (node, weight) pairs
GL3 = ((0.5 - 0.5 * np.sqrt(0.6), 5.0 / 18.0), (0.5, 8.0 / 18.0), (0.5 + 0.5 * np.sqrt(0.6), 5.0 / 18.0))
GAP_HEIGHTS = 5            # cross-sections through each gap, across the straight part of the wall


def _vertices(poly):
    """Return every vertex of a polygon or multipolygon, from the exterior and the holes, as an (n, 2) array."""
    pts = []
    for p in (list(poly.geoms) if hasattr(poly, 'geoms') else [poly]):
        if p.is_empty or p.geom_type != 'Polygon':
            continue
        pts.append(np.asarray(p.exterior.coords))
        pts.extend(np.asarray(r.coords) for r in p.interiors)
    return np.vstack(pts) if pts else np.zeros((0, 2))


RUN_MIN = 16               # the fewest consecutive posts from which the corners are found
MID_TOL = 1e-3             # how far a post's middle may lie from its interval's middle, as a fraction of the interval


def wall_corners(ks, frame, found, tip_in, mid, side):
    """Return {(k, side): v_k}, the trunk edge's distance from the axis at gap k, on one side, found from the posts'
    tips alone.

    Within one row interval the trunk's edge is straight, and its corners lie at the gaps, where the wall is open. A
    post sits at the middle of its interval, and its tip touches the edge there, so tip j = (v_j + v_{j+1}) / 2 for
    the corners v_j and v_{j+1} of its interval. Along a run of consecutive posts this fixes every corner from the
    first one: v_{j+1} = 2 tip_j - v_j. The first corner is not seen, and an error in it alternates in sign from gap
    to gap, so it is chosen to make the corners' second differences smallest (least squares): with p the corners from
    a first corner of 0 and s_i = (-1)^i, v = p + c s and c = sum(s_i d2p_i) / (4 m) over the m inner corners. A
    straight taper and a change of slope at one gap (where the taper meets its floor) are then both read exactly, to
    within the change of slope divided by 4m, whereas the mean of the two neighboring tips missed a change of slope
    by a quarter of it. Runs shorter than RUN_MIN posts, or whose posts are not at the middles of their intervals
    (MID_TOL), keep the mean of the two tips."""
    out = {}
    run = []

    def close(run):
        if len(run) < RUN_MIN:
            return
        tips = np.array([tip_in[(j, side)] for j in run], float)
        p = np.zeros(len(run) + 1)
        for i, t in enumerate(tips):
            p[i + 1] = 2.0 * t - p[i]
        sgn = np.array([1.0 if i % 2 == 0 else -1.0 for i in range(len(p))])
        d2 = p[:-2] - 2.0 * p[1:-1] + p[2:]
        c = float(np.sum(sgn[1:-1] * d2)) / (4.0 * len(d2))
        v = p + c * sgn
        for i, j in enumerate(run):
            out[(j, side)] = float(v[i])
        out[(run[-1] + 1, side)] = float(v[-1])

    prev = None
    for k in ks:
        ok = (k, side) in found
        if ok and prev is not None and (k - 1, side) in found and (k + 1, side) in found:
            # the post's middle against the middle of its interval, whose ends are the midpoints to its neighbors
            x0 = 0.5 * (mid[(k - 1, side)] + mid[(k, side)])
            x1 = 0.5 * (mid[(k, side)] + mid[(k + 1, side)])
            ok = abs(mid[(k, side)] - 0.5 * (x0 + x1)) <= MID_TOL * (x1 - x0)
        if ok and prev is not None and k == prev + 1 and frame[k][0] is frame[prev][0]:
            run.append(k)
        else:
            close(run)
            run = [k] if ok else []
        prev = k if ok else None
    close(run)
    return out


def measure_straight_rows(pm, sc, info, posts=None):
    """Measure every ordinary straight row on the polygon read from the photomask, with the auditor's own code (see
    the module's docstring).

    The function returns dict(Rc={k: R}, RsA, RsB, RgA={k: R}, RgB, wc_gap={k: w}, gaps={(k, side): (x, G, l_g)},
    x_dev_max_um, missing, posts), in which k is 0-based, Rc and Rs belong to segment k (gap k to gap k+1), and Rg
    belongs to gap k."""
    spec = sc.spec
    h, c_ent = float(spec.h), float(spec.c_ent)
    rr = lambda w: S.rr(w, h)
    if posts is None:
        posts = find_posts(pm, sc, info)
    found, missing, frame = posts
    ks = sorted({k for k, _ in found})
    tip_in, tip_out, mid, half = {}, {}, {}, {}
    for (k, side), (xa, xb, ya, yb) in found.items():
        # Side A lies at positive y, so its inner tip is the post's lowest point and its outer tip the highest. Side B
        # is mirrored.
        tip_in[(k, side)] = ya if side == 'A' else -yb
        tip_out[(k, side)] = yb if side == 'A' else -ya
        mid[(k, side)] = 0.5 * (xa + xb)
        half[(k, side)] = 0.5 * (xb - xa)
    ws_design = np.asarray(sc.ws, float)
    # The fluid near each block of 32 posts of one side of one leg is clipped once, because sections on the clip are
    # much cheaper than on the whole mask. Its vertices are kept in the leg's frame, with the side coordinate
    # positive away from the center.
    clips = {}

    def clip(k, side):
        leg, gx0 = frame[k]
        key = (id(leg), side, k // 32)
        if key not in clips:
            sg = 1.0 if side == 'A' else -1.0
            blk = [kk for kk in ks if kk // 32 == k // 32 and frame[kk][0] is leg and (kk, side) in found]
            reach = 2.0 * float(np.nanmax(ws_design[max(0, blk[0] - 2):blk[-1] + 3])) + 60.0      # where to look
            x_lo = min(mid[(kk, side)] for kk in blk) - 200.0
            x_hi = max(mid[(kk, side)] for kk in blk) + 200.0
            y_lo = min(tip_in[(kk, side)] for kk in blk) - 10.0
            y_hi = max(tip_out[(kk, side)] for kk in blk) + reach + 400.0
            q = [_to_mask(leg, gx0, xx, sg * yy) for xx in (x_lo, x_hi) for yy in (y_lo, y_hi)]
            cl = pm.intersection(box(min(p[0] for p in q), min(p[1] for p in q), max(p[0] for p in q),
                                     max(p[1] for p in q)))
            v = _vertices(cl)
            if leg['dirn'] > 0:
                vx, vy = v[:, 0] - gx0, v[:, 1] - leg['Y']
            else:
                vx, vy = gx0 - v[:, 0], leg['Y'] - v[:, 1]
            clips[key] = (cl, vx, sg * vy, reach)
        return clips[key]

    out = dict(Rc={}, RsA={}, RsB={}, RgA={}, RgB={}, wc_gap={}, gaps={}, missing=list(missing))
    # Gap k is the opening between post k-1 and post k. It is read on GAP_HEIGHTS cross-sections across the straight
    # part of the wall (between the posts' rounded ends), together with the wall it crosses (the posts' extent across
    # the flow).
    for k in ks:
        for side in ('A', 'B'):
            a, b = found.get((k - 1, side)), found.get((k, side))
            if a is None or b is None or frame[k - 1][0] is not frame[k][0]:
                continue
            leg, gx0 = frame[k]
            sg = 1.0 if side == 'A' else -1.0
            cl = clip(k, side)[0]
            r_end = max(half[(k - 1, side)], half[(k, side)])          # the rounded ends' extent, a post's half length
            y0 = max(tip_in[(k - 1, side)], tip_in[(k, side)]) + r_end
            y1 = min(tip_out[(k - 1, side)], tip_out[(k, side)]) - r_end
            ys = np.linspace(y0, y1, GAP_HEIGHTS) if y1 > y0 else np.array([0.5 * (y0 + y1)])
            Gs = np.array([cl.intersection(LineString([_to_mask(leg, gx0, mid[(k - 1, side)], sg * yy),
                                                       _to_mask(leg, gx0, mid[(k, side)], sg * yy)])).length
                           for yy in ys])
            l_g = 0.5 * ((tip_out[(k - 1, side)] - tip_in[(k - 1, side)]) + (tip_out[(k, side)] - tip_in[(k, side)]))
            if np.all(Gs > 0):
                G = float(S.rinv(float(np.mean([rr(g) for g in Gs])), h))      # the width of the mean resistance
                out['Rg' + side][k] = S.r_gap(G, l_g, h, c_ent)
            else:
                G = 0.0
                out['Rg' + side][k] = float('nan')
            out['gaps'][(k, side)] = (0.5 * (a[1] + b[0]), G, l_g)
    # The center width is read at each post tip. At a gap the wall is open, so the width there is found from the
    # tips (wall_corners); the side channel's inner edge, l_g outside the trunk's edge, follows the same corners.
    wcp = {k: tip_in[(k, 'A')] + tip_in[(k, 'B')] for k in ks if (k, 'A') in found and (k, 'B') in found}
    corner = {}                                             # (k, side) -> the trunk edge's half-width at gap k
    for side in ('A', 'B'):
        corner.update(wall_corners(ks, frame, found, tip_in, mid, side))
    for k in ks:
        if k in wcp and k - 1 in wcp and frame[k - 1][0] is frame[k][0]:
            if (k, 'A') in corner and (k, 'B') in corner:
                out['wc_gap'][k] = corner[(k, 'A')] + corner[(k, 'B')]
            else:
                out['wc_gap'][k] = 0.5 * (wcp[k - 1] + wcp[k])
    # Segment k runs from gap k to gap k+1, and its measurement needs both gaps and the post between them, on both
    # sides.
    x_design = np.asarray(sc.x, float)
    xdev = 0.0
    for k in ks:
        if not all((k, s) in out['gaps'] and (k + 1, s) in out['gaps'] for s in ('A', 'B')) or k not in wcp:
            continue
        leg, gx0 = frame[k]
        off = leg['p0'] if leg['p0'] is not None else 0.0
        xg = {s: (out['gaps'][(k, s)][0], out['gaps'][(k + 1, s)][0]) for s in ('A', 'B')}
        lam = 0.5 * ((xg['A'][1] - xg['A'][0]) + (xg['B'][1] - xg['B'][0]))
        xdev = max(xdev, abs(xg['A'][0] - (x_design[k] - off)), abs(xg['B'][0] - (x_design[k] - off)))
        out['Rc'][k] = (rr(out['wc_gap'][k]) + 4.0 * rr(wcp[k]) + rr(out['wc_gap'][k + 1])) / 6.0 * lam * UM
        for side, sg in (('A', 1.0), ('B', -1.0)):
            # The inner edge of the side channel is the trunk's edge l_g further out: within the segment, the line
            # from the corner at gap k to the corner at gap k+1, both carried out by the post's own extent across the
            # flow (where the corners are not known, the line through the posts' outer tips at their middles). The
            # side width is the fluid from that edge to the outer wall. Both edges are straight between the
            # polygon's vertices, so the integral of r(w) is taken piecewise between them, with the 3-point
            # Gauss-Legendre rule on each piece.
            x0, x1, xm = xg[side][0], xg[side][1], mid[(k, side)]
            if (k, side) in corner and (k + 1, side) in corner:
                lgk = tip_out[(k, side)] - tip_in[(k, side)]
                xs_tip = [x0 - 1.0, x0, x1, x1 + 1.0]
                y0_, y1_ = corner[(k, side)] + lgk, corner[(k + 1, side)] + lgk
                d_ = (y1_ - y0_) / (x1 - x0)
                ys_tip = [y0_ - d_, y0_, y1_, y1_ + d_]
            else:
                xs_tip = [mid[(k - 1, side)], mid[(k, side)], mid[(k + 1, side)]]
                ys_tip = [tip_out[(k - 1, side)], tip_out[(k, side)], tip_out[(k + 1, side)]]
            cl, vx, vy, reach = clip(k, side)
            y_in_v = np.interp(vx, xs_tip, ys_tip)
            sel = (vx > x0 + 1e-6) & (vx < x1 - 1e-6) & (vy > y_in_v - 1e-3) & (vy < y_in_v + reach)
            bps = np.unique(np.concatenate([[x0, xm, x1], vx[sel]]))
            bps = bps[np.concatenate([[True], np.diff(bps) > 1e-6])]
            R = 0.0
            for u0, u1 in zip(bps[:-1], bps[1:]):
                for t, wgt in GL3:
                    xx = u0 + t * (u1 - u0)
                    y_in = float(np.interp(xx, xs_tip, ys_tip)) + 1e-6
                    w = _first_run(cl, _to_mask(leg, gx0, xx, sg * y_in), _to_mask(leg, gx0, xx, sg * (y_in + reach)))
                    if w <= 0.0:
                        R = float('nan')
                        break
                    R += wgt * (u1 - u0) * rr(w)
                if not np.isfinite(R):
                    break
            out['Rs' + side][k] = R * UM
    out['x_dev_max_um'] = xdev
    out['posts'] = len(found)
    return out
