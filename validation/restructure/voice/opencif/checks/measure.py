"""opencif.checks.measure: measure a whole drawn device (any layout) and solve its as-drawn network, with both sides
explicit.

Every part is measured on the polygon, and the positions of the cuts and the corridors come from the construction:

  head units ........ the centerline length is the channel area in a corridor (with the junction squares excluded)
                      divided by the measured width, plus the excluded junction lengths (the 1-D convention), and the
                      widths come from cross-sections
  transition drop ... 201 cross-sections and the integral of r(w) dx
  linear segments ... the side widths on both sides, with r(w) integrated over each row's drawn profile and the width
                      through every post center read as well, and the center width at the post tip near mid-segment
  gaps .............. the width at mid-wall (both sides), with R_gap = r(G)[l_g + c_ent D_h] (the spec 14.4 model)
  loops ............. on each side, cross-sections perpendicular to the path every ~1 um of center path and the
                      integral of r(w) ds along that side's own centerline
  outlet ............ the stretch widths and each lead's width along its drawn centerline, with no trust in the design

The network is the center rail with side A and side B (no symmetry is assumed), an inlet flow of 1, and both outlets
at 0."""
import math
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl
from shapely.geometry import LineString, box, Point
from opencif.physics import schedule as S
from opencif.geometry import shapes as Gm
from opencif.geometry import layout as L
from shapely.ops import unary_union


def xsec(fluid, p0, p1):
    """Return the length of fluid on the segment from p0 to p1, which is the width that a cross-section reads."""
    return fluid.intersection(LineString([p0, p1])).length


def leg_of(legs, k):
    """Return the leg whose rows include row k."""
    for lg in legs:
        if lg['rows'][0] <= k <= lg['rows'][1]:
            return lg
    raise ValueError(k)


def l2g(leg, gx0, x, y):
    """Map a point (x, y) in a leg's local frame to the global frame, in which the leg starts at x = gx0. A leg that
    runs in the negative direction is turned by 180 degrees."""
    if leg['dirn'] > 0:
        return (gx0 + x, leg['Y'] + y)
    return (gx0 - x, leg['Y'] - y)


def solve3(Rc, RsA, RsB, RgA, RgB, Rc_out, RsA_out, RsB_out):
    """Solve the explicit two-sided ladder for a unit inlet flow with both outlets at 0. The unknown pressures are the
    center P_c(i) and the sides P_A(i) and P_B(i), for i = 0 to N-1."""
    N = len(RgA)
    n = 3 * N
    I, J, V = [], [], []
    b = np.zeros(n)

    def cond(a, c, g):
        I.extend([a, a, c, c]); J.extend([a, c, c, a]); V.extend([g, -g, g, -g])
    # the node indices
    cI = lambda i: i
    aI = lambda i: N + i
    bI = lambda i: 2 * N + i
    diag_out = np.zeros(n)
    for i in range(N):
        cond(cI(i), aI(i), 1.0 / RgA[i])
        cond(cI(i), bI(i), 1.0 / RgB[i])
        if i < N - 1:
            cond(cI(i), cI(i + 1), 1.0 / Rc[i])
            cond(aI(i), aI(i + 1), 1.0 / RsA[i])
            cond(bI(i), bI(i + 1), 1.0 / RsB[i])
    diag_out[cI(N - 1)] += 1.0 / Rc_out
    diag_out[aI(N - 1)] += 1.0 / RsA_out
    diag_out[bI(N - 1)] += 1.0 / RsB_out
    A = sp.csr_matrix((V, (I, J)), shape=(n, n)) + sp.diags(diag_out)
    b[cI(0)] = 1.0
    x = spl.spsolve(A.tocsc(), b)
    Pc, PA, PB = x[:N], x[N:2 * N], x[2 * N:]
    QA = (Pc - PA) / RgA
    QB = (Pc - PB) / RgB
    Qc_in = np.concatenate([[1.0], (Pc[:-1] - Pc[1:]) / np.asarray(Rc)[:N - 1]])
    fA, fB = QA / Qc_in, QB / Qc_in
    F = 1.0 - Pc[-1] / Rc_out
    return dict(fA=fA, fB=fB, F=F, Pc=Pc, PA=PA, PB=PB)


_GLX, _GLW = np.polynomial.legendre.leggauss(6)


TRANS_NSUB = 10


def trans_nodes(xa, xb, nsub=TRANS_NSUB):
    """Return the nodes and weights of a composite 6-point Gauss-Legendre rule over [xa, xb] with nsub
    sub-intervals."""
    e = np.linspace(xa, xb, nsub + 1)
    xs = np.concatenate([e[i] + 0.5 * (_GLX + 1.0) * (e[i + 1] - e[i]) for i in range(nsub)])
    ws = np.concatenate([0.5 * _GLW * (e[i + 1] - e[i]) for i in range(nsub)])
    return xs, ws


def measure_device(fluid, sc, info, delta=4.0):
    """Measure the drawn fluid polygon of the device that sc plans and info describes, and return the measured
    resistances (Rc, RsA, RsB, RgA, RgB, and the outlet's), the widths behind them, and the planned features that the
    drawing lacks (not_found). A head unit's corridor and cross-sections reach delta (um) beyond its half-width."""
    spec = sc.spec
    h, lg, Lam = spec.h, spec.l_g, spec.Lam_L
    N, NM = sc.N, sc.NM
    legs, gx = info['legs'], info['gx']
    _clips = {}
    ws_max = float(np.nanmax(sc.ws)) if np.any(np.isfinite(sc.ws)) else spec.GS
    y_clip = float(np.max(sc.wc)) / 2.0 + lg + max(ws_max, spec.GS) + 80.0

    def clip_for(leg, j, k, off, chunk=32):
        """Return the fluid clipped to a block of rows of one leg, on which cross-sections are ~50x cheaper than on
        the whole mask."""
        key = (j, k // chunk)
        if key not in _clips:
            k0 = (k // chunk) * chunk
            k1 = min(k0 + chunk, N - 1)
            xa, xb = sc.x[k0] - off - 2.0 * Lam, sc.x[k1] - off + 2.0 * Lam
            pa, pb = l2g(leg, gx[j], xa, -y_clip), l2g(leg, gx[j], xb, y_clip)
            _clips[key] = fluid.intersection(box(min(pa[0], pb[0]), min(pa[1], pb[1]), max(pa[0], pb[0]), max(pa[1], pb[1])))
        return _clips[key]
    out = dict(units=[], lin=[], gaps=[], wc=[], loops=[], not_found=[])
    # A measured width that enters r(w) must be positive. Where a section finds no fluid at the place where the plan
    # puts a feature (a head unit, gap, or side channel that the drawing does not contain), the feature is recorded
    # once per (row, kind, side) in out['not_found'] and its design width stands in, so that the measurement completes.
    # pipeline.assess then fails the build by name ('every planned feature found on the drawing'). A positive width is
    # returned unchanged, so this substitution leaves the measurement of a good drawing as it is (E5, Reviewer A).
    _nf_seen = {}

    def found(w, w_design, k, kind, side):
        w = float(w)
        if w > 0.0:
            return w
        key = (int(k), kind, side)
        if key in _nf_seen:
            _nf_seen[key]['sections_missed'] += 1
        else:
            _nf_seen[key] = dict(row=int(k) + 1, kind=kind, side=side, measured_width_um=w,
                                 design_width_um=float(w_design), sections_missed=1)
            out['not_found'].append(_nf_seen[key])
        return float(w_design)
    RsA = np.array(sc.Rs, float).copy(); RsB = RsA.copy()
    RgA = np.array(sc.Rg, float).copy(); RgB = RgA.copy()
    Rc = np.array(sc.Rc, float).copy()
    # ---- head units (leg 1 is the identity frame; side B is measured too, at the mirror image of side A's coordinates)
    units_ = list(info['units'])
    for iu, u in enumerate(units_):
        k, w = u['k'], u['w']
        p = Gm.unit_path(spec, u, from_J=True)
        yJ = spec.yJ(w)
        m = 7.0
        cut = [box(u['X'] - w / 2 - m, yJ - w / 2 - m, u['X'] + w / 2 + m, yJ + w / 2 + m)]
        if u['kind'] in ('T', 'A'):
            Xn = u['X'] + u['ell']
            # The junction with the next unit is excluded too. Where the width changes (from narrow to full units),
            # that junction's stem is wider, and its fillets must fall inside the excluded square. In review, row 20 of
            # gifford_targets_meander_head_19um read +6.9e-4 long without this.
            wn = float(units_[iu + 1]['w']) if iu + 1 < len(units_) else w
            wq = max(w, wn)
            yJn = spec.yJ(wn)
            cut.append(box(Xn - wq / 2 - m - spec.r_fillet, min(yJ, yJn) - wq / 2 - m,
                           Xn + wq / 2 + m, max(yJ, yJn) + wq / 2 + m))
        # The centerline length inside the excluded squares is added back exactly.
        excl = LineString(p.p).intersection(unary_union(cut)).length
        # A square corner's roundings are counted by its corner term (c_corner, in bend_excess), so their net area is
        # taken off the drawn area before the area becomes a length; corners inside the excluded squares are not in the
        # area. An arc path has no corners, so its measurement is unaffected (DINH-1).
        a_round = sum(1 for ci in getattr(p, 'corners', []) if not any(c.contains(Point(p.p[ci])) for c in cut)) * \
            Gm.corner_round_area(spec)
        res = {}
        for side, sg in (('A', 1.0), ('B', -1.0)):
            corr = LineString([(x, sg * y) for x, y in p.p]).buffer(w / 2.0 + delta, cap_style='flat',
                                                                    join_style='round', quad_segs=8)
            for c in cut:
                cc = c if sg > 0 else box(c.bounds[0], -c.bounds[3], c.bounds[2], -c.bounds[1])
                corr = corr.difference(cc)
            area = fluid.intersection(corr).area
            P = np.asarray(p.p)
            ws = []
            for a, bb in zip(P[:-1], P[1:]):
                Ls = math.hypot(*(bb - a))
                if Ls > 3 * w:
                    for t in (0.25, 0.5, 0.75):
                        c = a + t * (bb - a)
                        if any(cc.contains(Point(c)) for cc in cut):
                            continue
                        d = (bb - a) / Ls
                        nn = np.array([-d[1], d[0]])
                        q0, q1 = c - nn * (w / 2 + delta), c + nn * (w / 2 + delta)
                        ws.append(xsec(fluid, (q0[0], sg * q0[1]), (q1[0], sg * q1[1])))
            wm = float(np.median(ws)) if ws else w
            wm = found(wm, w, k, 'head unit ' + u['kind'], side)
            Lm = (area - a_round) / wm + excl
            # the drawn length plus the bend correction of the unit's arcs (at the design radii and the measured width)
            bx = Gm.bend_excess(spec, p, wm) + Gm.junction_excess(spec, u, wm)
            R = S.rr(wm, h) * (Lm + bx) * S.UM
            if u['kind'] == 'X':
                s_, rq = spec.s, spec.s / 2.0
                x0, x1 = u['X'] + s_ + rq, u['X'] + u['ell']
                y_in = spec.wc0 / 2.0 + lg
                xs = np.linspace(x0, x1, 201)
                wd = np.array([found(xsec(fluid, (x, sg * y_in), (x, sg * (y_in + spec.GS + 10))), spec.GS, k,
                                     'head unit X stem', side) for x in xs])
                R += float(np.trapezoid([S.rr(ww, h) for ww in wd], xs) * S.UM)
            res[side] = dict(L=Lm, L_eff=Lm + bx, w=wm, wmin=min(ws) if ws else None, wmax=max(ws) if ws else None, R=R)
        RsA[k], RsB[k] = res['A']['R'], res['B']['R']
        out['units'].append(dict(k=k, kind=u['kind'], L_req=u.get('L'),    # for the transition unit, its loop length
                                 L_drawn=u.get('L_drawn'), A=res['A'], B=res['B'],
                                 R_target=sc.Rs[k]))
    # ---- linear segments, gaps, and center widths
    turns = set(sc.info.get('turns', []))
    paper_, wR_, kL_ = L.transition_params(sc)
    pulled_ = L.pulled_rows(sc)
    wdraw_ = lambda k: L.drawn_side_width(sc, k, paper_, wR_, pulled_)
    if not paper_:
        wR_ = 0.0
    for k in range(N):
        leg = leg_of(legs, k)
        j = leg['j']
        off = leg['p0'] if leg['p0'] is not None else 0.0
        # gap k
        xl = sc.x[k] - off
        hw = L.wc_at(sc, sc.x[k]) / 2.0
        wdes = sc.wL[k]
        for side, sg in (('A', 1.0), ('B', -1.0)):
            ym = sg * (hw + lg / 2.0)
            p0 = l2g(leg, gx[j], xl - wdes / 2 - 6, ym)
            p1 = l2g(leg, gx[j], xl + wdes / 2 + 6, ym)
            wg = found(xsec(fluid, p0, p1), wdes, k, 'gap', side)
            Rgk = S.r_gap(wg, lg, h, spec.c_ent)
            if side == 'A':
                RgA[k] = Rgk
            else:
                RgB[k] = Rgk
            out['gaps'].append((k, side, wg, wdes))
        if k >= N - 1 or k < NM or k in turns:
            continue
        # Side segment k is straight and lies within one leg, and r(w) is integrated over its drawn profile. An ordinary
        # row takes 6 Gauss-Legendre sections on each half of the row. The head-to-linear transition rows (the
        # rectilinear rows and the first linear row), where the papers' transition makes the width vary within a row,
        # take TRANS_NSUB x 6 sections. The width at mid-segment is read on every row (the drawn wall is piecewise
        # linear through it), and it is the width recorded for a transition row.
        xm = 0.5 * (sc.x[k] + sc.x[k + 1])
        hwm = L.wc_at(sc, xm) / 2.0
        integ = (NM <= k <= kL_) and NM > 0
        for side, sg in (('A', 1.0), ('B', -1.0)):
            y0, y1 = sg * (hwm + lg - 2.0), sg * (hwm + lg + max(sc.ws[k], wR_) + 30.0)
            wsm = xsec(fluid, l2g(leg, gx[j], xm - off, y0), l2g(leg, gx[j], xm - off, y1))
            if not integ:
                # On every ordinary straight segment, r(w) is integrated over the drawn profile, because the legacy
                # check read the midpoint width only and could not see the wall rule's resistance error (opencif fix
                # F1).
                cl = clip_for(leg, j, k, off)
                xa_, xb_ = sc.x[k], sc.x[k + 1]
                Rk = 0.0
                for u0, u1 in ((xa_, xm), (xm, xb_)):
                    for ti, wi in zip(0.5 * (_GLX + 1.0), _GLW):
                        xx = u0 + ti * (u1 - u0)
                        hw_ = L.wc_at(sc, xx) / 2.0
                        w_ = found(xsec(cl, l2g(leg, gx[j], xx - off, sg * (hw_ + lg + 1e-6)),
                                        l2g(leg, gx[j], xx - off, sg * (hw_ + lg + max(sc.ws[k], wR_) + 30.0))),
                                   sc.ws[k], k, 'side channel', side)
                        Rk += 0.5 * (u1 - u0) * wi * S.rr(w_, h)
                Rk *= S.UM
                w_eq = S.rinv(Rk / (sc.ell[k] * S.UM), h)
            elif integ:
                # A composite Gauss-Legendre rule over the whole row (TRANS_NSUB x 6 sections) uses the same nodes and
                # weights as the intended network's integration of the drawing code's profile
                # (build.intended_network), so a perfect drawing gives the intended resistance (review R3-2).
                xs_, wq_ = trans_nodes(sc.x[k], sc.x[k + 1])
                ws_ = []
                for xx in xs_:
                    hw_ = L.wc_at(sc, xx) / 2.0
                    # The section starts on the side channel's inner edge, because across a gap or an opening the band
                    # is fluid too, and a section starting 2 um inside it would add the opening to the side width.
                    ws_.append(found(xsec(fluid, l2g(leg, gx[j], xx - off, sg * (hw_ + lg + 1e-6)),
                                          l2g(leg, gx[j], xx - off, sg * (hw_ + lg + max(sc.ws[k], wR_) + 30.0))),
                                     wdraw_(k), k, 'side channel (transition row)', side))
                Rk = float(np.sum(wq_ * np.array([S.rr(w_, h) for w_ in ws_]))) * S.UM
            if integ:
                w_eq = wsm
            if side == 'A':
                RsA[k] = Rk
            else:
                RsB[k] = Rk
            # Each entry is (row, side, measured width, drawn intent). An ordinary row reports the
            # resistance-equivalent width of the whole drawn segment, and a transition row of the papers reports the
            # mid-segment width (compared with ws_rect).
            out['lin'].append((k, side, w_eq, wdraw_(k)))
        # The center width is read at the post tip, as the narrowest section (relative to the design taper) within
        # +-0.3 ell of the row's middle. An ordinary row has its post centered there, whereas the papers' transition
        # puts the first pill off-center (by 2.5 um in Dinh), where a mid-row section read the rounded post flank
        # 0.21 um short (Rc -1.4%).
        w0, w1 = sc.wc[k], sc.wc[k + 1]
        devs = []
        clc = clip_for(leg, j, k, off)                      # the clipped fluid, on which sections are ~50x cheaper
        for xx in xm + np.linspace(-0.3, 0.3, 13) * sc.ell[k]:
            hwx = L.wc_at(sc, xx) / 2.0
            devs.append(xsec(clc, l2g(leg, gx[j], xx - off, -hwx - 3), l2g(leg, gx[j], xx - off, hwx + 3)) - 2.0 * hwx)
        dw = float(min(devs)) + (hwm - 0.25 * (w0 + w1)) * 2.0       # the deviation, referred to the row's mean width
        if min(w0, w1) + dw <= 0.0:                          # no center channel at the post tip (E5)
            found(0.0, 0.5 * (w0 + w1), k, 'centre', 'C')
            dw = 0.0
        wcm = 0.5 * (w0 + w1) + dw
        Rc[k] = (S.rr(w0 + dw, h) + 4 * S.rr(wcm, h) + S.rr(w1 + dw, h)) / 6.0 * sc.ell[k] * S.UM
        out['wc'].append((k, wcm, 0.5 * (w0 + w1)))
    # The center segments of the head rows (constant wc0) are measured at the wall center between the stubs.
    for k in range(NM):
        xm = 0.5 * (sc.x[k] + sc.x[k + 1])
        wcm = found(xsec(fluid, (xm, -spec.wc0 / 2 - 3), (xm, spec.wc0 / 2 + 3)), spec.wc0, k, 'centre (head)', 'C')
        Rc[k] = S.rr(wcm, h) * sc.ell[k] * S.UM
        out['wc'].append((k, wcm, spec.wc0))
    # ---- loops (r(w) ds integrated along each side's centerline)
    for lp in info['loops']:
        k = lp['k']
        R_t = lp.get('R', sc.info.get('R_turn', 0.0))
        th0_ = lp.get('th0', -math.pi / 2.0)
        leg = leg_of(legs, k)
        j = leg['j']
        nxt = legs[j + 1]
        off = leg['p0'] if leg['p0'] is not None else 0.0
        a_, b_ = Lam / 2.0, Lam / 2.0 + math.pi * R_t
        rec = {}
        for side, sg in (('A', 1.0), ('B', -1.0)):
            f = lp['fA'] if side == 'A' else lp['fB']
            sA = -1.0 if lp['ccw'] else 1.0
            s_rad = sA * sg                                  # +1 when this side lies outward of the axis radius
            inner = s_rad < 0

            def wmeas(ps, f=f, sg=sg, s_rad=s_rad):
                out_w = []
                for pp in np.atleast_1d(ps):
                    hwp = L.wc_at(sc, sc.x[k] + pp) / 2.0
                    wdes = float(f(pp))
                    if pp <= a_ + 1e-9:
                        xl = sc.x[k] + pp - off
                        q0 = l2g(leg, gx[j], xl, sg * (hwp + lg))
                        q1 = l2g(leg, gx[j], xl, sg * (hwp + lg + wdes + 15))
                    elif pp >= b_ - 1e-9:
                        xl = (sc.x[k] + pp) - nxt['p0']
                        q0 = l2g(nxt, gx[j + 1], xl, sg * (hwp + lg))
                        q1 = l2g(nxt, gx[j + 1], xl, sg * (hwp + lg + wdes + 15))
                    else:
                        u = (pp - a_) / (math.pi * R_t)
                        th = th0_ + (math.pi * u if lp['ccw'] else -math.pi * u)
                        r_in = R_t + s_rad * (hwp + lg)
                        r_out = R_t + s_rad * (hwp + lg + wdes + 15)
                        q0 = (lp['cx'] + r_in * math.cos(th), lp['cy'] + r_in * math.sin(th))
                        q1 = (lp['cx'] + r_out * math.cos(th), lp['cy'] + r_out * math.sin(th))
                    out_w.append(found(xsec(fluid, q0, q1), wdes, k, 'loop side', 'A' if sg > 0 else 'B'))
                return np.array(out_w)
            Rk = L.loop_R(sc, k, inner, wmeas, n_str=60, n_arc=600)
            wsample = wmeas(np.linspace(0, Lam + math.pi * R_t, 41))
            rec[side] = dict(R=Rk, w_min=float(wsample.min()), w_mid=float(wsample[20]))
            if side == 'A':
                RsA[k] = Rk
            else:
                RsB[k] = Rk
        # the center channel through the loop, which is measured on the drawing too rather than taken at its design
        # value (review CR-8)
        cx_, cy_ = lp['cx'], lp['cy']

        def wmeas_c(ps):
            ws_ = []
            for pp in np.atleast_1d(ps):
                hwp = L.wc_at(sc, sc.x[k] + pp) / 2.0
                if pp <= a_ + 1e-9:
                    xl = sc.x[k] + pp - off
                    q0, q1 = l2g(leg, gx[j], xl, -(hwp + 3.0)), l2g(leg, gx[j], xl, hwp + 3.0)
                elif pp >= b_ - 1e-9:
                    xl = (sc.x[k] + pp) - nxt['p0']
                    q0, q1 = l2g(nxt, gx[j + 1], xl, -(hwp + 3.0)), l2g(nxt, gx[j + 1], xl, hwp + 3.0)
                else:
                    u = (pp - a_) / (math.pi * R_t)
                    th = th0_ + (math.pi * u if lp['ccw'] else -math.pi * u)
                    q0 = (cx_ + (R_t - hwp - 3.0) * math.cos(th), cy_ + (R_t - hwp - 3.0) * math.sin(th))
                    q1 = (cx_ + (R_t + hwp + 3.0) * math.cos(th), cy_ + (R_t + hwp + 3.0) * math.sin(th))
                ws_.append(found(xsec(fluid, q0, q1), 2.0 * hwp, k, 'loop centre', 'C'))
            return np.array(ws_)
        from scipy.integrate import simpson as _simpson
        # The arc is integrated section by section. The two straight half rows at its ends follow the straight rows'
        # convention, in which the center width is read at the post center (p = a_, b_), because elsewhere a section
        # would add the gap opening or the rounded post tips, which the 1-D model leaves out everywhere.
        PLc = Lam + math.pi * R_t
        ps_ = np.linspace(a_, b_, 301)
        wcs_ = wmeas_c(ps_)
        isub = np.arange(0, len(ps_), 10)
        bfc = np.interp(ps_, ps_[isub], [S.bend_factor(spec, R_t - wcs_[i] / 2.0, R_t + wcs_[i] / 2.0) for i in isub])
        Rck = _simpson(np.array([S.rr(w_, h) for w_ in wcs_]) * bfc, x=ps_)
        if getattr(spec, 'bend_model', 'stokes') != 'none':          # the arc's two junctions (S.KAPPA_J)
            wmc = float(np.median(wcs_))
            Rck -= min(S.kappa_j(spec, wmc) * wmc, math.pi * R_t) * (
                S.bend_factor(spec, R_t - wmc / 2.0, R_t + wmc / 2.0) - 1.0) * S.rr(wmc, h)
        w_ends = wmeas_c(np.array([a_ - 1e-6, b_ + 1e-6]))
        Rck += S.rr(w_ends[0], h) * a_ + S.rr(w_ends[1], h) * (PLc - b_)
        Rc[k] = Rck * S.UM
        out['loops'].append(dict(k=k, target=sc.Rs[k], A=rec['A'], B=rec['B'], Rc=Rc[k], Rc_design=sc.Rc[k]))
    # ---- outlet (the stretch widths at its midpoint, then the three leads measured along their drawn centerlines)
    o = info['outlet']
    leg = legs[-1]
    j = leg['j']
    xm = o['x_e_local'] - spec.L_out / 2.0
    hwe = sc.wc[-1] / 2.0
    wce_m = found(xsec(fluid, l2g(leg, gx[j], xm, -hwe - 3), l2g(leg, gx[j], xm, hwe + 3)), sc.wc[-1], N - 1,
                  'outlet stretch centre', 'C')
    Rc_str = S.rr(wce_m, h) * spec.L_out * S.UM
    rec = dict(stretch_wc=wce_m)
    R_side = {}
    for side, sg in (('A', 1.0), ('B', -1.0)):
        y0, y1 = sg * (hwe + lg - 2.0), sg * (hwe + lg + sc.ws[N - 2] + 30.0)
        wsm = found(xsec(fluid, l2g(leg, gx[j], xm, y0), l2g(leg, gx[j], xm, y1)), sc.ws[N - 2], N - 1,
                    'outlet stretch side', side)
        Rl, wl_s = measure_path(fluid, o['paths'][side], o['w_lead'], h)
        if not np.all(wl_s > 0.0):                          # a lead section that found no fluid (E5)
            wl_s = np.array([found(x, o['w_lead'], N - 1, 'side lead', side) for x in wl_s])
            Rl = measure_path_R(wl_s, o['paths'][side], o['w_lead'], h)
        # The lead's arcs (at the design radii and the measured width) carry the curved-duct correction, as in the
        # design.
        Rb = S.rr(float(np.median(wl_s)), h) * L.arcs_excess(spec, o.get('bend_arcs', {}).get(side, []),
                                                             float(np.median(wl_s))) * S.UM
        if L.flare_on(spec):
            # The stretch widens here (the flare). Its side width is measured on the drawing at 601 stations, from
            # the end of the long wall's rounded tip (G/2 + l_g/2 + 2 um after the last gap node) to 15 um before the
            # lead path, on a cross-section that reaches past the lead width, and r(w) is integrated (the exact wall
            # rule, as in the design). The tip itself is counted at the first station's width, as RC4 counted the
            # whole stretch. The last 15 um (the junction with the lead ribbon, where the fillet smoothing shaves
            # ~0.1 um) is counted at the design profile, as measure_path counts a lead's end zones. The cross-sections
            # are taken on the fluid clipped to a box around the stretch, which gives the same lengths, whereas the
            # whole mask is slow.
            EZ = 15.0
            x0s, x1s = o['x_e_local'] - spec.L_out + spec.G / 2.0 + spec.l_g / 2.0 + 2.0, o['x_e_local'] - EZ
            xst = np.linspace(x0s, x1s, 601)
            y1f = sg * (hwe + lg + max(o['w_lead'], sc.ws[N - 2]) + 30.0)
            cq = [l2g(leg, gx[j], xq, yq) for xq in (x0s, x1s) for yq in (y0, y1f)]
            loc = fluid.intersection(box(min(q[0] for q in cq) - 5.0, min(q[1] for q in cq) - 5.0,
                                         max(q[0] for q in cq) + 5.0, max(q[1] for q in cq) + 5.0))
            prof = L.flare_profile(spec, sc.ws[N - 2], o['w_lead'], spec.L_out)
            ps_, pw_ = np.array([q[0] for q in prof]), np.array([q[1] for q in prof])
            s_st = o['x_e_local'] - spec.L_out                  # the stretch's start (the profile's s = 0)
            wst = np.array([found(xsec(loc, l2g(leg, gx[j], xq, y0), l2g(leg, gx[j], xq, y1f)),
                                  float(np.interp(xq - s_st, ps_, pw_)), N - 1, 'outlet flare', side) for xq in xst])
            s_ez = np.linspace(spec.L_out - EZ, spec.L_out, 61)
            R_ez = float(np.trapezoid([S.rr(x, h) for x in np.interp(s_ez, ps_, pw_)], s_ez))
            R_str = (np.trapezoid([S.rr(x, h) for x in wst], xst) + S.rr(wst[0], h) * (x0s - (o['x_e_local'] - spec.L_out))
                     + R_ez) * S.UM
            wsm = float(wst[0])
        else:
            R_str = S.rr(wsm, h) * spec.L_out * S.UM
        # The stretch end's excess (schedule.outlet_exp_c, which is the abrupt width change into the lead or the
        # flare's excess over its integral) takes its coefficient at the design widths, as the arcs' bend terms use the
        # design radii, so that drawing noise cannot invent a step, and it is scaled by the measured narrower width.
        wl_d, ws_d = float(o['w_lead']), float(sc.ws[N - 2])
        wn_m = wsm if wl_d >= ws_d else float(np.median(wl_s))
        c_end = S.outlet_exp_c(spec, ws_d, wl_d)
        R_end = S.rr(wn_m, h) * c_end * wn_m * S.UM
        R_str += R_end
        R_side[side] = R_str + Rl + Rb
        rec[side] = dict(stretch_ws=wsm, R_stretch=R_str, R_end=R_end, c_end=c_end, R_lead=Rl, R_bend=Rb,
                         w_lead_min=float(wl_s.min()), w_lead_max=float(wl_s.max()), w_lead_med=float(np.median(wl_s)))
    Rcl, wc_s = measure_path(fluid, o['paths']['C'], o['w_centre_lead'], h)
    if not np.all(wc_s > 0.0):                              # a center-lead section that found no fluid (E5)
        wc_s = np.array([found(x, o['w_centre_lead'], N - 1, 'centre lead', 'C') for x in wc_s])
        Rcl = measure_path_R(wc_s, o['paths']['C'], o['w_centre_lead'], h)
    rec['C'] = dict(R_lead=Rcl, w_lead_min=float(wc_s.min()), w_lead_max=float(wc_s.max()))
    Rxc, Rxs = o.get('R_ext_c', 0.0), o.get('R_ext_side', 0.0)      # off-chip tubing, at its design values (unmeasured)
    Rc_m = Rc_str + Rcl + Rxc
    RsA_m, RsB_m = R_side['A'] + Rxs, R_side['B'] + Rxs
    out['outlet'] = dict(Rc_out=Rc_m, RsA_out=RsA_m, RsB_out=RsB_m,
                         Rc_out_rel=Rc_m / sc.Rc_out - 1, RsA_out_rel=RsA_m / sc.Rs_out - 1,
                         RsB_out_rel=RsB_m / sc.Rs_out - 1, R_ext_c=Rxc, R_ext_side=Rxs, detail=rec)
    # The drawn center width at every gap (for the realized cutoff D_i) comes from the measured mid-segment widths,
    # taken as offsets from the design and interpolated to the gap positions.
    if out['wc']:
        km = np.array([t[0] for t in out['wc']])
        off_ = np.array([t[1] - t[2] for t in out['wc']])
        xm_ = 0.5 * (sc.x[km] + sc.x[np.minimum(km + 1, N - 1)])
        o_ = np.argsort(xm_)
        out['wc_gap'] = np.asarray(sc.wc, float) + np.interp(sc.x, xm_[o_], off_[o_])
    out.update(Rc=Rc, RsA=RsA, RsB=RsB, RgA=RgA, RgB=RgB)
    return out


def measure_path(fluid, pts, w, h, delta=4.0, n=41):
    """Return the resistance of a drawn lead along its centerline polyline, with the section widths. The polyline's
    first point is 1 um inside the upstream channel and its last point is 1 um inside the port. The lead is counted
    from 1 um after its first point to 1 um before its last, that is, from the upstream channel's end to the port's
    edge, because the port itself carries no resistance in the 1-D model. Widths are measured by normal cross-sections
    at n stations from 1 + 1.5 w after the first point to 1 + 1.5 w before the last, and the two end zones (the
    junction and port necks, where the fillets widen the lead) are counted at the design width. If a section finds no
    fluid (width <= 0), the resistance is nan, and the caller substitutes the design width there and recomputes with
    measure_path_R (E5)."""
    P = np.asarray(pts, float)
    seg = np.hypot(*np.diff(P, axis=0).T)
    s = np.concatenate([[0.0], np.cumsum(seg)])
    Lt = s[-1]
    ez = min(1.5 * w, 0.3 * (Lt - 2.0))          # end zones (junction, port neck); a short lead keeps its middle 40%
    s0, s1 = 1.0 + ez, Lt - 1.0 - ez
    st = np.linspace(s0, s1, n)
    wd = []
    for q in st:
        i = min(int(np.searchsorted(s, q, side='right')) - 1, len(seg) - 1)
        while seg[i] <= 1e-12:
            i -= 1
        t = (q - s[i]) / seg[i]
        c = P[i] + t * (P[i + 1] - P[i])
        d = (P[i + 1] - P[i]) / seg[i]
        nn = np.array([-d[1], d[0]])
        wd.append(xsec(fluid, tuple(c - nn * (w / 2 + delta)), tuple(c + nn * (w / 2 + delta))))
    wd = np.array(wd)
    if not np.all(wd > 0.0):
        return float('nan'), wd
    return measure_path_R(wd, pts, w, h, n=n), wd


def measure_path_R(wd, pts, w, h, n=41):
    """Return measure_path's resistance from given section widths, with the same stations, end zones, and formula."""
    P = np.asarray(pts, float)
    seg = np.hypot(*np.diff(P, axis=0).T)
    s = np.concatenate([[0.0], np.cumsum(seg)])
    Lt = s[-1]
    ez = min(1.5 * w, 0.3 * (Lt - 2.0))
    s0, s1 = 1.0 + ez, Lt - 1.0 - ez
    st = np.linspace(s0, s1, n)
    R = (np.trapezoid([S.rr(x, h) for x in wd], st) + S.rr(w, h) * ((s0 - 1.0) + (Lt - 1.0 - s1))) * S.UM
    return float(R)


def network(sc, m):
    """Solve the measured network m of the device that sc plans, and return the emergent F, each side's fractions fA
    and fB with their ratios to the planned f, and each gap's realized cutoff D_i = kappa_g x(f_i) w_c(i) on sides A
    and B (dA, dB) at the drawn center width. The measured outlet is used when m has one, and the design outlet
    otherwise."""
    o = m.get('outlet')
    if o is not None:                                   # as drawn, outlet stretch and leads measured on the drawing
        r = solve3(m['Rc'][:-1], m['RsA'], m['RsB'], m['RgA'], m['RgB'], o['Rc_out'], o['RsA_out'], o['RsB_out'])
    else:
        r = solve3(m['Rc'][:-1], m['RsA'], m['RsB'], m['RgA'], m['RgB'], sc.Rc_out, sc.Rs_out, sc.Rs_out)
    kap = S.KAPPA_G
    f = sc.f
    res = dict(F=r['F'], ratioA=r['fA'] / f, ratioB=r['fB'] / f)
    wcg = m.get('wc_gap', sc.wc)                       # the drawn center width at each gap
    res['dA'] = np.array([kap * S.x_of(fe) * wc if 0 < fe < 0.5 else float('nan') for fe, wc in zip(r['fA'], wcg)])
    res['dB'] = np.array([kap * S.x_of(fe) * wc if 0 < fe < 0.5 else float('nan') for fe, wc in zip(r['fB'], wcg)])
    res['wc_gap'] = np.asarray(wcg, float)
    res['fA'], res['fB'] = r['fA'], r['fB']
    return res
