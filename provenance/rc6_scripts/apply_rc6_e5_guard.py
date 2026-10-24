"""rc6 (Reviewer A's E5): measure2.measure_device's width guard, in place, each edit by exact string replacement (it
stops if a target is not found exactly once).  Every measured width that enters r(w) passes through found(): a positive
width is returned unchanged (a good drawing is measured exactly as before); a section that finds no fluid where the
plan puts the feature is recorded once per (row, kind, side) and its design width stands in, so the measurement
completes and pipeline.assess fails the build by name.  measure_path returns nan when a lead section is empty, and the
caller recomputes with the substituted widths (measure_path_R, the same stations and formula).

Re-applied on 2026-10-17 after the first rc6 working tree was lost (the guard as recorded in the rc6 reconstruction
record).  Usage: python3 apply_rc6_e5_guard.py OPENCIF_ROOT"""
import os
import sys

R = sys.argv[1]
P = os.path.join(R, 'opencif', 'engine', 'measure2.py')
s = open(P).read()


def edit(old, new, count=1):
    global s
    assert s.count(old) == count, (old[:80], s.count(old))
    s = s.replace(old, new)


edit("""    out = dict(units=[], lin=[], gaps=[], wc=[], loops=[])
""", """    out = dict(units=[], lin=[], gaps=[], wc=[], loops=[], not_found=[])
    # rc6 (E5, Reviewer A): a measured width that enters r(w) must be positive.  Where a section finds no fluid at the
    # place the plan puts a feature (a head unit, gap or side channel the drawing does not contain), the feature is
    # recorded once per (row, kind, side) in out['not_found'] and its design width stands in, so that the measurement
    # completes; pipeline.assess then fails the build by name ('every planned feature found on the drawing').  A
    # positive width is returned unchanged, so a good drawing is measured exactly as before.
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
""")
# head units: the median section width; the X unit's stem
edit("""            wm = float(np.median(ws)) if ws else w
""", """            wm = float(np.median(ws)) if ws else w
            wm = found(wm, w, k, 'head unit ' + u['kind'], side)
""")
edit("""                wd = np.array([xsec(fluid, (x, sg * y_in), (x, sg * (y_in + spec.GS + 10))) for x in xs])""",
     """                wd = np.array([found(xsec(fluid, (x, sg * y_in), (x, sg * (y_in + spec.GS + 10))), spec.GS, k,
                                     'head unit X stem', side) for x in xs])""")
# gaps
edit("""            wg = xsec(fluid, p0, p1)
            Rgk = S.r_gap(wg, lg, h, spec.c_ent)""", """            wg = found(xsec(fluid, p0, p1), wdes, k, 'gap', side)
            Rgk = S.r_gap(wg, lg, h, spec.c_ent)""")
# side segments: ordinary rows (Gauss-Legendre) and the transition rows
edit("""                        w_ = xsec(cl, l2g(leg, gx[j], xx - off, sg * (hw_ + lg + 1e-6)),
                                  l2g(leg, gx[j], xx - off, sg * (hw_ + lg + max(sc.ws[k], wR_) + 30.0)))""",
     """                        w_ = found(xsec(cl, l2g(leg, gx[j], xx - off, sg * (hw_ + lg + 1e-6)),
                                        l2g(leg, gx[j], xx - off, sg * (hw_ + lg + max(sc.ws[k], wR_) + 30.0))),
                                   sc.ws[k], k, 'side channel', side)""")
edit("""                    ws_.append(xsec(fluid, l2g(leg, gx[j], xx - off, sg * (hw_ + lg + 1e-6)),
                                    l2g(leg, gx[j], xx - off, sg * (hw_ + lg + max(sc.ws[k], wR_) + 30.0))))""",
     """                    ws_.append(found(xsec(fluid, l2g(leg, gx[j], xx - off, sg * (hw_ + lg + 1e-6)),
                                          l2g(leg, gx[j], xx - off, sg * (hw_ + lg + max(sc.ws[k], wR_) + 30.0))),
                                     wdraw_(k), k, 'side channel (transition row)', side))""")
# the centre at the post tip; the head's centre rows
edit("""        dw = float(min(devs)) + (hwm - 0.25 * (w0 + w1)) * 2.0       # deviation, referred to the row's mean width
""", """        dw = float(min(devs)) + (hwm - 0.25 * (w0 + w1)) * 2.0       # deviation, referred to the row's mean width
        if min(w0, w1) + dw <= 0.0:                          # rc6 (E5): no centre channel at the post tip
            found(0.0, 0.5 * (w0 + w1), k, 'centre', 'C')
            dw = 0.0
""")
edit("""        wcm = xsec(fluid, (xm, -spec.wc0 / 2 - 3), (xm, spec.wc0 / 2 + 3))""",
     """        wcm = found(xsec(fluid, (xm, -spec.wc0 / 2 - 3), (xm, spec.wc0 / 2 + 3)), spec.wc0, k, 'centre (head)', 'C')""")
# loops: the sides and the centre
edit("""                    out_w.append(xsec(fluid, q0, q1))""",
     """                    out_w.append(found(xsec(fluid, q0, q1), wdes, k, 'loop side', 'A' if sg > 0 else 'B'))""")
edit("""                ws_.append(xsec(fluid, q0, q1))
            return np.array(ws_)""", """                ws_.append(found(xsec(fluid, q0, q1), 2.0 * hwp, k, 'loop centre', 'C'))
            return np.array(ws_)""")
# outlet: the stretch (centre, side, flare) and the leads
edit("""    wce_m = xsec(fluid, l2g(leg, gx[j], xm, -hwe - 3), l2g(leg, gx[j], xm, hwe + 3))""",
     """    wce_m = found(xsec(fluid, l2g(leg, gx[j], xm, -hwe - 3), l2g(leg, gx[j], xm, hwe + 3)), sc.wc[-1], N - 1,
                  'outlet stretch centre', 'C')""")
edit("""        wsm = xsec(fluid, l2g(leg, gx[j], xm, y0), l2g(leg, gx[j], xm, y1))
        Rl, wl_s = measure_path(fluid, o['paths'][side], o['w_lead'], h)""",
     """        wsm = found(xsec(fluid, l2g(leg, gx[j], xm, y0), l2g(leg, gx[j], xm, y1)), sc.ws[N - 2], N - 1,
                    'outlet stretch side', side)
        Rl, wl_s = measure_path(fluid, o['paths'][side], o['w_lead'], h)
        if not np.all(wl_s > 0.0):                          # rc6 (E5): a lead section found no fluid
            wl_s = np.array([found(x, o['w_lead'], N - 1, 'side lead', side) for x in wl_s])
            Rl = measure_path_R(wl_s, o['paths'][side], o['w_lead'], h)""")
edit("""            wst = np.array([xsec(loc, l2g(leg, gx[j], xq, y0), l2g(leg, gx[j], xq, y1f)) for xq in xst])
            prof = L.flare_profile(spec, sc.ws[N - 2], o['w_lead'], spec.L_out)
            ps_, pw_ = np.array([q[0] for q in prof]), np.array([q[1] for q in prof])""",
     """            prof = L.flare_profile(spec, sc.ws[N - 2], o['w_lead'], spec.L_out)
            ps_, pw_ = np.array([q[0] for q in prof]), np.array([q[1] for q in prof])
            s_st = o['x_e_local'] - spec.L_out                  # the stretch's start (the profile's s = 0)
            wst = np.array([found(xsec(loc, l2g(leg, gx[j], xq, y0), l2g(leg, gx[j], xq, y1f)),
                                  float(np.interp(xq - s_st, ps_, pw_)), N - 1, 'outlet flare', side) for xq in xst])""")
edit("""    Rcl, wc_s = measure_path(fluid, o['paths']['C'], o['w_centre_lead'], h)
""", """    Rcl, wc_s = measure_path(fluid, o['paths']['C'], o['w_centre_lead'], h)
    if not np.all(wc_s > 0.0):                              # rc6 (E5): a centre-lead section found no fluid
        wc_s = np.array([found(x, o['w_centre_lead'], N - 1, 'centre lead', 'C') for x in wc_s])
        Rcl = measure_path_R(wc_s, o['paths']['C'], o['w_centre_lead'], h)
""")
# measure_path: nan on an empty section; measure_path_R: the same stations and formula from given widths
edit("""    those two end zones are counted at the design width.\"\"\"""",
     """    those two end zones are counted at the design width.  rc6 (E5): if a section finds no fluid (width <= 0) the
    resistance is nan; the caller substitutes the design width there and recomputes with measure_path_R.\"\"\"""")
edit("""    wd = np.array(wd)
    R = (np.trapezoid([S.rr(x, h) for x in wd], st) + S.rr(w, h) * ((s0 - 1.0) + (Lt - 1.0 - s1))) * S.UM
    return float(R), wd
""", """    wd = np.array(wd)
    if not np.all(wd > 0.0):
        return float('nan'), wd
    return measure_path_R(wd, pts, w, h, n=n), wd


def measure_path_R(wd, pts, w, h, n=41):
    \"\"\"measure_path's resistance from given section widths (the same stations, end zones and formula)\"\"\"
    P = np.asarray(pts, float)
    seg = np.hypot(*np.diff(P, axis=0).T)
    s = np.concatenate([[0.0], np.cumsum(seg)])
    Lt = s[-1]
    ez = min(1.5 * w, 0.3 * (Lt - 2.0))
    s0, s1 = 1.0 + ez, Lt - 1.0 - ez
    st = np.linspace(s0, s1, n)
    R = (np.trapezoid([S.rr(x, h) for x in wd], st) + S.rr(w, h) * ((s0 - 1.0) + (Lt - 1.0 - s1))) * S.UM
    return float(R)
""")
open(P, 'w').write(s)
print('applied to', P)
