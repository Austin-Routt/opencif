"""opencif.build: plan, draw, and validate one device from an engine Spec. This is the pipeline of the replica
generator's device.run, parameterized, with validators V1 to V8 as there and the per-row assessment of the realized
cutoff in addition.

Every check measures the drawn mask, that is, the fluid polygon written to the DXF. The measured network (V5) solves
the 3-rail Kirchhoff network on resistances integrated over the drawn polygons (the straight segments, the head units,
the loops, and the outlet leads) and returns each gap's emergent fraction f_i, and hence the realized cutoff
D_i = kappa_g x(f_i) w_c(i).
"""
import json
import math
import os
import pickle
import subprocess
import sys
import time

import numpy as np
from shapely.geometry import LineString, Polygon, box
from shapely.ops import unary_union

from opencif.physics import schedule as S  # noqa: E402
from opencif.geometry import shapes as Gm  # noqa: E402
from opencif.geometry import layout as L  # noqa: E402
from opencif.checks import measure as M2  # noqa: E402
from opencif.checks import auditor as AU  # noqa: E402
from opencif.checks import design_code as DCR  # noqa: E402
from opencif.physics import scenarios as RB  # noqa: E402
from opencif.geometry import dxf as io_dxf  # noqa: E402

DRC = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'checks', 'drc.py')
# The validity envelope. Its limits are those of the byte-frozen reference implementation's Tolerances (CIF_v2's M13
# gap dominance and M14 channel and Stokes gates), read from opencif.physics.kernels, which returns the design code's
# numbers bit for bit (tests/test_kernels.py). App. B's thesis cites them ("Re_stokes = 10^-2 on the gap and lamina and
# Re_max = 10 on the channel ... refuses a design that violates either (its M14)"). The kappa plateau edge u >= 0.035
# and the foot 0.012 are schedule.U_LO and U_FF, from the frozen Pairing (glossary \plateauUwc, \footUwc).
_TOL = S.C.Tolerances()
M13_MAX = _TOL.gap_dom          # 0.05
RE_MAX = _TOL.Re_max            # 10
RE_STOKES = _TOL.Re_stokes      # 1e-2
RHO_F = 1000.0                  # water (CIF_v2 Process.rho_f); MU is schedule.MU = 1e-3 Pa s
# Two supplied constants (glossary: external) come from spec v5 sec. 2.3 and CIF_v2's Pairing. The first is gamma_GS,
# the admissible fold width over gap GS/G, which is a pre-flight gate in spec v5 and CIF_v2 design(); here it is
# checked when a meander head is drawn, and a design outside it is refused in strict mode and reported otherwise. The
# second is f_soft, a soft upper bound on f (spec v5: "if max(f) > f_soft: WARN"), which is reported in every record
# and never gated.
_PAIR = S.C.Pairing()
GAMMA_GS = (float(_PAIR.gamma_GS[0]), float(_PAIR.gamma_GS[1]))   # (1.4, 1.5)
F_SOFT = float(_PAIR.f_soft)                                       # 5e-4
# Two gates come from CIF_v2's march() (spec v5 E11 and M16) with CIF_v2's own tolerances. They are the thin-lamina
# limit delta/h <= 0.10 and the quasi-developed taper |w_c(i) - w_c(i-1)| / Lambda <= 1e-2, with Lambda = G + g_p
# (CIF_v2's uniform pitch).
EPS_LAM = float(_TOL.eps_lam)       # 0.10
EPS_TAPER = float(_TOL.eps_taper)   # 1e-2
# The aspect-ratio floor is that of spec v5 sec. 18 and CIF_v2, Process.w_min = max(F_min, h/AR_max). AR_max is a
# supplied process limit (glossary: external), and opencif's process.aspect_ratio_max defaults to CIF_v2's value.
AR_MAX_DEFAULT = float(S.C.Process().AR_max)   # 10
# The limits of CIF_v2's embedded physics auditor (App. B "P8"; CIF_v2 _PHYSICS_SRC, PhysicsSpec.Re_max, .Q_in_uLmin,
# and check P8_stokes_gates) are its own inlet channel limit, its default operating point, and its entrance-length
# criterion 0.05 Re D_h <= Lambda/2. The port (opencif.checks.auditor, check P8) evaluates them on the drawn mask, and
# gates() gives their design-stage values. They are reported rather than gated (docs/LIMITATIONS.md 6), and
# tests/test_appendix_trace.py pins the values to the CIF_v2 source text.
P8_RE_MAX = AU.STRICT_RE_MAX
P8_Q_UL_MIN = AU.STRICT_Q_UL_MIN
# Dean (1928) gives the leading inertial correction to the flux of a curved pipe, 1 - 0.03058 (K/576)^2 with
# K = 2 Re^2 a/R, which is used here as an order-of-magnitude estimate for the rectangular ducts
# (opencif/physics/schedule.py, bend_ratio).
DEAN_C2 = 0.03058


# The kernel pair that every record declares (App. A Results 2 and 5-6; App. B C1-C3; glossary kappa_g) is the cubic
# skim kernel with its own prefactor, which is never mixed with kappa_F, together with the plateau edge and foot and the
# Hill curve used below the knee.
KERNEL_PAIR = dict(kernel='cubic', f_of_x='f = 3x^2 - 2x^3, x = delta/w_c (delta/w_c < 1/2)', cutoff='D = kappa_g delta*',
                   kappa_g=S.KAPPA_G, kappa_ci=[float(S.KAPPA_CI[0]), float(S.KAPPA_CI[1])],
                   plateau_u=S.U_LO, foot_u=S.U_FF, hill_below_knee=dict(kmax=S.HILL[0], u_half=S.HILL[1], n=S.HILL[2]),
                   source='opencif.physics.kernels.Pairing (= CIF_v2 v2.1.0 Pairing, reference/CIF_v2.py); Hill fit: glossary macros')


def hill_kappa(u):
    r"""Return the kappa(u) collapse curve below the knee, which is the frozen Hill fit of App. B P-C7.1,
    kappa(u) = kmax u^n / (u^n + u_half^n), with the glossary macros \hillKmax, \hillUhalf, and \hillN
    (schedule.HILL)."""
    kmax, uh, n = S.HILL
    u = np.asarray(u, float)
    return kmax * u ** n / (u ** n + uh ** n)


def ring_radii(xy, ds=1.0):
    """Return the local radius at every vertex of a closed ring, which is the radius of the circle through the vertex
    and the points ds before and after it along the ring."""
    P = xy[:-1]
    n = len(P)
    seg = np.hypot(*np.diff(np.vstack([P, P[:1]]), axis=0).T)
    s = np.concatenate([[0.0], np.cumsum(seg)])
    Lr = s[-1]
    if n < 3 or Lr < 4 * ds:
        return np.array([])

    def at(q):
        q = np.mod(q, Lr)
        i = np.clip(np.searchsorted(s, q, side='right') - 1, 0, n - 1)
        t = (q - s[i]) / np.where(seg[i] > 0, seg[i], 1.0)
        j = (i + 1) % n
        return P[i] + t[:, None] * (P[j] - P[i])
    A, B, C = at(s[:-1] - ds), P, at(s[:-1] + ds)
    a = np.hypot(*(B - C).T); b = np.hypot(*(A - C).T); c = np.hypot(*(A - B).T)
    cross = np.abs((B - A)[:, 0] * (C - A)[:, 1] - (B - A)[:, 1] * (C - A)[:, 0])
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.where(cross > 1e-12, a * b * c / (2 * cross), np.inf)


def min_curvature_radius(fluid, ds=1.0):
    """Return the smallest local radius over every ring of the fluid, the vertex where it occurs, and the number of
    vertices whose radius is below 1, 3, and 4.5 um."""
    rmin, where, n_lt = 1e9, None, {1.0: 0, 3.0: 0, 4.5: 0}
    for g in (list(fluid.geoms) if hasattr(fluid, 'geoms') else [fluid]):
        for ring in [g.exterior] + list(g.interiors):
            xy = np.asarray(ring.coords)
            R = ring_radii(xy, ds)
            if not len(R):
                continue
            i = int(np.argmin(R))
            if R[i] < rmin:
                rmin, where = float(R[i]), xy[i].tolist()
            for t in n_lt:
                n_lt[t] += int(np.sum(R < t))
    return rmin, where, {str(k): v for k, v in n_lt.items()}


def rounding_check(fluid, spec):
    """Compare the smallest local radius of curvature of the drawn boundary with the spec's criterion (V2)."""
    rmin, where, nlt = min_curvature_radius(fluid)
    return dict(min_local_radius_um=rmin, at=where, vertices_with_radius_below=nlt, criterion_um=spec.rmin_ok,
                ok=bool(rmin >= spec.rmin_ok))


def outlet_and_loops_measured(V, m, spec):
    """Write the outlet and U-turn blocks of the record from the measurement, and add the loops' narrowest width to the
    side-channel continuity block."""
    V['outlet_measured'] = {k: v for k, v in m['outlet'].items() if k != 'detail'}
    det_ = m['outlet'].get('detail') or {}
    # the leads' narrowest measured widths (the aspect-ratio inventory) and each side's stretch-end term
    V['outlet_measured']['lead_width_min_um'] = {k_: float(det_[k_]['w_lead_min']) for k_ in ('A', 'B', 'C')
                                                 if isinstance(det_.get(k_), dict) and 'w_lead_min' in det_[k_]}
    V['outlet_measured']['stretch_end'] = {k_: dict(c=float(det_[k_].get('c_end', 0.0)), R=float(det_[k_].get('R_end', 0.0)),
                                                    stretch_width_um=float(det_[k_]['stretch_ws']))
                                           for k_ in ('A', 'B') if isinstance(det_.get(k_), dict) and 'stretch_ws' in det_[k_]}
    V['loops_measured'] = [dict(row=int(lp['k']) + 1, A=dict(w_min=lp['A']['w_min'], R_rel=lp['A']['R'] / lp['target'] - 1),
                                B=dict(w_min=lp['B']['w_min'], R_rel=lp['B']['R'] / lp['target'] - 1),
                                Rc_rel_to_design=lp['Rc'] / lp['Rc_design'] - 1) for lp in m['loops']]
    if m['loops'] and 'side_continuity' in V:
        lwm = min(min(lp['A']['w_min'], lp['B']['w_min']) for lp in m['loops'])
        sc_v = V['side_continuity']
        sc_v['loop_narrowest_um'] = float(lwm)
        sc_v['narrowest_um'] = float(min(sc_v['narrowest_um'], lwm))
        sc_v['above_floor'] = bool(sc_v['narrowest_um'] >= spec.F_min - 0.05)


def drawn_parameters(m, sc_f, fluid, info, kL):
    """Return the quantities that the papers print, read back from the drawn mask as generated values for comparison
    with the printed ones in docs/REPRODUCTIONS.md (V8). RC3 called this block 'published'."""
    spec = sc_f.spec
    first_lin = [a for k, s_, a, b in m['lin'] if k == kL and s_ == 'A']
    last_lin = [a for k, s_, a, b in m['lin'] if k == sc_f.N - 2 and s_ == 'A']
    wc_list = sorted(m['wc'])
    return dict(depth_um=spec.h, G_meas=float(np.median([a for k, s_, a, b in m['gaps'] if k >= sc_f.NM])),
                G_meas_min=float(np.min([a for k, s_, a, b in m['gaps'] if k >= sc_f.NM])),
                wc_in_meas=float(wc_list[0][1]), wc_out_meas=float(wc_list[-1][1]),
                ws_first_linear_eq=float(first_lin[0]) if first_lin else None,
                ws_out_eq=float(last_lin[0]) if last_lin else None,
                footprint_mm=[(fluid.bounds[2] - fluid.bounds[0]) / 1000, (fluid.bounds[3] - fluid.bounds[1]) / 1000],
                array_length_mm=float(sc_f.x[-1] / 1000),
                head_length_mm=float((sc_f.x[sc_f.NM] - info['x_top']) / 1000) if sc_f.NM else 0.0)


def intended_and_declared(sc_f, info, m, spec, D_star, cut_tol, mode, seal):
    """Return the intended network (the closure plus the declared deliberate deviations) and the declared rows as
    (inten, nets, decl), both from the design and the drawing code's own profiles, before the drawn polygon is looked
    at."""
    inten = intended_network(sc_f, info, m)
    f_star_rows = np.array([S.f_of(spec.delta_star, w) for w in sc_f.wc])
    turns_ = sc_f.info.get('turns', [])
    nets = dict(full=inten)
    opts_ = active_options(spec, sc_f.f, f_star_rows, turns_, mode, seal)
    if OPT_PAPER in opts_:
        nets['no_paper'] = intended_network(sc_f, info, m, drop=('paper',))
    if OPT_PULL in opts_:
        nets['no_pull'] = intended_network(sc_f, info, m, drop=('pulled',))
    if OPT_PAPER in opts_ and OPT_PULL in opts_:
        nets['none'] = intended_network(sc_f, info, m, drop=('paper', 'pulled'))
    elif 'no_paper' in nets or 'no_pull' in nets:
        nets['none'] = nets.get('no_paper', nets.get('no_pull'))
    decl = declared_rows(spec, sc_f.f, f_star_rows, turns_, D_star, cut_tol, nets, mode, seal)
    return inten, nets, decl


def run_drc(dxf, outdir, name, spec):
    """Run the DRC (opencif/checks/drc.py) on the photomask file in its own process, at 0.1 um under the process
    floor, and return its result or the error."""
    d_um = spec.F_min - 0.1
    aj = os.path.join(outdir, name + '_audit.json')
    if os.path.exists(aj):
        os.remove(aj)                                 # never read a previous run's result
    r = subprocess.run([sys.executable, DRC, aj, dxf, '%.2f' % d_um, '--islands'], capture_output=True, text=True)
    try:
        if r.returncode != 0:
            raise RuntimeError('DRC exited %d' % r.returncode)
        return json.load(open(aj))
    except Exception as e:
        return dict(error='%s: %s' % (e, r.stderr[-2000:]))


def components_drawn(fluid, sc, info):
    """Return what was drawn. The posts in the inlet stretch (holes of the fluid between the inlet port and the head)
    and the outlet ports present (discs inside the fluid) are read from the drawing, and the head (its class, units,
    arms and their lane pairs, unit widths, and whether gap 1 is sealed), the outlet style, the legs, and the U-turns
    come from the drawn units and legs. `pipeline.assess` compares these with the request (the spec) rather than with
    the plan."""
    spec = sc.spec
    units = info['units']
    arms = [u for u in units if u['kind'] == 'A']
    inl = info['inlet']
    # The posts in the inlet stretch are counted as small solid islands (holes of the drawn fluid under 10^4 um^2)
    # whatever the inlet style. They are counted from the inlet port to 5 um before the head's top (or the first post),
    # or, with arms, to 5 um before gap 1, and within the inlet's own half-width (+10 um) and between the innermost arm
    # lanes.
    n_posts = 0
    if fluid.geom_type == 'Polygon':
        x_in = float(inl['x_port_c']) + float(inl.get('port_r', spec.port_r_in))
        y_lim = max(spec.wc0 / 2.0, float(spec.stem_half), float(inl.get('Hf') or 0.0)) + 10.0
        if arms:
            x_end, y_lim = -5.0, min(y_lim, spec.R0 - spec.w_arm / 2.0 - 1.0)
        else:
            x_end = float(info['x_top']) - 5.0
        for r in fluid.interiors:
            P_ = Polygon(r)
            c = P_.centroid
            if x_in < c.x < x_end and abs(c.y) < y_lim and P_.area < 1e4:
                n_posts += 1
    o = info['outlet']
    style = 'ring' if o.get('style') == 'ring' else ('teardrop' if o.get('kind') == 'combined' else 'separate')
    # Each distinct outlet port must be on the drawing, as a disc of 90% of its radius inside the fluid.
    ports = None
    if o.get('ports'):
        from shapely.geometry import Point
        pr = o.get('port_r') or {}
        seen, ports = set(), 0
        for k, (x, y) in o['ports'].items():
            key = (round(float(x), 3), round(float(y), 3))
            if key in seen:
                continue
            seen.add(key)
            if fluid.contains(Point(float(x), float(y)).buffer(0.9 * float(pr.get(k, spec.port_r)), quad_segs=16)):
                ports += 1
    f1 = float(sc.f_seed / sc.info['fstar']) if sc.NM > 0 else None
    # The head is a geometric description rather than a style name, because the request is compared separately
    # (pipeline.requested_vs_drawn).
    head = 'none'
    if sc.NM:
        ws_ = sorted({round(float(u['w']), 3) for u in units})
        narrow = getattr(spec, 'w_narrow', None) is not None and min(ws_) < float(spec.GS) - 1e-6
        head = ('%s%s%s' % ('narrow head (units narrower than GS %.1f um): ' % float(spec.GS) if narrow else '',
                            '%d arms (lane pairs %s) + hairpins' % (len(arms), [int(u['n']) for u in arms]) if arms else 'hairpins',
                            (', unit widths %s um' % ws_) if len(ws_) > 1 else ''))
    return dict(head_drawn=head, head_units=int(sc.NM), head_arms=len(arms), arm_lane_pairs=[int(u['n']) for u in arms],
                head_unit_widths_um=sorted({round(float(u['w']), 3) for u in units}),
                gap1_f_over_fstar=f1, gap1_sealed=(abs(f1 - 1.0) < 1e-9) if f1 is not None else True,
                filter_posts_counted=int(n_posts), outlet_style=style, outlet_ports=ports, legs=len(info['legs']),
                u_turns=len(info['loops']))


def topology_counts(fluid, sc, info):
    """Return whether the fluid is one polygon, its holes, and the number of holes that the plan expects, which is two
    per head unit, the inlet's posts, two per straight-leg row, two per loop, and the outlet's own."""
    single = fluid.geom_type == 'Polygon'
    nposts = sum((lg['rows'][1] - max(lg['rows'][0], sc.NM)) for lg in info['legs'])
    exp_holes = 2 * sc.NM + info['inlet'].get('n_posts', 0) + 2 * nposts + 2 * len(info['loops']) + info['outlet']['holes']
    holes = [Polygon(r) for r in fluid.interiors] if single else []
    return single, holes, exp_holes


def port_separation(fluid, sc, info):
    """Cut the outlet stretch and check that the product (center) port shares a component with neither filtrate port
    nor the inlet, and that no filtrate port shares one with the inlet. A combined outlet has one filtrate port
    (A = B)."""
    from shapely.geometry import Point
    spec = sc.spec
    o = info['outlet']
    leg, gx0 = info['legs'][-1], info['gx'][-1]
    xl = o['x_e_local'] - spec.L_out / 2.0
    # The cut reaches past the lead width too, because a flared stretch is wider.
    hw = sc.wc[-1] / 2.0 + spec.l_g + max(sc.ws[sc.N - 2], float(o.get('w_lead', 0.0) or 0.0)) + 30.0
    cut = L.to_global(box(xl - 5.0, -hw, xl + 5.0, hw), leg, gx0)
    parts = fluid.difference(cut)
    comps = list(parts.geoms) if hasattr(parts, 'geoms') else [parts]

    def comp_of(pt):
        P = Point(pt)
        return next((i for i, c in enumerate(comps) if c.contains(P)), None)
    ic, ia, ib = comp_of(o['ports']['C']), comp_of(o['ports']['A']), comp_of(o['ports']['B'])
    iin = comp_of((info['inlet']['x_port_c'], 0.0))
    ok = None not in (ic, ia, ib, iin) and ic not in (ia, ib, iin) and iin not in (ia, ib)
    if o.get('kind') == 'combined':
        ok = ok and ia == ib
    return dict(n_components_after_cut=len(comps), outlet_kind=o.get('kind', 'separate'), separated=bool(ok))


def side_continuity(fluid, sc, info, pad=15.0):
    """Check that every straight side channel is one connected strip of fluid on every leg and both sides, which the
    DRC's single-component test cannot show because the center channel joins everything through the gaps. The strip is
    the fluid beyond the post band, up to the drawn width + pad, from the first straight row of the leg to its last,
    and its narrowest cross-section, sampled every 0.5 um over the first 3 rows of the leg, every 2 um where the drawn
    width is within 2x the floor, and every 10 um elsewhere, is compared with the floor."""
    from opencif.checks import measure as M2_
    from shapely.geometry import Polygon as P_, LineString as LS_
    spec = sc.spec
    N, NM, Lam, lg = sc.N, sc.NM, spec.Lam_L, spec.l_g
    paper, wR, kL = L.transition_params(sc)
    pulled = L.pulled_rows(sc)
    wdr = np.array([L.drawn_side_width(sc, k, paper, wR, pulled) if (NM <= k < N - 1) else np.nan for k in range(N)])
    res = []
    for leg in info['legs']:
        j = leg['j']
        ka, kb = leg['rows']
        k0 = max(ka, NM)
        k1 = kb if leg['turn_after'] is not None else N - 1
        if k1 <= k0:
            continue
        off = leg['p0'] if leg['p0'] is not None else 0.0
        # The strip starts at the upstream edge of the first flow-carrying gap, because a sealed start's rounded dead
        # end upstream of gap 1 carries no flow and narrows to a point by construction, so it is not a channel.
        xa = sc.x[k0] - (spec.G / 2.0 if (k0 == 0 and spec.head == 'none') else 0.0)
        xb = sc.x[k1] + (0.0 if leg['turn_after'] is not None else 0.25 * spec.L_out)
        xs = np.linspace(xa, xb, max(200, int((xb - xa) / 2.0)))
        wloc = np.interp(xs, sc.x[k0:k1 + 1], np.nan_to_num(wdr[k0:k1 + 1], nan=float(np.nanmax(wdr[k0:k1 + 1]))))
        yin = L.wc_at(sc, xs) / 2.0 + lg + 0.3
        yout = yin + np.maximum(wloc, spec.F_min) + pad
        for sg, side in ((1.0, 'A'), (-1.0, 'B')):
            pts = [M2_.l2g(leg, info['gx'][j], xx - off, sg * yy) for xx, yy in zip(xs, yin)] + \
                  [M2_.l2g(leg, info['gx'][j], xx - off, sg * yy) for xx, yy in zip(xs[::-1], yout[::-1])]
            strip = P_(pts).buffer(0)
            piece = fluid.intersection(strip)
            comps = [c for c in (piece.geoms if hasattr(piece, 'geoms') else [piece]) if c.area > 1.0]
            # The narrowest section is found chunk by chunk, so that each cut works on a small polygon.
            wmin, at = 1e9, None
            nch = max(1, int((xb - xa) / (8 * Lam)))
            edges = np.linspace(xa, xb, nch + 1)
            for c0, c1 in zip(edges[:-1], edges[1:]):
                # steps of 0.5 um over the first 3 rows, 2 um where the drawn width is within 2x the floor, and 10 um
                # elsewhere
                w_here = float(np.min(np.interp([c0, c1], xs, wloc)))
                step = 0.5 if c0 < sc.x[min(k0 + 3, N - 1)] else (2.0 if w_here < 2 * spec.F_min else 10.0)
                xx = np.arange(c0 + 0.25, c1, step)
                if not len(xx):
                    continue
                ia, ib = M2_.l2g(leg, info['gx'][j], c0 - 1 - off, 0.0), M2_.l2g(leg, info['gx'][j], c1 + 1 - off, 0.0)
                yy_max = float(np.max(yout)) + 5
                bx = box(min(ia[0], ib[0]) - (yy_max if abs(ia[1] - ib[1]) > 1 else 0), min(ia[1], ib[1]) - yy_max,
                         max(ia[0], ib[0]) + (yy_max if abs(ia[1] - ib[1]) > 1 else 0), max(ia[1], ib[1]) + yy_max)
                loc = piece.intersection(bx)
                for x_ in xx:
                    yi = float(L.wc_at(sc, x_)) / 2.0 + lg + 0.3
                    yo = yi + float(np.interp(x_, xs, yout - yin))
                    w_ = loc.intersection(LS_([M2_.l2g(leg, info['gx'][j], x_ - off, sg * yi),
                                               M2_.l2g(leg, info['gx'][j], x_ - off, sg * yo)])).length + 0.3
                    if w_ < wmin:
                        wmin, at = w_, (float(x_), M2_.l2g(leg, info['gx'][j], x_ - off, sg * (yi + 0.5 * w_)))
            res.append(dict(leg=j + 1, side=side, rows=[k0 + 1, k1 + 1], pieces=len(comps), connected=len(comps) == 1,
                            narrowest_um=float(wmin), narrowest_row=int(np.searchsorted(sc.x, at[0])) if at else None,
                            narrowest_xy=at[1] if at else None))
    ok = all(r['connected'] for r in res)
    wn = min(r['narrowest_um'] for r in res) if res else None
    # The 0.05 um allowance is the cross-section sampling's own resolution.
    return dict(ok=bool(ok), narrowest_um=wn, floor_um=spec.F_min, above_floor=bool(wn is None or wn >= spec.F_min - 0.05),
                legs=res)


def intended_network(sc, info, m, drop=()):
    """Solve the network that the design intends to draw, for the sealed device. The network takes the closure's
    resistances everywhere except at the declared deliberate deviations, which take their intended drawn values (for
    the pull-in, the pulled widths and the loops' drawn width profiles; for the papers' rectilinear transition rows,
    the rows as drawn, because they have no separate intent). Comparing the drawn mask's network with this one tells,
    row by row, whether the mask draws its design, and comparing this one with the targets tells whether the design
    realizes D and F. The latter comparison also predicts, before the drawing is looked at, where a declared departure
    such as a floor start or a seed will move the cutoff."""
    from opencif.checks import measure as M2_
    spec = sc.spec
    N, NM = sc.N, sc.NM
    paper, wR, kL = L.transition_params(sc)
    pulled = L.pulled_rows(sc) if 'pulled' not in drop else set()
    if 'paper' in drop:
        paper = False                                     # the closure at the transition rows instead
    turns = set(sc.info.get('turns', []))
    RsA = np.array(sc.Rs, float).copy()
    RsB = RsA.copy()
    for k in range(NM, N - 1):
        if k in turns:
            continue
        if paper and NM <= k <= kL + 1:
            continue                                      # the drawn transition profile, below
        elif k in pulled:
            w = L.drawn_side_width(sc, k, False, wR, pulled)
            RsA[k] = RsB[k] = S.rr(w, spec.h) * sc.ell[k] * S.UM
    loops_pulled = (getattr(spec, 'loop_pull_in', 0.0) or 0.0) > 0.0
    for lp in info['loops']:
        k = lp['k']
        if 'pulled' in drop and loops_pulled:
            # Xia's pull-in narrows the loops too (layout.loop_profile subtracts it), so without the option the loop
            # rows take their closure targets, which the unpulled loop rule meets exactly.
            RsA[k] = RsB[k] = float(sc.Rs[k])
            continue
        RsA[k] = L.loop_R(sc, k, bool(lp['ccw']), lp['fA'], n_str=60, n_arc=600)
        RsB[k] = L.loop_R(sc, k, not bool(lp['ccw']), lp['fB'], n_str=60, n_arc=600)
    RgA = np.array(sc.Rg, float)
    RgB = RgA.copy()
    Rc = np.array(sc.Rc, float)
    if paper:
        # The side widths of the papers' rectilinear transition (reproduction-only) are the profile that the layout
        # draws (ws_rect beside the wall blocks, and the smoothstep narrowing to the first post) rather than the
        # closure's. The intent is that profile as the drawing code laid it down, before any union or smoothing
        # (layout.SIDE_PROFILES), integrated row by row, and the gaps and the center keep their design values. Nothing
        # here is read from the drawn polygon.
        leg0 = info['legs'][0]
        for iw, Rs_ in ((0, RsA), (1, RsB)):
            prof = np.array(L.SIDE_PROFILES[(leg0['j'], iw)], float)
            for k in range(NM, min(N - 1, kL + 2)):
                if k <= kL:                               # the measurement's own quadrature (measure.trans_nodes)
                    xs, wq = M2_.trans_nodes(sc.x[k], sc.x[k + 1])
                else:                                     # row kL+1 is an ordinary row, with 2 x 6 Gauss-Legendre nodes
                    xs, wq = M2_.trans_nodes(sc.x[k], sc.x[k + 1], 2)
                w = np.interp(xs, prof[:, 0], prof[:, 1])
                Rs_[k] = float(np.sum(wq * np.array([S.rr(ww, spec.h) for ww in w]))) * S.UM
    r = M2_.solve3(Rc[:-1], RsA, RsB, RgA, RgB, sc.Rc_out, sc.Rs_out, sc.Rs_out)
    wc = np.asarray(sc.wc, float)
    dA = np.array([S.KAPPA_G * S.x_of(fe) * w if 0 < fe < 0.5 else np.nan for fe, w in zip(r['fA'], wc)])
    dB = np.array([S.KAPPA_G * S.x_of(fe) * w if 0 < fe < 0.5 else np.nan for fe, w in zip(r['fB'], wc)])
    return dict(F=float(r['F']), fA=np.asarray(r['fA'], float), fB=np.asarray(r['fB'], float), dA=dA, dB=dB)


OPT_SEED = 'seed: gap 1 above f*'
OPT_DEPICTED = 'head drawn as depicted'
OPT_PAPER = "papers' rectilinear transition"
OPT_PULL = 'loop pull-in'
OPT_BOTH = "papers' transition + loop pull-in"
OPT_FLOOR = 'floor start'


def active_options(spec, f, fstar_cut, turns, mode, seal):
    """Return the declared options in force, each of which deliberately moves some gaps off D. They apply in reproduce
    mode only, except for a seed explicitly allowed in evaluate mode (strict mode refuses `seal: allow` at the spec)."""
    opts = []
    if abs(f[0] / fstar_cut[0] - 1.0) > 1e-9 and (mode == 'reproduce' or seal == 'allow'):
        opts.append(OPT_SEED)
    if mode != 'reproduce':
        return opts
    if spec.f_head is not None:
        opts.append(OPT_DEPICTED)
    if getattr(spec, 'transition', 'exact') == 'paper' and spec.head != 'none':
        opts.append(OPT_PAPER)
    if (getattr(spec, 'loop_pull_in', 0.0) or 0.0) > 0 and len(turns):
        opts.append(OPT_PULL)
    if (getattr(spec, 'S0', 0.0) or 0.0) > 0:
        opts.append(OPT_FLOOR)
    return opts


def _outside(net, D_target, tol):
    """Return a mask of the rows of net that have no cutoff on either side, or whose defined cutoffs differ from
    D_target by more than the relative tolerance tol. A side without a cutoff is ignored when the other side has one."""
    hi = np.fmax(net['dA'], net['dB']) / D_target
    lo = np.fmin(net['dA'], net['dB']) / D_target
    return ~(np.isfinite(hi) & np.isfinite(lo) & (np.abs(hi - 1) <= tol) & (np.abs(lo - 1) <= tol))


def declared_rows(spec, f, fstar_cut, turns, D_target, tol, nets, mode, seal):
    """Return the rows that each declared option moves off D +- tol, attributed row by row before the drawing is
    measured. The attribution uses intended networks built from the design and the drawing code's own profiles, never
    from the drawn polygon:

      nets['full']      every deliberate deviation in force
      nets['no_paper']  the same with the papers' transition replaced by the closure (when that option is in force)
      nets['no_pull']   the same with the pull-in replaced by the closure
      nets['none']      both replaced

    A row outside D +- tol in the full network is declared under the papers' transition or the pull-in when removing
    that option alone brings it inside, and under both together when only removing both does. Otherwise, when it is
    still outside without them, it is declared under the seed (gap 1 only), a depicted head (its own rows only), or a
    floor start (the leading run of departing rows from gap 1, which ends at the first 10 consecutive rows inside).
    Else it is declared under nothing and fails Q2. Thus a fault anywhere, of any size, that no option explains is
    never exempt."""
    opts = active_options(spec, f, fstar_cut, turns, mode, seal)
    if not opts:
        return {}
    N = len(f)
    full = _outside(nets['full'], D_target, tol)
    no_p = _outside(nets['no_paper'], D_target, tol) if 'no_paper' in nets else full
    no_l = _outside(nets['no_pull'], D_target, tol) if 'no_pull' in nets else full
    none = _outside(nets['none'], D_target, tol) if 'none' in nets else full
    run_end = 0
    if OPT_FLOOR in opts:
        ins = ~none
        for i in range(N - 10):
            if ins[i:i + 10].all():
                run_end = i
                break
    out = {o: [] for o in opts}
    if OPT_PAPER in opts and OPT_PULL in opts:
        out[OPT_BOTH] = []
    for r in np.where(full)[0]:
        r = int(r)
        if OPT_PAPER in opts and not no_p[r]:
            out[OPT_PAPER].append(r)
        elif OPT_PULL in opts and not no_l[r]:
            out[OPT_PULL].append(r)
        elif OPT_PAPER in opts and OPT_PULL in opts and not none[r]:
            out[OPT_BOTH].append(r)
        elif OPT_SEED in opts and r == 0:
            out[OPT_SEED].append(r)
        elif OPT_DEPICTED in opts and r < len(spec.f_head):
            out[OPT_DEPICTED].append(r)
        elif OPT_FLOOR in opts and r < run_end:
            out[OPT_FLOOR].append(r)
    return {k: v for k, v in out.items() if v or k in opts}


def drawing_fidelity(nw, inten, fstar_cut, N, D_target=None):
    """Compare the drawn network with the intended network row by row, judging each row on its own scale (Q1). Where
    both define a cutoff, the deviation is

        |D_drawn - D_intended| / max(D_intended, D)

    so that a row at or below the target is judged in units of D (and a row far below D is not magnified), whereas a
    deliberate departure far above D (a seed, or the papers' transition, as in Dinh's row 258 at 2.9 D) is judged in
    units of its own intended cutoff. Where the intended row has no cutoff (reverse flow, f >= 1/2, or near-zero f),
    the deviation is

        |f_drawn - f_intended| / max(|f_int|, f*)

    The deviation of each row's own cutoff and the raw f deviation in units of f* are reported too."""
    dd = np.concatenate([np.asarray(nw['dA']), np.asarray(nw['dB'])])
    di = np.concatenate([np.asarray(inten['dA']), np.asarray(inten['dB'])])
    fs = np.concatenate([fstar_cut, fstar_cut])
    fi = np.concatenate([inten['fA'], inten['fB']])
    fd = np.concatenate([np.asarray(nw['fA']), np.asarray(nw['fB'])])
    both = np.isfinite(dd) & np.isfinite(di) & (np.abs(fi) > 0.05 * fs)
    Dt = D_target if D_target else np.nanmedian(di)
    di_ = np.where(both, di, Dt)
    dev_D = np.where(both, np.abs(dd - di_) / np.maximum(di_, Dt), 0.0)
    dev_own = np.where(both, np.abs(dd / di_ - 1), 0.0)
    dev_f = np.where(~both, np.abs(fd - fi) / np.maximum(np.abs(fi), fs), 0.0)
    dev_f_star = np.where(~both, np.abs(fd - fi) / fs, 0.0)
    iD, iF, iO, iS = (int(np.argmax(v)) for v in (dev_D, dev_f, dev_own, dev_f_star))
    return dict(cutoff_dev_max=float(dev_D[iD]), cutoff_dev_row=iD % N + 1,
                own_cutoff_rel_dev_max=float(dev_own[iO]), own_cutoff_rel_dev_row=iO % N + 1,
                f_dev_max_where_no_cutoff=float(dev_f[iF]), f_dev_row=iF % N + 1,
                f_dev_over_fstar_max_where_no_cutoff=float(dev_f_star[iS]), f_dev_over_fstar_row=iS % N + 1,
                F_drawn=float(nw['F']), F_intended=float(inten['F']), F_abs_dev=float(abs(nw['F'] - inten['F'])))


def realised_cutoff(sc, nw, D_target, tol=0.01, declared=None, extra_tols=(0.01, 0.02, 0.05)):
    """Return the realized cutoff on the drawn mask row by row, from the measured network's emergent f_i (both sides)
    and the drawn center width. Every row counts, so a row with reverse flow (f < 0) or with f >= 1/2 has no cutoff
    and is a departure in its own right. Shares are of the net filtrate (extraction minus return), and the rows in
    `declared` are reported apart."""
    fA, fB = np.asarray(nw['fA'], float), np.asarray(nw['fB'], float)
    dA, dB = np.asarray(nw['dA'], float), np.asarray(nw['dB'], float)
    N = len(fA)
    keep = np.cumprod(1.0 - fA - fB)
    Pin = np.concatenate([[1.0], keep[:-1]])
    q = (fA + fB) * Pin                                                # the net flow leaving the center at each row
    Fnet = float(np.sum(q))
    hi = np.fmax(dA, dB) / D_target
    lo = np.fmin(dA, dB) / D_target
    undefined = ~np.isfinite(dA) | ~np.isfinite(dB)                   # f < 0 (reverse) or f >= 1/2 on a side
    reverse = (fA < 0) | (fB < 0)
    dec = np.zeros(N, bool)
    for rows in (declared or {}).values():
        dec[[r for r in rows if 0 <= r < N]] = True
    out = dict(D_target_um=float(D_target), tolerance=tol, F_net=Fnet,
               reverse_flow_rows=[int(i) + 1 for i in np.where(reverse)[0]][:50],
               reverse_flow_share_of_inlet=float(-np.sum(q[q < 0])),
               declared_rows={k: [min(v) + 1, max(v) + 1, len(v)] if v else [] for k, v in (declared or {}).items()})
    fin = ~undefined
    out['D_realised_max_um'] = float(np.nanmax(np.where(fin, np.fmax(dA, dB), np.nan)))
    out['D_realised_min_um'] = float(np.nanmin(np.where(fin, np.fmin(dA, dB), np.nan)))
    out['row1_cutoff_over_D'] = float(hi[0]) if np.isfinite(hi[0]) else None
    # Rows that extract f >= 1/2 on a side have no cutoff and pass every particle up to the gap (codereview2). The
    # count is of f >= 1/2 itself rather than of every row without a cutoff, because a side at f = 0 extracts nothing
    # (codereview3 B2).
    out['rows_f_half_or_more'] = int(np.sum((fA >= 0.5) | (fB >= 0.5)))
    for t in sorted(set(extra_tols) | {tol}):
        inside = fin & (np.abs(hi - 1) <= t) & (np.abs(lo - 1) <= t)
        key = '%g%%' % (100 * t)
        out['share_of_F_within_' + key] = float(np.sum(q[inside]) / Fnet)
        out['share_of_F_oversize_' + key] = float(np.sum(np.clip(q[fin & (hi > 1 + t)], 0, None)) / Fnet)
        out['rows_outside_' + key] = int(np.sum(~inside))
    bad = ~(fin & (np.abs(hi - 1) <= tol) & (np.abs(lo - 1) <= tol))
    und = bad & ~dec
    out['undeclared_rows_outside'] = int(np.sum(und))
    out['undeclared_rows_list'] = [int(i) + 1 for i in np.where(und)[0]][:60]
    out['undeclared_worst_over'] = float(np.nanmax(np.where(und & fin, hi, np.nan))) if np.any(und & fin) else None
    out['undeclared_worst_under'] = float(np.nanmin(np.where(und & fin, lo, np.nan))) if np.any(und & fin) else None
    out['undeclared_undefined_rows'] = int(np.sum(und & undefined))
    per = {}
    for name_, rows_ in (declared or {}).items():
        mk = np.zeros(N, bool)
        mk[[r for r in rows_ if 0 <= r < N]] = True
        mk &= bad
        per[name_] = dict(rows_outside=int(np.sum(mk)),
                          first_last=[int(np.where(mk)[0][0]) + 1, int(np.where(mk)[0][-1]) + 1] if mk.any() else [],
                          worst_over=float(np.nanmax(np.where(mk & fin, hi, np.nan))) if np.any(mk & fin) else None,
                          worst_under=float(np.nanmin(np.where(mk & fin, lo, np.nan))) if np.any(mk & fin) else None,
                          undefined_rows=int(np.sum(mk & undefined)),
                          share_of_F_oversize=float(np.sum(np.clip(q[mk & fin & (hi > 1 + tol)], 0, None)) / Fnet))
    out['declared_by_option'] = per
    out['declared_rows_outside'] = int(np.sum(bad & dec))
    out['declared_worst_over'] = float(np.nanmax(np.where(bad & dec & fin, hi, np.nan))) if np.any(bad & dec & fin) else None
    out['declared_share_of_F_oversize'] = float(np.sum(np.clip(q[bad & dec & fin & (hi > 1 + tol)], 0, None)) / Fnet)
    return out


def model_scenarios(sc, info, m, nw, D_target, tol):
    """Solve scenarios on the measured network for what the 1-D model cannot resolve (reviews R2-4, R2-7, and R3-1).

    The first scenarios concern the junction physics. Every arc (the U-turn center and sides, the head hairpins and
    arms, the transition unit, and the outlet leads) carries the exact fully developed Stokes curved-duct correction
    less its two junctions (kappa_j w), every head unit ends in a T-junction counted as c_T w, and every outlet stretch
    whose end changes width carries the excess c_exp w_narrow (schedule.outlet_exp_c). In the schedule, each of these
    terms is interpolated between its deep and shallow 2-D values by a depth weight lam(h/w). The band re-solves the
    drawn device with lam moved to both ends of its range (all terms together, with c_T also widened and c_exp also
    scaled +-30%) and, for information, also with the junction terms ignored, with the T-junction term alone ignored,
    and with no bend correction. Dean (inertial) terms are not modeled, because Dean's leading term is ~1.5e-4 of a
    bend's resistance at Re = 10 and a/R = 0.2, and 3.7e-3 at a/R = 1, the most that any drawable arc reaches. gates()
    reports the bound for every bend of the device (dean_leading_term_bound, from a Reynolds bound that holds in every
    channel) where that Reynolds bound is <= 10.

    The off-chip outlet resistances are design values (often assumed), and the last rows can depend on them strongly,
    so the tubing length (x0.5, x2) and the center I.D. (+-0.025 mm) are varied. Finally, the outlet balance moves
    both side outlets' resistance by +-1%, for what the 1-D outlet model does not count."""
    import copy
    from opencif.checks import measure as M2_
    from opencif.geometry import shapes as Gm_
    spec = sc.spec
    o = m['outlet']
    N = sc.N
    wcg = np.asarray(m.get('wc_gap', sc.wc), float)
    base = dict(Rc=np.array(m['Rc'], float), RsA=np.array(m['RsA'], float), RsB=np.array(m['RsB'], float),
                RgA=np.array(m['RgA'], float), RgB=np.array(m['RgB'], float),
                Rco=float(o['Rc_out']), RAo=float(o['RsA_out']), RBo=float(o['RsB_out']))

    def solve(e):
        r = M2_.solve3(e['Rc'][:-1], e['RsA'], e['RsB'], e['RgA'], e['RgB'], e['Rco'], e['RAo'], e['RBo'])
        d = lambda ff: np.array([S.KAPPA_G * S.x_of(fe) * w if 0 < fe < 0.5 else np.nan for fe, w in zip(ff, wcg)])
        return d(r['fA']) / D_target, d(r['fB']) / D_target, float(r['F'])

    a0, b0 = np.asarray(nw['dA']) / D_target, np.asarray(nw['dB']) / D_target
    in0 = np.isfinite(a0) & np.isfinite(b0) & (np.abs(np.fmax(a0, b0) - 1) <= tol) & (np.abs(np.fmin(a0, b0) - 1) <= tol)
    fA0, fB0 = np.asarray(nw['fA'], float), np.asarray(nw['fB'], float)
    q0 = (fA0 + fB0) * np.concatenate([[1.0], np.cumprod(1.0 - fA0 - fB0)[:-1]])
    F0 = float(np.sum(q0))

    def summary(a, b, F_):
        hi, lo = np.fmax(a, b), np.fmin(a, b)
        fin = np.isfinite(a) & np.isfinite(b)
        outside = ~(fin & (np.abs(hi - 1) <= tol) & (np.abs(lo - 1) <= tol))
        # Changes are counted only on rows inside D +- tol on the drawn mask, because declared departures, rows with no
        # cutoff, and near-zero fractions would magnify a relative change (review R3 minor 4).
        mv = np.where(in0, np.fmax(np.abs(a - a0), np.abs(b - b0)), np.nan)      # in units of D
        lost = in0 & ~fin                                # an in-tolerance row that loses its cutoff (reverse, f >= 1/2)
        newly = outside & in0
        fin_mv = np.where(in0 & fin, mv, -1.0)
        i = int(np.argmax(fin_mv))
        return dict(F=F_, rows_outside_tol=int(outside.sum()), rows_newly_outside_tol=int(newly.sum()),
                    # every row, with no cap (RC4 kept 40, and RC5's band has 47)
                    newly_outside_rows=[int(r) + 1 for r in np.where(newly)[0]],
                    rows_losing_their_cutoff=int(lost.sum()),
                    share_of_F_in_those_rows=float(np.sum(np.abs(q0[newly])) / F0) if F0 else None,
                    rows_moved_over_half_tol=int(np.sum((fin_mv > tol / 2)) + lost.sum()),
                    max_cutoff_change=float(fin_mv[i]) if fin_mv[i] >= 0 else None, max_change_row=i + 1,
                    cutoff_over_D_range=[float(np.nanmin(lo)), float(np.nanmax(hi))])

    def variant(**kw):
        sp_ = copy.copy(spec)
        for k_, v_ in kw.items():
            setattr(sp_, k_, v_)
        sc_ = copy.copy(sc)
        sc_.spec = sp_
        return sc_

    def elements_for(scv):
        """Return the measured elements with every arc's bend term re-evaluated under scv.spec, applied on top of the
        measurement (as a ratio for the loops and as an added excess for the head units and the outlet leads)."""
        spv = scv.spec
        e = {k_: (v.copy() if isinstance(v, np.ndarray) else v) for k_, v in base.items()}
        h = spec.h
        for lp in info.get('loops', []):
            k = int(lp['k'])
            R = S.R_of(sc, k)
            wm = 0.5 * (float(sc.wc[k]) + float(sc.wc[min(k + 1, N - 1)]))
            Lb = spec.Lam_L + math.pi * R + S.arc_excess_len(spec, R, math.pi, wm)
            Lv = spec.Lam_L + math.pi * R + S.arc_excess_len(spv, R, math.pi, wm)
            e['Rc'][k] *= Lv / Lb
            for side, fw in (('RsA', lp['fA']), ('RsB', lp['fB'])):
                inner = bool(lp['ccw']) if side == 'RsA' else not bool(lp['ccw'])
                e[side][k] *= L.loop_R(scv, k, inner, fw, n_str=60, n_arc=600) / L.loop_R(sc, k, inner, fw, n_str=60, n_arc=600)
        units = {int(u['k']): u for u in info.get('units', [])}
        for rec in m.get('units', []):
            u = units.get(int(rec['k']))
            if u is None:
                continue
            p = Gm_.unit_path(spec, u, from_J=True)
            for side, key in (('A', 'RsA'), ('B', 'RsB')):
                wm = float(rec[side]['w'])
                dl = (Gm_.bend_excess(spv, p, wm) + Gm_.junction_excess(spv, u, wm)
                      - Gm_.bend_excess(spec, p, wm) - Gm_.junction_excess(spec, u, wm))
                e[key][int(rec['k'])] += S.rr(wm, h) * dl * S.UM
        # The outlet leads' arcs are the design outlet's (info['outlet']), and their width is the measured median lead
        # width, or the design width when that is missing. RC4 read both from the measured outlet record (m['outlet']),
        # which has neither, so its junction band never moved the leads' arc terms.
        oi = info.get('outlet') or {}
        arcs = oi.get('bend_arcs') or {}
        det = o.get('detail') or {}
        for side, key in (('A', 'RAo'), ('B', 'RBo')):
            if arcs.get(side) and side in det:
                wl = float(det[side].get('w_lead_med', oi.get('w_lead', 0.0)) or 0.0)
                if wl > 0:
                    e[key] += S.rr(wl, h) * (L.arcs_excess(spv, arcs[side], wl) - L.arcs_excess(spec, arcs[side], wl)) * S.UM
            if side in det:                                       # the stretch end's term moves with the band
                ws_d, wl_d = float(sc.ws[N - 2]), float(oi.get('w_lead', 0.0) or 0.0)
                wn_m = float(det[side].get('stretch_ws', 0.0) or 0.0) if wl_d >= ws_d else \
                    float(det[side].get('w_lead_med', 0.0) or 0.0)
                if wl_d > 0 and wn_m > 0:
                    e[key] += S.rr(wn_m, h) * (S.outlet_exp_c(spv, ws_d, wl_d) - S.outlet_exp_c(spec, ws_d, wl_d)) * wn_m * S.UM
        return e

    out = {}
    on = getattr(spec, 'bend_model', 'stokes') != 'none'
    wref = float(spec.GS)
    out['bend_model'] = dict(model=getattr(spec, 'bend_model', 'stokes'),
                             junction_lambda_at_GS=S.depth_lambda(spec.h, wref), kappa_j_at_GS=S.kappa_j(spec, wref),
                             c_T_at_GS=S.c_t(spec, wref),
                             limits=dict(kappa_j=[S.KAPPA_J_DEEP, S.KAPPA_J_SHALLOW], c_T=[S.C_T_DEEP, S.C_T_SHALLOW]))
    if getattr(spec, 'turns', 'arc') == 'square':          # the square corners' term, as the record's (DINH-1)
        c_, dp_, sh_, lam_, held_ = S.c_corner(spec, wref, parts=True)
        spl_, sph_ = copy.copy(spec), copy.copy(spec)
        spl_.junction_band, sph_.junction_band = 'lo', 'hi'
        out['bend_model']['square_corner'] = dict(
            c_corner_at_GS=c_, deep=dp_, shallow=sh_, depth_lambda=lam_, band=[S.c_corner(sph_, wref), S.c_corner(spl_, wref)],
            r_out_um=float(spec.r_corner_out), r_in_um=float(spec.r_corner_in), r_in_over_w=float(spec.r_corner_in) / wref,
            held_at_table_end=bool(held_), table_r_in_over_w=list(S.CORNER_RIN),
            r_out_over_w=float(spec.r_corner_out) / wref, r_out_solved_to_over_w=S.CORNER_ROUT_SOLVED,
            round_area_um2=float(Gm_.corner_round_area(spec)),
            radii_note='rc7 (code review R1): the radii drawn and modelled, max(spec radius, head fillet)',
            junctions_ignored_note=('rc7 (code review R14): the junctions_ignored scenario (junction_band none) drops the '
                                    'square corners\' term too, while an arc head keeps its developed bend factor there: '
                                    'for a square-turn head that scenario is not like-for-like with an arc head\'s'),
            source='provenance/rc7_scripts/corner/corner_table.json (2-D Hele-Shaw and Stokes solves); the excess over the '
                   'centreline through the corner point, per corner, in units of w')
    if on:
        for lab, kw in (('junction_band_lo', dict(junction_band='lo')), ('junction_band_hi', dict(junction_band='hi')),
                        ('t_junction_ignored', dict(t_junction=False)), ('junctions_ignored', dict(junction_band='none')),
                        ('bend_correction_off', dict(bend_model='none'))):
            out[lab] = summary(*solve(elements_for(variant(**kw))))
        lp_info = []
        for lp in info.get('loops', [])[:8]:
            k = int(lp['k'])
            R = S.R_of(sc, k)
            wm = 0.5 * (float(sc.wc[k]) + float(sc.wc[min(k + 1, N - 1)]))
            lp_info.append(dict(row=k + 1, R_um=float(R), centre_conductance_ratio=float(S.bend_ratio(R - wm / 2, R + wm / 2, spec.h)),
                                centre_excess_um=float(S.arc_excess_len(spec, R, math.pi, wm))))
        out['bend_correction_per_loop'] = lp_info
    # The outlet balance. The last rows' gap flows are a small difference of the outlet pressures, so they move with any
    # error in the side outlet's resistance that the 1-D model does not count (the step from the side channel into a
    # wider lead, the outlet bifurcation's 3-D flow, and the port). Both side outlets are therefore moved by +-1% (an
    # illustrative size, as in the other scenarios), and the last row's sensitivity d(D_N/D)/d ln R_s,out is recorded.
    def side_out(ks):
        e = dict(base)
        e['RAo'], e['RBo'] = base['RAo'] * ks, base['RBo'] * ks
        return e
    out['outlet_side_resistance_+1%'] = summary(*solve(side_out(1.01)))
    out['outlet_side_resistance_-1%'] = summary(*solve(side_out(0.99)))
    a1, b1, _ = solve(side_out(1.0 + 1e-3))
    out['dDlast_dlnRside_out'] = float((np.fmax(a1[N - 1], b1[N - 1]) - np.fmax(a0[N - 1], b0[N - 1])) / 1e-3)
    Rxc, Rxs = float(o.get('R_ext_c', 0.0) or 0.0), float(o.get('R_ext_side', 0.0) or 0.0)
    if Rxc > 0 or Rxs > 0:
        def tub(kc, ks):
            e = dict(base)
            e['Rco'] = base['Rco'] + (kc - 1) * Rxc
            e['RAo'] = base['RAo'] + (ks - 1) * Rxs
            e['RBo'] = base['RBo'] + (ks - 1) * Rxs
            return summary(*solve(e))
        out['tubing_length_x0.5'] = tub(0.5, 0.5)
        out['tubing_length_x2'] = tub(2.0, 2.0)
        cid = getattr(spec, 'tubing_centre_id_mm', None)
        if cid:
            out['tubing_centre_id_minus_0.025mm'] = tub((cid / (cid - 0.025)) ** 4, 1.0)
            out['tubing_centre_id_plus_0.025mm'] = tub((cid / (cid + 0.025)) ** 4, 1.0)
        # the sensitivity of the last row's cutoff to the center tubing, d(D_N/D)/d ln R_ext_c
        e = dict(base)
        e['Rco'] = base['Rco'] + 1e-3 * Rxc
        a1, b1, _ = solve(e)
        out['dDlast_dlnRextc'] = float((np.fmax(a1[N - 1], b1[N - 1]) - np.fmax(a0[N - 1], b0[N - 1])) / 1e-3)
    return out


def gates(spec, sc_f, D_star):
    r"""Evaluate the design gates for every design, whatever the mode. Strict mode refuses on them, and evaluate and
    reproduce modes report them:

      M13      gap dominance, max R_gap f / R_c < 0.05 (spec v5 sec. 19; CIF_v2 M13)
      M14      channel Reynolds number, Re <= 10 on the center channel (the flow entering each row, CIF_v2's Qc), and
               the Stokes gates, Re_gap <= 1e-2 and Re_delta <= 1e-2 (CIF_v2 M14; App. B thesis scope)
      plateau  u_min = D/max(w_c) >= 0.035 (glossary \plateauUwc; App. B C8 CQ7)
      E11      thin lamina, delta*/h <= eps_lam = 0.10 (CIF_v2 march(), Tolerances.eps_lam; spec v5 E11)
      M16      quasi-developed taper, max_i |w_c(i) - w_c(i-1)| / Lambda <= eps_taper = 1e-2 (CIF_v2 march(),
               Tolerances.eps_taper; spec v5 M16), with CIF_v2's discrete predecessor w_c(0) = W(0) = w_c0 and its
               uniform pitch Lambda = G + g_p (CIF_v2 Process.Lambda), over every row

    CIF_v2 evaluates both Stokes gates at row 1 at the full inlet flow, with Re_gap = reynolds(f* Q_in, G, h) and
    Re_delta = rho C_shear Q_in/(a b^2) delta*^2/mu. Here both are taken at every row at the design values that set the
    cutoff (f*(i) = g(delta*/w_c(i)) and the lamina delta*, with each gap taken through its own opening, G or a head
    mouth), and CIF_v2's literal row-1 Re_gap is included, so the gate is never weaker than CIF_v2's. Re_gap is the
    Reynolds number of the gap flow in the gap's own section (CIF_v2 `reynolds`), which is (w_c + h)/(G + h) times
    (1.2 to 2.8 on the presets) the flux product f Re_channel that App. B's thesis section quotes, and thus the
    stricter of the two. The Stokes gates are cutoff conditions, which CIF_v2 does not gate where no cutoff is claimed,
    so a row planned away from f* (an unsealed gap 1, which is a declared departure) is reported in Re_*_planned_max
    rather than gated.

    Below the plateau, the cutoff is reported by App. B C7's method (P-C7.1). Kappa comes from the frozen Hill fit at
    u = kappa_g delta*/w_c (the widest center channel gives the lowest kappa), and the band is
    [kappa(u) delta*, kappa_g delta*], which is the band App. B gives for a device whose lamina is delta* at every gap
    (its Xia-scaled reading). The band is plateau-anchored rather than the implicit fixed point (App. B C7 CQ3). Below
    the foot, u < 0.012 (glossary \footUwc), no band is given.

    CIF_v2 applies E11 to the rows whose f it computes from delta* (its S1 path, where the lamina is delta*) rather
    than to a prescribed schedule (S2). Here every cutoff row's lamina is delta*, which is gated, and the largest
    planned lamina, which includes rows planned away from f* (declared departures), is reported in
    lamina_over_depth_planned_max rather than gated.

    Four further quantities are reported and not gated. The first is the largest design f against kappa_g's
    calibration cap 5.76e-4 (glossary \fgapDesignCap) and against spec v5's soft bound f_soft = 5e-4 (CIF_v2 Pairing;
    spec v5 sec. 18 and validator D3 of sec. 19, "WARN"). The second is P8, CIF_v2's physics auditor, which in App. B's
    words "checks the drawn device's inlet against its own, stricter channel limit". From its own definition, the
    inlet Re = rho Q D_h / (mu w_c0 h), with D_h = 2 w_c0 h / (w_c0 + h), is compared with 1.0, and the naive entrance
    length 0.05 Re D_h with Lambda/2, at the auditor's default operating point (5 uL/min) and at the design's flow. P8
    is not gated (docs/LIMITATIONS.md 6), because 0.05 Re D_h is the form that spec v5 (its L_e row), the glossary, and
    CIF_v2's own entrance_length (S16-c) call wrong in this regime, and because the Re limit is the auditor's own
    rather than the design envelope's (App. B gives Re_max = 10 on the channel, gated here as M14). The third is a
    Reynolds bound for any channel of the device. A duct of depth h and width w carrying Q_x has
    Re = 2 rho Q_x / (mu (w + h)), and every channel carries at most the inlet flow and is at least the DRC's tested
    width F_min - 0.1 um wide, so Re < 2 rho Q_in / (mu (h + F_min - 0.1 um)), which is a bound conditional on the
    drawing's DRC check. The fourth, which follows from that bound, is Dean's leading inertial term 0.03058 (K/576)^2,
    with K = 2 Re^2 a/R, at a/R = 1, where a is the half-width in the plane of the bend (any drawable arc has a
    centerline radius R >= a). It is an order-of-magnitude bound for every bend (a circular-pipe result, and the
    per-bend a/R is not computed), and it is given where that Reynolds bound is <= 10, the range the estimate is
    quoted for.

    The fluid is the spec's (fluid.density_kg_m3 and viscosity_Pa_s, with water, rho 1000 kg/m^3 and mu 1e-3 Pa s, by
    default; ADJ-C5), and the Reynolds numbers scale with rho/mu. The schedule does not depend on the fluid, because
    its resistance ratios are independent of mu for a Newtonian fluid."""
    p = spec
    Q = p.Q_in_uL_min * 1e-9 / 60.0
    mu, h = float(getattr(p, 'mu_f', S.MU)), p.h
    RHO_F = float(getattr(p, 'rho_f', 1000.0))                              # local, the spec's fluid (ADJ-C5)
    wc = np.asarray(sc_f.wc, float)
    Pin = np.concatenate([[1.0], np.asarray(sc_f.P[:-1], float)])          # the center flow entering row i (CIF_v2 Qc)
    Re = RHO_F * Q * Pin * S.Dh(wc, h) * S.UM / (mu * h * S.UM * wc * S.UM)
    wg = np.asarray(sc_f.wL, float)
    fst = np.array([S.f_of(spec.delta_star, w) for w in wc])                # the design fraction, which sets the cutoff

    def re_duct(qg, w):
        return RHO_F * qg * S.Dh(w, h) * S.UM / (mu * h * S.UM * w * S.UM)
    Re_gap = re_duct(fst * Pin * Q, wg)                                      # one gap's flow (per side)
    Re_gap_v2 = float(re_duct(S.f_of(spec.delta_star, p.wc0) * Q, p.G))     # CIF_v2's literal row-1 term
    Re_gap_pl = re_duct(np.abs(np.asarray(sc_f.Qg, float)) * Q, wg)
    a = np.maximum(wc, h) * S.UM
    b = np.minimum(wc, h) * S.UM
    cs = {}
    C = []
    for x in a / b:
        k = round(float(x), 9)
        if k not in cs:
            cs[k] = S.C.C_shear(float(x))
        C.append(cs[k])
    gam = np.array(C) * Q * Pin / (a * b ** 2)                               # the maximum wall shear (CIF_v2 C_shear)
    # The record reports the peak wall shear along the center channel via C(alpha), as spec v5 S:606 asks ("REPORT ...
    # max gamma_dot via C(alpha) [NEVER 6]"), and, for the parallel-plate taper (reproduce mode only), how far
    # 6 Q/(a b^2) understates it (foundation check 2).
    C_ = np.array(C)
    shear_peak = dict(formula='C(alpha) Q_c / (a b^2), the centre flow entering each row, on the planned centre widths at '
                              'the operating flow (spec v5 S16-b, S:449)',
                      gamma_max_per_s=float(np.max(gam)), gamma_min_per_s=float(np.min(gam)),
                      row1_per_s=float(gam[0]), last_row_per_s=float(gam[-1]),
                      taper=getattr(p, 'taper', 'const'),
                      taper_shear_formula=(getattr(p, 'shear_formula', 'exact') if getattr(p, 'taper', 'const') == 'shear' else None),
                      note_rows='every row\'s peak falls with its centre flow once the centre is at its floor (or constant)')
    if getattr(p, 'taper', 'const') == 'shear':
        # The taper is the linear rows above the planner's floor; the head rows sit at w_c0, outside it
        # (codereview3 B1).
        floor_ = max(float(getattr(p, 'wc_min', 0.0) or 0.0), 0.3 * float(p.wc0))
        tp = wc > floor_ + 1e-6
        kind_ = getattr(sc_f, 'kind', None)                                 # absent on reassess's and tests' namespaces
        if kind_ is not None:
            tp &= np.asarray(kind_) == 'L'
        if np.any(tp):
            shear_peak['over_the_taper'] = dict(rows=int(np.sum(tp)), max_per_s=float(np.max(gam[tp])),
                                                min_per_s=float(np.min(gam[tp])),
                                                max_over_min=float(np.max(gam[tp]) / np.min(gam[tp])))
    if getattr(p, 'taper', 'const') == 'shear' and getattr(p, 'shear_formula', 'exact') == 'parallel_plates':
        shear_peak['true_over_parallel_plates_minus_1'] = [float(np.min(C_ / 6.0 - 1.0)), float(np.max(C_ / 6.0 - 1.0))]
        shear_peak['note'] = ('the parallel-plate taper holds 6 Q/(a b^2) constant (a reproduce-mode departure from spec '
                              'v5 S16-b); the true peak shear is C(alpha)/6 times it and drifts along the taper')
    Re_d = RHO_F * gam * (spec.delta_star * S.UM) ** 2 / mu
    Re_d_pl = RHO_F * gam * (np.asarray(sc_f.delta, float) * S.UM) ** 2 / mu
    gap_max = max(float(np.max(Re_gap)), Re_gap_v2)
    gapdom = float(np.max(sc_f.Rg * np.abs(sc_f.f) / sc_f.Rc))
    Le = S.C.entrance_length(float(Re[0]), S.Dh(p.wc0, p.h) * S.UM) / S.UM
    u = D_star / wc
    u_min = float(np.min(u))
    # The smallest planned row cutoff over the center width, over every row that carries flow with the declared rows
    # included (a depicted head or a two-cutoff reading can plan rows below the device's D), is reported rather than
    # gated, because the gate is the design's D. Rows that extract less than 5% of the constant-cutoff fraction at their
    # width carry no flow to speak of (Xia's depicted rectilinear rows have f/f* ~ 2e-7) and are left out, and the share
    # of the extracted flow at rows below the plateau is reported beside it (code review R6).
    f_pl = np.asarray(sc_f.f, float)
    q = 2 * np.clip(f_pl, 0, None) * Pin
    fs_row = np.array([S.f_of(spec.delta_star, float(w_)) for w_ in np.asarray(wc, float)])
    carry = (f_pl > 0.0) & (f_pl < 0.5) & (f_pl >= 0.05 * fs_row)
    u_all = np.array([S.KAPPA_G * S.x_of(float(fe)) if c_ else np.nan for fe, c_ in zip(f_pl, carry)])
    u_min_all = float(np.nanmin(u_all)) if np.any(carry) else None
    share_below_all = float(np.sum(q[carry & (np.nan_to_num(u_all, nan=1.0) < S.U_LO)]) / np.sum(q)) if np.sum(q) > 0 else 0.0
    knee_rows = u < S.U_LO
    below = u_min < S.U_LO
    foot = u_min >= S.U_FF
    ci = [S.KAPPA_CI[0] * spec.delta_star, S.KAPPA_CI[1] * spec.delta_star]
    if below and foot:
        k_pass = float(hill_kappa(S.KAPPA_G * spec.delta_star / float(np.max(wc))))
        c7 = [k_pass * spec.delta_star, S.KAPPA_G * spec.delta_star]
        band = dict(band_um=c7, basis='App. B C7 method below the plateau (Hill kappa %.3f at the widest centre channel)'
                    % k_pass)
    elif below:
        k_pass, c7 = None, None
        band = dict(band_um=None, basis='below the foot of the kappa curve (u < %.3f): no band' % S.U_FF)
    else:
        k_pass, c7 = None, None
        band = dict(band_um=ci, basis='on the plateau: the kappa_g interval x delta*')
    band.update(retention_side=ci, passage_edge_below_knee=(c7[0] if c7 else None), c7_band_constant_lamina=c7,
                kappa_at_passage_edge=k_pass,
                method=('on the plateau (u_min >= 0.035): the kappa_g interval x delta*.  Below it the record uses App. B '
                        'P-C7.1 (the frozen Hill fit; glossary \\hillKmax, \\hillUhalf, \\hillN): band [kappa(u) delta*, '
                        'kappa_g delta*], plateau-anchored (App. B C7 CQ3)' if not below else
                        'below the plateau (u_min < 0.035): App. B P-C7.1, the frozen Hill fit kappa(u) (glossary \\hillKmax, '
                        '\\hillUhalf, \\hillN) at u = kappa_g delta*/w_c on the widest centre channel; band [kappa(u) delta*, '
                        'kappa_g delta*], App. B C7\'s band for a lamina held at delta* (its Xia-scaled reading; for Dinh V1 '
                        '[\\dinhCollapseLo, \\dinhBandHi] = [1.80, 2.17] um, inside the wider [1.25, 2.17] it records because '
                        'Dinh\'s schedule is not printed).  Plateau-anchored, not the implicit fixed point d = kappa(d/w_c) '
                        'delta, which can lie far lower (App. B C7 CQ3).  retention_side: the kappa_g interval x delta*, '
                        'the band on the plateau'))
    lam = float(spec.delta_star / h)                                          # E11 on the cutoff rows (lamina delta*)
    dl_pl = np.asarray(sc_f.delta, float)
    lam_pl = float(np.nanmax(dl_pl)) / h if np.any(np.isfinite(dl_pl)) else float('nan')
    tap = np.abs(np.diff(np.concatenate([[float(p.wc0)], wc]))) / (p.G + p.g_p)   # M16, CIF_v2's discrete form
    wc0m, hm = p.wc0 * S.UM, h * S.UM                                        # P8, from its own definition

    def p8(q_ul_min):
        q_ = q_ul_min * 1e-9 / 60.0
        dh0 = 2.0 * wc0m * hm / (wc0m + hm)
        re_ = RHO_F * (q_ / (wc0m * hm)) * dh0 / mu
        le_ = 0.05 * re_ * dh0 / S.UM
        return dict(flow_uL_min=float(q_ul_min), inlet_Re=float(re_), inlet_Re_within_limit=bool(re_ <= P8_RE_MAX),
                    naive_Le_um=float(le_), naive_Le_over_half_pitch=float(le_ / (0.5 * p.Lam_L)),
                    passes=bool(re_ <= P8_RE_MAX and le_ <= 0.5 * p.Lam_L))
    w_floor = max(float(getattr(p, 'F_min', 0.0)) - 0.1, 0.0) * S.UM         # the DRC's tested width
    Re_any = 2.0 * RHO_F * Q / (mu * (hm + w_floor))                          # any channel of the device
    dean = (DEAN_C2 * (2.0 * Re_any ** 2 * 1.0 / 576.0) ** 2
            if Re_any <= RE_MAX * (1.0 + 1e-9) else None)                     # quoted for Re <= 10 only
    return dict(M13_gap_dominance_max=gapdom, M13_ok=gapdom < M13_MAX, Re_max=float(np.max(Re)), Re_ok=float(np.max(Re)) <= RE_MAX,
                Re_gap_max=gap_max, Re_gap_row=(int(np.argmax(Re_gap)) + 1 if np.max(Re_gap) >= Re_gap_v2 else 1),
                Re_gap_row1_cifv2=Re_gap_v2,
                Re_delta_max=float(np.max(Re_d)), Re_delta_row=int(np.argmax(Re_d)) + 1, Re_delta_row1=float(Re_d[0]),
                stokes_ok=bool(gap_max <= RE_STOKES and np.max(Re_d) <= RE_STOKES),
                Re_gap_planned_max=float(np.nanmax(Re_gap_pl)), Re_delta_planned_max=float(np.nanmax(Re_d_pl)),
                f_design_max=float(np.max(fst)), f_planned_max=float(np.max(np.asarray(sc_f.f, float))),
                f_calibration_cap=S.FGAP_CAP, f_beyond_calibration=bool(np.max(fst) > S.FGAP_CAP),
                f_soft=F_SOFT, f_above_soft=bool(np.max(fst) > F_SOFT),
                lamina_over_depth=lam, lamina_limit=EPS_LAM, lamina_ok=bool(lam <= EPS_LAM),
                lamina_over_depth_planned_max=lam_pl,
                lamina_all_rows_ok=bool(not (lam_pl > EPS_LAM) and lam <= EPS_LAM),
                taper_rate_max=float(np.max(tap)), taper_rate_row=int(np.argmax(tap)) + 1, taper_rate_limit=EPS_TAPER,
                taper_ok=bool(float(np.max(tap)) <= EPS_TAPER), taper_pitch_um=float(p.G + p.g_p),
                P8=dict(inlet_Re_limit=P8_RE_MAX, at_auditor_flow=p8(P8_Q_UL_MIN), at_design_flow=p8(p.Q_in_uL_min),
                        status='design-stage values (planned widths) of the physics auditor\'s own inlet criteria, '
                               'check P8; the audit evaluates them on the drawn mask (validation.physics_audit, P8) '
                               'and gates the design envelope there (check M14: Re <= 10 on the channel, Re_gap and Re_delta '
                               '<= 1e-2); P8 is reported, not gated (docs/LIMITATIONS.md 6): its 0.05 Re D_h is the '
                               'entrance length spec v5, the glossary and CIF_v2 S16-c call wrong in this regime, and '
                               'its Re limit 1.0 is CIF_v2\'s settable default (App. B prints no value), stricter than '
                               'the design envelope\'s (Re_max = 10, M14)'),
                Re_any_channel_bound=float(Re_any), dean_leading_term_bound=dean,
                dean_note=('Dean\'s leading term 0.03058 (K/576)^2, K = 2 Re^2 a/R, at Re = Re_any_channel_bound and a/R = 1 '
                           '(a: half-width in the plane of the bend): an order-of-magnitude screening estimate for every '
                           'bend, not a validated bound (Dean\'s result is for a circular pipe; no per-bend a/R or 3-D '
                           'solve); Re_any_channel_bound = 2 rho Q/(mu (h + F_min - 0.1 um)), every channel '
                           'being at least the DRC\'s tested width' if dean is not None else
                           'no bound: the channel Reynolds bound 2 rho Q/(mu (h + F_min - 0.1 um)) = %.2f exceeds %g (the '
                           'gated centre-channel Re is Re_max)' % (Re_any, RE_MAX)),
                shear_peak=shear_peak,
                fluid=dict(rho_kg_m3=RHO_F, mu_Pa_s=mu, rheology='newtonian',
                           note=('water (the default); the Reynolds gates scale with rho/mu' if (RHO_F, mu) == (1000.0, 1.0e-3)
                                 else 'the spec\'s fluid section; the Reynolds gates scale with rho/mu') +
                                '; the schedule does not depend on the fluid (resistance ratios, Newtonian); the absolute '
                                'resistances in the record (Pa s/m^3) are at water\'s viscosity'),
                Le_over_Lambda_L=Le / p.Lam_L,
                plateau_u_min=u_min, plateau_ok=bool(u_min >= S.U_LO), plateau_u_min_all_planned_rows=u_min_all,
                plateau_F_share_below_all_planned_rows=share_below_all,
                plateau_all_planned_rows_note=('rc7: planned row cutoffs over the centre width, declared rows included, over '
                                               'the rows extracting at least 5 % of the constant-cutoff fraction at their '
                                               'width; the share is of the flow all rows extract'),
                plateau_rows_in_knee=int(np.sum(knee_rows)),
                plateau_F_share_in_knee=float(np.sum(q[knee_rows]) / np.sum(q)),
                cutoff_band_um=band,
                kappa_band_um=ci)


class LayoutCollision(Exception):
    """The exception for a drawn outline whose parts cross or come closer than the process floor, or in which a port
    crowds another (the clearance guarantee). It is raised before any photomask is written, pipeline.run records it as
    a refusal (exit 5), and it carries the collision checks (`collisions`), the outline, and the plan, for the
    refusal's numbers and its figure."""

    def __init__(self, coll, fluid=None, sc=None, info=None):
        Exception.__init__(self, '; '.join(coll['why']))
        self.coll, self.fluid, self.sc, self.info = coll, fluid, sc, info


def collisions(fluid, sc, info, spec):
    """Check whether a drawn outline's parts collide (the clearance guarantee). Every part (the inlet with the head,
    each leg, each U-turn, and the outlet) must stay at least the process floor from every part it is not joined to
    (layout.part_clearance), the outlet ports must be separated from each other and from the inlet, and the inlet port
    must be clear of every outlet port (layout.port_spacing). These are the layout's own conditions, whose remedy is a
    layout change, so an outline that violates them is refused before the photomask is written. The outline's
    integrity (one polygon with the planned holes, continuous side channels, every feature found, and the DRC) stays a
    check of the drawn mask, so a failure there is a FAIL, and its photomask is written as <name>_FAILED.dxf
    (pipeline.name_mask). The topology counts are computed here once for the record. The result is dict(ok, why,
    single, holes, exp_holes, separation, part_clearance, port_spacing, pairs_below_floor)."""
    single, holes, exp_holes = topology_counts(fluid, sc, info)
    sep = port_separation(fluid, sc, info)
    pc = L.part_clearance(info['parts'], spec.F_min)
    ps = L.port_spacing(info, spec)
    below = []                                        # every unrelated pair below the floor (part_clearance lists 3)
    if not pc['ok']:
        prt = info['parts']
        for i in range(len(prt)):
            for j in range(i + 2, len(prt)):
                (na, a), (nb, b) = prt[i], prt[j]
                if a.distance(b) < spec.F_min - 0.1:
                    below.append((na, nb))
    why = []
    if not sep['separated']:
        why.append('the outlet ports are not separated from each other and from the inlet')
    if not pc['ok']:
        c0 = pc['closest'][0]
        why.append("'%s' and '%s' are %.1f um apart (%s; they must keep the %.1f um process floor)"
                   % (c0['parts'][0], c0['parts'][1], c0['gap_um'], 'they overlap' if c0['gap_um'] <= 0.0 else
                      'closer than the floor', pc['floor_um']))
    if not ps['ok']:
        why.append('outlet port %s is %.0f um from the inlet port (at least %.0f)'
                   % (ps['port'], ps['inlet_to_nearest_outlet_port_um'], ps['required_um']))
    return dict(ok=not why, why=why, single=single, holes=holes, exp_holes=exp_holes, separation=sep,
                part_clearance=pc, port_spacing=ps, pairs_below_floor=below)


OUTLET_EXT_STEP = 500.0      # the smallest step of the outlet-stretch extension, in um (plan_and_draw)
OUTLET_EXT_MAX = 12          # the extensions tried before the outline is left to the collision gate


def _outlet_meets_inlet(coll):
    """Return True for the collisions that an outlet-stretch extension can cure, which are an outlet port too close to
    the inlet port, outlet ports that are not separated, or the outlet part too close to another part, provided that no
    pair without the outlet is below the floor, because then no extension can clear the gate (code review 4, F8)."""
    if any('outlet' not in pr for pr in coll.get('pairs_below_floor') or []):
        return False
    pc = coll['part_clearance']
    return (not coll['port_spacing']['ok'] or not coll['separation']['separated'] or
            (not pc['ok'] and 'outlet' in pc['closest'][0]['parts']))


def plan_and_draw(spec, kind, n_legs, R_turn, failsafe, say, name, t0):
    """Plan the design on its layout (with the floor and envelope re-plans) and draw its outline in memory. This is
    the part of `build` before anything is written, and it is shared with the plan-stage check (pipeline.layout_check),
    so that the outline checked there is the outline that a build would write. The function mutates spec as the
    re-plans ask (cap_frac, ws_in, and L_out_extra), so the caller passes a copy when the spec must stay as it is.

    When the outline's outlet meets the inlet end (a folded layout whose last leg ends beside its inlet with nothing
    there to clear, for example a paperclip of two nested legs), the outlet stretch is lengthened and the design is
    re-planned, so that its closure takes the longer stretch (the clearance guarantee). Each extension is at least
    OUTLET_EXT_STEP and is sized by the ports' shortfall, up to OUTLET_EXT_MAX times, and only while no pair without
    the outlet collides. The extensions stop when a try does not help or when the longer stretch no longer closes the
    outlet, and the best outline drawn then goes to the gate, with its tries in coll['extensions'] (code review 4, F5,
    F6, F8). A design whose outline is clear at once is drawn without an extension. A paperclip's ring outlet is not
    extended, because the ring sets its stretch, whereas a linear device's stretch is lengthened too (F7). The function
    returns (fluid, sc_f, info, exact, plan_kinds, final_kinds, coll)."""
    ring = getattr(spec, 'outlet', 'separate') == 'combined' and getattr(spec, 'outlet_style', 'teardrop') == 'ring'
    ring_sets_stretch = ring and kind == 'paperclip' and n_legs > 1      # only the paperclip's ring sets the stretch
    best, tries = None, []
    base_ = (spec.cap_frac, spec.ws_in)                  # each try re-plans from these (third-party review 1 of rc7)
    for ext_try in range(OUTLET_EXT_MAX + 1):
        if ext_try:
            # The floor and envelope re-plans (cap_frac, ws_in) start again from the spec as given, so that a try's
            # outline depends on its stretch alone rather than on the tries before it. RC7 carried them over, so that
            # an extended build was a function of its tries.
            spec.cap_frac, spec.ws_in = base_
        try:
            got = _plan_and_draw_once(spec, kind, n_legs, R_turn, failsafe, say, name, t0)
        except (L.OutletInfeasible, L.OutletOutsideTables, RuntimeError) as e:
            if ext_try == 0 or (isinstance(e, RuntimeError) and 'floor not met' not in str(e)):
                raise                                     # a later try that cannot meet the floor ends like F6
            # When the longer stretch no longer closes the outlet, the best outline drawn goes to the gate
            # (code review 4, F6).
            say('[%s] %s with the stretch %.0f um longer (%s): the gate gets the best outline drawn'
                % (name, 'the floor is not met' if isinstance(e, RuntimeError) else 'the outlet does not close',
                   spec.L_out_extra, e))
            tries.append(dict(extension_um=float(spec.L_out_extra), outlet_refused=str(e)[:200]))
            break
        fluid, sc_f, info, exact, plan_kinds, final_kinds = got
        coll = collisions(fluid, sc_f, info, spec)
        # How far the outlet is from clearing (lower is better) is the sum of the ports' shortfall, the closest pair's
        # shortfall and their overlap's area (in 10^3 um^2), and the ports' separation.
        pc_ = coll['part_clearance']
        ov_ = 0.0
        if not pc_['ok']:
            parts_ = dict(info['parts'])
            a_, b_ = (parts_.get(n_) for n_ in pc_['closest'][0]['parts'])
            ov_ = float(a_.intersection(b_).area) / 1e3 if a_ is not None and b_ is not None else 0.0
        short = (max(0.0, coll['port_spacing']['required_um'] - coll['port_spacing']['inlet_to_nearest_outlet_port_um']) +
                 (0.0 if pc_['ok'] else max(0.0, pc_['floor_um'] - pc_['closest'][0]['gap_um']) + ov_) +
                 (0.0 if coll['separation']['separated'] else 1.0))
        ext_now = float(getattr(spec, 'L_out_extra', 0.0) or 0.0)
        tries.append(dict(extension_um=ext_now, shortfall=round(short, 3), why=coll['why']))
        # When the last extension did not help, the tries stop and the gate refuses.
        stuck = best is not None and short >= best[0] - 1e-6
        if best is None or short < best[0]:
            best = (short, (fluid, sc_f, info, exact, plan_kinds, final_kinds, coll), ext_now)
            best_state = (spec.cap_frac, spec.ws_in)      # the re-plans of the outline kept
        if coll['ok'] or ring_sets_stretch or not _outlet_meets_inlet(coll) or stuck or ext_try == OUTLET_EXT_MAX:
            break
        ps = coll['port_spacing']
        step = max(OUTLET_EXT_STEP, (ps['required_um'] - ps['inlet_to_nearest_outlet_port_um'] + 200.0)
                   if not ps['ok'] else 0.0)
        spec.L_out_extra = ext_now + step
        say('[%s] the outlet meets the inlet end (%s): outlet stretch lengthened by %.0f um (%.0f in all), re-plan' % (
            name, '; '.join(coll['why']), step, spec.L_out_extra))
    short_b, out, ext_b = best
    if len(tries) > 1:                                    # the extensions are part of the record (code review 4, F5)
        out[6]['extensions'] = dict(tried=tries, kept_extension_um=ext_b)
        if ext_b != float(getattr(spec, 'L_out_extra', 0.0) or 0.0):
            say('[%s] kept the outline with the stretch %.0f um longer (the best of the tries)' % (name, ext_b))
        spec.L_out_extra = ext_b                          # the spec describes the outline returned
        spec.cap_frac, spec.ws_in = best_state            # with the re-plans that drew it
    return out


def _plan_and_draw_once(spec, kind, n_legs, R_turn, failsafe, say, name, t0):
    """Plan the design on its layout and draw its outline once, planning at most 8 times, re-planning while the
    narrowest side channel or loop is below the floor or a taper unit outgrows its capacity, and raising RuntimeError
    when the floor is still not met. It returns (fluid, sc_f, info, exact, plan_kinds, final_kinds)."""
    for attempt in range(8):
        sc = S.plan(spec)
        say('[%s] planned: A=%d (sealed would need %s, lateral cap %s), f1/f*=%.3f, N=%d, NM=%d, w_s(in) target %.2f (%.0fs)' % (
            name, sc.A, sc.info.get('A_need_sealed'), sc.info.get('A_cap'), sc.f_seed / sc.info['fstar'], sc.N, sc.NM,
            spec.ws_in, time.time() - t0))
        sc_l, legs = L.plan_layout(sc, kind, n_legs, R_turn)
        loop_w = {}
        for k in sc_l.info['turns']:
            wi, _ = L.loop_profile(sc_l, k, True)
            wo, _ = L.loop_profile(sc_l, k, False)
            loop_w[k] = 0.5 * (wi + wo)
        sc_l.info['loop_w'] = loop_w
        units = Gm.plan_head(sc_l)
        if units:
            Gm.solve_transition(sc_l, units[-1])
        exact = dict(F=sc_l.F, N=sc_l.N)
        if failsafe:
            sc_f, units = RB.failsafe(sc_l.spec, sc_l, units)
            sc_f.info['Rc_lead_len'] = sc_l.info['Rc_lead_len']
            sc_f.info['loop_w'] = loop_w
            say('[%s] fail-safe margin on %d rows; F %.5f -> %.5f, N %d -> %d' % (
                name, len(sc_f.info['failsafe_mult']), sc_l.F, sc_f.F, sc_l.N, sc_f.N))
        else:
            sc_f = sc_l
        paper_, wR_, kL_ = L.transition_params(sc_f)
        ws_lin = np.array([wR_ if (paper_ and sc_f.kind[k] == 'R') else sc_f.ws[k] for k in range(sc_f.NM, sc_f.N - 1)])
        ws_min = float(np.nanmin(ws_lin))
        lw_min = min([min(L.loop_profile(sc_f, k, True)[0], L.loop_profile(sc_f, k, False)[0]) for k in sc_f.info['turns']]
                     or [1e9])
        short = spec.F_min - min(ws_min, lw_min)
        if spec.head == 'none':
            short = -1.0  # classic sealed start, whose narrowest width is the design (gated by assess's floor check)
        over = [k for k in range(sc_f.NM - 1) if sc_f.kind[k] in ('T', 'N')
                and sc_f.L_req[k] > sc_f.spec.L_cap_w(sc_f.wM[k]) / sc_f.spec.cap_frac]
        if over:
            ratio = max(sc_f.L_req[k] / (sc_f.spec.L_cap_w(sc_f.wM[k]) / sc_f.spec.cap_frac) for k in over)
            new_frac = spec.cap_frac / ratio * 0.995
            say('[%s] envelope: %d taper units outgrow their capacity (worst x%.4f): cap_frac %.4f -> %.4f, re-plan' % (
                name, len(over), ratio, spec.cap_frac, new_frac))
            spec.cap_frac = new_frac
            continue
        if short <= 0.005:
            break
        say('[%s] floor: narrowest straight side %.3f um, loop %.3f um < F_min %.1f: w_s(in) target %.2f -> %.2f, re-plan' % (
            name, ws_min, lw_min, spec.F_min, spec.ws_in, spec.ws_in + short + 0.1))
        spec.ws_in += short + 0.1
    else:
        raise RuntimeError('floor not met after re-planning')
    plan_kinds = {k: int(np.sum(sc_l.kind == k)) for k in 'ANTRL'}
    final_kinds = {k: int(np.sum(sc_f.kind == k)) for k in 'ANTRL'}
    legs = L.make_legs(sc_f)
    fluid, sc_f, info = L.assemble(sc_f, legs)
    say('[%s] drawn: bounds %s (%.0fs)' % (name, [round(v, 1) for v in fluid.bounds], time.time() - t0))
    return fluid, sc_f, info, exact, plan_kinds, final_kinds


def build(spec, kind='linear', n_legs=1, outdir='.', name='device', failsafe=False, do_audit=True, verbose=True,
          R_turn=None, cut_tol=0.01, mode='strict', seal='require', ar_max=AR_MAX_DEFAULT):
    """Plan, draw, and validate one device from an engine Spec, and return (rep, fluid, sc_f, info). The photomask is
    written to outdir as <name>_UNASSESSED.dxf, with the record (<name>.json) and the arrays and outline (<name>.pkl)
    beside it, and an outline whose parts collide raises LayoutCollision before anything is written."""
    t0 = time.time()
    os.makedirs(outdir, exist_ok=True)
    log = []

    def say(msg):
        if verbose:
            print(msg, flush=True)
        log.append(msg)
    fluid, sc_f, info, exact, plan_kinds, final_kinds, coll = plan_and_draw(spec, kind, n_legs, R_turn, failsafe, say,
                                                                            name, t0)
    # An outline whose parts collide is never written, because it is refused before the photomask exists (the
    # clearance guarantee).
    if not coll['ok']:
        say('[%s] parts collide, nothing written: %s' % (name, '; '.join(coll['why'])))
        raise LayoutCollision(coll, fluid, sc_f, info)
    # The mask is written as <name>_UNASSESSED.dxf, and pipeline.name_mask gives it its name once the verdict exists
    # (<name>.dxf or <name>_FAILED.dxf), so a run that stops before its verdict leaves no deliverable name
    # (code review 4, F4).
    dxf = os.path.join(outdir, name + '_UNASSESSED.dxf')
    io_dxf.write_dxf(fluid, dxf)
    V = {}
    # ---- V3 topology
    single, holes, exp_holes = coll['single'], coll['holes'], coll['exp_holes']
    V['topology'] = dict(single_component=single, valid=bool(fluid.is_valid), holes=len(holes),
                         holes_expected=int(exp_holes), holes_ok=bool(len(holes) == exp_holes),
                         rectilinear_rows=int(np.sum(sc_f.kind == 'R')), arms_per_side=int(sc_f.A),
                         meander_units_per_side=int(sc_f.NM), row_kinds_planned=plan_kinds,
                         row_kinds_drawn=final_kinds,
                         row_kinds_survive=bool(all(plan_kinds[k] == final_kinds[k] for k in 'ANTR')))
    # ---- V3c head envelope
    hz = max(L.head_halfwidth(sc_f), spec.W / 2.0) + 1000.0
    if sc_f.NM > 0:
        head = fluid.intersection(box(info['x_top'] - 10.0, -hz, sc_f.x[sc_f.NM], hz))
        axis_ = LineString([(info['x_top'], 0.0), (sc_f.x[sc_f.NM], 0.0)])
        parts_ = [g for g in (head.geoms if hasattr(head, 'geoms') else [head]) if g.intersects(axis_)]
        head = unary_union(parts_) if parts_ else head
        hb = head.bounds
    else:
        hb = (0.0, -spec.wc0 / 2.0, 0.0, spec.wc0 / 2.0)
    V['envelope'] = dict(W_um=spec.W, head_halfwidth_um=float(max(abs(hb[1]), abs(hb[3]))),
                         ok=bool(max(abs(hb[1]), abs(hb[3])) <= spec.W / 2.0 + 0.5))
    # ---- V2 curvature
    V['rounding'] = rounding_check(fluid, spec)
    rmin = V['rounding']['min_local_radius_um']
    say('[%s] topology holes %d/%d; envelope %s; r_min %.2f um (%.0fs)' % (
        name, len(holes), exp_holes, V['envelope']['ok'], rmin, time.time() - t0))
    V['outlet_separation'] = coll['separation']
    V['part_clearance'] = coll['part_clearance']
    V['port_spacing'] = coll['port_spacing']
    V['components_drawn'] = components_drawn(fluid, sc_f, info)
    say('[%s] part clearance: %s (closest %s); drawn components %s' % (
        name, 'ok' if V['part_clearance']['ok'] else 'PARTS TOO CLOSE', V['part_clearance']['closest'][:1],
        {k: v for k, v in V['components_drawn'].items() if k != 'head_unit_widths_um'}))
    V['side_continuity'] = side_continuity(fluid, sc_f, info)
    say('[%s] side-channel continuity: %s; narrowest %.3f um (floor %.1f) (%.0fs)' % (
        name, 'every leg and side one strip' if V['side_continuity']['ok'] else 'BROKEN %s' % [
            (r['leg'], r['side'], r['pieces']) for r in V['side_continuity']['legs'] if not r['connected']],
        V['side_continuity']['narrowest_um'], spec.F_min, time.time() - t0))
    # ---- V4/V5 measured network
    m = M2.measure_device(fluid, sc_f, info)
    # the features that the measurement probes did not find where the plan puts them (measure.found; E5)
    V['features_found'] = dict(ok=not m['not_found'], not_found=m['not_found'][:20], n_not_found=len(m['not_found']),
                               by_kind={kd: sum(1 for q in m['not_found'] if q['kind'] == kd)
                                        for kd in sorted({q['kind'] for q in m['not_found']})})
    if m['not_found']:
        f0 = m['not_found'][0]
        say('[%s] DRAWING DEFECT: %d planned feature(s) not found on the drawing (first: row %d, %s, side %s)'
            % (name, len(m['not_found']), f0['row'], f0['kind'], f0['side']))
    nw = M2.network(sc_f, m)
    rA, rB = nw['ratioA'], nw['ratioB']
    dev = np.maximum(np.abs(rA - 1), np.abs(rB - 1))
    units_rel = [max(abs(rec['A'].get('L_eff', rec['A']['L']) / rec['L_req'] - 1),
                     abs(rec['B'].get('L_eff', rec['B']['L']) / rec['L_req'] - 1))
                 for rec in m['units'] if rec['L_req']]               # resistance lengths (drawn plus bend correction)
    paper_ = L.transition_params(sc_f)[0]
    kL = int(next(k for k in range(sc_f.NM, sc_f.N) if sc_f.kind[k] == 'L'))
    lin = np.array([(a, b) for k, s_, a, b in m['lin'] if not (paper_ and sc_f.NM <= k <= kL)]) \
        if m['lin'] else np.zeros((0, 2))
    gaps = np.array([(a, b) for k, s_, a, b in m['gaps']])
    V['fidelity'] = dict(unit_length_rel_max=float(max(units_rel)) if units_rel else 0.0,
                         side_eqwidth_dev_max_um=float(np.max(np.abs(lin[:, 0] - lin[:, 1]))) if len(lin) else 0.0,
                         gap_width_dev_max_um=float(np.max(np.abs(gaps[:, 0] - gaps[:, 1]))),
                         loop_R_rel_max=float(max([max(abs(lp['A']['R'] / lp['target'] - 1), abs(lp['B']['R'] / lp['target'] - 1))
                                                    for lp in m['loops']] or [0.0])),
                         outlet_R_rel_max=float(max(abs(m['outlet']['Rc_out_rel']), abs(m['outlet']['RsA_out_rel']),
                                                    abs(m['outlet']['RsB_out_rel']))))
    D_star = S.KAPPA_G * spec.delta_star
    outlet_and_loops_measured(V, m, spec)
    V['network'] = dict(F_design=sc_f.F, F_target=spec.F_target, F_as_drawn=nw['F'], F_rel_to_design=nw['F'] / sc_f.F - 1,
                        F_rel_to_target=nw['F'] / spec.F_target - 1,
                        row_f_ratio_max_dev=float(np.nanmax(dev)), row_f_ratio_max_dev_row=int(np.nanargmax(dev)) + 1,
                        cutoff_design_um=D_star)
    inten, nets, decl = intended_and_declared(sc_f, info, m, spec, D_star, cut_tol, mode, seal)
    f_star_rows = np.array([S.f_of(spec.delta_star, w) for w in sc_f.wc])
    V['drawing_fidelity'] = drawing_fidelity(nw, inten, f_star_rows, sc_f.N, D_star)
    paper_, wR_, kL_ = L.transition_params(sc_f)
    V['realised_cutoff'] = realised_cutoff(sc_f, nw, D_star, cut_tol, decl)
    V['intended_cutoff'] = realised_cutoff(sc_f, inten, D_star, cut_tol, decl)
    V['network']['F_intended'] = float(inten['F'])
    try:
        V['model_scenarios'] = model_scenarios(sc_f, info, m, nw, D_star, cut_tol)
    except Exception as e:                                # reported, never silently dropped
        V['model_scenarios'] = dict(error='%s: %s' % (type(e).__name__, e))
    say('[%s] fidelity %s' % (name, V['fidelity']))
    say('[%s] network F %.6f (design %.6f, target %.6f); realised cutoff %s (%.0fs)' % (
        name, nw['F'], sc_f.F, spec.F_target, {k: round(v, 5) if isinstance(v, float) else v
                                               for k, v in V['realised_cutoff'].items()}, time.time() - t0))
    V['gates'] = gates(spec, sc_f, D_star)
    # ---- V7 (`robustness_design_network`), perturbations of the design network
    res = RB.scenarios(sc_f, info['units'])
    f_star_cut = np.array([S.f_of(spec.delta_star, w) for w in sc_f.wc])
    rob = {}
    for nm, r in res.items():
        fa = r * sc_f.f
        dc = np.array([S.x_of(v) * w if 0 < v < 0.5 else np.nan for v, w in zip(fa, sc_f.wc)]) / spec.delta_star
        inner = dc[1:-1]
        rob[nm] = dict(cutoff_max_over_Dstar=float(np.nanmax(inner)), cutoff_min_over_Dstar=float(np.nanmin(inner)),
                       rows_without_cutoff=int(np.sum(~np.isfinite(dc))), rows_outside_1pct=int(np.sum(~(np.abs(dc - 1) <= 0.01))),
                       cutoff_row1_over_Dstar=float(dc[0]) if np.isfinite(dc[0]) else None,
                       cutoff_lastrow_over_Dstar=float(dc[-1]) if np.isfinite(dc[-1]) else None)
    V['robustness_design_network'] = rob     # perturbations of the design network rather than of the drawing
    # ---- V8 parameters read back from the drawing
    V['drawn_parameters'] = drawn_parameters(m, sc_f, fluid, info, kL)
    # Spec v5 G9 asks for a "Length budget: sum over leg k of Lambda <= l_leg for all k", so the record gives each
    # leg's path length along its rows. Leg 1 is measured from the head's top, so the unfolded head counts; an
    # intermediate leg runs from half a pitch before its first gap to half a pitch after its last; and the last leg
    # runs to the end of the outlet stretch, x_N + L_out. U-turns and outlet leads are not counted, so the length is
    # at least G9's sum of Lambda and thus conservative, and the record carries this definition as
    # length_budget.definition. pipeline.assess sets the length beside the chip's leg length. It is reported and never
    # gated (in App. B, the design code "reports the result without failing").
    V['length_budget'] = dict(leg_path_lengths_mm=[float(((lg['p1'] if lg['p1'] is not None else sc_f.x[-1]) -
                                                          (lg['p0'] if lg['p0'] is not None else info['x_top'])) / 1000)
                                                   for lg in info['legs']],
                              total_mm=None,
                              definition='each leg along its rows: leg 1 from the head\'s top, an intermediate leg from '
                                         'half a pitch before its first gap to half a pitch after its last, the last '
                                         'leg to the end of the outlet stretch; U-turns and outlet leads not counted '
                                         '(at least spec v5 G9\'s sum of Lambda: conservative)')
    V['length_budget']['total_mm'] = float(sum(V['length_budget']['leg_path_lengths_mm']))
    design = dict(name=name, kernel_pair=KERNEL_PAIR, spec={k: (v if isinstance(v, (int, float, str)) or v is None else str(v))
                                   for k, v in spec.__dict__.items()},
                  layout=kind, n_legs=n_legs, A=int(sc_f.A), f1_over_fstar=float(sc_f.f_seed / sc_f.info['fstar']),
                  gap1_cutoff_um=float(S.KAPPA_G * S.x_of(sc_f.f_seed) * spec.wc0), D_star_um=D_star,
                  N=int(sc_f.N), NM=int(sc_f.NM), F=float(sc_f.F), exact_before_margin=exact,
                  failsafe_rows=len(sc_f.info.get('failsafe_mult', {})),
                  K_narrow=int(sc_f.info.get('K_narrow', 0) or 0),
                  turns=[int(k) + 1 for k in sc_f.info.get('turns', [])],
                  leg_Y=sc_f.info.get('leg_Y'), leg_d=sc_f.info.get('leg_d'),
                  ws_min_um=float(np.nanmin(sc_f.ws[sc_f.NM:sc_f.N - 1])),
                  wall_rule=L.WALL_RULE, x_top=info['x_top'],
                  # the element's hydraulic resistance from gap 1 to the outlets (all outlets at the same pressure),
                  # per unit inlet flow, in SI units (Pa s m^-3), which is used to split the flow in a parallel array
                  R_element_SI=float(-(sc_f.Pc[sc_f.N - 1] - sc_f.Rc_out * sc_f.P[sc_f.N - 1])),
                  R_element_SI_note=('rc7 (code review R15): at water\'s viscosity (schedule.MU = 1e-3 Pa s), as every '
                                     'absolute resistance in the record; it scales with the fluid\'s viscosity, the '
                                     'schedule (resistance ratios) does not'),
                  # The outlet as designed is recorded with the stretch's end (schedule.outlet_exp_c, which gives the
                  # kind, the width ratio, c with its deep and shallow values and the depth weight, the excess as a
                  # length and a resistance, and for a flare its angle, length, and straight run) and with the side
                  # lead's design width and length. L_out_um is the drawn outlet stretch, the layout's (sc_f.spec,
                  # which includes a ring's loop clearance, a serpentine's or paperclip's turn, or the flare's straight
                  # run plus the flare), rather than the spec's base value.
                  outlet=dict(style=('separate' if getattr(spec, 'outlet', 'separate') == 'separate'
                                     else str(getattr(spec, 'outlet_style', 'teardrop'))),
                              w_lead_um=float(info['outlet'].get('w_lead', 0.0) or 0.0),
                              L_lead_um=float(info['outlet'].get('L_lead', 0.0) or 0.0), L_out_um=float(sc_f.spec.L_out),
                              R_stretch_design_SI=info['outlet'].get('R_stretch_design'),
                              stretch_end=info['outlet'].get('expansion_term')))
    if getattr(spec, 'L_out_extra', 0.0):                 # only when it was needed (the clearance guarantee)
        design['outlet_stretch_extended_um'] = float(spec.L_out_extra)
    if do_audit:
        V['drc'] = run_drc(dxf, outdir, name, spec)
        say('[%s] DRC %s' % (name, {k: V['drc'].get(k) for k in ('n_components', 'invalid_rings', 'n_lost_open',
                                                                  'n_lost_solid', 'd_um')}))
    # The physics auditor (P0 to P10) runs on the photomask file read back. P0 requires the file to equal the drawn
    # polygon vertex for vertex, so the measurement above is the measurement of the file.
    V['physics_audit'] = AU.audit(dxf, sc_f, info, m, V, decl, mode, ar_max, fluid=fluid, nw=nw, intended=inten)
    say('[%s] %s (%.0fs)' % (name, AU.summary_line(V['physics_audit']), time.time() - t0))
    # the design code's mandatory reports and acceptance validators, each by its own formula (reported, as there)
    V['design_code_reports'] = DCR.reports(sc_f, info, m, spec, float(nw['F']))
    rep = dict(design=design, validation=V, log=log, seconds=time.time() - t0)
    arrays = dict(x=sc_f.x, wc=sc_f.wc, ws=sc_f.ws, f=sc_f.f, P=sc_f.P, S=sc_f.S, kind=sc_f.kind, NM=sc_f.NM,
                  N=sc_f.N, ratioA=rA, ratioB=rB, dA=nw['dA'], dB=nw['dB'], fA=nw['fA'], fB=nw['fB'],
                  fA_int=inten['fA'], fB_int=inten['fB'], dA_int=inten['dA'], dB_int=inten['dB'], wc_gap=nw.get('wc_gap'),
                  nets={k_: dict(dA=v_['dA'], dB=v_['dB'], fA=v_['fA'], fB=v_['fB'], F=v_['F']) for k_, v_ in nets.items()},
                  declared=decl, kL=kL_,
                  fstar_cut=f_star_cut, turns=sc_f.info.get('turns', []), legs=info['legs'], gx=info['gx'],
                  loops=[{k: v for k, v in lp.items() if k not in ('fA', 'fB')} for lp in info['loops']],
                  units=[{k: v for k, v in u.items() if k != 'drop_poly'} for u in info['units']],
                  inlet=info['inlet'], outlet={k: v for k, v in info['outlet'].items()}, x_top=info['x_top'],
                  meas_lin=m['lin'], meas_wc=m['wc'], spec=spec)
    with open(os.path.join(outdir, name + '.pkl'), 'wb') as fh:
        pickle.dump(dict(fluid=fluid, arrays=arrays), fh)
    rep['seconds'] = time.time() - t0
    with open(os.path.join(outdir, name + '.json'), 'w') as fh:
        json.dump(rep, fh, indent=1, default=float)
    say('[%s] done in %.0fs' % (name, time.time() - t0))
    return rep, fluid, sc_f, info
