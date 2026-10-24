"""cfd_kit.py: the CFD validation kit, a window of a built device ready for an independent 3-D Stokes solve (future
work, docs/LIMITATIONS.md items 8-10).

The 1-D network's weakest rows are those right after a U-turn and at a meander head's seam. Their gap flows are set by
a small pressure difference that is left after each channel has lost a large pressure through the turn or along the
head, so small errors of the resistance model there (bends, the head's T-junctions) move their cutoffs most (the
record's model band). This tool cuts such a window out of the drawn mask and writes what a Stokes solver needs to test
the network there:

  window.dxf      the fluid polygon of the window (layer CHANNELS), in um. Every open boundary is a straight cut across
                  one channel, perpendicular to its flow, for the center channel and the two side channels. In the array
                  rows a cut lies at a post's middle, and in a head window it crosses a meander unit at the apex of its
                  U-turn, so that the window holds that unit's return line, drop, and T-junction and every later unit
                  whole.
  boundary.csv    every cut, with its channel, its end points (um), the drawn width there, the flow through it as a
                  fraction of the device inlet flow and in uL/min at the spec's flow rate, and its role (in or out)
  gaps.csv        every gap in the window, with its row, its side, a section across the gap (end points, um) over which
                  to integrate the flux, and the network's prediction for the drawn mask, which is f (of the center flow
                  arriving at that row), the gap flow (as a fraction of the device inlet and in uL/min), and the cutoff
                  D_i / D
  window.png      the window with its cuts and gap sections marked
  README.md       the boundary conditions and the comparison to make

The network numbers are those of the measured (drawn) network in the record, so a difference is a model error (bends,
junctions, gap mouths) rather than a drawing error. Two setups answer different questions.
  (a) With all cut flows prescribed except the downstream center cut (p = 0), the solve tests how the window
      distributes its gap flows, since each side's total gap flow over the window is then fixed by the cuts.
  (b) Prescribing the inflow cuts and holding every outflow cut at the network's relative pressure is not possible
      from the record, because the record stores fractions rather than pressures. To test the window's totals as well,
      solve a longer window and compare only its inner rows.

Usage: python3 tools/cfd_kit.py RUNS_DIR NAME uturn J [ROWS]       (J-th U-turn, 0-based; ROWS each side, default 12)
       python3 tools/cfd_kit.py RUNS_DIR NAME seam [UNITS] [ROWS]  (last UNITS head units (default 3) + ROWS array rows)
       python3 tools/cfd_kit.py RUNS_DIR NAME rows K0 K1           (straight rows K0..K1, 1-based, within one leg)
"""
import csv
import json
import os
import pickle
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from opencif.physics import schedule as S  # noqa: E402
from opencif.geometry import dxf as io_dxf  # noqa: E402
from shapely.geometry import box, Point  # noqa: E402
from shapely.ops import unary_union  # noqa: E402


def leg_of(legs, k):
    """Return the leg that holds row k."""
    for lg in legs:
        if lg['rows'][0] <= k <= lg['rows'][1]:
            return lg
    raise ValueError(k)


def l2g(leg, gx0, x, y):
    """Map a point (x, y) of a leg's own frame, measured from gx0 along the leg and from the leg's axis Y across it, to
    the device's coordinates. A leg that runs toward -x is turned through 180 degrees."""
    if leg['dirn'] > 0:
        return (gx0 + x, leg['Y'] + y)
    return (gx0 - x, leg['Y'] - y)


def main(argv):
    """Cut the window that argv names out of a finished build, write its kit, and return a summary of it."""
    runs, name, kind = argv[0], argv[1], argv[2]
    d = os.path.join(runs, name)
    P = pickle.load(open(os.path.join(d, name + '.pkl'), 'rb'))
    rec = json.load(open(os.path.join(d, name + '.json')))
    arr, fluid = P['arrays'], P['fluid']
    spec = arr['spec']
    N, NM = int(arr['N']), int(arr['NM'])
    x = np.asarray(arr['x'], float)
    legs, gx = arr['legs'], arr['gx']
    lg, Lam = spec.l_g, spec.Lam_L
    fA, fB = np.asarray(arr['fA'], float), np.asarray(arr['fB'], float)
    keep = np.cumprod(1.0 - fA - fB)                      # the center flow after row k (fraction of the device inlet)
    Pin = np.concatenate([[1.0], keep[:-1]])
    qA, qB = fA * Pin, fB * Pin
    SA, SB = np.cumsum(qA), np.cumsum(qB)                 # side flows after row k
    Q_uL = float(getattr(spec, 'Q_in_uL_min', 30.0))
    D = S.KAPPA_G * spec.delta_star
    wc_at = lambda xp: float(np.interp(xp, x, arr['wc']))
    meas = {}
    for k, side, w_eq, w_int in arr.get('meas_lin', []):
        meas[(int(k), side)] = float(w_eq)

    def ws_drawn(k, side):
        if (k, side) in meas:
            return meas[(k, side)]
        w = arr['ws'][k] if k < len(arr['ws']) else np.nan
        return float(np.nan_to_num(w, nan=spec.GS))
    units = {int(u['k']): u for u in arr.get('units', [])}
    turns = list(arr.get('turns', []))
    head_cut = None
    if kind == 'uturn':
        j = int(argv[3])
        m = int(argv[4]) if len(argv) > 4 else 12
        kt = int(turns[j])
        k0, k1 = kt - m, kt + 1 + m
        tag = 'uturn%d' % j
    elif kind == 'seam':
        if NM == 0:
            raise SystemExit('no meander head: no seam')
        nu = int(argv[3]) if len(argv) > 3 else 3
        m = int(argv[4]) if len(argv) > 4 else 16
        u0 = NM - 1 - nu                                  # the window starts on the drop of head unit u0
        if u0 < 0 or units.get(u0, {}).get('kind') != 'T':
            raise SystemExit('the seam window needs taper units before the transition unit')
        head_cut = units[u0]
        k0, k1 = u0 + 1, min(NM + m, N - 2)
        tag = 'seam'
    else:
        k0, k1 = int(argv[3]) - 1, int(argv[4]) - 1
        tag = 'rows%d-%d' % (k0 + 1, k1 + 1)
    out = os.path.join(d, 'cfd_kit_' + tag)
    os.makedirs(out, exist_ok=True)
    ku, kd = k0 - 1, k1
    region, cuts = [], []
    reach = 40.0 + max([ws_drawn(k, 'A') for k in range(max(ku, NM), kd + 1)] + [spec.GS])

    def add_cuts(xm, leg, kseg, role, sides):
        jj = leg['j']
        off = leg['p0'] if leg['p0'] is not None else 0.0
        hw = wc_at(xm) / 2.0
        flows = dict(centre=keep[kseg], A=SA[kseg], B=SB[kseg])
        for ch, (ya, yb) in [('centre', (-hw, hw))] + sides:
            pa, pb = l2g(leg, gx[jj], xm - off, ya), l2g(leg, gx[jj], xm - off, yb)
            cuts.append(dict(channel=ch, row_after=kseg + 1, role=role, x0=pa[0], y0=pa[1], x1=pb[0], y1=pb[1],
                             width_um=abs(yb - ya), flow_fraction_of_inlet=float(flows[ch]),
                             flow_uL_min=float(flows[ch] * Q_uL)))
    # upstream cuts
    if head_cut is not None:
        # The cut crosses head unit u0 at the apex of its U-turn (x = X + s/2, where the turn's channel runs along +x),
        # so the window keeps the unit's return line, quarter arc, drop, and T-junction, and every later unit whole.
        u = head_cut
        w = float(u['w'])
        s_ = w + spec.wall_M
        xc = u['X'] + s_ / 2.0
        yc = float(u['R']) + s_ / 2.0
        add_cuts(xc, legs[0], u0, 'in', [('A', (yc - w / 2.0, yc + w / 2.0)), ('B', (-(yc + w / 2.0), -(yc - w / 2.0)))])
        Rmax = max(float(units[k]['R'] or 0.0) + spec.GS + 60.0 for k in range(u0, NM) if k in units)
        xb = 0.5 * (x[kd] + x[kd + 1])
        Y = max(Rmax, wc_at(xb) / 2 + lg + reach)
        region.append(box(xc, -Y, xb, Y))
    else:
        leg = leg_of(legs, ku)
        xm = 0.5 * (x[ku] + x[ku + 1])
        hw = wc_at(xm) / 2.0
        wa, wb = ws_drawn(ku, 'A'), ws_drawn(ku, 'B')
        add_cuts(xm, leg, ku, 'in', [('A', (hw + lg, hw + lg + wa)), ('B', (-(hw + lg + wb), -(hw + lg)))])
    # downstream cuts
    leg = leg_of(legs, kd)
    xm = 0.5 * (x[kd] + x[kd + 1])
    hw = wc_at(xm) / 2.0
    wa, wb = ws_drawn(kd, 'A'), ws_drawn(kd, 'B')
    add_cuts(xm, leg, kd, 'out', [('A', (hw + lg, hw + lg + wa)), ('B', (-(hw + lg + wb), -(hw + lg)))])
    # the window region, which is the band of each leg's own channels between the cuts, plus the loop annulus
    if head_cut is None:
        for leg in legs:
            a, b = max(leg['rows'][0], ku), min(leg['rows'][1], kd)
            if a > b:
                continue
            jj = leg['j']
            off = leg['p0'] if leg['p0'] is not None else 0.0
            xa = 0.5 * (x[a] + x[a + 1]) if a == ku else x[a] - Lam / 2.0 - 1.0
            xb = 0.5 * (x[b] + x[b + 1]) if b == kd else x[b] + Lam / 2.0 + 1.0
            Y = max(wc_at(xa), wc_at(xb)) / 2.0 + lg + reach
            pa, pb = l2g(leg, gx[jj], xa - off, -Y), l2g(leg, gx[jj], xb - off, Y)
            region.append(box(min(pa[0], pb[0]), min(pa[1], pb[1]), max(pa[0], pb[0]), max(pa[1], pb[1])))
        for lp in arr.get('loops', []):
            k = int(lp['k'])
            if ku < k < kd:
                R = float(lp['R'])
                rr_ = wc_at(x[k]) / 2.0 + lg + reach
                ring = Point(lp['cx'], lp['cy']).buffer(R + rr_, quad_segs=128).difference(
                    Point(lp['cx'], lp['cy']).buffer(max(R - rr_, 0.0), quad_segs=128))
                sgn = 1.0 if leg_of(legs, k)['dirn'] > 0 else -1.0
                half = box(min(lp['cx'], lp['cx'] + sgn * (R + rr_ + 10)), lp['cy'] - R - rr_ - 10,
                           max(lp['cx'], lp['cx'] + sgn * (R + rr_ + 10)), lp['cy'] + R + rr_ + 10)
                region.append(ring.intersection(half))
    reg = unary_union(region)
    win = fluid.intersection(reg)
    open_len = win.boundary.intersection(reg.boundary.buffer(1e-3)).length
    cut_len = sum(c['width_um'] for c in cuts)
    check = dict(open_boundary_um=open_len, cuts_um=cut_len, ok=abs(open_len - cut_len) <= 0.02 * cut_len + 2.0)
    if not check['ok']:
        print('WARNING: the window is open %.1f um along its boundary against %.1f um of cuts: see window.png' % (open_len, cut_len))
    io_dxf.write_dxf(win, os.path.join(out, 'window.dxf'))
    rows_out = []
    for k in range(k0, k1 + 1):
        leg = leg_of(legs, k)
        jj = leg['j']
        off = leg['p0'] if leg['p0'] is not None else 0.0
        if k < NM:                                        # a head gap, the stub from the trunk up to the junction J_k
            u = units.get(k)
            w = float(u['w']) if u else spec.GS
            ys = spec.wc0 / 2.0 + lg / 2.0
            xs0, xs1 = (u['X'] if u else x[k]) - w / 2.0 - 2.0, (u['X'] if u else x[k]) + w / 2.0 + 2.0
        else:
            hw = wc_at(x[k]) / 2.0
            ys = hw + lg / 2.0
            xs0, xs1 = x[k] - off - spec.G / 2.0 - 2.0, x[k] - off + spec.G / 2.0 + 2.0
        for side, sg, f_, q_ in (('A', 1.0, fA[k], qA[k]), ('B', -1.0, fB[k], qB[k])):
            pa, pb = l2g(leg, gx[jj], xs0, sg * ys), l2g(leg, gx[jj], xs1, sg * ys)
            dd = np.asarray(arr['dA' if side == 'A' else 'dB'], float)[k] / D
            rows_out.append(dict(row=k + 1, side=side, kind=str(arr['kind'][k]), x0=pa[0], y0=pa[1], x1=pb[0], y1=pb[1],
                                 f_of_arriving_centre_flow=float(f_), gap_flow_fraction_of_inlet=float(q_),
                                 gap_flow_uL_min=float(q_ * Q_uL), cutoff_over_D=float(dd) if np.isfinite(dd) else ''))
    for fn, rows in (('gaps.csv', rows_out), ('boundary.csv', cuts)):
        with open(os.path.join(out, fn), 'w', newline='') as fh:
            w = csv.DictWriter(fh, list(rows[0].keys()))
            w.writeheader()
            [w.writerow(r) for r in rows]
    qmin = min(abs(r['gap_flow_fraction_of_inlet']) for r in rows_out if r['gap_flow_fraction_of_inlet'] != 0)
    cmax = max(c['flow_fraction_of_inlet'] for c in cuts)
    readme = """# CFD validation kit: %s, %s (rows %d-%d)

Device record: `%s.json` (code %s).  Depth h = %.1f um; mu = 1e-3 Pa s; device inlet flow %.1f uL/min; D = %.4f um.
Window check: open boundary %.1f um = cuts %.1f um: %s.

**Boundary conditions.** No-slip on every wall (side walls, posts, top, bottom).  Prescribe the flows of `boundary.csv`
(fully developed rectangular-duct profiles over the drawn width) on every cut except the downstream centre cut, which
is the pressure reference (p = 0).  The flows balance (inflow = outflow).

**Comparison.** Integrate the solver's flux across each section of `gaps.csv` (it crosses the gap, or a head unit's
stub, at mid-depth of the post zone with 2 um of solid at each end) and compare with `gap_flow_fraction_of_inlet`.
Report per gap (q_solver - q_network)/q_network and the implied cutoff change (D ~ f^(1/2) at small f).  With all cut
flows prescribed, each side's TOTAL gap flow over the window is fixed: the test is how the window distributes it
(which is what the seam and turn rows get wrong when a junction term is wrong).

**Precision.** The smallest gap flow is %.2e of the device inlet flow, against channel flows up to %.2e: a 1 %% check of
a gap flow needs fluxes conserved to ~%.0e of the channel flow.  Check the window's own mass balance and mesh-converge
the gap fluxes, not only the pressure drop.

**What a difference calibrates.** The record's model band (`validation.model_scenarios`) shows how far the junction
terms' range moves these rows.  A 3-D result outside it calls for re-fitting the junction constants
(schedule: KAPPA_J_DEEP/SHALLOW for bends, C_T_DEEP/SHALLOW for the head's T-junctions, and the depth weight).
""" % (name, tag, k0 + 1, k1 + 1, name, rec.get('code', {}).get('sha256', '?'), spec.h, Q_uL, D, check['open_boundary_um'],
       check['cuts_um'], 'ok' if check['ok'] else 'NOT CLEAN', qmin, cmax, 0.01 * qmin / cmax)
    open(os.path.join(out, 'README.md'), 'w').write(readme)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon as MP
    fig, ax = plt.subplots(figsize=(10, 8))
    for g in (list(win.geoms) if hasattr(win, 'geoms') else [win]):
        if g.is_empty:
            continue
        ax.add_patch(MP(np.asarray(g.exterior.coords), closed=True, fc='#c9dcf0', ec='#1f4e79', lw=0.4))
        for hole in g.interiors:
            ax.add_patch(MP(np.asarray(hole.coords), closed=True, fc='white', ec='#1f4e79', lw=0.4))
    for c in cuts:
        ax.plot([c['x0'], c['x1']], [c['y0'], c['y1']], color='#c0392b' if c['role'] == 'in' else '#27ae60', lw=2)
    for r in rows_out:
        ax.plot([r['x0'], r['x1']], [r['y0'], r['y1']], color='#7f8c8d', lw=0.8)
    b = win.bounds
    ax.set_xlim(b[0] - 20, b[2] + 20)
    ax.set_ylim(b[1] - 20, b[3] + 20)
    ax.set_aspect('equal')
    ax.set_title('%s %s: rows %d-%d; red = inflow cuts, green = outflow cuts, grey = gap sections' % (name, tag, k0 + 1, k1 + 1),
                 fontsize=9, loc='left')
    fig.tight_layout()
    fig.savefig(os.path.join(out, 'window.png'), dpi=160)
    plt.close(fig)
    res = dict(out=out, rows=[k0 + 1, k1 + 1], gaps=len(rows_out), cuts=len(cuts), window_check=check,
               inflow=sum(c['flow_fraction_of_inlet'] for c in cuts if c['role'] == 'in'),
               outflow=sum(c['flow_fraction_of_inlet'] for c in cuts if c['role'] == 'out'))
    print(json.dumps(res, indent=1, default=float))
    return res


if __name__ == '__main__':
    main(sys.argv[1:])
