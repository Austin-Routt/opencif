"""opencif.geometry.depicted: the heads as the papers' figures draw them, and what they would do if built that way.

For each depicted row, the figure gives a unit geometry (the reach of a hairpin, or the lane and top of an arm), and
hence a drawn centerline length L_dep(k) in this work's own unit geometry and a side resistance
R_s,dep(k) = r(w) L_dep(k). The gap fractions f(k) at which the ladder's exact closure (schedule.build, the closure of
spec v5, section 13; internal patch RP-B5) produces exactly those resistances follow from the exact closure of row k,

        S(k) R_s(k) = [P_c(k) - P_c(k+1)] - R_gap [Q_gap(k) - Q_gap(k+1)],

which is linear in the cumulative side flows S. A tridiagonal solve, repeated three times so that f* of the first
linear row follows the taper, gives S(k) (f_head_from_lengths), and then
f(k) = (S(k) - S(k-1)) / P(k-1). An earlier fixed-point iteration, which relied on the gap term being small (M13),
diverged on Xia's 97-row head. S(1) = f(1) is the gap-1 extraction that the depicted first unit implies. The result,
f_head, is what a device built to the figure would extract at each of those gaps. Dividing it by the constant-cutoff
f*(k) gives the over-extraction (> 1) or under-extraction (< 1) of the depicted head, and kappa x(f) w_c gives the
cutoff that each row would have.

The sources are the Dinh 2024 Fig. 2A inset (r4/meas_dinh_inset.py, with 18 rows, 3 arms, and 9.89 um/px) and the
Xia 2016 Fig. 3b micrograph (r4/meas_xia_head2.py, with 10 rows). The micrograph is read at 1.798 um/px, the scale
that the printed widths and loop curvature of panels b-g give (data/r6/rc6_values.json), whereas rc5 read it at the
bar of panel i, 1.953 um/px."""
import json
import math
import os
import numpy as np
from opencif.physics import schedule as S
from opencif.geometry import shapes as Gm
from opencif import published as D

HERE = os.path.dirname(os.path.abspath(__file__))
R4 = os.path.join(HERE, '..', 'data', 'r4')

# ---- Dinh 2024 Fig. 2A inset (native px, axis 192.6, row-1 mouth at y 599.5, 9.89 um/px) ---------------------
DINH_SC, DINH_AX = 9.89, 192.6
# the white finger extents (left x_start, right x_stop) of the taper rows 4 to 18 (meas_dinh_inset.py, threshold 160)
DINH_FIN_L = [36, 55, 70, 82, 92, 100, 107, 113, 118, 122, 126, 130, 133, 136, 138]
DINH_FIN_R = [350, 331, 316, 304, 294, 286, 279, 273, 268, 263, 259, 256, 253, 250, 248]
# the arms (rows 1 to 3), with the top of each black U-turn at the inner lane columns, and the row-1 mouth
DINH_ARM_TOP_PX = [453.0, 534.0, 594.0]
DINH_ROW1_PX = 599.5


def dinh_depicted_geometry(spec):
    """Return each depicted row's kind and its reach R (taper) or its lane R and top t (arm), in this work's unit
    geometry."""
    sa, w = spec.s_a, spec.w_arm
    rows = []
    for k, t_px in enumerate(DINH_ARM_TOP_PX):
        x_top_outer = (t_px - DINH_ROW1_PX) * DINH_SC - 3.0       # the black edge (the raster thins fluid ~3 um/side)
        t = x_top_outer + sa / 2.0 + w / 2.0                        # U-turn center = outer edge + s_a/2 + w/2
        rows.append(dict(k=k, kind='A', R=spec.R0 + 2 * k * sa, t=t, n=1))
    for j, (xl, xr) in enumerate(zip(DINH_FIN_L, DINH_FIN_R)):
        reach_white = 0.5 * ((xr - 0.5 - DINH_AX) + (DINH_AX - (xl - 0.5))) * DINH_SC - 3.0
        R = reach_white - spec.wall_M / 2.0                          # finger end = U-turn center + wall/2
        rows.append(dict(k=3 + j, kind='T', R=R))
    return rows


# ---- Xia 2016 Fig. 3b micrograph (native px, axis 242, 1.798 um/px from data/r6/rc6_values.json) --------------
_R6 = json.load(open(os.path.join(HERE, '..', 'data', 'r6', 'rc6_values.json')))
XIA_SC, XIA_AX = _R6['xia']['scale_um_per_px'], 242.0
XIA_TIP_PX = [723, 677, 644, 618, 596, 577, 559, 546, 533, 523]     # outer edge of each U-turn (rows 1 to 10)


# On Fig. 3a, read at its 5 mm bar (r4/meas_xia_fig3a.py), the head's black band starts 1.39 mm below the center of
# the inlet port. Its rendered edge lies half a drawn line (2 px, 56 um) outside the true edge, so the first unit's
# top is 1.45 mm below the port's center. The first pill lies 13.75 +- 0.2 mm below the port's center (box c, located
# on the Fig. 3c micrograph).
_M3A = json.load(open(os.path.join(R4, 'meas_xia_fig3a.json')))
# Fig. 3a's bar is read with its anti-aliased ends (178.32 px, where r4 used 179), which scales every Fig. 3a length
# by 1.003827.
XIA_HEAD_TOP_MM = _R6['xia']['layout']['head_top_below_inlet_mm']          # 1.453 (1.447 in rc5)
XIA_FIRST_PILL_MM = _R6['xia']['layout']['first_pill_below_inlet_mm']      # 13.786 (13.734 in rc5)


def xia_depicted_geometry(spec, full=True, n_rows=12):
    """Return each depicted row of Xia's head, together with the tip law. Units 1 to 10 come from Fig. 3b (the
    measured tips). The measured tips follow tip(k) = a - b ln k (rms 12 um = 6 px over k = 1 to 10,
    val/xia_tip_sensitivity.py), and the further units continue this law, shifted to meet unit 10, until the smallest
    drawable unit (U-turn lane R = yJ + s). The law reaches that unit at about unit 33, which is where Fig. 3a's band,
    too, has shrunk to the smallest unit. The rows from there to Fig. 3a's first pill are rectilinear (blocks with a
    GS-wide side channel; Figs. 1c and 3c). With full=False, the function returns the earlier 12-unit head
    (Fig. 1b(i)), which is kept for comparison."""
    meas = [(x - XIA_AX) * XIA_SC for x in XIA_TIP_PX]
    k = np.arange(1, len(meas) + 1)
    A_ = np.vstack([np.ones_like(k, dtype=float), -np.log(k)]).T
    a, b = np.linalg.lstsq(A_, np.array(meas), rcond=None)[0]
    tip_min = spec.yJ(spec.GS) + spec.s + spec.s / 2.0 + spec.GS / 2.0
    tips = list(meas)
    if full:
        off = (a - b * math.log(len(meas))) - meas[-1]          # continuous at the last measured unit
        while True:
            t = a - b * math.log(len(tips) + 1) - off
            if t < tip_min:
                break
            tips.append(t)
    else:
        d = np.diff(tips)
        q = (d[-1] / d[-2] + d[-2] / d[-3]) / 2.0
        while len(tips) < n_rows:
            d = np.append(d, d[-1] * q)
            tips.append(tips[-1] + d[-1])
    rows = []
    for j, tip in enumerate(tips):
        R = tip - spec.s / 2.0 - spec.GS / 2.0                       # outer edge = U-turn center + s/2 + GS/2
        rows.append(dict(k=j, kind='T', R=R, tip_um=tip, extrapolated=(j >= len(XIA_TIP_PX))))
    if full:
        n_tot = int(round((XIA_FIRST_PILL_MM - XIA_HEAD_TOP_MM) * 1000.0 / (2.0 * spec.s)))
        for j in range(len(tips), n_tot):
            rows.append(dict(k=j, kind='R', R=None, extrapolated=False))
    return rows, dict(a=float(a), b=float(b), offset=float(off) if full else 0.0, tip_min=float(tip_min))


def unit_length(spec, row, X, ell):
    """Return the drawn centerline length of a depicted unit, an arm's by arm_length and a taper unit's as
    2 (R - yJ) + taper_fixed."""
    if row['kind'] == 'A':
        return Gm.arm_length(spec, X, ell, row['R'], row['n'], row['t'], spec.w_arm)
    yJ = spec.yJ(spec.GS)
    return 2.0 * (row['R'] - yJ) + Gm.taper_fixed(spec, ell, spec.GS)


def f_head_from_lengths(spec, rows, iters=60):
    """Return f(1) to f(n), at which the exact closure gives R_s(k) = r(w) L_dep(k) for k = 1 to n-1 and the n-th
    row's unit is a taper unit of the depicted reach as well (its closure involves f(n+1), the constant-cutoff value).
    The schedule, the closure and depicted lengths, and the iteration count are returned with it."""
    h = spec.h
    n = len(rows)
    A = sum(1 for r in rows if r['kind'] == 'A')
    nR = sum(1 for r in rows if r['kind'] == 'R')
    over = {'NM_force': None} if nR == 0 else {'NM_force': n - nR, 'n_rect': nR}
    wR = spec.ws_rect or spec.GS
    # The head's geometry (pitches, R_c, R_gap) does not depend on f, and the exact closure of each row,
    #     S_k R_s,k = R_c,k (1 - 2 S_k) - R_g,k (S_k - S_k-1) + R_g,k+1 (S_k+1 - S_k),
    # is linear in the cumulative side flows S. Therefore one tridiagonal solve gives the S (and hence the f) that
    # make the closure draw exactly the depicted resistances. The row after the head is at the constant cutoff,
    # S_n+1 = S_n + f*(1 - 2 S_n). A damped fixed-point iteration, which was used before, diverged on Xia's 97-row
    # head.
    fs0 = S.f_of(spec.delta_star, spec.wc0)
    f = np.full(n, fs0)
    sc = S.build(spec.__class__(**{**spec.__dict__, 'f_head': tuple(f), **over}), f[0], A)
    n_head = n - nR

    def bend_x(r):          # the depicted unit's arcs and end T-junction carry the model's corrections too
        if r['kind'] not in ('A', 'T'):
            return 0.0
        w_ = spec.GS if r['kind'] == 'T' else spec.w_arm
        u = dict(kind=r['kind'], X=float(sc.x[r['k']]), ell=float(sc.ell[r['k']]), w=w_, R=r['R'], n=r.get('n'), t=r.get('t'))
        # The last head row is the transition unit only when the head ends at the depicted rows, as in Xia, where
        # rectilinear rows follow. Otherwise, as in Dinh, the closure continues the head and every depicted row is an
        # ordinary unit.
        last_is_transition = nR > 0 and r['k'] == n_head - 1
        tj = 0.0 if last_is_transition else Gm.junction_excess(spec, u, w_)
        return Gm.bend_excess(spec, Gm.unit_path(spec, u, from_J=True), w_) + tj
    Rs_dep = np.array([(S.rr(wR, h) * sc.ell[r['k']] if r['kind'] == 'R' else
                        S.rr(spec.GS if r['kind'] == 'T' else spec.w_arm, h) *
                        (unit_length(spec, r, sc.x[r['k']], sc.ell[r['k']]) + bend_x(r))) * S.UM for r in rows])
    it = 0
    for it in range(1, 4):                     # f* of the first linear row depends weakly on the taper, so refine
        Rc, Rg = np.asarray(sc.Rc, float), np.asarray(sc.Rg, float)
        fs1 = S.f_of(spec.delta_star, sc.wc[n])
        M = np.zeros((n, n))
        bvec = Rc[:n].copy()
        for k in range(n):
            M[k, k] = Rs_dep[k] + 2.0 * Rc[k] + Rg[k] + Rg[k + 1]
            if k > 0:
                M[k, k - 1] = -Rg[k]
            if k < n - 1:
                M[k, k + 1] = -Rg[k + 1]
            else:
                M[k, k] -= Rg[k + 1] * (1.0 - 2.0 * fs1)
                bvec[k] += Rg[k + 1] * fs1
        Ssol = np.linalg.solve(M, bvec)
        Pp = np.concatenate([[1.0], 1.0 - 2.0 * Ssol[:-1]])
        f = np.diff(np.concatenate([[0.0], Ssol])) / Pp
        sc = S.build(spec.__class__(**{**spec.__dict__, 'f_head': tuple(f), **over}), f[0], A)
    sc = S.build(spec.__class__(**{**spec.__dict__, 'f_head': tuple(f), **over}), f[0], A)
    wd = lambda r: wR if r['kind'] == 'R' else (spec.GS if r['kind'] == 'T' else spec.w_arm)
    L_req = np.array([sc.Rs[r['k']] / (S.rr(wd(r), h) * S.UM) for r in rows])
    L_dep = np.array([sc.ell[r['k']] if r['kind'] == 'R' else unit_length(spec, r, sc.x[r['k']], sc.ell[r['k']])
                      for r in rows])
    return f, sc, L_req, L_dep, it


def report(name, spec, rows):
    """Return the record of a depicted head, with each row's depicted and closure lengths, f, f/f*, and cutoff."""
    f, sc, L_req, L_dep, it = f_head_from_lengths(spec, rows)
    fstar = np.array([S.f_of(spec.delta_star, w) for w in sc.wc[:len(rows)]])
    Dstar = S.KAPPA_G * spec.delta_star
    cut = np.array([S.KAPPA_G * S.x_of(v) * w if 0 < v < 0.5 else float('nan') for v, w in zip(f, sc.wc)])
    out = dict(name=name, iterations=it, rows=[], D_star_um=Dstar, f_star=float(fstar[0]))
    for r, fi, fs, c, lr, ld in zip(rows, f, fstar, cut, L_req, L_dep):
        out['rows'].append(dict(row=r['k'] + 1, kind=r['kind'], R_um=r['R'], t_um=r.get('t'),
                                extrapolated=r.get('extrapolated', False), L_dep_um=float(ld),
                                L_closure_um=float(lr), f=float(fi), f_over_fstar=float(fi / fs),
                                cutoff_um=float(c), cutoff_over_Dstar=float(c / Dstar)))
    out['f_head'] = [float(v) for v in f]
    out['S_after_head'] = float(sc.S[len(rows) - 1])
    out['NM'] = int(sum(1 for r in rows if r['kind'] != 'R'))
    out['n_rect'] = int(sum(1 for r in rows if r['kind'] == 'R'))
    out['ws_first_pill_um'] = float(sc.ws[len(rows)]) if len(rows) < sc.N - 1 else None
    return out


if __name__ == '__main__':
    # Run as a script, the module compares its result with the shipped data/depicted.json and writes the file only
    # with --write (codereview2 N3). The caches are cleared before each entry (BUILD-1), the dinh entry carries the
    # turns metadata of the regeneration script, and the xia entry is kept from the file byte for byte unless --xia is
    # given. provenance/rc7_scripts/regen_depicted_dinh_square.py wrote the shipped file this way.
    import sys
    path_ = os.path.join(os.path.dirname(HERE), 'data', 'depicted.json')
    with open(path_) as fh_:
        shipped = json.load(fh_)
    res = {}
    S.clear_caches()
    sp = D.dinh_fig2('V1')
    # The head is drawn with square quarter turns (outer 7 um, inner 10 um), so its depicted fractions are derived for
    # them (DINH-1, code review R2). The shipped entry is written by
    # provenance/rc7_scripts/regen_depicted_dinh_square.py, which also keeps the 'xia' entry byte for byte (a run of
    # this file recomputes it in one process).
    sp = sp.__class__(**{**sp.__dict__, 'turns': 'square', 'r_corner_out': 7.0, 'r_corner_in': 10.0})
    res['dinh'] = report('Dinh 2024 Fig. 2A inset (V1), square quarter-turns (outer 7 um, inner 10 um; rc7 DINH-1, code '
                         'review R2)', sp, dinh_depicted_geometry(sp))
    res['dinh']['turns'] = dict(turns='square', r_corner_out_um=7.0, r_corner_in_um=10.0,
                                note='rc7 (code review R2): regenerated for the square-turn head (rc6 derived it for arc '
                                     'turns; see provenance/rc7_scripts/regen_depicted_dinh_square.py)')
    if '--xia' in sys.argv:
        S.clear_caches()
        sx = D.xia_fig()
        rows_x, law = xia_depicted_geometry(sx)
        res['xia'] = report('Xia 2016 Figs. 3b + 3a (RS-7): measured units, extrapolated units, rectilinear rows', sx, rows_x)
        res['xia']['tip_law'] = law
    else:
        res['xia'] = shipped['xia']
    differs = []
    for key_ in ('dinh', 'xia'):
        if key_ == 'xia' and '--xia' not in sys.argv:
            print('xia: kept from the file (not recomputed; --xia recomputes it)')
            continue
        same_ = json.dumps(res[key_], sort_keys=True) == json.dumps(shipped[key_], sort_keys=True)
        print('%s: %s the shipped entry' % (key_, 'equals' if same_ else 'DIFFERS FROM'))
        if not same_:
            differs.append(key_)
    # A differing entry moves every mask that uses it, so it is written only on purpose (codereview3 m3).
    if '--write' in sys.argv and differs and '--force' not in sys.argv:
        print('not written: %s differ%s from the shipped file (every mask drawn from it would move); add --force to '
              'write anyway' % (', '.join(differs), 's' if len(differs) == 1 else ''))
        sys.exit(1)
    if '--write' in sys.argv:
        with open(path_, 'w') as fh_:
            json.dump(res, fh_, indent=1)                     # read by spec.py
        print('written', path_)
    for key, r in res.items():
        print(r['name'], 'iterations', r['iterations'], 'D*', round(r['D_star_um'], 3))
        for row in r['rows']:
            print('  row %2d %s R %7s L_dep %8.1f L_closure %8.1f  f/f* %7.3f  cutoff/D* %6.3f%s' % (
                row['row'], row['kind'], row['R_um'], row['L_dep_um'], row['L_closure_um'], row['f_over_fstar'],
                row['cutoff_over_Dstar'], '  (extrapolated)' if row['extrapolated'] else ''))
