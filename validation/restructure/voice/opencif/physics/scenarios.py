"""opencif.physics.scenarios: what the 1-D design cannot know, that is, the model and fabrication errors, propagated
through the ladder.

The scenarios perturb the design network (a frozen_ref-style solve, with both sides equal) in four ways:

    junc    the 1-D junction convention, with each meander unit's length moved by +- w/2 (one junction's half-width)
    bias    a uniform fabrication bias dw on every lateral width (channels, gaps, and posts shrink accordingly)
    rgap    the meander-mouth R_gap model off by +-20 % relative to the linear gaps
    noise   random per-segment side-width noise (sigma, um), which stands for the mask address grid and edge roughness

The report gives f_em/f at the seam rows (and the loops) and at gap 1, the cutoff, and F."""
import math
import numpy as np
from opencif.physics import schedule as S
from opencif.physics import network as NW


def arc_length_unit(spec, u):
    """Return the centerline arc length (um) inside one meander unit (its U-turns and quarter arcs), from the
    construction."""
    w = u['w']
    s = w + spec.wall_M
    if u['kind'] in ('T', 'X'):
        return math.pi * s / 2.0 + math.pi * (s / 2.0) / 2.0
    n = u['n']
    return math.pi * (s / 2.0) / 2.0 * 3 + math.pi * s / 2.0 * n + math.pi * s / 2.0 * (n - 1)


def solve(sc, Rc=None, Rs=None, Rg=None):
    """Solve the design network of the schedule sc, with any of Rc, Rs, and Rg replaced by the values given."""
    Rc = sc.Rc[:-1] if Rc is None else Rc
    Rs = sc.Rs if Rs is None else Rs
    Rg = sc.Rg if Rg is None else Rg
    return NW.solve(Rc, Rs, Rg, sc.Rc_out, sc.Rs_out)


def scenarios(sc, units, eps_b=0.032, dw=0.5, rgap_rel=0.2, noise=0.01, seed=1):
    """Return, for each error scenario by name, the ratio f_em/f of the perturbed design network's gap fractions to
    the unperturbed ones, row by row."""
    spec = sc.spec
    h = spec.h
    NM = sc.NM
    out = {}
    base = solve(sc)['f']
    # The old 'bend +3.2 %' scenario is gone, because the design now carries the exact curved-duct correction, and its
    # remaining uncertainty is reported on the measured network by build.model_scenarios.
    # the junction convention
    for sgn in (+1, -1):
        Rs = sc.Rs.copy()
        for u in units:
            Rs[u['k']] += sgn * S.rr(u['w'], h) * (u['w'] / 2.0) * S.UM
        out['junction %+d w/2' % sgn] = solve(sc, Rs=Rs)['f'] / base
    # a uniform fabrication bias on every lateral width
    for d in (+dw, -dw):
        Rc = np.array([S.rr(0.5 * (sc.wc[k] + sc.wc[min(k + 1, sc.N - 1)]) + d, h) * sc.ell[k] * S.UM
                       for k in range(sc.N - 1)])
        Rs = sc.Rs.copy()
        for u in units:
            Rs[u['k']] *= S.rr(u['w'] + d, h) / S.rr(u['w'], h)
        for k in range(NM, sc.N - 1):
            if np.isfinite(sc.ws[k]):
                Rs[k] = S.rr(sc.ws[k] + d, h) * sc.ell[k] * S.UM
        Rg = np.array([S.r_gap(sc.wL[k] + d, spec.l_g - d, h, spec.c_ent) for k in range(sc.N)])
        for k in sc.info.get('turns', []):                  # loops, scaled by the mean of their two side widths
            wl = sc.info.get('loop_w', {}).get(k, sc.ws[k - 1])
            Rs[k] = sc.Rs[k] * S.rr(wl + d, h) / S.rr(wl, h)
        wce, wse = sc.wc[-1], sc.ws[sc.N - 2]
        r = NW.solve(Rc, Rs, Rg, sc.Rc_out * S.rr(wce + d, h) / S.rr(wce, h),
                     sc.Rs_out * S.rr(wse + d, h) / S.rr(wse, h))
        out['bias %+.1f um' % d] = r['f'] / base
    # the meander-mouth R_gap model
    for sgn in (+1, -1):
        Rg = sc.Rg.copy()
        mouth = np.array([k < NM or sc.kind[k] == 'R' for k in range(sc.N)])
        Rg[mouth] *= (1 + sgn * rgap_rel)
        out['mouth R_gap %+d%%' % (sgn * 100 * rgap_rel)] = solve(sc, Rg=Rg)['f'] / base
    # the random per-segment side-width noise
    rng = np.random.default_rng(seed)
    Rs = sc.Rs.copy()
    for k in range(NM, sc.N - 1):
        if np.isfinite(sc.ws[k]):
            Rs[k] = S.rr(sc.ws[k] + rng.normal(0, noise), h) * sc.ell[k] * S.UM
    out['noise %.3f um' % noise] = solve(sc, Rs=Rs)['f'] / base
    return out


def summarize(sc, res, rows_seam=6):
    """Tabulate each scenario's f_em/f at the seam (rows_seam rows on either side of NM), at gap 1, and in the bulk
    beyond row NM + 20 (the rms and the maximum deviation from 1)."""
    NM = sc.NM
    lo, hi = NM - rows_seam, NM + rows_seam
    tab = []
    for name, r in res.items():
        seam = r[lo:hi]
        tab.append(dict(name=name, seam_max=float(np.max(seam)), seam_min=float(np.min(seam)),
                        row1=float(r[0]), bulk_rms=float(np.sqrt(np.mean((r[NM + 20:] - 1) ** 2))),
                        bulk_max=float(np.max(np.abs(r[NM + 20:] - 1)))))
    return tab


def failsafe_rows(sc, window=None):
    """Return the rows on which the margin may act. By default (window None), these are every row 2 .. N-1, plus gap 1
    when the inlet is sealed (f1 = f*).

    A fixed window around the seam and the loops (window = n rows) was tried first and is not enough. A uniform
    +0.5 um bias still raised the cutoff 1.3 % just outside a 20-row window on the Gifford bead device, because the
    disturbance that a seam releases decays over tens of rows at f ~ 3e-4."""
    rows = set(range(1, sc.N - 1))
    if window is not None:
        rows = set(range(1, min(sc.N - 1, sc.NM + window)))
        for k in sc.info.get('turns', []):
            rows |= set(range(k - window, min(sc.N - 1, k + window)))
    if abs(sc.f_seed / sc.info['fstar'] - 1.0) < 1e-9:
        rows.add(0)
    return sorted(rows)


def failsafe(spec, sc, units, window=None, iters=8, tol=1e-3, **kw):
    """Make the head, the seam, and the loops fail-safe, which is Xia 2016's 'pull-in' made quantitative.

    For every row in failsafe_rows, the design f is lowered by the largest over-extraction that any (non-random) error
    scenario produces there, relative to the constant-cutoff f*(i) of that row, and the device is re-planned with the
    seam fixed. This repeats until no scenario exceeds f*(i) by more than tol (relative) on those rows. Rows may end
    with f < 0 by design, which returns fluid to the center channel and is safe for retention. Random per-segment
    noise is reported rather than margined, because it has no sign."""
    from opencif.geometry import shapes as Gm
    mult = {}
    sc_cur, units_cur = sc, units
    hist = []
    for it in range(iters):
        res = scenarios(sc_cur, units_cur, **kw)
        rows = failsafe_rows(sc_cur, window)
        # The constant-cutoff f*(i) is that of the current plan. With a constant-shear taper, w_c(i) follows the
        # margined march, so f*(i) must be recomputed every iteration (a stale f* left +0.37 % f on Xia under the
        # +0.5 um bias).
        fstar_rows = np.array([S.f_of(spec.delta_star, w) for w in sc_cur.wc])
        worst_rel = 0.0
        for i in rows:
            worst = max(r[i] * sc_cur.f[i] for name, r in res.items() if not name.startswith('noise'))
            excess = worst - fstar_rows[i]
            worst_rel = max(worst_rel, excess / fstar_rows[i])
            if excess > 0:
                mult[i] = (sc_cur.f[i] - excess) / fstar_rows[i] if i not in mult else \
                    mult[i] - excess / fstar_rows[i]
        hist.append(worst_rel)
        if worst_rel <= tol and it > 0:
            break
        st = sc.info.get('st')
        # The re-plan keeps the narrow rows, because a dropped K_narrow once turned them back into 28 um units.
        sc_cur = S.build(spec, sc.f_seed, sc.A, st, f_mult=mult, NM_force=sc.NM,
                         K_narrow=sc.info.get('K_narrow', 0))
        if 'R_turn' in sc.info:                       # keep the layout, that is, the loops and the outlet lead
            Rt = sc.info.get('R_turns') or sc.info['R_turn']
            sc_cur = S.reclose(sc_cur, sc.info.get('turns', []), Rt, sc.info.get('Rc_lead', 0.0) + spec.R_ext_c)
            sc_cur.info['Rc_lead'] = sc.info.get('Rc_lead', 0.0)
            sc_cur.info['loop_w'] = sc.info.get('loop_w', {})
            # The layout keys are copied, so that the flare survives the fail-safe (rc5 review).
            for key in ('Rc_lead_len', 'leg_Y', 'leg_d', 'layout', 'ring',
                        'L_out_base', 'w_lead', 'flare_len', 'flare_kind'):
                if key in sc.info:
                    sc_cur.info[key] = sc.info[key]
        units_cur = Gm.plan_head(sc_cur)
        if units_cur:
            Gm.solve_transition(sc_cur, units_cur[-1])
    sc_cur.info['failsafe_mult'] = mult
    sc_cur.info['failsafe_rows'] = rows
    sc_cur.info['failsafe_history'] = hist
    sc_cur.info['A_need_sealed'] = sc.info.get('A_need_sealed')
    sc_cur.info['A_cap'] = sc.info.get('A_cap')
    return sc_cur, units_cur
