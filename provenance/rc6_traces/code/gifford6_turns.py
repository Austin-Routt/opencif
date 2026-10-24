"""Gifford et al. 2014, Fig. 6 micrograph: the turns beyond the legs' ends, and the leg axes they imply.

Beyond the divider walls the channels of legs 5, 3 and 1 curve away (their U-turns) while leg 7's run straight on (the
outlet).  Each wall line of the turns is followed on the image (the darkness ridge, step by step from just beyond the
walls' straight starts, in the rows' frame of panels_gifford6.py) and all of them are fitted with circles about ONE
common centre (each wall its own radius).  If the fit holds to a fraction of a pixel, the turns are concentric arcs
about that centre, as opencif's paperclip draws the loops at one end (layout.py: "every loop turning the same way ...
so loops at one end are concentric"); a semicircular turn about it then takes leg 2j - 1 at Y to leg 2j at
2 Y_c - Y.  The legs' axes 7, 5, 3, 1 are panels_gifford6's (the post rows' tip lines); legs 2, 4, 6 follow from the
centre.  Gifford prints no layout dimension (the schematic is not to scale, p. 4502): these are readings of the
micrograph; rc6's Fig. 6 presets take them (provenance/rc6_scripts/gifford_layout_measurement.json,
apply_rc6_layout_presets.py).

usage: python3 -I gifford6_turns.py ROOT RASTERS OUT.json [OUT.png]; panels_gifford6.py imports follow_walls and
fit_common for the turn comparison against a mask."""
import json
import os
import sys

import numpy as np
from scipy import optimize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import panels_gifford6 as P
import trace_lib as T

DU = P.DU
U_START, U_ARC = 70.0, 60.0          # px: just beyond the walls' starts (u ~ 83); the arcs are fitted from u <= 60


def ridge(D, us, vs, u, vguess, half=2.5):
    j = int(round((u - us[0]) / DU))
    i0, i1 = int(round((vguess - half - vs[0]) / DU)), int(round((vguess + half - vs[0]) / DU))
    if j < 0 or i0 < 0 or i1 >= len(vs):
        return None
    p, vv = D[i0:i1, j], vs[i0:i1]
    if len(p) < 5 or np.any(np.isnan(p)):
        return None
    i = int(np.nanargmax(p))
    if p[i] < 15 or i in (0, len(p) - 1):
        return None
    a, b, c = p[i - 1], p[i], p[i + 1]
    den = a - 2 * b + c
    return float(vv[i] + (0.5 * (a - c) / den * DU if den < 0 else 0.0))


def follow_walls(D, us, vs, v_top_leg):
    """every wall line crossing u = U_START, followed step by step along -u on its darkness ridge; returns the walls of
    the top leg (straight on: the outlet) and those of the turns, each dict(v_start, points (u, v) px, turned_deg)"""
    j0 = int(round((U_START - us[0]) / DU))
    col = D[:, j0]
    starts = [float(vs[i]) for i in range(1, len(vs) - 1)
              if not np.isnan(col[i]) and col[i] >= col[i - 1] and col[i] > col[i + 1] and col[i] > 25]
    walls = []
    for v0 in starts:
        pts, v, dv = [(U_START, v0)], v0, 0.0
        for u in np.arange(U_START - 1.0, -300.0, -1.0):
            vn = ridge(D, us, vs, u, v + dv)
            if vn is None:
                break
            dv = 0.7 * dv + 0.3 * (vn - v)
            v = vn
            pts.append((float(u), v))
        C = np.array(pts)
        if len(C) <= 60:
            continue
        turned = float(np.degrees(np.arctan(abs(C[-1, 1] - C[-6, 1]) / 5.0)))
        walls.append(dict(v_start=v0, points=C, turned_deg=turned))
    straight = [w for w in walls if w['v_start'] < v_top_leg + 70.0]          # the top leg (7) and its outer walls
    arcs = [w for w in walls if w['v_start'] >= v_top_leg + 70.0]
    return straight, arcs


def fit_common(Ps, c0=(80.0, 1200.0)):
    """circles about ONE common centre, each point set (u, v px) its own radius: returns centre, its SD, radii, residuals"""
    def resid(x):
        return np.concatenate([np.hypot(Q[:, 0] - x[0], Q[:, 1] - x[1]) - x[2 + i] for i, Q in enumerate(Ps)])
    x0 = list(c0) + [float(np.hypot(Q[0, 0] - c0[0], Q[0, 1] - c0[1])) for Q in Ps]
    sol = optimize.least_squares(resid, x0)
    cov = np.linalg.inv(sol.jac.T @ sol.jac) * np.var(sol.fun)
    return (float(sol.x[0]), float(sol.x[1])), np.sqrt(np.diag(cov))[:2], np.asarray(sol.x[2:], float), sol.fun


def main():
    ROOT, RAS, OUT = sys.argv[1:4]
    PNG = sys.argv[4] if len(sys.argv) > 4 else None
    M = json.load(open(os.path.join(ROOT, 'provenance', 'rc6_scripts', 'gifford_fig6_measurement.json')))
    r = T.Raster(os.path.join(RAS, 'gif7-000.png'), P.FIG6_SHA)
    th, s_um = float(M['rows_direction_deg']), float(M['um_per_px'])
    D, us, vs = P.resample(r.im, th)
    img = P.image_legs(D, us, vs, M)
    straight, arcs = follow_walls(D, us, vs, img[0]['centre_v'])
    walls = straight + arcs
    Ps = [w['points'][w['points'][:, 0] <= U_ARC] for w in arcs]
    (cu, cv), sd, radii, fun = fit_common(Ps)
    # each wall alone (its own centre), as a check of the common centre
    own = []
    for Q in Ps:
        s1 = optimize.least_squares(lambda x: np.hypot(Q[:, 0] - x[0], Q[:, 1] - x[1]) - x[2],
                                    [cu, cv, float(np.hypot(Q[0, 0] - cu, Q[0, 1] - cv))])
        own.append(dict(centre_v=float(s1.x[1]), centre_u=float(s1.x[0]), radius_um=float(s1.x[2] * s_um)))
    names = {q['name']: q['centre_v'] for q in img}
    Y = {7: -(names['bottom'] - names['top']) * s_um, 5: -(names['bottom'] - names['second']) * s_um,
         3: -(names['bottom'] - names['third']) * s_um, 1: 0.0}
    Yc = (cv - names['bottom']) * s_um
    leg_Y = [Y[1], 2 * Yc - Y[1], Y[3], 2 * Yc - Y[3], Y[5], 2 * Yc - Y[5], Y[7]]
    # where the turn's centre lies along the flow: against the legs' last pills (their mean u) and the walls' starts
    last_pill_u = float(np.mean([q['last_pill_u'] for q in img[1:]]))
    out = dict(
        method=__doc__.split('\n\n')[1].replace('\n', ' '),
        raster_sha256=r.sha256, um_per_px=s_um, rows_direction_deg=th,
        walls_followed=len(walls), walls_straight_top_leg=len(straight), walls_in_the_turns=len(arcs),
        points=int(sum(len(Q) for Q in Ps)), residual_sd_px=float(np.std(fun)), residual_sd_um=float(np.std(fun) * s_um),
        centre_uv_px=[cu, cv], centre_sd_px=[float(sd[0]), float(sd[1])],
        wall_radii_um=[float(x * s_um) for x in radii],
        visible_turn_deg=[float(np.ptp(np.degrees(np.unwrap(np.arctan2(Q[:, 1] - cv, Q[:, 0] - cu))))) for Q in Ps],
        each_wall_own_centre_v_px=dict(min=min(o['centre_v'] for o in own), max=max(o['centre_v'] for o in own),
                                       sd=float(np.std([o['centre_v'] for o in own]))),
        straight_top_leg_turned_deg=[w['turned_deg'] for w in straight],
        centre_from_leg1_axis_um=Yc, centre_beyond_last_pills_um=float((last_pill_u - cu) * s_um),
        legs_axes_um={'7': Y[7], '5': Y[5], '3': Y[3], '1': Y[1]},
        leg_Y_um=[round(v, 1) for v in leg_Y],
        note=('leg_Y in opencif\'s convention (leg 1 at 0, +Y towards leg 2): legs 1, 3, 5, 7 from the post rows; 2, 4, 6 '
              'as 2 Y_c - Y of the leg before them (a semicircular turn about the measured centre); the mask y axis along the '
              'image\'s +v (the non-mirrored map of panels_gifford6.py)'))
    json.dump(out, open(OUT, 'w'), indent=1)
    print(json.dumps({k: v for k, v in out.items() if k not in ('method', 'wall_radii_um')}, indent=1))
    if PNG:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        x0_, y0_, x1_, y1_ = P.CROP
        fig = plt.figure(figsize=(6.6, 5.2), dpi=220)
        ax = fig.add_axes([0, 0.08, 1, 0.92])
        ax.imshow(r.im[y0_:y1_, x0_:x1_], cmap='gray', vmin=0, vmax=255, extent=(x0_ - .5, x1_ - .5, y1_ - .5, y0_ - .5))
        for Q in Ps:
            xy = P.to_raster(Q[:, 0], Q[:, 1], th)
            ax.plot(xy[:, 0], xy[:, 1], '-', color='#00d5ff', lw=0.6)
        tt = np.linspace(0, 2 * np.pi, 720)
        for R_ in radii:
            ux, vy = cu + R_ * np.cos(tt), cv + R_ * np.sin(tt)
            xy = P.to_raster(ux, vy, th)
            ax.plot(xy[:, 0], xy[:, 1], '--', color='#ffb000', lw=0.4, alpha=0.8)
        ax.set_xlim(x0_ - .5, x1_ - .5)
        ax.set_ylim(y1_ - .5, y0_ - .5)
        ax.axis('off')
        fig.text(0.01, 0.01, 'Gifford et al. 2014 Fig. 6: the %d turn walls followed on the image (cyan) and the circles fitted about one '
                 'common centre (amber, dashed): residual SD %.2f px (%.2f µm); centre %.0f µm beyond leg 1\'s axis.'
                 % (len(Ps), np.std(fun), np.std(fun) * s_um, Yc), fontsize=5.5, wrap=True)
        fig.savefig(PNG)


if __name__ == '__main__':
    main()
