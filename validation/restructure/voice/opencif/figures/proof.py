"""opencif.figures.proof: the multipanel proof figures of a built device.

  <name>_layout.png   the whole mask (in strips for long devices) inside the chip outline, with the zoom windows marked
  <name>_zooms.png    zooms with the same style of scale bar, at the inlet, the head, the seam, the first rows, the
                      side channel at the first rows, the middle of the array, the loop, the last rows, and the outlet
  <name>_physics.png  the row-by-row proof from the drawn mask, with the widths, the gap fractions, the realized
                      cutoff against D, and F
"""
import math
import os

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from opencif.figures import drawing as RD  # noqa: E402
from opencif.checks import measure as M2  # noqa: E402
from opencif.geometry import layout as L  # noqa: E402

INK = '#1f1f1f'
DESIGN = '#8a8f98'
DRAWN_A = '#1f5fa8'
DRAWN_B = '#c2410c'
TARGET = '#127a43'
BAND = '#d7efe0'


def _scalebar(ax, frac=0.2, fs=7):
    """Draw a red scale bar whose length is 1, 2, or 5 times a power of ten and at most frac of the axis width."""
    x0, x1 = ax.get_xlim()
    span = (x1 - x0) * frac
    nice = 10 ** math.floor(math.log10(span))
    for m in (5, 2, 1):
        if m * nice <= span:
            L = m * nice
            break
    lab = ('%g mm' % (L / 1000)) if L >= 1000 else ('%g µm' % L)
    RD.scalebar(ax, L, lab, loc=(0.05, 0.06), fs=fs, color='#d62828', lw=2.5)


def _glob(info, k_local_x, y, leg):
    """Map a point in the leg's local coordinates to the device frame (measure.l2g)."""
    return M2.l2g(leg, info['gx'][leg['j']], k_local_x, y)


def _leg_of(info, k):
    """Return the leg that holds row k, or the last leg when none does."""
    for lg in info['legs']:
        if lg['rows'][0] <= k <= lg['rows'][1]:
            return lg
    return info['legs'][-1]


def _row_window(sc, info, k0, k1, pad_y=None):
    """Return the global window that covers rows k0 to k1 of the array, both side channels included."""
    leg = _leg_of(info, k0)
    off = leg['p0'] if leg['p0'] is not None else 0.0
    hw = float(np.nanmax(sc.wc[k0:k1 + 1])) / 2 + sc.spec.l_g + float(np.nanmax(np.nan_to_num(sc.ws[k0:k1 + 1], nan=sc.spec.GS))) + 25
    hw = pad_y or hw
    a = _glob(info, sc.x[k0] - off - sc.spec.Lam_L, -hw, leg)
    b = _glob(info, sc.x[k1] - off + sc.spec.Lam_L, hw, leg)
    return (min(a[0], b[0]), min(a[1], b[1]), max(a[0], b[0]), max(a[1], b[1]))


def _square(win, aspect=1.6):
    """Grow a window about its center until its width over its height equals aspect."""
    x0, y0, x1, y1 = win
    w, h = x1 - x0, y1 - y0
    cx, cy = 0.5 * (x0 + x1), 0.5 * (y0 + y1)
    if w / h < aspect:
        w = h * aspect
    else:
        h = w / aspect
    return (cx - w / 2, cy - h / 2, cx + w / 2, cy + h / 2)


def zoom_windows(sc, info, fluid=None):
    """Return the zoom windows of a device as (label, window) pairs in the device frame."""
    spec = sc.spec
    N, NM = sc.N, sc.NM
    wins = []
    inl = info['inlet']
    r_in = spec.port_r_in
    xp = inl.get('x_port_c', info['x_top'] - 500)
    wins.append(('inlet', _square((xp - 1.3 * r_in, -1.3 * r_in, max(info['x_top'], xp + 1.3 * r_in), 1.3 * r_in))))
    if NM > 0:
        from shapely.geometry import box as _box
        W = float(getattr(spec, 'W', 4000.0))
        hd = fluid.intersection(_box(info['x_top'] - 40, -W / 2 - 3000, sc.x[NM] + 40, W / 2 + 3000)) if fluid is not None else None
        hh = (max(abs(hd.bounds[1]), abs(hd.bounds[3])) if hd is not None and not hd.is_empty else W / 2) + 40
        wins.append(('meander head', _square((info['x_top'] - 40, -hh, sc.x[NM] + 40, hh))))
        kL = next((k for k in range(NM, N) if sc.kind[k] == 'L'), NM)
        wins.append(('head → array seam', _square(_row_window(sc, info, max(NM - 3, 0), min(kL + 5, N - 2)))))
    k0 = 0 if NM == 0 else NM
    wins.append(('first straight rows %d-%d' % (k0 + 1, k0 + 4), _square(_row_window(sc, info, k0, min(k0 + 3, N - 2)))))
    # The next window shows the side channel at the first straight rows, on side A only and magnified about 5x across
    # the channel, because these are the narrowest drawn channels.
    leg = _leg_of(info, k0)
    off = leg['p0'] if leg['p0'] is not None else 0.0
    yi = float(sc.wc[k0]) / 2 + spec.l_g
    try:                                                  # the drawn width (the papers' transition draws ws_rect)
        wsk = float(L.drawn_side_width(sc, k0))
    except Exception:
        wsk = float('nan')
    wsk = float(np.nan_to_num(wsk, nan=spec.GS))
    a = _glob(info, sc.x[k0] - off - spec.G, yi - 0.6 * spec.l_g, leg)
    b = _glob(info, sc.x[min(k0 + 2, N - 1)] - off, yi + wsk + 20, leg)
    wins.append(('side channel at row %d (drawn w_s %.2f µm)' % (k0 + 1, wsk),
                 _square((min(a[0], b[0]), min(a[1], b[1]), max(a[0], b[0]), max(a[1], b[1])))))
    km = (NM + N) // 2
    wins.append(('mid-array (row %d)' % (km + 1), _square(_row_window(sc, info, km, min(km + 6, N - 2)))))
    if info['loops']:
        lp = info['loops'][0]
        R = lp.get('R', sc.info.get('R_turn', 300.0))
        hw = float(sc.wc[lp['k']]) / 2 + spec.l_g + float(np.nan_to_num(sc.ws[lp['k']], nan=50)) + 40
        wins.append(('first U-turn', _square((lp['cx'] - R - hw, lp['cy'] - R - hw, lp['cx'] + R + hw, lp['cy'] + R + hw))))
    wins.append(('last rows (row %d)' % N, _square(_row_window(sc, info, N - 7, N - 2))))
    o = info['outlet']
    pts = np.array([o['ports'][k] for k in ('A', 'B', 'C')])
    rr = max(spec.port_r, spec.port_r_c)
    wins.append(('outlet ports', _square((pts[:, 0].min() - rr * 1.3, pts[:, 1].min() - rr * 1.3,
                                          pts[:, 0].max() + rr * 1.3, pts[:, 1].max() + rr * 1.3))))
    return wins


def layout_figure(fluid, sc, info, path, title, chip=None, wins=None):
    """Draw the whole mask with the chip outline and the lettered zoom windows, and save it to path. A mask six or
    more times longer than it is high is cut into up to six strips."""
    b = fluid.bounds
    W, H = b[2] - b[0], b[3] - b[1]
    n_strips = 1 if W / max(H, 1) < 6 else int(min(6, math.ceil(W / H / 5.0)))
    fig, axes = plt.subplots(n_strips, 1, figsize=(12, 1.2 + 2.6 * n_strips * (1 if n_strips > 1 else max(0.6, min(1.8, 4.5 * H / W)))))
    axes = np.atleast_1d(axes)
    seg = W / n_strips
    for i, ax in enumerate(axes):
        win = (b[0] + i * seg - 50, b[1] - 50, b[0] + (i + 1) * seg + 50, b[3] + 50)
        if chip and n_strips == 1:
            cw, ch, (cx0, cy0) = chip
            win = (min(win[0], cx0) - 200, min(win[1], cy0) - 200, max(win[2], cx0 + cw) + 200, max(win[3], cy0 + ch) + 200)
        RD.draw(ax, fluid, win, style='layout')
        if chip:
            cw, ch, (cx0, cy0) = chip
            ax.add_patch(Rectangle((cx0, cy0), cw, ch, fill=False, ec=TARGET, lw=1.2, ls='--'))
        for j, (lab, w) in enumerate(wins or []):
            if w[2] < win[0] or w[0] > win[2]:
                continue
            ax.add_patch(Rectangle((w[0], w[1]), w[2] - w[0], w[3] - w[1], fill=False, ec='#d62828', lw=0.9))
            ax.text(max(w[0], win[0]), min(w[3], win[3]), ' %s' % chr(65 + j), color='#d62828', fontsize=8, va='top',
                    fontweight='bold', clip_on=True)
        ax.set_xticks([]); ax.set_yticks([])
        _scalebar(ax, 0.12)
        if n_strips > 1:
            ax.set_ylabel('strip %d' % (i + 1), fontsize=8)
    import textwrap
    fig.suptitle('\n'.join(textwrap.wrap(title, 150)), fontsize=9, x=0.01, ha='left')
    fig.tight_layout()
    fig.savefig(path, dpi=170)
    plt.close(fig)


def zooms_figure(fluid, sc, info, path, title, wins):
    """Draw each zoom window in the micrograph style, four to a row and lettered as on the layout figure, and save the
    figure to path."""
    n = len(wins)
    cols = 4
    rows = int(math.ceil(n / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(18, 3.9 * rows))
    axes = np.atleast_1d(axes).ravel()
    for j, ax in enumerate(axes):
        if j >= n:
            ax.axis('off')
            continue
        lab, w = wins[j]
        RD.draw(ax, fluid, w, style='micro')
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title('%s  %s' % (chr(65 + j), lab), fontsize=9, loc='left')
        _scalebar(ax, 0.25)
    import textwrap
    fig.suptitle('\n'.join(textwrap.wrap(title, 170)) + '\n(fluid light grey, solid dark grey, every edge outlined)', fontsize=9,
                 x=0.01, ha='left')
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)


def physics_figure(arr, rep, path, title):
    """Draw the six panels of the row-by-row proof on the drawn mask and save them to path. The panels show the
    widths, the gap fractions, the realized cutoff, the filtrate fraction, and two close views of the realized
    cutoff."""
    sc_N = int(arr['N'])
    k = np.arange(1, sc_N + 1)
    D = rep['validation']['network']['cutoff_design_um']
    tol = rep.get('assessment', {}).get('cutoff_tolerance', 0.01)
    fig, ax = plt.subplots(3, 2, figsize=(14, 11.4))
    a = ax[0, 0]
    a.plot(k, arr['wc'], color=INK, lw=1.2, label='centre width $w_c$ (design)')
    ws = np.asarray(arr['ws'], float)
    a.plot(k, ws, color=DESIGN, lw=1.2, label='side width $w_s$ (closure)')
    lin = arr.get('meas_lin') or []
    if lin:
        kk = np.array([t[0] for t in lin if t[1] == 'A']) + 1
        we = np.array([t[2] for t in lin if t[1] == 'A'])
        a.plot(kk, we, '.', ms=1.6, color=DRAWN_A, label='side width, drawn (resistance-equivalent)')
    NM = int(arr['NM'])
    if NM:
        a.axvspan(1, NM, color='#f1f1f1', zorder=0, label='meander head')
    a.set_xlabel('row (gap) i'); a.set_ylabel('µm'); a.set_title('Widths', fontsize=9, loc='left')
    a.legend(fontsize=7, frameon=False)
    b = ax[0, 1]
    f = np.asarray(arr['f'], float)
    with np.errstate(divide='ignore', invalid='ignore'):
        rA = 100 * (np.asarray(arr['fA']) / f - 1)
        rB = 100 * (np.asarray(arr['fB']) / f - 1)
    b.axhline(0, color=DESIGN, lw=1.0)
    b.plot(k, rA, color=DRAWN_A, lw=0.8, label='side A')
    b.plot(k, rB, color=DRAWN_B, lw=0.8, ls='--', label='side B')
    lim = np.nanpercentile(np.abs(np.concatenate([rA[1:-1], rB[1:-1]])), 99.9) if sc_N > 3 else 1.0
    lim = max(lim * 1.3, 0.05)
    b.set_ylim(-lim, lim)
    nclip = int(np.sum((np.abs(rA) > lim) | (np.abs(rB) > lim)))
    b.set_xlabel('row (gap) i'); b.set_ylabel('drawn f / design f − 1  (%)')
    b.set_title('Gap fractions on the drawn mask (measured network)%s' % (
        ('; %d row(s) beyond the axis' % nclip) if nclip else ''), fontsize=9, loc='left')
    if NM:
        b.axvspan(1, NM, color='#f1f1f1', zorder=0)
    b.legend(fontsize=7, frameon=False)
    c = ax[1, 0]
    c.axhspan(D * (1 - tol), D * (1 + tol), color=BAND, zorder=0, label='D ± %g%%' % (100 * tol))
    g = rep['validation']['gates']
    band = g.get('cutoff_band_um', {})
    if band.get('passage_edge_below_knee'):
        c.axhline(band['passage_edge_below_knee'], color=TARGET, lw=0.8, ls=':', label='passage edge below the κ knee')
    c.axhline(D, color=TARGET, lw=1.2, label='target D = %.3f µm' % D)
    if arr.get('dA_int') is not None:
        c.plot(k, arr['dA_int'], color=DESIGN, lw=2.2, alpha=0.6, label='intended design (network, sealed)')
    c.plot(k, arr['dA'], color=DRAWN_A, lw=0.8, label='realised cutoff, side A')
    c.plot(k, arr['dB'], color=DRAWN_B, lw=0.8, ls='--', label='realised cutoff, side B')
    rc = rep['validation']['realised_cutoff']
    lo, hi = min(rc['D_realised_min_um'], D * (1 - 3 * tol)), max(rc['D_realised_max_um'], D * (1 + 3 * tol))
    c.set_ylim(lo - 0.05 * (hi - lo), hi + 0.05 * (hi - lo))
    tr = rc.get('trimmed_last_gap')
    if tr:
        c.annotate('last gap trimmed so F = target exactly\n(cutoff %.2f D, %.1e of F: retentive)' % (tr['cutoff_over_D'], tr['share_of_F']),
                   xy=(sc_N, max(lo, D * tr['cutoff_over_D'])), xytext=(0.62, 0.12), textcoords='axes fraction', fontsize=7,
                   arrowprops=dict(arrowstyle='->', lw=0.6, color=INK))
    if NM:
        c.axvspan(1, NM, color='#f1f1f1', zorder=0)
    c.set_xlabel('row (gap) i'); c.set_ylabel('µm')
    c.set_title('Realised cutoff $D_i = κ_g\\,x(f_i)\\,w_c(i)$ on the drawn mask: %.1f%% of F within ±%g%%' % (
        100 * rc['share_of_F_within_%g%%' % (100 * tol)], 100 * tol), fontsize=9, loc='left')
    c.legend(fontsize=7, frameon=False, ncol=2)
    d = ax[1, 1]
    P = np.asarray(arr['P'], float)
    d.plot(k, 1 - P, color=DESIGN, lw=1.2, label='design $1 - P(i)$')
    net = rep['validation']['network']
    d.axhline(net['F_target'], color=TARGET, lw=1.2, label='target F = %.5f' % net['F_target'])
    d.plot([sc_N], [net['F_as_drawn']], 'o', ms=6, mfc='none', mec=DRAWN_A, label='F on the drawn mask = %.6f' % net['F_as_drawn'])
    d.set_xlabel('row (gap) i'); d.set_ylabel('cumulative filtrate fraction')
    d.set_title('Filtrate fraction: decoupled from the cutoff (N set by F alone)', fontsize=9, loc='left')
    d.legend(fontsize=7, frameon=False, loc='lower right')
    import textwrap
    # Panels (e) and (f) show the entrance and the first U-turn (or the last rows) row by row, because departures
    # concentrate there.
    k0 = max(0, NM - 5)
    turns = [int(t) for t in (arr.get('turns') or [])]
    windows = [('entrance: rows %d-%d' % (k0 + 1, min(k0 + 80, sc_N)), k0, min(k0 + 80, sc_N))]
    if turns:
        t0 = turns[0]
        windows.append(('around the first U-turn (after row %d)' % (t0 + 1), max(0, t0 - 40), min(sc_N, t0 + 41)))
    else:
        windows.append(('last rows %d-%d' % (max(1, sc_N - 80), sc_N), max(0, sc_N - 80), sc_N))
    for axx, (lab, a0, a1) in zip((ax[2, 0], ax[2, 1]), windows):
        kk = k[a0:a1]
        axx.axhspan(D * (1 - tol), D * (1 + tol), color=BAND, zorder=0)
        axx.axhline(D, color=TARGET, lw=1.0)
        if arr.get('dA_int') is not None:
            axx.plot(kk, np.asarray(arr['dA_int'])[a0:a1], '-', lw=2.4, alpha=0.6, color=DESIGN, label='intended design')
        axx.plot(kk, np.asarray(arr['dA'])[a0:a1], 'o-', ms=2.2, lw=0.8, color=DRAWN_A, label='side A')
        axx.plot(kk, np.asarray(arr['dB'])[a0:a1], 's--', ms=2.0, lw=0.8, color=DRAWN_B, label='side B')
        if NM and a0 < NM:
            axx.axvspan(a0 + 1, NM, color='#f1f1f1', zorder=0, label='meander head')
        for t in turns:
            if a0 <= t < a1:
                axx.axvline(t + 1.5, color=INK, lw=0.6, ls=':')
        axx.set_xlabel('row (gap) i'); axx.set_ylabel('realised cutoff (µm)')
        axx.set_title('Realised cutoff, %s' % lab, fontsize=9, loc='left')
        axx.legend(fontsize=7, frameon=False)
    fig.suptitle('\n'.join(textwrap.wrap(title, 170)), fontsize=9, x=0.01, ha='left')
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def all_figures(fluid, sc, info, arr, rep, outdir, name, title, chip=None):
    """Draw the layout, zoom, and physics figures of a build into outdir and return their three paths."""
    wins = zoom_windows(sc, info, fluid)
    p1 = os.path.join(outdir, name + '_layout.png')
    p2 = os.path.join(outdir, name + '_zooms.png')
    p3 = os.path.join(outdir, name + '_physics.png')
    layout_figure(fluid, sc, info, p1, title, chip, wins)
    zooms_figure(fluid, sc, info, p2, title, wins)
    physics_figure(arr, rep, p3, title)
    return [p1, p2, p3]


def refused_layout(fluid, info, coll, outdir, name, verdict):
    """Draw the figure of a layout refused because its parts collide, for which no photomask is written (the
    clearance guarantee). The figure shows the outline made in memory, the two closest unrelated parts outlined in red,
    and a zoom where they meet or overlap, and the function returns its file name, <name>_refused_layout.png."""
    import textwrap
    from shapely.ops import nearest_points
    path = os.path.join(outdir, name + '_refused_layout.png')
    pc = coll.get('part_clearance') or {}
    parts = dict(info.get('parts') or [])
    pair = (pc.get('closest') or [{}])[0].get('parts') if not pc.get('ok', True) else None
    b = fluid.bounds
    W, H = b[2] - b[0], b[3] - b[1]
    h_top = max(1.5, min(6.0, 11.5 * (H + 200.0) / (W + 200.0)))
    if pair:
        fig, axes = plt.subplots(2, 1, figsize=(12, h_top + 6.4), gridspec_kw=dict(height_ratios=[h_top, 6.0]))
    else:
        fig, axes = plt.subplots(1, 1, figsize=(12, h_top + 1.0))
    axes = np.atleast_1d(axes)
    RD.draw(axes[0], fluid, (b[0] - 100, b[1] - 100, b[2] + 100, b[3] + 100), style='layout')
    zoom = None
    if pair and pair[0] in parts and pair[1] in parts:
        A, B = parts[pair[0]], parts[pair[1]]
        for g_, ax_ in ((A, axes[0]), (B, axes[0])):
            for poly in (g_.geoms if hasattr(g_, 'geoms') else [g_]):
                if hasattr(poly, 'exterior'):
                    x_, y_ = poly.exterior.xy
                    ax_.plot(x_, y_, color='#d62828', lw=0.8)
        ov = A.intersection(B)
        if not ov.is_empty:
            ob = ov.bounds
            cx, cy, r = 0.5 * (ob[0] + ob[2]), 0.5 * (ob[1] + ob[3]), max(400.0, 0.6 * max(ob[2] - ob[0], ob[3] - ob[1]))
        else:
            pa, pb = nearest_points(A, B)
            cx, cy, r = 0.5 * (pa.x + pb.x), 0.5 * (pa.y + pb.y), 400.0
        zoom = (cx - 1.6 * r, cy - r, cx + 1.6 * r, cy + r)
        axes[0].add_patch(Rectangle((zoom[0], zoom[1]), zoom[2] - zoom[0], zoom[3] - zoom[1], fill=False, ec='#d62828',
                                    lw=1.2))
        RD.draw(axes[1], fluid, zoom, style='layout')
        for g_ in (A, B):
            for poly in (g_.geoms if hasattr(g_, 'geoms') else [g_]):
                if hasattr(poly, 'exterior'):
                    x_, y_ = poly.exterior.xy
                    axes[1].plot(x_, y_, color='#d62828', lw=1.0)
        if not ov.is_empty:
            for poly in (ov.geoms if hasattr(ov, 'geoms') else [ov]):
                if hasattr(poly, 'exterior'):
                    axes[1].fill(*poly.exterior.xy, color='#d62828', alpha=0.45, lw=0)
        axes[1].set_xlim(zoom[0], zoom[2]); axes[1].set_ylim(zoom[1], zoom[3])
        axes[1].set_title("'%s' and '%s': %.1f um apart (floor %.1f um); red: their outlines and overlap"
                          % (pair[0], pair[1], pc['closest'][0]['gap_um'], pc['floor_um']), fontsize=9, loc='left')
    for ax in axes:
        ax.set_xticks([]); ax.set_yticks([])
        _scalebar(ax, 0.12)
    head_ = verdict.split('  Remedies:')[0].rstrip('.')   # the remedies are in the record
    fig.suptitle('\n'.join(textwrap.wrap('%s -- %s (no photomask written; the record names the remedies)' % (name, head_),
                                          150)), fontsize=8, x=0.01, ha='left')
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return os.path.basename(path)


def rerender(outdir, name):
    """Redraw the three figures of a finished build from its record and its pickle, without rebuilding the device."""
    import json
    import pickle
    from types import SimpleNamespace
    rep = json.load(open(os.path.join(outdir, name + '.json')))
    pk = pickle.load(open(os.path.join(outdir, name + '.pkl'), 'rb'))
    arr, fluid = pk['arrays'], pk['fluid']
    sc = SimpleNamespace(x=np.asarray(arr['x']), wc=np.asarray(arr['wc']), ws=np.asarray(arr['ws']),
                         kind=np.asarray(arr['kind']), N=int(arr['N']), NM=int(arr['NM']), spec=arr['spec'],
                         info=dict(turns=list(arr['turns']), R_turn=rep['design'].get('R_turn') or 300.0))
    info = dict(legs=arr['legs'], gx=arr['gx'], loops=arr['loops'], inlet=arr['inlet'], outlet=arr['outlet'],
                x_top=arr['x_top'])
    chip = None
    ch = rep.get('input', {}).get('chip') or {}
    if ch.get('width_mm') and ch.get('length_mm'):
        b = fluid.bounds
        cw, chh = ch['length_mm'] * 1000.0, ch['width_mm'] * 1000.0
        chip = (cw, chh, (0.5 * (b[0] + b[2]) - cw / 2, 0.5 * (b[1] + b[3]) - chh / 2))
    title = rep.get('title') or name
    return all_figures(fluid, sc, info, arr, rep, outdir, name, title, chip)
