"""dissertation_figures.py: the main-text figures of the dissertation, drawn from a finished build.

The script lies outside the code stamp, because it draws what a build already holds and computes nothing that a record
does not. The per-build QA figures (`<name>_layout.png`, `_zooms.png`, and `_physics.png`, drawn by
`opencif/figures/proof.py`, which is stamped) stay the validation figures, and these are their main-text counterparts
(adjudication of the rc7 reviews, section 10). Each figure has one quantity per axis, two validated categorical colors
(dataviz `validate_palette.js`, light surface: #2a78d6 / #eb6834), text in ink, the Q2 band and the declared rows
marked, a vector PDF for the linework, and a PNG at 300 dpi.

  python3 tools/dissertation_figures.py build RUN_DIR [OUT_DIR]
      RUN_DIR is runs/<name>/ of a finished build, with <name>.json, <name>.pkl, and <name>.dxf (or <name>_FAILED.dxf).
      <name>_overview.{pdf,png}   the photomask as written (the DXF, with CHANNELS filled and ISLANDS cut out), a scale
                                  bar, the inlet and outlet ports labeled, and four zoom windows (A-D) marked and drawn
                                  beneath, each with its own scale bar
      <name>_physics.{pdf,png}    (a) the realized cutoff on the drawn mask, D_i / D - 1 (%), per gap row for both side
                                  channels, with the +-tolerance band, the declared rows hatched and named, and the rows
                                  off the axis counted; (b) the cumulative filtrate fraction 1 - P(i), the target F, and
                                  F on the drawn mask; (c) the 1-D model's scenarios (`validation.model_scenarios`),
                                  each scenario's largest change of any row's cutoff against the tolerance, with the
                                  rows that it moves beyond
      <name>_provenance.md        the input's hash, the code stamp, the tests' hash, the photomask's sha256 (recomputed
                                  and compared with the record's), D, F, N, NM, the floor, the statuses, and the verdict
  python3 tools/dissertation_figures.py v13 OUT_STEM V13_JSON[=TITLE] [V13_JSON[=TITLE] ...]
      the taper-direction sign test (V13, `provenance/rc7_traces/code/taper_sign.py`), with, for each reading, the
      outer wall-to-wall width along the leg on the image and on the mask, each with its least-squares line, the
      verdict, the slopes with the image's standard error, and the coverage. No paper raster is drawn (widths only), so
      the figure may be public.

The pickle is the build's own (written by `opencif/pipeline.py`), and only a pickle that one has built oneself should be
loaded."""
import hashlib
import json
import math
import os
import pickle
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt                       # noqa: E402
from matplotlib.patches import Rectangle              # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import opencif                                        # noqa: E402,F401  (the package of the pickles' classes)

INK, INK2, GRID, SURF = '#0b0b0b', '#52514e', '#e4e3df', '#fcfcfb'
C1, C2 = '#2a78d6', '#eb6834'                         # categorical slots 1 and 2 (validated, light)
BAND = '#efeeea'
FLUID = '#2a78d6'
MARK = '#0b0b0b'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8.5, 'axes.edgecolor': GRID, 'axes.linewidth': 0.8,
                     'xtick.color': INK2, 'ytick.color': INK2, 'axes.labelcolor': INK, 'text.color': INK,
                     'xtick.major.size': 0, 'ytick.major.size': 0, 'figure.facecolor': SURF, 'axes.facecolor': SURF,
                     'pdf.fonttype': 42})


def sha256(path):
    """Return the sha256 of a file, read in blocks of 1 MiB."""
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for b in iter(lambda: fh.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def _clean(ax):
    for s_ in ('top', 'right'):
        ax.spines[s_].set_visible(False)
    ax.grid(True, axis='y', color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)


# ---------------------------------------------------------------------------------------------------------- the mask
def read_mask(dxf):
    """Return the photomask as written, its CHANNELS outlines and its ISLANDS (posts), as point arrays."""
    import ezdxf
    doc = ezdxf.readfile(dxf)
    ch, isl = [], []
    for e in doc.modelspace():
        if e.dxftype() != 'LWPOLYLINE':
            continue
        p = np.array([(float(q[0]), float(q[1])) for q in e.get_points()])
        (ch if e.dxf.layer == 'CHANNELS' else isl if e.dxf.layer == 'ISLANDS' else []).append(p)
    return ch, isl


def draw_mask(ax, ch, isl, win=None, lw=0.0):
    from matplotlib.collections import PolyCollection
    ax.add_collection(PolyCollection(ch, facecolors=FLUID, edgecolors=FLUID if lw else 'none', linewidths=lw, zorder=1))
    if win is not None:
        x0, y0, x1, y1 = win
        isl = [p for p in isl if p[:, 0].max() >= x0 and p[:, 0].min() <= x1 and p[:, 1].max() >= y0 and p[:, 1].min() <= y1]
    ax.add_collection(PolyCollection(isl, facecolors=SURF, edgecolors='none', zorder=2))
    ax.set_aspect('equal')
    ax.set_xticks([]); ax.set_yticks([])
    for s_ in ax.spines.values():
        s_.set_visible(False)


def scalebar(ax, frac=0.2, loc=(0.04, 0.06), fs=7.5, color=INK):
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    span = (x1 - x0) * frac
    nice = 10 ** math.floor(math.log10(span))
    L = next(m * nice for m in (5, 2, 1) if m * nice <= span)
    lab = ('%g mm' % (L / 1000)) if L >= 1000 else ('%g µm' % L)
    xa, ya = x0 + loc[0] * (x1 - x0), y0 + loc[1] * (y1 - y0)
    ax.plot([xa, xa + L], [ya, ya], color=color, lw=2.2, solid_capstyle='butt', zorder=5)
    ax.text(xa + L / 2, ya + 0.025 * (y1 - y0), lab, ha='center', va='bottom', fontsize=fs, color=color, zorder=5,
            bbox=dict(boxstyle='round,pad=0.12', fc=SURF, ec='none', alpha=0.85))


def _windows(rec, arr, fluid):
    """Pick four zooms from the stamped renderer's own windows (opencif.figures.proof.zoom_windows), which are the inlet
    or the head, the side channel of the first straight rows (the narrowest drawn channel), the first U-turn or the
    middle of the array, and the outlet ports."""
    from types import SimpleNamespace
    from opencif.figures import proof as RR
    sc = SimpleNamespace(x=np.asarray(arr['x']), wc=np.asarray(arr['wc']), ws=np.asarray(arr['ws']),
                         kind=np.asarray(arr['kind']), N=int(arr['N']), NM=int(arr['NM']), spec=arr['spec'],
                         info=dict(turns=list(arr['turns']), R_turn=rec['design'].get('R_turn') or 300.0))
    info = dict(legs=arr['legs'], gx=arr['gx'], loops=arr['loops'], inlet=arr['inlet'], outlet=arr['outlet'],
                x_top=arr['x_top'])
    w = RR.zoom_windows(sc, info, fluid)
    pick = []
    for keys in (('meander head', 'inlet'), ('side channel at row',), ('first U-turn', 'mid-array'), ('outlet ports',)):
        for k in keys:
            hit = [x for x in w if x[0].startswith(k)]
            if hit:
                pick.append(hit[0])
                break
    return pick, info


def overview(run_dir, out_dir):
    import textwrap
    name = os.path.basename(os.path.normpath(run_dir))
    rec = json.load(open(os.path.join(run_dir, name + '.json')))
    pk = pickle.load(open(os.path.join(run_dir, name + '.pkl'), 'rb'))
    arr, fluid = pk['arrays'], pk['fluid']
    dxf = next(os.path.join(run_dir, f) for f in (name + '.dxf', name + '_FAILED.dxf')
               if os.path.exists(os.path.join(run_dir, f)))
    ch, isl = read_mask(dxf)
    wins, info = _windows(rec, arr, fluid)
    allp = np.vstack(ch)
    bx0, by0 = allp.min(axis=0)
    bx1, by1 = allp.max(axis=0)
    W, H = bx1 - bx0, by1 - by0
    pad = 0.015 * W
    xl, xr = bx0 - pad, bx1 + pad
    yb, yt = by0 - pad - 0.10 * H, by1 + pad
    FW = 12.0                                             # figure width (in); layout in inches, then fractions
    top_w = FW - 0.3
    top_h = min(max(top_w * (yt - yb) / (xr - xl), 0.9), 5.5)
    z_w = (FW - 0.3 - 3 * 0.25) / 4
    z_h = z_w / 1.6                                       # the renderer's zoom windows have aspect 1.6
    t_h, gap, zt_h, foot = 0.45, 0.30, 0.42, 0.35
    FH = t_h + top_h + gap + zt_h + z_h + foot
    fig = plt.figure(figsize=(FW, FH))
    ax = fig.add_axes([0.15 / FW, (foot + z_h + zt_h + gap) / FH, top_w / FW, top_h / FH])
    draw_mask(ax, ch, isl)
    ax.set_xlim(xl, xr); ax.set_ylim(yb, yt)
    ax.set_anchor('W')
    boxes = []
    for j, (lab, w) in enumerate(wins):          # a window taller than the overview is drawn to the overview's edges
        y0_, y1_ = max(w[1], yb), min(w[3], yt)
        boxes.append(((w[0] - xl) / (xr - xl), (w[2] - xl) / (xr - xl), (y0_ - yb) / (yt - yb), (y1_ - yb) / (yt - yb)))
        ax.add_patch(Rectangle((w[0], y0_), w[2] - w[0], y1_ - y0_, fill=False, ec=C2, lw=1.0, zorder=5))
        up = (j % 2 == 0) and w[3] <= yt              # the letter above an unclamped box's top edge, or inside a clamped one
        down = (j % 2 == 1) and w[1] >= yb
        ax.annotate('ABCD'[j], (w[2], y1_ if j % 2 == 0 else y0_), xytext=(2, 1 if up else -1 if (down or j % 2 == 0) else 1),
                    textcoords='offset points', color=C2, fontsize=9, fontweight='bold', zorder=6, ha='left',
                    va='bottom' if (up or (j % 2 == 1 and not down)) else 'top')
    for lx in (0.02, 0.30, 0.55, 0.80):          # the scale bar goes where no zoom box crosses it
        if not any(bx0_ < lx + 0.15 and bx1_ > lx - 0.01 and by0_ < 0.2 for bx0_, bx1_, by0_, _ in boxes):
            break
    scalebar(ax, 0.12, loc=(lx, 0.05))
    d = rec['design']
    rc = rec['validation']['realised_cutoff']
    fig.text(0.15 / FW, 1 - 0.12 / FH, '%s: the photomask' % name, fontsize=10, fontweight='bold', va='top')
    fig.text(0.15 / FW, 1 - 0.32 / FH, 'N %d gap rows (%d in the head) · D %.3f µm · F %.5f on the drawn mask · '
             'outline %.1f × %.2f mm' % (d['N'], d['NM'], rc['D_target_um'], rec['validation']['network']['F_as_drawn'],
                                         W / 1000, H / 1000), fontsize=8.5, va='top', color=INK2)
    ports = info['outlet'].get('ports') or {}
    names = {'A': 'side A', 'B': 'side B', 'C': 'centre'}
    for j, (lab, w) in enumerate(wins):
        a = fig.add_axes([(0.15 + j * (z_w + 0.25)) / FW, foot / FH, z_w / FW, z_h / FH])
        draw_mask(a, ch, isl, win=w)
        a.set_xlim(w[0], w[2]); a.set_ylim(w[1], w[3])
        fig.text((0.15 + j * (z_w + 0.25)) / FW, (foot + z_h + 0.04) / FH,
                 '\n'.join(textwrap.wrap('%s  %s' % ('ABCD'[j], lab), 52)), fontsize=8, va='bottom', color=INK)
        scalebar(a, 0.25)
        for s_ in a.spines.values():
            s_.set_visible(True); s_.set_edgecolor(C2); s_.set_linewidth(1.0)
        if lab.startswith('outlet ports'):
            seen = {}                                     # a combined (ring) outlet has A and B at one port
            for k_, p in ports.items():
                key = (round(float(p[0]), 3), round(float(p[1]), 3))
                seen.setdefault(key, []).append(k_)
            for key, ks in seen.items():
                p = key
                rr = float(arr['spec'].port_r_c if ks == ['C'] else arr['spec'].port_r)
                lab_ = names['C'] if ks == ['C'] else ('sides A + B' if len(ks) > 1 else names.get(ks[0], ks[0]))
                a.annotate(lab_, (p[0] - rr, p[1]), xytext=(-4, 0), textcoords='offset points',
                           ha='right', va='center', fontsize=7.5, color=INK,
                           bbox=dict(boxstyle='round,pad=0.15', fc=SURF, ec='none', alpha=0.85))
        if lab.startswith('inlet') and info['inlet'].get('x_port_c') is not None:
            xi = info['inlet']['x_port_c']
            ri = float(info['inlet'].get('port_r') or 0)
            a.annotate('inlet port', (xi, -ri), xytext=(0, -3), textcoords='offset points', ha='center', va='top',
                       fontsize=7.5, color=INK, bbox=dict(boxstyle='round,pad=0.15', fc=SURF, ec='none', alpha=0.85))
    fig.text(0.15 / FW, 0.08 / FH, 'Fluid blue, solid (posts and the chip) white; drawn from the photomask %s (sha256 '
             '%s…); zooms A–D marked on the device.' % (os.path.basename(dxf), sha256(dxf)[:16]), fontsize=7,
             color=INK2)
    stem = os.path.join(out_dir, name + '_overview')
    fig.savefig(stem + '.pdf')
    fig.savefig(stem + '.png', dpi=300)
    plt.close(fig)
    return stem


# ------------------------------------------------------------------------------------------------------- the physics
def physics(run_dir, out_dir):
    import textwrap
    name = os.path.basename(os.path.normpath(run_dir))
    rec = json.load(open(os.path.join(run_dir, name + '.json')))
    arr = pickle.load(open(os.path.join(run_dir, name + '.pkl'), 'rb'))['arrays']
    N, NM = int(arr['N']), int(arr['NM'])
    k = np.arange(1, N + 1)
    rc = rec['validation']['realised_cutoff']
    net = rec['validation']['network']
    D, tol = float(rc['D_target_um']), float(rc.get('tolerance', 0.01))
    rA = 100 * (np.asarray(arr['dA'], float) / D - 1)
    rB = 100 * (np.asarray(arr['dB'], float) / D - 1)
    dec = rc.get('declared_by_option') or {}
    rows_of = {k_: sorted(int(i) + 1 for i in (arr.get('declared') or {}).get(k_, []))   # the pickle's rows are 0-based
               or list(range(v_['first_last'][0], v_['first_last'][1] + 1)) for k_, v_ in dec.items()}
    undeclared = np.ones(N, bool)
    for rr in rows_of.values():
        undeclared[np.asarray(rr, int) - 1] = False
    worst = float(np.nanmax(np.abs(np.concatenate([rA[undeclared], rB[undeclared]])))) if undeclared.any() else 0.0
    YL = max(1.5 * 100 * tol, min(1.15 * worst, 5 * 100 * tol))   # the band in view, undeclared rows up to 5 x tol
    fig = plt.figure(figsize=(12.0, 7.8))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.3, 1], width_ratios=[1, 1.15], hspace=0.45, wspace=0.42,
                          top=0.86, bottom=0.08, left=0.07, right=0.98)
    a = fig.add_subplot(gs[0, :])
    _clean(a)
    a.axhspan(-100 * tol, 100 * tol, color=BAND, zorder=0, linewidth=0)
    a.text(N * 0.995, -100 * tol, '±%g %% (Q2)' % (100 * tol), fontsize=7.5, color=INK2, va='top', ha='right')
    if NM:
        a.axvspan(0.5, NM + 0.5, color=GRID, alpha=0.45, zorder=0, linewidth=0)
        a.text(N * 0.004, -YL * 0.93, 'meander head (design.NM): rows 1–%d' % NM, fontsize=7.5, color=INK2, ha='left', va='bottom', zorder=6,
               bbox=dict(boxstyle='round,pad=0.2', fc=SURF, ec='none', alpha=0.9))
    lines = []
    for opt, v in dec.items():
        rr = rows_of[opt]
        runs, st_ = [], rr[0]
        for i_ in range(1, len(rr) + 1):                 # contiguous runs of the declared rows, each hatched
            if i_ == len(rr) or rr[i_] != rr[i_ - 1] + 1:
                runs.append((st_, rr[i_ - 1]))
                if i_ < len(rr):
                    st_ = rr[i_]
        for r0, r1 in runs:
            a.axvspan(r0 - 0.5, r1 + 0.5, facecolor='none', edgecolor='#b5b3ad', hatch='//', linewidth=0, zorder=1)
        r0, r1 = v['first_last']
        lines.append((('row %d: %s (up to %.2f D)' % (r0, opt, v['worst_over'])) if r0 == r1 else
                      ('rows %d–%d: %s (up to %.2f D)' % (r0, r1, opt, v['worst_over']))) if len(rr) == r1 - r0 + 1 else
                     ('%d rows in %d–%d: %s (up to %.2f D)' % (len(rr), r0, r1, opt, v['worst_over'])))
    if lines:
        a.text(0.995, 0.975, 'declared rows (hatched; off the axis where beyond it)\n' + '\n'.join(lines),
               transform=a.transAxes, ha='right', va='top', fontsize=7.5, color=INK2, zorder=6,
               bbox=dict(boxstyle='round,pad=0.3', fc=SURF, ec=GRID, alpha=0.95))
    a.plot(k, np.clip(rA, -YL, YL), color=C1, linewidth=1.0, label='side channel A', zorder=3)
    a.plot(k, np.clip(rB, -YL, YL), color=C2, linewidth=0.8, label='side channel B', zorder=3)
    off = int(np.sum((np.abs(rA) > YL) | (np.abs(rB) > YL)))
    a.set_xlim(0, N + 1); a.set_ylim(-YL, YL)
    a.set_xlabel('gap row i'); a.set_ylabel('realised cutoff, D_i / D − 1 (%)')
    a.set_title('(a) The realised cutoff at every gap of the drawn mask (D = %.3f µm): %.1f %% of F within ±%g %%%s'
                % (D, 100 * rc['share_of_F_within_%g%%' % (100 * tol)], 100 * tol,
                   ('; %d row(s) beyond the axis, clipped' % off) if off else ''), fontsize=9, loc='left')
    a.legend(loc='lower right', frameon=False, ncol=2, fontsize=8)
    b = fig.add_subplot(gs[1, 0])
    _clean(b)
    P = np.asarray(arr['P'], float)
    S0 = float(rec['design']['spec'].get('S0') or 0.0)
    # the network's cumulative filtrate fraction on the drawn mask
    real = 1 - np.cumprod(1 - np.asarray(arr['fA'], float) - np.asarray(arr['fB'], float))
    if S0 == 0.0:                                          # for a sealed start, the design march is the device's own
        b.plot(k, 1 - P, color=C1, linewidth=2.4, alpha=0.45, label='design 1 − P(i)')
    b.plot(k, real, color=C1, linewidth=1.2, label='on the drawn mask, 1 − Π(1 − f_A − f_B)')
    b.axhline(net['F_target'], color=INK2, lw=0.8, ls='--', label='target F = %.5f' % net['F_target'])
    b.plot([N], [net['F_as_drawn']], 'o', ms=7, mfc='none', mec=C2, mew=1.4,
           label='F on the drawn mask = %.6f' % net['F_as_drawn'])
    note_b = ('\nfloor start: the design march begins with a fictitious side flow (2σ0 = %.4f of the inlet, App. A);\n'
              'the drawn mask has none, so the curve is the network\'s on the drawn mask' % (2 * S0)) if S0 > 0.0 else ''
    b.set_xlim(0, N * 1.03); b.set_ylim(0, 1.0)
    b.set_xlabel('gap row i' + note_b, fontsize=8.5 if not note_b else 7.5); b.set_ylabel('cumulative filtrate fraction')
    b.set_title('(b) The cumulative filtrate fraction (N set by F alone)', fontsize=9, loc='left')
    b.legend(loc='lower right', frameon=False, fontsize=7.5)
    c = fig.add_subplot(gs[1, 1])
    ms = rec['validation'].get('model_scenarios') or {}
    rows = [(kk, v) for kk, v in ms.items() if isinstance(v, dict) and 'max_cutoff_change' in v]
    FLOOR = 1e-3                                          # % of D, the log axis' floor (smaller changes drawn at it)
    vals = [max(100 * float(v['max_cutoff_change']), FLOOR) for _, v in rows]
    y = np.arange(len(rows))[::-1]
    c.grid(True, axis='x', color=GRID, linewidth=0.6); c.set_axisbelow(True)
    for s_ in ('top', 'right'):
        c.spines[s_].set_visible(False)
    c.set_xscale('log')
    c.barh(y, [v_ - FLOOR * 0.5 for v_ in vals], left=FLOOR * 0.5, height=0.55, color=C1, zorder=2)
    c.axvline(100 * tol, color=INK2, lw=0.8, ls='--', zorder=3)
    c.text(100 * tol, -0.9, '±%g %%' % (100 * tol), fontsize=7.5, color=INK2, ha='center', va='top')
    for yy, (kk, v), val in zip(y, rows, vals):
        n_ = int(v.get('rows_newly_outside_tol', 0))
        raw = 100 * float(v['max_cutoff_change'])
        txt = '%s%s' % ('%.3g %%' % raw if raw >= FLOOR else '< %g %%' % FLOOR,
                        ('; %d row(s) newly beyond, %.2g %% of F' % (
                            n_, 100 * float(v.get('share_of_F_in_those_rows', 0)))) if n_ else '')
        if val >= 5.0 or (n_ and val >= 1.0):             # a long bar carries its text inside, at its end
            c.text(val / 1.15, yy, txt, fontsize=7, color='white', va='center', ha='right', zorder=4)
        else:
            c.text(val * 1.15, yy, txt, fontsize=7, color=INK2, va='center', zorder=4)
    c.set_yticks(y); c.set_yticklabels([kk.replace('_', ' ') for kk, _ in rows], fontsize=7.5)
    c.set_xlim(FLOOR * 0.5, max(max(vals + [100 * tol]) * 2.0, 100 * tol * 3000))
    c.set_ylim(-1.3, len(rows) - 0.4)
    c.set_xlabel('largest change of any row\'s cutoff under the scenario (% of D, log scale)')
    c.set_title('(c) The 1-D model\'s scenarios (record: validation.model_scenarios)', fontsize=9, loc='left')
    fig.text(0.01, 0.985, '%s: the physics on the drawn mask' % name, fontsize=10, va='top', fontweight='bold')
    fig.text(0.01, 0.955, '\n'.join(textwrap.wrap('verdict (exit %s): %s' % (rec['assessment'].get('exit_code'),
                                                                            rec['assessment']['verdict']), 175)[:3]),
             fontsize=7.5, va='top', color=INK2)
    stem = os.path.join(out_dir, name + '_physics')
    fig.savefig(stem + '.pdf')
    fig.savefig(stem + '.png', dpi=300)
    plt.close(fig)
    return stem


# ---------------------------------------------------------------------------------------------------- the provenance
def provenance(run_dir, out_dir):
    name = os.path.basename(os.path.normpath(run_dir))
    rec = json.load(open(os.path.join(run_dir, name + '.json')))
    dxf = next(os.path.join(run_dir, f) for f in (name + '.dxf', name + '_FAILED.dxf')
               if os.path.exists(os.path.join(run_dir, f)))
    h = sha256(dxf)
    rec_h = (rec.get('files') or {}).get('dxf_sha256')
    inp = hashlib.sha256(json.dumps(rec.get('input', {}), sort_keys=True).encode()).hexdigest()
    a, d, v = rec['assessment'], rec['design'], rec['validation']
    rc, net, pr = v['realised_cutoff'], v['network'], (rec.get('input') or {}).get('process') or {}
    st = a.get('status_codes') or a.get('status') or {}
    L = ['# %s: provenance' % name, '',
         '| field | value |', '|---|---|',
         '| input (the record\'s `input` block, canonical JSON) sha256 | `%s` |' % inp,
         '| code | opencif %s, stamp `%s` (%s files) |' % (rec['code'].get('version'), rec['code'].get('sha256'),
                                                          rec['code'].get('files')),
         '| tests | `%s` (%s files) |' % ((rec.get('environment') or {}).get('tests_sha256'),
                                         (rec.get('environment') or {}).get('tests_files')),
         '| photomask | `%s`, sha256 `%s` (%s the record\'s) |' % (os.path.basename(dxf), h,
                                                                   'equal to' if h == rec_h else 'NOT EQUAL TO'),
         '| mode | %s |' % a.get('mode'),
         '| target cutoff D | %.4f µm (realised on the drawn mask %.4f–%.4f µm) |' % (
             rc['D_target_um'], rc['D_realised_min_um'], rc['D_realised_max_um']),
         '| filtrate fraction F | target %.6f, on the drawn mask %.6f |' % (net['F_target'], net['F_as_drawn']),
         '| share of F within ±%g %% of D | %.4f |' % (100 * rc.get('tolerance', 0.01),
                                                    rc['share_of_F_within_%g%%' % (100 * rc.get('tolerance', 0.01))]),
         '| declared rows | %s |' % ('; '.join(('%s: %d row(s) in rows %d–%d' % (k_, x[2], x[0], x[1])) if x[2] != x[1] - x[0] + 1
                                             else ('%s: row %d' % (k_, x[0])) if x[0] == x[1]
                                             else ('%s: rows %d–%d' % (k_, x[0], x[1]))
                                             for k_, x in (rc.get('declared_rows') or {}).items()) or 'none'),
         '| gap rows N (head NM) | %d (%d) |' % (d['N'], d['NM']),
         '| layout | %s, %s leg(s) |' % (d.get('layout'), d.get('n_legs')),
         '| the process floor | min feature %s µm, depth %s µm, gap %s µm |' % (pr.get('min_feature_um'),
                                                                               pr.get('depth_um'), pr.get('gap_um')),
         '| narrowest side channel | planned %.3f µm (`design.ws_min_um`); drawn %s (`validation.passage`, straight side channel) |' % (
             float(d.get('ws_min_um') or float('nan')),
             ('%.3f µm' % v['passage']['classes']['straight side channel']['width_um'])
             if ((v.get('passage') or {}).get('classes') or {}).get('straight side channel') else 'not measured'),
         '| statuses | %s |' % ', '.join('%s %s' % (k_, x) for k_, x in st.items()),
         '| verdict (exit) | %s (%s) |' % (a.get('verdict'), a.get('exit_code')), '']
    p = os.path.join(out_dir, name + '_provenance.md')
    open(p, 'w', encoding='utf-8').write('\n'.join(L))
    return p, h == rec_h


# --------------------------------------------------------------------------------------------------------------- V13
def v13(out_stem, items):
    """Draw one panel per reading, with the image's and the mask's outer wall-to-wall width along the leg, each about
    its own mean, and their least-squares lines. The image's channel can be wider than the mask's by a constant, which
    the test does not judge, so the panel states it."""
    n = len(items)
    fig, axes = plt.subplots(1, n, figsize=(4.9 * n, 4.3), squeeze=False)
    for ax, (path, title) in zip(axes[0], items):
        r = json.load(open(path))
        _clean(ax)
        rows = r.get('rows') or []
        if not rows:
            ax.text(0.5, 0.5, 'no positions found\n(%s)' % r.get('verdict'), ha='center', va='center',
                    transform=ax.transAxes, color=INK2)
            ax.set_title(title, fontsize=8.5, loc='left')
            continue
        X = np.array([q['X'] for q in rows]) / 1000.0
        wi = np.array([q['width_image'] for q in rows])
        wm = np.array([q['width_mask'] for q in rows])
        Xm = X.mean()
        ax.plot(X, wi - wi.mean(), 'o', ms=2.6, color=C1, label='image (walls found)', zorder=3)
        ax.plot(X, wm - wm.mean(), '-', color=C2, lw=1.4, label='mask', zorder=2)
        for sl, col in ((r['slope_image'], C1), (r['slope_mask'], C2)):
            ax.plot(X, sl * 1000.0 * (X - Xm), '--', color=col, lw=0.9, zorder=4)
        ax.set_xlabel('X along the leg (mm)')
        ax.set_ylabel('outer wall-to-wall width − its mean (µm)')
        ax.set_title('%s\n%s: image %+.5f ± %.5f µm/µm (%.1f SE), mask %+.5f' % (
            title, r['verdict'], r['slope_image'], r['se_image'], abs(r['z_image']), r['slope_mask']), fontsize=8,
            loc='left')
        ax.text(0.02, 0.03, 'window ±%g µm; %d of %d eligible positions (coverage %.2f)\nmean widths: image %.1f µm, '
                'mask %.1f µm' % (r['window_um'], r['positions'], r['eligible'], r['coverage'] or 0, wi.mean(), wm.mean()),
                transform=ax.transAxes, fontsize=7, color=INK2, va='bottom',
                bbox=dict(boxstyle='round,pad=0.2', fc=SURF, ec='none', alpha=0.85))
    h, l = axes[0][0].get_legend_handles_labels()
    fig.legend(h, l, loc='upper center', ncol=2, frameon=False, fontsize=8.5, bbox_to_anchor=(0.5, 0.0))
    fig.text(0.01, 1.0, 'V13, the taper-direction sign test: the sign of the image\'s width slope against the mask\'s '
             '(SE: the larger of Newey-West, lag 100 µm, and the independent-sample one; the magnitude is reported, not judged)', fontsize=9.5, fontweight='bold', va='bottom')
    fig.tight_layout()
    fig.savefig(out_stem + '.pdf', bbox_inches='tight')
    fig.savefig(out_stem + '.png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    return out_stem


def main(argv):
    if len(argv) >= 2 and argv[0] == 'build':
        run_dir = argv[1]
        out = argv[2] if len(argv) > 2 else run_dir
        os.makedirs(out, exist_ok=True)
        s1 = overview(run_dir, out)
        s2 = physics(run_dir, out)
        p, same = provenance(run_dir, out)
        print('wrote %s.{pdf,png}, %s.{pdf,png}, %s (photomask sha256 %s the record\'s)'
              % (s1, s2, p, 'equal to' if same else 'NOT EQUAL TO'))
        return 0 if same else 1
    if len(argv) >= 3 and argv[0] == 'v13':
        items = []
        for a in argv[2:]:
            path, _, title = a.partition('=')
            items.append((path, title or os.path.basename(path)))
        print('wrote %s.{pdf,png}' % v13(argv[1], items))
        return 0
    print(__doc__)
    return 2


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
