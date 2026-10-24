"""rc6: the figure-derived values re-read with every scale bar measured with its anti-aliased ends, and Xia's Fig. 3
panels b-g at the scale their printed widths and loop curvature give.  Writes the measurements and the values the
engine reads (opencif/data/r6/rc6_values.json, hashed).

Re-implemented on 2026-10-17 after the first rc6 working tree was lost in a container restart; the methods are those
recorded in the rc6 traces report (provenance/rc6_traces) and the rc6 reconstruction record, and the outputs are
checked against the values the lost run recorded (printed at the end; every rounded value the engine reads must
agree).

Sources: the publisher rasters (sha256 checked); the rc6 traces' registered per-column pill fits on Dinh's Fig. 2B
panels 2-4 (results/dinh/*.json, measured.pitch_um_by_column: image-only measurements, the same for any mask); the
rc6 traces' Xia arc fits (measurements/xia_fig3d_arc_scale.json); rc5's Dinh flare onset readings
(opencif/data/r5/figure_anchors.json, read at the 46 px bar); r4's Fig. 3a layout (opencif/data/r4/meas_xia_fig3a.json,
read at 179 px).

usage: python3 -I measure_rc6.py OPENCIF_ROOT RASTERS OUT.json VALUES.json"""
import json
import os
import sys

import numpy as np

R, RAS, OUT, VALS = sys.argv[1:5]
CODE = os.path.join(R, 'provenance', 'rc6_traces', 'code')
sys.path.insert(0, CODE)
import trace_lib as T                     # noqa: E402
import panels_xia as PX                   # noqa: E402
import panels_gifford as PG               # noqa: E402

SHA = {'dinh5-000.png': 'a21e2e3ce0710e32', 'gif5-001.png': 'e70a8b38b87f9639', 'xia5-000.png': 'e9fb7fbb348567d3'}
im = {f: T.Raster(os.path.join(RAS, f), s).im for f, s in SHA.items()}
out = {}

# ---- scale bars, anti-aliased ends (trace_lib.bar_length_aa; the windows of measure_scales.py)
bars = {}
for i, r0 in enumerate((470, 952, 1436, 1920)):
    b = T.bar_length_aa(im['dinh5-000.png'], r0, r0 + 20, 870, 975)
    bars['dinh_fig2B_panel%d' % (i + 1)] = dict(um=50, px_aa=b[0], full_px=b[3])
b = T.bar_length_aa(im['gif5-001.png'], 550, 600, 690, 829)
bars['gifford_fig4_inset'] = dict(um=300, px_aa=b[0], full_px=b[3])
b = T.bar_length_aa(im['xia5-000.png'], 3085, 3120, 1860, 2025, dark=False)
bars['xia_fig3_panel_i'] = dict(um=250, px_aa=b[0], full_px=b[3])
b = T.bar_length_aa(im['xia5-000.png'], 2925, 2960, 1000, 1230)
bars['xia_fig3a'] = dict(um=5000, px_aa=b[0], full_px=b[3])
for v in bars.values():
    v['um_per_px'] = v['um'] / v['px_aa']
out['bars'] = bars

# ---- Dinh Fig. 2B: the row pitch on panels 2-4 (two post columns each), the flare onset at the corrected bar
dres = json.load(open(os.path.join(R, 'provenance', 'rc6_traces', 'results', 'dinh', 'dinh_results_gp26.json')))
read = {p: [float(v) for v in dres[p]['measured']['pitch_um_by_column']] for p in ('p2', 'p3', 'p4')}
for p in read:                              # the traces read them at the same bar (46.906 px): check
    assert abs(dres[p]['measured']['bar_px_aa'] - bars['dinh_fig2B_panel2']['px_aa']) < 0.05
allr = np.array([v for p in read for v in read[p]])
a5 = json.load(open(os.path.join(R, 'opencif', 'data', 'r5', 'figure_anchors.json')))['dinh2024_fig2b_panel4']
f_bar = a5['scale_bar_px'] / bars['dinh_fig2B_panel4']['px_aa']          # 46 / 46.906: every rc5 length x this
on5 = a5['flare_onset_from_last_gap_centre_um']
po5 = a5['flare_onset_from_last_post_um']
onset = dict(left=on5['left'] * f_bar, right=on5['right'] * f_bar)
out['dinh'] = dict(pitch_readings_um=read, row_pitch_um=float(allr.mean()), row_pitch_sd_um=float(allr.std(ddof=1)),
                   G_printed_um=19.0, g_p_um=float(allr.mean()) - 19.0, rc5_to_rc6_length_factor=f_bar,
                   flare_onset_from_last_gap_centre_um=onset, flare_onset_um=0.5 * (onset['left'] + onset['right']),
                   flare_onset_from_last_post_um=dict(left=po5['left'] * f_bar, right=po5['right'] * f_bar))

# ---- Gifford Fig. 4 inset: the row pitch on all four post columns (each column's pills indexed along the leg and
# fitted), and l_g = each leg's post columns' centre-to-centre separation minus the printed 100 um centre
gi = im['gif5-001.png']
kg = bars['gifford_fig4_inset']['px_aa'] / 300.0
box = T.micrograph_box(gi, (190, 40, 829, 560))
C = PG.pill_centroids(gi, box)
u0 = np.array([1.0, 1.0]) / np.sqrt(2)
cols, v = PG.columns(C, u0)
u = PG.direction(cols, u0)
cols, v = PG.columns(C, u)
vimg = np.array([u[1], -u[0]])
cols = [cols[i] for i in np.argsort([np.mean(P[:, :2] @ vimg) for P in cols])]
assert len(cols) == 4
pitch_guess_px = 46.1 * kg                  # only indexes the pills (rc5's pitch); the fit sets the value
gp = []
for P in cols:
    pa = np.sort(P[:, :2] @ u)
    idx = np.round((pa - pa[0]) / pitch_guess_px).astype(int)
    gp.append(float(np.polyfit(idx, pa, 1)[0]) / kg)
cent = [float(np.mean(P[:, :2] @ vimg)) for P in cols]
sep = dict(leg2=(cent[1] - cent[0]) / kg, leg4=(cent[3] - cent[2]) / kg)
out['gifford'] = dict(scale_um_per_px=1.0 / kg, leg_direction_deg=float(np.degrees(np.arctan2(u[1], u[0]))),
                      row_pitch_readings_um=gp, row_pitch_um=float(np.mean(gp)), G_record_um=20.5,
                      g_p_um=float(np.mean(gp)) - 20.5, post_column_separation_um=sep,
                      l_g_um=float(np.mean(list(sep.values()))) - 100.0, posts_per_column=[len(P) for P in cols])

# ---- Xia Fig. 3: the scale of panels b-g
xi = im['xia5-000.png']
ev = PX.scale_evidence(xi)                  # printed centre widths 300 (panel c) / 150 (panel g)
c_rows = PX.centre_rows(xi, 1519.5, 1725.5, 450, 745, 21, 1400, 1850)
g_rows = PX.centre_rows(xi, 1769.4, 1891.6, 1600, 1975, 21, 1700, 1950)
Cc = float(np.mean((c_rows[:, 3] + c_rows[:, 4] - c_rows[:, 1] - c_rows[:, 2]) / 2.0))     # post columns, centre to centre
Cg = float(np.mean((g_rows[:, 3] + g_rows[:, 4] - g_rows[:, 1] - g_rows[:, 2]) / 2.0))
f3a = 179.0 / bars['xia_fig3a']['px_aa']    # r4 read Fig. 3a at 179 px
m3a = json.load(open(os.path.join(R, 'opencif', 'data', 'r4', 'meas_xia_fig3a.json')))
legY = {int(k_): v_ * f3a for k_, v_ in m3a['leg_Y_um'].items()}
legY_r = {k_: float(int(round(v_))) for k_, v_ in legY.items()}
arcs = json.load(open(os.path.join(R, 'provenance', 'rc6_traces', 'measurements', 'xia_fig3d_arc_scale.json')))
rad = arcs['left_set']['radii_px']          # loop 4's arcs on panel d (concentric fit), inner to outer
axis_px = 0.5 * (rad[2] + rad[3])           # its centre channel's two walls: the loop's axis radius
Rt4 = 0.5 * (legY_r[4] - legY_r[5])         # loop 4 joins legs 4 and 5 (Fig. 3a's axes at its bar, as the layout)
scales = dict(printed_centre_widths=ev['um_per_px_printed'], post_columns_centre_to_centre=150.0 / (Cc - Cg),
              loop4_axis_curvature=Rt4 / axis_px,
              fig3a_leg_spacing=abs(legY[3] - legY[5]) / ev['panel_d_axis_separation_px'],   # legs 3 and 5 at Fig. 3a's bar
              bar_panel_i=bars['xia_fig3_panel_i']['um_per_px'])
s_x = round(float(np.mean([scales['printed_centre_widths'], scales['post_columns_centre_to_centre'],
                           scales['loop4_axis_curvature']])), 3)


def band_period(c0, c1, r0, r1, guess, m_max=10):
    """the period (px) of a band's row profile along the image rows: the autocorrelation's peaks near m x guess
    (parabolic refinement), m = 1..m_max, fitted through the origin"""
    a = xi[r0:r1, c0:c1].mean(axis=1)
    a = a - a.mean()
    n = a.shape[0]
    pk = []
    for m in range(1, m_max + 1):
        lags = np.arange(int(m * guess) - 5, int(m * guess) + 6)
        if lags[-1] + 2 >= n:
            break
        ac = np.array([np.sum(a[:n - l] * a[l:]) / (n - l) for l in lags])
        j = int(np.argmax(ac[1:-1])) + 1
        den = ac[j - 1] - 2 * ac[j] + ac[j + 1]
        pk.append((m, lags[j] + (0.5 * (ac[j - 1] - ac[j + 1]) / den if den != 0 else 0.0)))
    m_, q_ = np.array([v[0] for v in pk], float), np.array([v[1] for v in pk])
    return float(np.sum(m_ * q_) / np.sum(m_ * m_)), [float(v) for v in q_]


def pill_pitch(col, r0, r1, post_um=38.0, k=1.0 / 1.803):
    """the row pitch (px) along one post column: the pills' centres (pairs of dark outline lines on the column's
    profile, a pill post_um x k +-25 % long: a window, not a value), indexed by the median interval and fitted"""
    from scipy import ndimage as ndi
    p = ndi.gaussian_filter1d(T.col_profile(xi, col, r0, r1, 2), 0.8)
    L = T.dark_lines(p, 6, halfwin=5, off=r0)
    pc = np.array(T.pair_centres(L, 0.75 * post_um * k, 1.25 * post_um * k, skip=1))
    d = np.diff(pc)
    idx = np.round((pc - pc[0]) / np.median(d[(d > 20) & (d < 36)])).astype(int)
    return float(np.polyfit(idx, pc, 1)[0])


per = {lab: pill_pitch(col, r0, r1) for lab, (col, r0, r1) in {
    'c_L': (1519.5, 450, 745), 'c_R': (1725.5, 450, 745), 'g_L': (1769.4, 1600, 1975), 'g_R': (1891.6, 1600, 1975)}.items()}
Lam_px = float(np.mean(list(per.values())))
units = {lab: band_period(c0, c1, r0, r1, 65.0)[0] for lab, (c0, c1, r0, r1) in {
    'b_units_right': (330, 560, 90, 745), 'b_units_left': (5, 160, 90, 745),
    'c_blocks_L': (1500, 1545, 5, 440), 'c_blocks_R': (1705, 1750, 5, 440)}.items()}
u_px = float(np.mean(list(units.values())))
xia = dict(scale_um_per_px=scales, scale_used_um_per_px=s_x,
           scale_note='the mean of the three determinations precise to < 1 %% (printed centre widths, post columns '
                      'centre to centre, loop-4 axis curvature); the panel-i bar gives %.4f' % scales['bar_panel_i'],
           fig3a_length_factor_vs_rc5=f3a, post_columns_centre_to_centre_px=dict(panel_c=Cc, panel_g=Cg),
           loop4=dict(axis_radius_px=axis_px, axis_radius_um=Rt4),
           row_pitch_px=per, row_pitch_um=Lam_px * s_x, G_printed_um=16.0, g_p_um=Lam_px * s_x - 16.0,
           l_g_um=dict(panel_g=Cg * s_x - 150.0, panel_c=Cc * s_x - 300.0), unit_pitch_px=units,
           unit_pitch_2s_um=u_px * s_x, GS_printed_um=22.8, wall_M_um=u_px * s_x / 2.0 - 22.8,
           head_width_um=2000.0 * s_x / (250.0 / 128.0),
           fig3a_leg_straight_length_um=m3a['leg_straight_length_um'] * f3a)
xia['l_g_um_used'] = round(float(np.mean(list(xia['l_g_um'].values()))), 1)
xia['rows_per_leg_to_fill_fig3a'] = xia['fig3a_leg_straight_length_um'] / xia['row_pitch_um']
out['xia'] = xia
json.dump(out, open(OUT, 'w'), indent=1)

# ---- the values the engine reads
vals = dict(
    note='rc6 values re-read from the papers\' figures (provenance/rc6_scripts/measure_rc6.py; provenance/rc6_traces): '
         'scale bars measured with their anti-aliased ends; Xia panels b-g at the scale their printed widths and loop '
         'curvature give.  Rounded as the presets carry them.',
    dinh=dict(g_p_um=round(out['dinh']['g_p_um'], 1), flare_start_um=round(out['dinh']['flare_onset_um'], 1),
              row_pitch_um=round(out['dinh']['row_pitch_um'], 2)),
    gifford=dict(g_p_um=round(out['gifford']['g_p_um'], 1), l_g_um=round(out['gifford']['l_g_um'], 1),
                 row_pitch_um=round(out['gifford']['row_pitch_um'], 2)),
    xia=dict(scale_um_per_px=s_x, fig3a_length_factor=round(f3a, 6), g_p_um=round(xia['g_p_um'], 1),
             l_g_um=xia['l_g_um_used'], wall_M_um=round(xia['wall_M_um'], 1),
             head_width_um=round(xia['head_width_um'] / 10.0) * 10.0, row_pitch_um=round(xia['row_pitch_um'], 2),
             # Fig. 3a's layout (r4/meas_xia_fig3a.json, read at 179 px) at the bar's anti-aliased length: x f3a
             layout=dict(leg_Y_um=[legY_r[i] for i in range(1, 8)],
                         inlet_port_radius_um=int(round(m3a['inlet_port']['r_um'] * f3a)),
                         centre_port_radius_um=int(round(m3a['centre_port']['r_um'] * f3a)),
                         filtrate_port_radius_um=int(round(m3a['filtrate_port']['r_um'] * f3a)),
                         ring_radius_um=int(round(m3a['ring']['R_centreline_um'] * f3a)),
                         ring_centre_offset_um=int(round(m3a['ring']['centre_beyond_centre_port_um'] * f3a)),
                         ring_filtrate_port_offset_um=int(round(m3a['ring']['filtrate_port_beyond_centre_port_um'] * f3a)),
                         pc_inlet_offset_um=int(round(m3a['inlet_beyond_top_loops_um'] * f3a)),
                         pc_out_ext_um=int(round(m3a['ring']['top_beyond_bottom_loops_um'] * f3a)),
                         head_top_below_inlet_mm=round(m3a['head_top_below_inlet_mm'] * f3a, 3),
                         first_pill_below_inlet_mm=round(m3a['first_pill_below_inlet_mm'] * f3a, 3))))
vals['xia']['layout']['inlet_port_gap_um'] = int(round(vals['xia']['layout']['head_top_below_inlet_mm'] * 1000.0
                                                       - vals['xia']['layout']['inlet_port_radius_um']))
json.dump(vals, open(VALS, 'w'), indent=1)

# ---- the values the lost first rc6 run recorded (RECONSTRUCTION_PLAN.md A.11): every rounded value must agree
LOST = dict(dinh=dict(g_p_um=26.1, flare_start_um=65.4, row_pitch_um=45.11),
            gifford=dict(g_p_um=25.3, l_g_um=59.5, row_pitch_um=45.82),
            xia=dict(scale_um_per_px=1.798, fig3a_length_factor=1.003827, g_p_um=35.8, l_g_um=70.9, wall_M_um=35.8,
                     head_width_um=1840.0, row_pitch_um=51.78,
                     layout=dict(leg_Y_um=[0.0, 1332.0, -1745.0, 2229.0, -2650.0, 3267.0, -3975.0],
                                 inlet_port_radius_um=862, centre_port_radius_um=1723, filtrate_port_radius_um=1718,
                                 ring_radius_um=4195, ring_centre_offset_um=1191, ring_filtrate_port_offset_um=4781,
                                 inlet_port_gap_um=591, pc_inlet_offset_um=2553, pc_out_ext_um=3441,
                                 head_top_below_inlet_mm=1.453, first_pill_below_inlet_mm=13.786)))
bad = []
for dev, d in LOST.items():
    for k_, v_ in d.items():
        if isinstance(v_, dict):
            for k2, v2 in v_.items():
                if vals[dev][k_][k2] != v2:
                    bad.append((dev, k_, k2, vals[dev][k_][k2], v2))
        elif vals[dev][k_] != v_:
            bad.append((dev, k_, vals[dev][k_], v_))
print(json.dumps(vals, indent=1))
print('scales', json.dumps(scales))
print('against the lost run: %s' % ('all values agree' if not bad else 'DIFFERENT %s' % bad))
sys.exit(0 if not bad else 2)
