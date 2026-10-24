"""usage: run_dinh.py ROOT RASTERS OUT TAG [PRESET DXF RECORD LABEL] [--noscan] [--k K]"""
import json, sys, os, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import panels_dinh as PD
import panels_dinh_v2 as PV
a = [x for x in sys.argv[1:] if not x.startswith('--')]
scan = '--noscan' not in sys.argv
kov = None
if '--k' in sys.argv:
    kov = float(sys.argv[sys.argv.index('--k') + 1]); a = [x for x in a if x != sys.argv[sys.argv.index('--k') + 1]]
R, RAS, OUT, TAG = a[:4]
kw = dict(tag=TAG)
if len(a) > 4:
    kw.update(preset=a[4], dxf=a[5], record=a[6], label=a[7])
os.makedirs(OUT, exist_ok=True)
ctx = PD.setup(R, RAS, **kw)
if kov: ctx['k_override'] = kov
def conv(o):
    if isinstance(o, np.floating): return float(o)
    if isinstance(o, np.integer): return int(o)
    if isinstance(o, np.ndarray): return o.tolist()
    if isinstance(o, np.bool_): return bool(o)
    raise TypeError(type(o))
res = {}
t0 = time.time()
res['p4'] = PV.panel4(ctx, OUT, scan); print('p4 %.0fs' % (time.time() - t0), flush=True)
res['p3'] = PV.panel3(ctx, OUT, res['p4'], scan); print('p3 %.0fs' % (time.time() - t0), flush=True)
res['p2'] = PV.panel2(ctx, OUT, scan); print('p2 %.0fs' % (time.time() - t0), flush=True)
res['p1'] = PV.panel1(ctx, OUT, scan); print('p1 %.0fs' % (time.time() - t0), flush=True)
for k, v in res.items():
    m = v['measured']; V1 = v['V1_feature_vs_chamfer_px']; V2 = v['V2_chamfer']; mt = v['metric']; nl = v['metric_null']
    print('== %s %s' % (TAG, k))
    print('  k %.5f (bar %.2f, full px %d)  pitch img %s  mask %.1f' % (m['px_per_um'], m['bar_px_aa'], m['bar_px_full_pixels'], m.get('pitch_um'), v['pitch_um']))
    print('  V1 along %.2f across %.2f agree %s' % (V1['along'], V1['across'], V1['agree']))
    print('  V2 n %d best %.3f  1pitch %s  half %s  5um %s  median %.3f' % (V2['points'], V2['best_px'], np.round(V2['one_pitch_along_px'], 3), np.round(V2['half_pitch_along_px'], 3), np.round(V2['five_um_across_px'], 3), V2['landscape_median_px']))
    print('  metric median %.2f p90 %.2f signed %.2f n %d | null half %.2f  5um %.2f  1pitch %.2f' % (mt['median_abs_um'], mt['p90_abs_um'], mt['mean_signed_um'], mt['n'], nl['half_pitch_along']['median_abs_um'], nl['5um_across']['median_abs_um'], nl['one_pitch_along']['median_abs_um']), ' same-pts %.2f' % nl['registered_same_points']['median_abs_um'])
    sys_cls = {c: round(d['mean_signed_um'], 2) for c, d in v['V3_by_class'].items() if d['n'] >= 20}
    print('  V3', sys_cls)
    if 'V5_scale_scan' in v:
        s = v['V5_scale_scan']; print('  V5 along %.4f across %.4f  score best %.3f at bar %.3f edge %s' % (s['best_along'], s['best_across'], s['score_best'], s['score_at_bar'], s['on_grid_edge']))
    if 'image_widths_um' in m: print('  widths', {a: round(b, 2) for a, b in m['image_widths_um'].items()})
    for kk in ('mask_widths_um_last_row', 'mask_widths_um_at_row', 'mask_widths_um_third_post_row'):
        if kk in m: print('  mask  ', {a: round(b, 2) for a, b in m[kk].items()})
    for kk in ('row_from_post_columns', 'row_from_outer_walls', 'rows_traced', 'fraction_along_linear', 'first_post_row'):
        if kk in m: print('  ', kk, m[kk])
    print('  rotation %.3f deg (sd %.2f px)  walls %s' % (m['rotation_deg'], m.get('rotation_fit_sd_px', -1), m.get('wall_angles_deg')))
    if 'feature_fit_residual_um' in m:
        print('  feature residuals um along', np.round(m['feature_fit_residual_um']['along'], 2).tolist(), 'across', np.round(m['feature_fit_residual_um']['across'], 2).tolist())
    for kk in ('stem_width_um', 'first_junction_rows_px', 'mask_first_channel_edge_x_um'):
        if kk in m: print('  ', kk, m[kk])
json.dump(res, open(os.path.join(OUT, 'dinh_results_%s.json' % TAG), 'w'), default=conv, indent=1)
