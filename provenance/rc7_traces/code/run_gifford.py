import json, sys, os, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import panels_gifford as PG
a = [x for x in sys.argv[1:] if not x.startswith('--')]
R, RAS, OUT, TAG = a[:4]
kw = {}
if len(a) > 4:
    kw = dict(preset=a[4], dxf=a[5], record=a[6], label=a[7])
os.makedirs(OUT, exist_ok=True)
t0 = time.time()
v = PG.run(R, RAS, OUT, tag=TAG, scan='--noscan' not in sys.argv, **kw)
def conv(o):
    if isinstance(o, np.floating): return float(o)
    if isinstance(o, np.integer): return int(o)
    if isinstance(o, np.ndarray): return o.tolist()
    if isinstance(o, np.bool_): return bool(o)
    raise TypeError(type(o))
json.dump(v, open(os.path.join(OUT, 'gifford_results_%s.json' % TAG), 'w'), default=conv, indent=1)
m = v['measured']; V1 = v['V1_feature_vs_chamfer_px']; V2 = v['V2_chamfer']; mt = v['metric']; nl = v['metric_null']
print('time %.0fs' % (time.time() - t0))
print('k %.5f (bar %.2f, full %d) um/px %.4f  dir %.3f deg' % (m['px_per_um'], m['bar_px_aa'], m['bar_px_full_pixels'], m['um_per_px'], m['leg_direction_deg']))
print('columns', np.round(m['columns_across_px'], 1), m['posts_per_column'], 'leg4', m['leg4'])
print('feature fit', m['feature_fit']['rms_um'], np.round(m['feature_fit']['along_um'], 2), np.round(m['feature_fit']['across_um'], 2))
print('post col sep', m['post_column_separation_um'], 'leg sep', m['leg_centreline_separation_um'], 'pitch', m['pitch_um'])
print('V1', V1)
print('V2 n %d best %.3f 1p %s half %s 5um %s med %.3f' % (V2['points'], V2['best_px'], np.round(V2['one_pitch_along_px'], 3), np.round(V2['half_pitch_along_px'], 3), np.round(V2['five_um_across_px'], 3), V2['landscape_median_px']))
print('metric', mt, '\n null', {a: b.get('median_abs_um') for a, b in nl.items()})
print('V3', {c: (round(d['mean_signed_um'], 2) if d['mean_signed_um'] is not None else None, d['n']) for c, d in v['V3_by_class'].items()})
if 'V5_scale_scan' in v:
    s = v['V5_scale_scan']; print('V5 along %.4f across %.4f best %.3f at bar %.3f edge %s' % (s['best_along'], s['best_across'], s['score_best'], s['score_at_bar'], s['on_grid_edge']))
for L in m['leg4_profile_lines_um']: print('  prof', L)
