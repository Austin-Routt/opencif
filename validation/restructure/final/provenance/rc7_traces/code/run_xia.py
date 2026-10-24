"""usage: run_xia.py ROOT RASTERS OUT TAG SCALE(bar|printed|rc6) [PRESET DXF RECORD LABEL] [--noscan]"""
import json, sys, os, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
a = [x for x in sys.argv[1:] if not x.startswith('--')]
R, RAS, OUT, TAG, SCALE = a[:5]
os.environ['XIA_RASTERS'] = RAS
for p in (R,):
    sys.path.insert(0, p)
import trace_lib as T
import panels_xia as PX
from opencif import spec as SP
from opencif.physics import schedule as CS
name = 'xia2016_rs7_depicted'
preset = a[5] if len(a) > 5 else os.path.join(R, 'presets', name + '.yaml')
dxf = a[6] if len(a) > 6 else os.path.join(R, 'runs', name, name + '.dxf')
record = a[7] if len(a) > 7 else os.path.join(R, 'runs', name, name + '.json')
label = a[8] if len(a) > 8 else 'the shipped mask ' + name
rec = json.load(open(record))
mask = T.Mask(dxf)
assert mask.sha256 == rec['files']['dxf_sha256']
sp, op = SP.to_engine(SP.load(preset))
im = T.Raster(os.path.join(RAS, 'xia5-000.png'), PX.FIG3_SHA).im
kb, bar = PX.bar_scale(im)
ev = PX.scale_evidence(im)
if SCALE == 'bar':
    k, note = kb, 'the panel-i bar, 129.11 px per 250 µm'
elif SCALE == 'rc6':                 # the scale the rc6 presets are built at (data/r6/rc6_values.json)
    s6 = json.load(open(os.path.join(R, 'opencif', 'data', 'r6', 'rc6_values.json')))['xia']['scale_um_per_px']
    k, note = 1.0 / s6, ('the rc6 scale %.3f µm/px: the mean of the printed widths, the post columns\' separations and '
                         'the loop curvature (data/r6/rc6_values.json)' % s6)
else:
    k, note = ev['px_per_um_printed'], 'the printed widths w_c 300 / 150 µm (panels c, g)'
os.makedirs(OUT, exist_ok=True)
t0 = time.time()
v = PX.panel_g(mask, sp, rec, OUT, k, TAG, label, note, scan='--noscan' not in sys.argv)
v['scale_evidence'] = ev
v['bar'] = dict(px_aa=bar[0], full_pixels=bar[3])
def conv(o):
    if isinstance(o, np.floating): return float(o)
    if isinstance(o, np.integer): return int(o)
    if isinstance(o, np.ndarray): return o.tolist()
    if isinstance(o, np.bool_): return bool(o)
    raise TypeError(type(o))
json.dump(v, open(os.path.join(OUT, 'xia_results_%s.json' % TAG), 'w'), default=conv, indent=1)
m = v['measured']; V1 = v['V1_feature_vs_chamfer_px']; V2 = v['V2_chamfer']; mt = v['metric']; nl = v['metric_null']
print('time %.0fs  k %.5f (%.4f um/px)  rot %.3f' % (time.time() - t0, k, 1 / k, m['rotation_deg']))
print('meas', {a: b for a, b in m.items() if a not in ('feature_fit_residual_um',)})
print('feature resid along', np.round(m['feature_fit_residual_um']['along'], 2), 'across', np.round(m['feature_fit_residual_um']['across'], 2))
print('V1', V1)
print('V2 n %d best %.3f 1p %s half %s 5um %s med %.3f' % (V2['points'], V2['best_px'], np.round(V2['one_pitch_along_px'], 3), np.round(V2['half_pitch_along_px'], 3), np.round(V2['five_um_across_px'], 3), V2['landscape_median_px']))
print('metric', mt, '\n null', {a: b.get('median_abs_um') for a, b in nl.items()})
print('V3', {c: (round(d['mean_signed_um'], 2) if d['mean_signed_um'] is not None else None, d['n']) for c, d in v['V3_by_class'].items() if d['n'] >= 20})
if 'V5_scale_scan' in v:
    s = v['V5_scale_scan']; print('V5 along %.4f across %.4f best %.3f at scale %.3f edge %s' % (s['best_along'], s['best_across'], s['score_best'], s['score_at_bar'], s['on_grid_edge']))
