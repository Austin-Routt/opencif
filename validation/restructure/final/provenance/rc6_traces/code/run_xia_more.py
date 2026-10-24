"""usage: run_xia_more.py ROOT RASTERS OUT TAG PANELS(e.g. cdef) [PRESET DXF RECORD LABEL] [--noscan] [--fig3a JSON]
Xia Fig. 3 panels c, d, e, f at the rc6 scale (panels_xia_more.py); one results file per panel."""
import json, sys, os, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
F3A = None
argv = sys.argv[1:]
if '--fig3a' in argv:
    i_ = argv.index('--fig3a')
    F3A = json.load(open(argv[i_ + 1]))
    del argv[i_:i_ + 2]
a = [x for x in argv if not x.startswith('--')]
R, RAS, OUT, TAG, WHICH = a[:5]
for p in (R,):
    sys.path.insert(0, p)
import trace_lib as T
import panels_xia_more as PM
from opencif import spec as SP
name = 'xia2016_rs7_depicted'
preset = a[5] if len(a) > 5 else os.path.join(R, 'presets', name + '.yaml')
dxf = a[6] if len(a) > 6 else os.path.join(R, 'runs', name, name + '.dxf')
record = a[7] if len(a) > 7 else os.path.join(R, 'runs', name, name + '.json')
label = a[8] if len(a) > 8 else 'the shipped mask ' + name
rec = json.load(open(record))
mask = T.Mask(dxf)
assert mask.sha256 == rec['files']['dxf_sha256']
sp, op = SP.to_engine(SP.load(preset))
s6 = json.load(open(os.path.join(R, 'opencif', 'data', 'r6', 'rc6_values.json')))['xia']['scale_um_per_px']
k = 1.0 / s6
note = ('the rc6 scale %.3f µm/px: the mean of the printed widths, the post columns\' separations and the loop '
        'curvature (data/r6/rc6_values.json)' % s6)
os.makedirs(OUT, exist_ok=True)


def conv(o):
    if isinstance(o, np.floating): return float(o)
    if isinstance(o, np.integer): return int(o)
    if isinstance(o, np.ndarray): return o.tolist()
    if isinstance(o, np.bool_): return bool(o)
    raise TypeError(type(o))


for L in WHICH:
    t0 = time.time()
    v = PM.panel(L, mask, sp, rec, OUT, k, TAG, label, note, RAS, scan='--noscan' not in sys.argv, fig3a=F3A)
    json.dump(v, open(os.path.join(OUT, 'xia_results_fig3%s_%s.json' % (L, TAG)), 'w'), default=conv, indent=1)
    m = v['measured']; V1 = v['V1_feature_vs_chamfer_px']; mt = v['metric']; nl = v['metric_null']
    print('panel %s  %.0fs  rot %.3f  pitch %.2f (mask %.2f)' % (L, time.time() - t0, m['rotation_deg'], m['pitch_um']['image'], m['pitch_um']['mask']))
    for q in m['legs']:
        print('   leg %d centre %.1f um  located: %s' % (q['leg'], q['centre_um_image'], q['located']['by']))
    if 'leg_axis_separation_um' in m:
        print('   leg axes', m['leg_axis_separation_um'])
    print('   feature resid rms along %.2f across %.2f um' % (np.std(m['feature_fit_residual_um']['along']), np.std(m['feature_fit_residual_um']['across'])))
    print('   V1', V1)
    print('   metric median %.2f p90 %.2f  nulls %s' % (mt['median_abs_um'], mt['p90_abs_um'], {x: round(y.get('median_abs_um', np.nan), 2) for x, y in nl.items()}))
    if 'V5_scale_scan' in v:
        s = v['V5_scale_scan']; print('   V5 along %.4f across %.4f' % (s['best_along'], s['best_across']))
