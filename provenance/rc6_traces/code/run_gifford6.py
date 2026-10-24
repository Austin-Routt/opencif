"""usage: run_gifford6.py ROOT RASTERS OUT TAG [PRESET DXF RECORD LABEL] [--one-map [--turns]] [--noscan]
Gifford Fig. 6 traced with a Fig. 6 mask (panels_gifford6.py): one results file, a list of panel records --
gifford6_results_TAG.json (each leg on its own map), _onemap.json (--one-map --turns) or _onemap_noturns.json."""
import json, sys, os, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import panels_gifford6 as P6
a = [x for x in sys.argv[1:] if not x.startswith('--')]
R, RAS, OUT, TAG = a[:4]
kw = {}
if len(a) > 4:
    kw = dict(preset=a[4], dxf=a[5], record=a[6], label=a[7])
os.makedirs(OUT, exist_ok=True)
t0 = time.time()
recs = P6.run(R, RAS, OUT, tag=TAG, scan='--noscan' not in sys.argv, one_map='--one-map' in sys.argv,
             turns='--turns' in sys.argv, **kw)


def conv(o):
    if isinstance(o, np.floating): return float(o)
    if isinstance(o, np.integer): return int(o)
    if isinstance(o, np.ndarray): return o.tolist()
    if isinstance(o, np.bool_): return bool(o)
    raise TypeError(type(o))


# one record per panel, keyed by the panel's leg ('leg7', 'leg5', 'leg3', 'leg1') or 'four_legs' (full_device.py's "key")
keyed = {('four_legs' if 'four_legs' in v['name'] else v['name'].split('_')[2]): v for v in recs}
# the per-leg run and the one-map run of the same mask go to separate files (rc6: the shipped Fig. 6 preset has both)
# (a one-map run without --turns goes to its own file, so that it cannot replace the turns record summary_table reads)
suffix = ('_onemap' if '--turns' in sys.argv else '_onemap_noturns') if '--one-map' in sys.argv else ''
json.dump(keyed, open(os.path.join(OUT, 'gifford6_results_%s%s.json' % (TAG, suffix)), 'w'),
          default=conv, indent=1)
print('time %.0fs' % (time.time() - t0))
m = recs[0]['measured']
print('identification', m['identification']['best'])
for c in m['identification']['candidates']:
    if c['same_flow_direction']:
        print('   candidate', c)
print('leg spacing mask %.1f image %s; last leg end offset %.1f um' % (m['mask_leg_spacing_um'], np.round(m['image_leg_spacing_um'], 1), m['last_leg_end_offset_um']))
for t in m['side_channels_at_leg_ends']:
    print('   side', t)
for v in recs:
    V1 = v['V1_feature_vs_chamfer_px']; V2 = v['V2_chamfer']; mt = v['metric']; nl = v['metric_null']
    print(v['name'])
    print('   feature fit rms %.2f um (n %d)' % (v['measured']['feature_fit']['rms_um'], v['measured']['feature_fit']['n']))
    if 'offset_from_leg5_map_um' in v['measured']:
        print('   offset from leg 5 map', v['measured']['offset_from_leg5_map_um'])
    if 'turns' in v['measured']:
        tc = v['measured']['turns']
        print('   turns: image walls %d, mask walls %d; nearest mask edge median %.2f p90 %.2f max %.1f um; centre offset u %.1f v %.1f um; radii diff %.2f +- %.2f (max %.1f)' % (
            tc['image']['walls'], tc['mask']['walls'], tc['nearest_mask_edge_um']['median'], tc['nearest_mask_edge_um']['p90'], tc['nearest_mask_edge_um']['max'],
            tc['centre_offset_mask_minus_image_um']['u'], tc['centre_offset_mask_minus_image_um']['v'], tc['radii_diff_um']['mean'], tc['radii_diff_um']['sd'], tc['radii_diff_um']['max_abs']))
        rc = tc['radii_about_image_centre']
        print('   turns: radii about the image centre, mask - image %.2f +- %.2f (max %.1f) um' % (rc['mean_um'], rc['sd_um'], rc['max_abs_um']))
    if 'leg7_own_map_minus_common_um' in v['measured']:
        print('   leg 7 own map - common', v['measured']['leg7_own_map_minus_common_um'])
    print('   V1', V1)
    print('   V2 n %d best %.3f 1p %s half %s 5um %s med %.3f' % (V2['points'], V2['best_px'], np.round(V2['one_pitch_along_px'], 3), np.round(V2['half_pitch_along_px'], 3), np.round(V2['five_um_across_px'], 3), V2['landscape_median_px']))
    print('   metric', {x: round(y, 2) if isinstance(y, float) else y for x, y in mt.items()}, 'null', {x: round(y.get('median_abs_um', np.nan), 2) for x, y in nl.items()})
    print('   V3', {c: (round(d['mean_signed_um'], 2) if d['mean_signed_um'] is not None else None, d['n']) for c, d in v['V3_by_class'].items()})
    if 'V5_scale_scan' in v:
        s = v['V5_scale_scan']; print('   V5 along %.4f across %.4f best %.3f at 1: %.3f edge %s' % (s['best_along'], s['best_across'], s['score_best'], s['score_at_bar'], s['on_grid_edge']))
