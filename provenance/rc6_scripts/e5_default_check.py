"""Reviewer A's E5 device (test_negative 13's u05) at the rc6 default fillet: the full verdict and the failing checks
(LIMITATIONS 23).  usage: python3 e5_default_check.py OPENCIF_ROOT OUT_DIR  ->  OUT_DIR/e5_default_result.json (+ the build)
rc6 result (2026-10-17, code 4bbdaabef1225a24): e5_default_result.json beside this script."""
import copy, json, os, sys
R, out = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2])
os.makedirs(out, exist_ok=True)
sys.path[:0] = [R]
from opencif import api, pipeline as PL
u05 = api.build_spec(1.5, 0.50, process=dict(depth_um=60.0, gap_um=12.0, post_length_um=20.0, post_width_um=40.0,
                                              min_feature_um=10.0), chip=(25, 75), layout='linear', taper='constant_shear',
                     width_in_um=40.0, width_min_um=30.0, meander='narrow', inlet='plain', outlet='teardrop',
                     flow_uL_min=5.0, name='_t_e5')
r = PL.run(copy.deepcopy(u05), outdir=out, verbose=False, render=False)
a = r['assessment']
res = dict(code=r['code'], exit_code=a['exit_code'], verdict=a['verdict'],
           failing=[dict(check=c['check'], detail=str(c.get('detail', ''))[:400]) for c in a['checks'] if not c['ok']],
           features_found_ok=(r['validation'].get('features_found') or {}).get('ok'),
           fillet_um=r['design']['spec'].get('r_fillet'))
json.dump(res, open(os.path.join(out, 'e5_default_result.json'), 'w'), indent=1)
print(json.dumps(res, indent=1)[:3000])
