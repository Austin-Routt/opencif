"""Outlet-balance (+-1 % side outlet) numbers over the rc4 records."""
import json, glob, os, sys
R = sys.argv[1] if len(sys.argv) > 1 else 'opencif_rc4/runs'   # a tree's runs/
rows = []
for f in sorted(glob.glob(R + '/*/*.json')):
    n = os.path.basename(os.path.dirname(f))
    if os.path.basename(f) != n + '.json':
        continue
    r = json.load(open(f))
    if 'validation' not in r:
        continue
    ms = r['validation'].get('model_scenarios') or {}
    ob = [ms.get(k) for k in ('outlet_side_resistance_+1%', 'outlet_side_resistance_-1%')]
    if not all(ob):
        continue
    N = r['design']['N']
    nr = max(x['rows_newly_outside_tol'] for x in ob)
    mx = max(x['max_cutoff_change'] or 0 for x in ob)
    sh = max(x.get('share_of_F_in_those_rows') or 0 for x in ob)
    rws = sorted(set(ob[0]['newly_outside_rows']) | set(ob[1]['newly_outside_rows']))
    last = (N - min(rws) + 1) if rws else 0
    rows.append((n, N, nr, mx, sh, ms.get('dDlast_dlnRside_out'), last, r['assessment']['mode'], (r['input'].get('outlet') or {}).get('style', 'separate')))
for x in rows:
    print('%-48s N %5d rows %3d max %.3f D shareF %.4f dD/dln %8.1f within last %d rows %s %s' % x)
nr = [x[2] for x in rows]; mx = [x[3] for x in rows]; sh = [x[4] for x in rows]; last = [x[6] for x in rows]
print('designs', len(rows), 'rows range', min(nr), max(nr), 'max change range %.3f-%.3f' % (min(mx), max(mx)),
      'share F max %.4f' % max(sh), 'rows lie within the last', max(last), 'rows; designs with 0 rows:', sum(1 for v in nr if v == 0))
