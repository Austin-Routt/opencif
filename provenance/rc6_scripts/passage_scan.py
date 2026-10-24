"""rc6 (LIMITATIONS 26): the narrowest straight side channel of every shipped record against its cutoff D -- where a
particle the gaps pass could meet a side channel narrower than itself.  Reads runs/<name>/<name>.json (design.ws_min_um,
design.D_star_um); writes provenance/rc6_scripts/passage_scan.json and prints the sentence LIMITATIONS 26 carries.
opencif has no passage rule: this is a report over the shipped records, not a check the generator makes.
usage: python3 passage_scan.py OPENCIF_ROOT"""
import glob
import json
import os
import sys

R = sys.argv[1]
rows = []
for p in sorted(glob.glob(os.path.join(R, 'runs', '*', '*.json'))):
    n = os.path.basename(os.path.dirname(p))
    if os.path.basename(p) != n + '.json':
        continue
    try:
        d = json.load(open(p)).get('design') or {}
    except (ValueError, OSError):
        continue
    ws, D = d.get('ws_min_um'), d.get('D_star_um')
    if ws is None or not D:
        continue
    rows.append(dict(name=n, ws_min_um=float(ws), D_um=float(D), ratio=float(ws) / float(D)))
rows.sort(key=lambda r: r['ratio'])
below = [r for r in rows if r['ratio'] < 1.0]
within2 = [r for r in rows if 1.0 <= r['ratio'] < 2.0]
out = dict(note=__doc__.split('\n')[0], records=len(rows), below_D=below, between_1_and_2_D=within2, all=rows)
json.dump(out, open(os.path.join(R, 'provenance', 'rc6_scripts', 'passage_scan.json'), 'w'), indent=1)
fmt = lambda L: '; '.join('`%s` %.2f µm against D %.2f µm' % (r['name'], r['ws_min_um'], r['D_um']) for r in L)
print('Over the %d shipped records (`provenance/rc6_scripts/passage_scan.py`), %d draw a straight side channel narrower '
      'than D%s, and %d one narrower than 2 D%s.' % (
          len(rows), len(below), (' (' + fmt(below) + ')') if below else '', len(within2),
          (' (' + fmt(within2) + ')') if within2 else ''))
