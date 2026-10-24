"""rc7 (DINH-1): the square quarter-turn's excess length c (in w, against the centreline through the corner point) in
both limits, for outer-corner roundings r_out/w and inner fillets r_in/w; n = 40 and 80 cells per w, first-order
(staircase) extrapolation 2 c80 - c40, as validated on rc6's exact arc (validate: arc rm/w 1.0).
Usage: python3 -I corner_table.py OUT.json"""
import json, math, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import corner2d as C

out = dict(method=__doc__, cases=[])
cases = [('validate_arc_rm1.0', 'arc', dict(rm=1.0))]
for ro in (0.0, 0.214, 0.25, 0.307):
    for ri in (0.0, 0.18):
        cases.append(('square_rout%.3f_rin%.2f' % (ro, ri), 'square', dict(r_out=ro, r_in=ri)))
for lab, kind, kw in cases:
    t = time.time()
    row = dict(case=lab, kind=kind, **kw)
    for limit in ('shallow', 'deep'):
        vals = []
        for n in (40, 80):
            if limit == 'shallow':
                vals.append(C.hele_shaw(n, kind, **kw))
            else:
                vals.append(C.deep(n, kind, method='dissipation', **kw))
        row[limit] = dict(n40=vals[0], n80=vals[1], extrapolated=2 * vals[1] - vals[0])
    row['seconds'] = round(time.time() - t, 1)
    out['cases'].append(row)
    print(json.dumps(row), flush=True)
    json.dump(out, open(sys.argv[1], 'w'), indent=1)
