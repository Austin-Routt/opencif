"""rc5 -> rc6 record comparison.  For every spec: verdict, exit code, N, NM, F (design, drawn), realised cutoff range,
footprint, the failing checks; and whether the spec's input changed between the trees (process, meander, outlet,
paper notes) so that a difference is expected.  A spec whose engine-relevant input did not change must reproduce rc5's
numbers exactly, apart from the new check, timings, the code stamp and the defaults_used text.
usage: python3 compare_rc5_rc6.py RC5_ROOT RC6_ROOT OUT_STEM"""
import glob
import json
import math
import os
import sys

R5, R6, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
SKIP = ('seconds', 'code', 'dxf_sha256', 'title', 'files', 'time', 'log', 'defaults_used', 'notes')
NEW_CHECK = 'every planned feature found on the drawing'


def load(root, n):
    p = os.path.join(root, 'runs', n, n + '.json')
    return json.load(open(p)) if os.path.exists(p) else None


def num_diffs(a, b, path=''):
    out = []
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b)):
            # rc6's new fields: the features-found check (E5), the stretch-end table domain and the strict refusal flag
            # the spec now carries for it (the rc5 verification's deferred item (1))
            if any(s in k for s in SKIP) or k in ('features_found', 'table_domain', 'exp_domain_refuse'):
                continue
            # rc6 moved the unused default of spec.flare_start (67 -> 65.4); records without a flare carry it
            if k == 'flare_start' and path.endswith('.spec') and a.get('flare_deg') is None and b.get('flare_deg') is None:
                continue
            if k not in a or k not in b:
                out.append((path + '.' + k, 'only in ' + ('rc5' if k in a else 'rc6')))
            else:
                out += num_diffs(a[k], b[k], path + '.' + k)
    elif isinstance(a, list) and isinstance(b, list):
        if path.endswith('.checks'):
            a = [c for c in a if c.get('check') != NEW_CHECK]
            b = [c for c in b if c.get('check') != NEW_CHECK]
        if len(a) != len(b):
            return [(path, 'len %d vs %d' % (len(a), len(b)))]
        for i, (x, y) in enumerate(zip(a, b)):
            out += num_diffs(x, y, '%s[%d]' % (path, i))
    elif isinstance(a, bool) or isinstance(b, bool):
        if a != b:
            out.append((path, '%r vs %r' % (a, b)))
    elif isinstance(a, (int, float)) and isinstance(b, (int, float)):
        if not (a == b or (isinstance(a, float) and isinstance(b, float) and math.isnan(a) and math.isnan(b))):
            out.append((path, '%r vs %r' % (a, b)))
    elif a != b:
        out.append((path, '%s vs %s' % (str(a)[:80], str(b)[:80])))
    return out


def summary(r):
    if r is None:
        return None
    a = r.get('assessment', {})
    V = r.get('validation', {})
    rc = V.get('realised_cutoff', {}) or {}
    dp = V.get('drawn_parameters', {}) or {}
    d = r.get('design', {}) or {}
    return dict(verdict=a.get('verdict'), exit=a.get('exit_code'), N=d.get('N'), NM=d.get('NM'),
                F_design=d.get('F'), F_drawn=(V.get('network') or {}).get('F'),
                D_range=[rc.get('D_realised_min_um'), rc.get('D_realised_max_um')],
                footprint_mm=dp.get('footprint_mm'), array_length_mm=dp.get('array_length_mm'),
                failing=[c['check'] for c in a.get('checks', []) if not c.get('ok') and c.get('kind') == 'hard'],
                new_check=[c['ok'] for c in a.get('checks', []) if c.get('check') == NEW_CHECK])


names = sorted({os.path.basename(os.path.dirname(p)) for p in glob.glob(os.path.join(R6, 'runs', '*', '*.json'))
                if os.path.basename(p) == os.path.basename(os.path.dirname(p)) + '.json'})
rows = []
for n in names:
    a, b = load(R5, n), load(R6, n)
    ia = json.dumps({k: v for k, v in (a or {}).get('input', {}).items() if k != 'paper'}, sort_keys=True)
    ib = json.dumps({k: v for k, v in (b or {}).get('input', {}).items() if k != 'paper'}, sort_keys=True)
    input_changed = ia != ib
    paper_changed = (a or {}).get('input', {}).get('paper') != (b or {}).get('input', {}).get('paper')
    d = num_diffs(a, b) if (a and b) else [('record', 'missing in rc5' if b else 'missing in rc6')]
    # input differences are expected to show; strip the input block itself from the count
    d_out = [x for x in d if not x[0].startswith('.input')]
    rows.append(dict(name=n, input_changed=input_changed, paper_changed=paper_changed, rc5=summary(a), rc6=summary(b),
                     n_diffs=len(d_out), first_diffs=d_out[:12],
                     status=('identical' if not d_out else ('changed (input changed)' if input_changed else
                                                            'CHANGED WITHOUT AN INPUT CHANGE'))))
json.dump(rows, open(OUT + '.json', 'w'), indent=1, default=str)
with open(OUT + '.md', 'w') as fh:
    fh.write('| spec | input changed | status | rc5 verdict (exit) | rc6 verdict (exit) | N rc5 → rc6 | footprint mm rc5 → rc6 | new check |\n|---|---|---|---|---|---|---|---|\n')
    for r in rows:
        s5, s6 = r['rc5'] or {}, r['rc6'] or {}
        fp = lambda s: ('%.1f × %.1f' % tuple(s['footprint_mm'][:2])) if s.get('footprint_mm') else '—'
        fh.write('| %s | %s | %s | %s (%s) | %s (%s) | %s → %s | %s → %s | %s |\n' % (
            r['name'], 'yes' if r['input_changed'] else ('notes' if r['paper_changed'] else 'no'), r['status'],
            (s5.get('verdict') or '')[:60], s5.get('exit'), (s6.get('verdict') or '')[:60], s6.get('exit'),
            s5.get('N'), s6.get('N'), fp(s5), fp(s6), s6.get('new_check')))
print('%d specs: %d identical, %d changed with an input change, %d changed without' % (
    len(rows), sum(r['status'] == 'identical' for r in rows), sum(r['status'].startswith('changed') for r in rows),
    sum(r['status'].startswith('CHANGED') for r in rows)))
for r in rows:
    if r['status'].startswith('CHANGED'):
        print('UNEXPLAINED', r['name'], r['first_diffs'][:5])
