"""RC8: the full rebuild against the frozen rc7 tree (code 52ba91eddb3f07b8).  For every spec built in both trees: the
photomask's sha256, verdict, exit code, N, NM, F and every leaf of the record.

RC8 changes five things a record can show (the adjudication of the rc7 reviews, item 10; register GIF4-1):
  1. every full build records `validation.drc.file_final` (the name the verdict gave the photomask) beside
     `validation.drc.file` (the name the DRC audited): a key only in the RC8 record, on every spec that ran the DRC;
  2. a floor start outside reproduce mode carries one more note (K13: its rows are declared in reproduce mode only):
     the one shipped spec of that kind is m14_paperclip_const_nomeander_floorstart19 (evaluate, exit 2);
  3. `gifford2014_fig4_aligned` is the nested reading (D1/D6): its input changed, so its photomask and numbers differ;
  4. the two arrays (tools/make_arrays.py) carry the new check `inter_copy_port_spacing` (the 2 mm inlet rule between
     copies, the neighbouring outlet ports reported: adjudication D4);
  5. the code stamp, the version, the tests' hash, the timings, the log (skipped, as rc7's comparison did).
Anything else that differs is UNEXPLAINED.  usage: python3 -I compare_rc7_rc8.py RC7_ROOT RC8_ROOT OUT_STEM"""
import glob
import hashlib
import json
import math
import os
import sys

R0, R1, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
SKIP = ('seconds', 'seconds_total', 'code', 'element_code', 'dxf_sha256', 'title', 'files', 'log', 'defaults_used',
        'environment', 'reassessed')
K13_SPECS = ('m14_paperclip_const_nomeander_floorstart19',)
ARRAYS = ('m07_serpentine_const_pack_filter_separate_x3', 'm13_paperclip_const_nomeander_sealed_fine_x4')


def load(root, n):
    p = os.path.join(root, 'runs', n, n + '.json')
    return json.load(open(p)) if os.path.exists(p) else None


def mask(root, n):
    for f in (n + '.dxf', n + '_FAILED.dxf'):
        p = os.path.join(root, 'runs', n, f)
        if os.path.exists(p):
            return f, hashlib.sha256(open(p, 'rb').read()).hexdigest()
    return None, None


def leaf_diffs(a, b, path=''):
    out = []
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b)):
            if k in SKIP:
                continue
            if k not in a or k not in b:
                out.append((path + '.' + k, 'only in ' + ('rc7' if k in a else 'rc8'), None if k in a else b[k]))
            else:
                out += leaf_diffs(a[k], b[k], path + '.' + k)
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            return [(path, 'len %d vs %d' % (len(a), len(b)), (a, b))]
        for i, (x, y) in enumerate(zip(a, b)):
            out += leaf_diffs(x, y, '%s[%d]' % (path, i))
    elif isinstance(a, bool) or isinstance(b, bool):
        if a != b:
            out.append((path, '%r vs %r' % (a, b), None))
    elif isinstance(a, (int, float)) and isinstance(b, (int, float)):
        if not (a == b or (isinstance(a, float) and isinstance(b, float) and math.isnan(a) and math.isnan(b))):
            out.append((path, '%r vs %r' % (a, b), None))
    elif a != b:
        out.append((path, '%s vs %s' % (str(a)[:80], str(b)[:80]), None))
    return out


def expected(name, d):
    """which of the RC8 changes explains a difference (None: unexplained)"""
    p, what, extra = d
    if p == '.validation.drc.file_final' and what == 'only in rc8':
        return 'drc.file_final (new)'
    if name in ARRAYS and p.endswith('.inter_copy_port_spacing') and what == 'only in rc8':
        return 'inter_copy_port_spacing (new array check)'
    if name in K13_SPECS and what.startswith('len ') and isinstance(extra, tuple):
        a, b = extra
        added = [x for x in b if x not in a]
        if len(b) == len(a) + 1 and all(x in b for x in a) and len(added) == 1 and ('reproduce mode' in str(added[0])
                                                                                   or 'gates.mode: reproduce' in str(added[0])):
            return 'K13 note: ' + str(added[0])[:160]
    return None


def summary(r):
    if r is None:
        return None
    a = r.get('assessment', {})
    d = r.get('design', {}) or {}
    return dict(verdict=a.get('verdict', r.get('verdict')), exit=a.get('exit_code', r.get('exit_code')), N=d.get('N'),
                NM=d.get('NM'), F_design=d.get('F'))


def names(root):
    return {os.path.basename(os.path.dirname(p)) for p in glob.glob(os.path.join(root, 'runs', '*', '*.json'))
            if os.path.basename(p) == os.path.basename(os.path.dirname(p)) + '.json'}


rows = []
for n in sorted(names(R0) | names(R1)):
    a, b = load(R0, n), load(R1, n)
    if not (a and b):
        rows.append(dict(name=n, status='only in the RC8 tree' if b else 'only in the rc7 tree', rc7=summary(a),
                         rc8=summary(b), dxf_identical=None, input_changed=None, n_diffs=None, first_diffs=[]))
        continue
    ia = json.dumps({k: v for k, v in a.get('input', {}).items() if k != 'paper'}, sort_keys=True)
    ib = json.dumps({k: v for k, v in b.get('input', {}).items() if k != 'paper'}, sort_keys=True)
    input_changed = ia != ib
    (f0, s0), (f1, s1) = mask(R0, n), mask(R1, n)
    d_all = [x for x in leaf_diffs(a, b) if not x[0].startswith('.input')]
    explained = [(x[0], expected(n, x)) for x in d_all if expected(n, x)]
    d = [x[:2] for x in d_all if not expected(n, x)]
    same = (s0 == s1) if (s0 and s1) else None
    if not d and same in (True, None):
        st = 'identical' + (' apart from: ' + '; '.join(sorted({e for _, e in explained})) if explained else '') + (
            ' (input changed: yes)' if input_changed else '')
    elif input_changed:
        st = 'changed (input changed)'
    else:
        st = 'CHANGED WITHOUT AN INPUT CHANGE'
    rows.append(dict(name=n, input_changed=input_changed, dxf_identical=same, dxf_names=[f0, f1], dxf_sha256=[s0, s1],
                     rc7=summary(a), rc8=summary(b), n_diffs=len(d), first_diffs=d[:12], explained=explained,
                     status=st))
json.dump(rows, open(OUT + '.json', 'w'), indent=1, default=str)
SUM = []
SUM.append('%d specs: %d identical, %d changed with an input change, %d changed without, %d only in the RC8 tree, %d only '
      'in the rc7 tree' % (len(rows), sum(r['status'].startswith('identical') for r in rows),
                           sum(r['status'].startswith('changed') for r in rows),
                           sum(r['status'].startswith('CHANGED') for r in rows),
                           sum(r['status'] == 'only in the RC8 tree' for r in rows),
                           sum(r['status'] == 'only in the rc7 tree' for r in rows)))
SUM.append('photomasks: %d byte-identical, %d differ; renamed: %s' % (
    sum(r.get('dxf_identical') is True for r in rows), sum(r.get('dxf_identical') is False for r in rows),
    ', '.join('%s (%s -> %s)' % (r['name'], *r['dxf_names']) for r in rows
              if r.get('dxf_names') and r['dxf_names'][0] != r['dxf_names'][1]) or 'none'))
SUM.append('drc.file_final added on %d specs; K13 note on: %s' % (
    sum(any(e == 'drc.file_final (new)' for _, e in r.get('explained', [])) for r in rows),
    ', '.join(r['name'] for r in rows if any(e.startswith('K13') for _, e in r.get('explained', []))) or 'none'))
for r in rows:
    if r['status'].startswith('CHANGED') or r.get('dxf_identical') is False:
        SUM.append('%s %s: dxf %s; %s' % ('UNEXPLAINED' if r['status'].startswith('CHANGED') else 'expected', r['name'],
                                     r.get('dxf_identical'), r['first_diffs'][:5]))

for l_ in SUM:
    print(l_)
with open(OUT + '.md', 'w') as fh:
    fh.write('# RC8 records compared with rc7 final (%s)\n\n' % ' '.join(sys.argv[1:3]))
    for l_ in SUM:
        fh.write('- %s\n' % l_)
    fh.write('\n')
    fh.write('| spec | input changed | photomask | name (rc7 → RC8) | status | rc7 verdict (exit) | RC8 verdict (exit) '
             '| N |\n|---|---|---|---|---|---|---|---|\n')
    for r in rows:
        x0, x1 = r['rc7'] or {}, r['rc8'] or {}
        nm = r.get('dxf_names') or [None, None]
        st = r['status'].replace('drc.file_final (new)', '`drc.file_final` (new)')
        fh.write('| %s | %s | %s | %s | %s | %s (%s) | %s (%s) | %s → %s |\n' % (
            r['name'], {True: 'yes', False: 'no', None: '—'}[r.get('input_changed')],
            {True: 'byte-identical', False: 'differs', None: '—'}[r.get('dxf_identical')],
            '%s → %s' % (nm[0], nm[1]) if nm[0] != nm[1] else 'same', st,
            (x0.get('verdict') or '')[:50], x0.get('exit'), (x1.get('verdict') or '')[:50], x1.get('exit'),
            x0.get('N'), x1.get('N')))
