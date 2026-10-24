"""rc7: the final rebuild against the first rc7 freeze (code 4a624626f570918c, 2026-10-19 14:23), after the clearance
guarantee (register BUILD-2).  For every spec built in both trees: the photomask's sha256 (a FAIL's photomask is
<name>_FAILED.dxf in the final tree), verdict, exit code, N, NM, F and every numeric leaf of the record.  The guarantee
adds checks before writing and refusals; it changes no design that was clear, so a spec whose input did not change must
give a byte-identical photomask and the same numbers, apart from the code stamp, the timings, the log and the
defaults_used text.  A spec whose input changed is expected to differ only where its input does: the Dinh
figure-aligned refit (register INV-1, D 3.67 um) differs; the three Xia presets that depicted the head lost a value the
head replaced (meander.rectilinear_rows: 9), so their photomasks must not change.
usage: python3 compare_frozen_final.py FROZEN_ROOT FINAL_ROOT OUT_STEM"""
import glob
import hashlib
import json
import math
import os
import sys

R0, R1, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
SKIP = ('seconds', 'seconds_total', 'code', 'element_code', 'dxf_sha256', 'title', 'files', 'log', 'defaults_used',
        'environment', 'reassessed')


def drc_name(d):
    """the one expected record change of the guarantee on a design that was clear: the DRC reads the photomask before
    its verdict exists, so `validation.drc.file` names <name>_UNASSESSED.dxf, the file name_mask then renames (code
    review 4, F4; LIMITATIONS section E).  Not a change of the photomask or of any number."""
    p, v = d
    if p != '.validation.drc.file':
        return False
    a, _, b = v.partition(' vs ')
    return b == a[:-4] + '_UNASSESSED.dxf' or b == a[:-4].replace('_FAILED', '') + '_UNASSESSED.dxf'


def load(root, n):
    p = os.path.join(root, 'runs', n, n + '.json')
    return json.load(open(p)) if os.path.exists(p) else None


def mask(root, n):
    for f in (n + '.dxf', n + '_FAILED.dxf'):
        p = os.path.join(root, 'runs', n, f)
        if os.path.exists(p):
            return f, hashlib.sha256(open(p, 'rb').read()).hexdigest()
    return None, None


def num_diffs(a, b, path=''):
    out = []
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b)):
            if k in SKIP:
                continue
            if k not in a or k not in b:
                out.append((path + '.' + k, 'only in ' + ('frozen' if k in a else 'final')))
            else:
                out += num_diffs(a[k], b[k], path + '.' + k)
    elif isinstance(a, list) and isinstance(b, list):
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
        rows.append(dict(name=n, status='only in the final tree' if b else 'only in the frozen tree', frozen=summary(a),
                         final=summary(b), dxf_identical=None, input_changed=None, n_diffs=None, first_diffs=[]))
        continue
    ia = json.dumps({k: v for k, v in a.get('input', {}).items() if k != 'paper'}, sort_keys=True)
    ib = json.dumps({k: v for k, v in b.get('input', {}).items() if k != 'paper'}, sort_keys=True)
    input_changed = ia != ib
    (f0, s0), (f1, s1) = mask(R0, n), mask(R1, n)
    d_all = [x for x in num_diffs(a, b) if not x[0].startswith('.input')]
    d = [x for x in d_all if not drc_name(x)]
    drc_only = len(d_all) > len(d)
    same = (s0 == s1) if (s0 and s1) else None
    if not d and same in (True, None):
        st = 'identical' + (' apart from the DRC\'s file name' if drc_only else '') + (
            ' (input changed: yes)' if input_changed else '')
    elif input_changed:
        st = 'changed (input changed)'
    else:
        st = 'CHANGED WITHOUT AN INPUT CHANGE'
    rows.append(dict(name=n, input_changed=input_changed, dxf_identical=same, dxf_names=[f0, f1], frozen=summary(a),
                     final=summary(b), n_diffs=len(d), first_diffs=d[:12], status=st))
json.dump(rows, open(OUT + '.json', 'w'), indent=1, default=str)
with open(OUT + '.md', 'w') as fh:
    fh.write('| spec | input changed | photomask | name (frozen → final) | status | frozen verdict (exit) | final verdict (exit) '
             '| N |\n|---|---|---|---|---|---|---|---|\n')
    for r in rows:
        x0, x1 = r['frozen'] or {}, r['final'] or {}
        nm = r.get('dxf_names') or [None, None]
        fh.write('| %s | %s | %s | %s | %s | %s (%s) | %s (%s) | %s → %s |\n' % (
            r['name'], {True: 'yes', False: 'no', None: '—'}[r.get('input_changed')],
            {True: 'byte-identical', False: 'differs', None: '—'}[r.get('dxf_identical')],
            '%s → %s' % (nm[0], nm[1]) if nm[0] != nm[1] else 'same', r['status'],
            (x0.get('verdict') or '')[:50], x0.get('exit'), (x1.get('verdict') or '')[:50], x1.get('exit'),
            x0.get('N'), x1.get('N')))
print('%d specs: %d identical, %d changed with an input change, %d changed without, %d only in the final tree, %d only '
      'in the frozen tree' % (len(rows), sum(r['status'].startswith('identical') for r in rows),
                             sum(r['status'].startswith('changed') for r in rows),
                             sum(r['status'].startswith('CHANGED') for r in rows),
                             sum(r['status'] == 'only in the final tree' for r in rows),
                             sum(r['status'] == 'only in the frozen tree' for r in rows)))
print('photomasks: %d byte-identical, %d differ; renamed: %s' % (
    sum(r.get('dxf_identical') is True for r in rows), sum(r.get('dxf_identical') is False for r in rows),
    ', '.join('%s (%s -> %s)' % (r['name'], *r['dxf_names']) for r in rows
              if r.get('dxf_names') and r['dxf_names'][0] != r['dxf_names'][1]) or 'none'))
for r in rows:
    if r['status'].startswith('CHANGED') or r.get('dxf_identical') is False:
        print('%s %s: dxf %s; %s' % ('UNEXPLAINED' if r['status'].startswith('CHANGED') else 'expected', r['name'],
                                     r.get('dxf_identical'), r['first_diffs'][:5]))
