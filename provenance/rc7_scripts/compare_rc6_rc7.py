"""rc6 -> rc7 record comparison (adapted from provenance/rc6_scripts/compare_rc5_rc6.py).  For every spec built in both
trees: the DXF's sha256, verdict, exit code, N, NM, F, the realised cutoff range, the footprint, the failing checks, and
every numeric leaf of the record.  rc7's options are off by default, so a spec whose input did not change must give a
byte-identical photomask and the same numbers, apart from the fields rc7 adds (status codes, source trace, environment,
passage report, the fluid's rheology and note, the square-corner record), the code stamp, timings and the
defaults_used text (C8 adds origins and effective values).  A spec whose input changed (rc7's Dinh depicted preset:
square turns) is expected to differ; one built only in rc7 (the new presets) is listed.
usage: python3 compare_rc6_rc7.py RC6_ROOT RC7_ROOT OUT_STEM"""
import glob
import hashlib
import json
import math
import os
import sys

R6, R7, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
# codereview2 N6: exact keys (rc6's script matched substrings, which also dropped exit_code, element_code and topology)
SKIP = ('seconds', 'seconds_total', 'code', 'element_code', 'dxf_sha256', 'title', 'files', 'log', 'defaults_used', 'notes', 'environment')
RC7_NEW = ('status_codes', 'source_trace', 'passage', 'square_corner', 'rheology',
           # rc7 code review: R6's plateau fields, R15's resistance note, R16's reproduction status in every mode
           'plateau_u_min_all_planned_rows', 'plateau_F_share_below_all_planned_rows', 'plateau_all_planned_rows_note',
           'R_element_SI_note',
           # rc7's engine Spec fields (each at its rc6-equivalent default when the option is off) and the fields added
           # after the code reviews (foundation check 2: the peak shear; codereview2/3: the f >= 1/2 row count)
           'leg_order', 'mu_f', 'rho_f', 'r_corner_in', 'r_corner_out', 'r_junction', 'shear_formula', 'turns',
           'shear_peak', 'rows_f_half_or_more')


def load(root, n):
    p = os.path.join(root, 'runs', n, n + '.json')
    return json.load(open(p)) if os.path.exists(p) else None


def sha(root, n):
    # rc7 final (register BUILD-2): a FAIL's photomask is <name>_FAILED.dxf (m14)
    for f in (n + '.dxf', n + '_FAILED.dxf'):
        p = os.path.join(root, 'runs', n, f)
        if os.path.exists(p):
            return hashlib.sha256(open(p, 'rb').read()).hexdigest()
    return None


def num_diffs(a, b, path=''):
    out = []
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b)):
            if k in SKIP or k in RC7_NEW:
                continue
            if k == 'note' and path.endswith('.fluid'):          # rc7 (C5) rewords the fluid note
                continue
            if k == 'reproduction' and path.endswith('.status') and isinstance(b.get(k), str) and \
                    b[k].startswith('NOT CLAIMED') and k not in a:   # rc7 (R16): the status field in every mode
                continue
            if k not in a or k not in b:
                out.append((path + '.' + k, 'only in ' + ('rc6' if k in a else 'rc7')))
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
    V = r.get('validation', {}) or {}
    rc = V.get('realised_cutoff', {}) or {}
    dp = V.get('drawn_parameters', {}) or {}
    d = r.get('design', {}) or {}
    return dict(verdict=a.get('verdict', r.get('verdict')), exit=a.get('exit_code', r.get('exit_code')), N=d.get('N'),
                NM=d.get('NM'), F_design=d.get('F'), F_drawn=(V.get('network') or {}).get('F_as_drawn'),
                D_range=[rc.get('D_realised_min_um'), rc.get('D_realised_max_um')],
                footprint_mm=dp.get('footprint_mm'),
                failing=[c['check'] for c in a.get('checks', []) if not c.get('ok') and c.get('kind') == 'hard'])


def names(root):
    return {os.path.basename(os.path.dirname(p)) for p in glob.glob(os.path.join(root, 'runs', '*', '*.json'))
            if os.path.basename(p) == os.path.basename(os.path.dirname(p)) + '.json'}


rows = []
for n in sorted(names(R6) | names(R7)):
    a, b = load(R6, n), load(R7, n)
    ia = json.dumps({k: v for k, v in (a or {}).get('input', {}).items() if k != 'paper'}, sort_keys=True)
    ib = json.dumps({k: v for k, v in (b or {}).get('input', {}).items() if k != 'paper'}, sort_keys=True)
    input_changed = ia != ib
    s6, s7 = sha(R6, n), sha(R7, n)
    if not (a and b):
        rows.append(dict(name=n, status='only in rc7' if b else 'only in rc6', rc6=summary(a), rc7=summary(b),
                         dxf_identical=None, input_changed=None, n_diffs=None, first_diffs=[]))
        continue
    d_all = [x for x in num_diffs(a, b) if not x[0].startswith('.input')]
    # rc7 final (code review 4, F4; LIMITATIONS section E): the DRC reads the photomask before its verdict names it,
    # so validation.drc.file names <name>_UNASSESSED.dxf; not a change of the photomask or of any number
    d = [x for x in d_all if not (x[0] == '.validation.drc.file' and x[1].partition(' vs ')[2].endswith('_UNASSESSED.dxf')
                                  and x[1].partition(' vs ')[2].replace('_UNASSESSED.dxf', '.dxf')
                                  == x[1].partition(' vs ')[0])]
    same_dxf = (s6 == s7) if (s6 and s7) else None
    if not d and same_dxf in (True, None):
        st = 'identical' + (' apart from the DRC\'s file name' if len(d_all) > len(d) else '')
    elif input_changed:
        st = 'changed (input changed)'
    else:
        st = 'CHANGED WITHOUT AN INPUT CHANGE'
    rows.append(dict(name=n, input_changed=input_changed, dxf_identical=same_dxf, rc6=summary(a), rc7=summary(b),
                     n_diffs=len(d), first_diffs=d[:12], status=st))
json.dump(rows, open(OUT + '.json', 'w'), indent=1, default=str)
with open(OUT + '.md', 'w') as fh:
    fh.write('| spec | input changed | photomask | status | rc6 verdict (exit) | rc7 verdict (exit) | N rc6 → rc7 |\n'
             '|---|---|---|---|---|---|---|\n')
    for r in rows:
        x6, x7 = r['rc6'] or {}, r['rc7'] or {}
        fh.write('| %s | %s | %s | %s | %s (%s) | %s (%s) | %s → %s |\n' % (
            r['name'], {True: 'yes', False: 'no', None: '—'}[r.get('input_changed')],
            {True: 'byte-identical', False: 'differs', None: '—'}[r.get('dxf_identical')], r['status'],
            (x6.get('verdict') or '')[:60], x6.get('exit'), (x7.get('verdict') or '')[:60], x7.get('exit'),
            x6.get('N'), x7.get('N')))
print('%d specs: %d identical, %d changed with an input change, %d changed without, %d only in rc7, %d only in rc6' % (
    len(rows), sum(r['status'].startswith('identical') for r in rows), sum(r['status'].startswith('changed') for r in rows),
    sum(r['status'].startswith('CHANGED') for r in rows), sum(r['status'] == 'only in rc7' for r in rows),
    sum(r['status'] == 'only in rc6' for r in rows)))
print('photomasks: %d byte-identical, %d differ' % (sum(r.get('dxf_identical') is True for r in rows),
                                                    sum(r.get('dxf_identical') is False for r in rows)))
for r in rows:
    if r['status'].startswith('CHANGED') or r.get('dxf_identical') is False:
        print('%s %s: dxf %s; %s' % ('UNEXPLAINED' if r['status'].startswith('CHANGED') else 'expected', r['name'],
                                     r.get('dxf_identical'), r['first_diffs'][:5]))
