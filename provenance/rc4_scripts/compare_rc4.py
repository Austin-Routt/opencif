"""Compare every rc4 record (opencif_rc4/runs) with its rc3 record (opencif_next/runs): design, measured network, cutoff,
Q1, verdict, exit code, checks and notes.  Prints a per-record line and a summary."""
import glob, json, os, sys

import sys as _sys
A = _sys.argv[1] if len(_sys.argv) > 2 else 'opencif_next/runs'   # rc3 runs/
B = _sys.argv[2] if len(_sys.argv) > 2 else 'opencif_rc4/runs'    # rc4 runs/


def load(root, n):
    p = os.path.join(root, n, n + '.json')
    return json.load(open(p)) if os.path.exists(p) else None


names = sorted(os.path.basename(os.path.dirname(p)) for p in glob.glob(os.path.join(A, '*', '*.json'))
               if not p.endswith('_audit.json') and os.path.basename(os.path.dirname(p)) + '.json' == os.path.basename(p))
hashes, exits, issues = {}, {}, []
rows = []
for n in names:
    a, b = load(A, n), load(B, n)
    if b is None:
        issues.append((n, 'NO rc4 RECORD'))
        continue
    hashes[n] = (b.get('code') or {}).get('sha256')
    if 'element' in a:                                   # arrays
        diff = []
        for k in ('count', 'pitch_um', 'inlet'):
            if a.get(k) != b.get(k):
                diff.append((k, a.get(k), b.get(k)))
        ca = {k: v.get('ok') for k, v in a['checks'].items() if isinstance(v, dict)}
        cb = {k: v.get('ok') for k, v in b['checks'].items() if isinstance(v, dict)}
        if ca != cb:
            diff.append(('checks', ca, cb))
        for k in ('verdict', 'element_verdict', 'element_code'):     # rc4 check pass 2: verdicts compared too
            if a.get(k) != b.get(k):
                diff.append((k, str(a.get(k))[-140:], str(b.get(k))[-140:]))
        rows.append((n, a['exit_code'], b['exit_code'], diff, a['verdict'][:60], b['verdict'][:60]))
        exits[n] = (a['exit_code'], b['exit_code'])
        continue
    aa, ab = a['assessment'], b['assessment']
    exits[n] = (aa['exit_code'], ab['exit_code'])
    diff = []
    if 'refused' in a:
        for k in sorted(set(a['design']) | set(b['design'])):       # rc4 check pass 2: the whole record compared
            if a['design'].get(k) != b['design'].get(k):
                diff.append(('design.' + k, str(a['design'].get(k))[:80], str(b['design'].get(k))[:80]))
        if aa['verdict'] != ab['verdict']:
            diff.append(('verdict', aa['verdict'], ab['verdict']))
        ra_, rb_ = a.get('refused') or {}, b.get('refused') or {}
        for k in ('reasons', 'advice', 'gates_note'):
            if ra_.get(k) != rb_.get(k):
                diff.append(('refused.' + k, str(ra_.get(k))[:300], str(rb_.get(k))[:300]))
        ga_, gb_ = ra_.get('gates') or {}, rb_.get('gates') or {}
        for k in sorted(set(ga_) | set(gb_)):
            if ga_.get(k) != gb_.get(k):
                diff.append(('refused.gates.' + k, str(ga_.get(k))[:100], str(gb_.get(k))[:100]))
        for k in sorted(set(a) | set(b)):
            if k not in ('design', 'refused', 'assessment', 'code', 'created', 'timestamp') and a.get(k) != b.get(k):
                diff.append(('top.' + k, str(a.get(k))[:120], str(b.get(k))[:120]))
        rows.append((n, aa['exit_code'], ab['exit_code'], diff, aa['verdict'][:80], ab['verdict'][:80]))
        continue
    da, db = a['design'], b['design']
    for k in ('N', 'NM', 'layout', 'n_legs'):
        if da.get(k) != db.get(k):
            diff.append((k, da.get(k), db.get(k)))
    Va, Vb = a['validation'], b['validation']
    for k in ('F_as_drawn', 'F_intended'):
        x, y = Va['network'].get(k), Vb['network'].get(k)
        if x is None or y is None or abs(x - y) > 1e-12:
            diff.append((k, x, y))
    ra, rb = Va['realised_cutoff'], Vb['realised_cutoff']
    for k in ('D_realised_min_um', 'D_realised_max_um', 'share_of_F_within_1%'):
        if abs(ra[k] - rb[k]) > 1e-12:
            diff.append((k, ra[k], rb[k]))
    for k in ('undeclared_rows_outside', 'declared_rows'):
        if ra.get(k) != rb.get(k):
            diff.append((k, ra.get(k), rb.get(k)))
    if abs(Va['drawing_fidelity']['cutoff_dev_max'] - Vb['drawing_fidelity']['cutoff_dev_max']) > 1e-12:
        diff.append(('Q1', Va['drawing_fidelity']['cutoff_dev_max'], Vb['drawing_fidelity']['cutoff_dev_max']))
    if Va['published'] != Vb.get('drawn_parameters'):
        diff.append(('published->drawn_parameters', Va['published'], Vb.get('drawn_parameters')))
    if Va.get('robustness_design_network') != Vb.get('robustness_design_network'):
        diff.append(('robustness', 'differs', ''))
    msa, msb = Va.get('model_scenarios') or {}, Vb.get('model_scenarios') or {}
    if {k: msb.get(k) for k in msa} != msa:
        diff.append(('model_scenarios (rc3 keys)', 'differs', ''))
    ga, gb = Va['gates'], Vb['gates']
    for k in ga:
        if k in ('cutoff_band_um',):
            continue
        if ga[k] != gb.get(k):
            diff.append(('gates.' + k, ga[k], gb.get(k)))
    if ga['cutoff_band_um'].get('band_um') != gb['cutoff_band_um'].get('band_um'):
        diff.append(('band', ga['cutoff_band_um'].get('band_um'), gb['cutoff_band_um'].get('band_um')))
    ca = {c['check']: (c['ok'], c['kind']) for c in aa['checks']}
    cb = {c['check']: (c['ok'], c['kind']) for c in ab['checks']}
    mode_b = ab.get('mode')
    ar_b = ('aspect ratio of the narrowest CIF channel (gap openings, side channels, U-turn sides, head units; %s: '
            'LIMITATIONS 7)' % ('a pass condition in strict mode' if mode_b == 'strict' else
                                'reported, not a pass condition in %s mode' % mode_b))
    ren = {'aspect ratio of the narrowest CIF channel (gap openings, side channels, U-turn sides, head units; reported, '
           'not a pass condition: LIMITATIONS 7)': ar_b}
    ca2 = {ren.get(k, k): v for k, v in ca.items()}
    for k in set(ca2) | set(cb):
        va, vb = ca2.get(k), cb.get(k)
        if va is None or vb is None or va[0] != vb[0]:
            diff.append(('check ' + k[:40], va, vb))
    vals_a = {ren.get(c['check'], c['check']): c['value'] for c in aa['checks']}
    vals_b = {c['check']: c['value'] for c in ab['checks']}
    for k in vals_b:
        if k in vals_a and vals_a[k] != vals_b[k] and not k.startswith('validity envelope'):
            diff.append(('value ' + k[:40], str(vals_a[k])[:60], str(vals_b[k])[:60]))
    na, nb = set(aa['notes']), set(ab['notes'])
    new_notes = [x[:110] for x in nb - na]
    gone_notes = [x[:110] for x in na - nb]
    rows.append((n, aa['exit_code'], ab['exit_code'], diff, new_notes, gone_notes))
    if (a.get('files') or {}).get('dxf_sha256') != (b.get('files') or {}).get('dxf_sha256'):
        diff.append(('dxf_sha256', (a.get('files') or {}).get('dxf_sha256'), (b.get('files') or {}).get('dxf_sha256')))
    if aa['verdict'] != ab['verdict']:
        diff.append(('verdict', aa['verdict'][-90:], ab['verdict'][-120:]))
    if not (b.get('files') or {}).get('dxf_sha256'):
        issues.append((n, 'no files.dxf_sha256'))

for r in rows:
    print('%-50s exit %s -> %s  diffs %d' % (r[0], r[1], r[2], len(r[3])))
    for d in r[3]:
        print('      ', d)
    if len(r) > 4 and isinstance(r[4], list):
        for x in r[4]:
            print('      + note:', x)
        for x in r[5]:
            print('      - note:', x)
print('HASHES', sorted(set(hashes.values())))
print('EXIT CHANGES', {n: e for n, e in exits.items() if e[0] != e[1]})
print('ISSUES', issues)
