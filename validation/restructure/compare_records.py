"""Compare each 1.0.0 build record with its RC8 record: everything must be equal except what 1.0.0 adds or renames and
what differs between runs (timings, the code stamp, the environment, the log). What 1.0.0 adds or renames is named
below with the CHANGELOG item that makes it, and is set aside before the comparison; any other difference is reported.
One record changes for a named reason (CHANGED below): its difference is checked field by field and reported apart.
usage: python -I compare_records.py RC8_RUNS NEW_RUNS OUT.json"""
import json
import os
import re
import sys

RC8, NEW, OUT = sys.argv[1:4]
IGNORE_TOP = {'seconds', 'seconds_total', 'code', 'environment', 'log', 'title',
              'element_code'}                     # an array's element stamp: version, code stamp and file count
OLD_SRC = 'CIF_v2 v2.1.0 Pairing (byte-frozen, opencif/reference); Hill fit: glossary macros'
NEW_SRC = 'opencif.physics.kernels.Pairing (= CIF_v2 v2.1.0 Pairing, reference/CIF_v2.py); Hill fit: glossary macros'
# A record that 1.0.0 changes for a named reason: (record, path) -> (RC8 value, 1.0.0 value), with the reason.
CHANGED = {
    'm14_paperclip_const_nomeander_floorstart19': dict(
        reason='CHANGELOG 1.0.0 item 3: the physics audit\'s M14 (the design code\'s Reynolds and Stokes gates on the '
               'drawn mask) belongs to the validity envelope; the floor start\'s undeclared over-extracting rows exceed '
               'the Stokes gate (Re_gap 0.098 and Re_delta 0.18 against 0.01), so the record, which already FAILs Q2 '
               '(exit 2, unchanged), is also outside the envelope, and its verdict names the audit\'s failure (P5)',
        fields={'m14_paperclip_const_nomeander_floorstart19/assessment/status/validity_envelope': ('INSIDE', 'OUTSIDE'),
                'm14_paperclip_const_nomeander_floorstart19/assessment/status_codes/validity_envelope': ('INSIDE',
                                                                                                         'OUTSIDE'),
                'm14_paperclip_const_nomeander_floorstart19/assessment/verdict': (
                    'FAIL (Q2 cutoff D at every gap (undeclared rows))',
                    'FAIL (Q2 cutoff D at every gap (undeclared rows); physics audit of the photomask file (P0-P10))')})}
RENAMES = [(r'\bcifsched\b', 'schedule'), (r'\bmeasure2\b', 'measure')]


def norm_str(s):
    for a, b in RENAMES:
        s = re.sub(a, b, s)
    return s


def walk(a, b, path, out):
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b)):
            p = path + '/' + str(k)
            if k in ('seconds', 'time_s', 'elapsed_s', 't_s'):
                continue
            if k not in a:
                out.append((p, 'only in 1.0.0'))
            elif k not in b:
                out.append((p, 'only in RC8'))
            else:
                walk(a[k], b[k], p, out)
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            out.append((path, 'length %d vs %d' % (len(a), len(b))))
            return
        for i, (x, y) in enumerate(zip(a, b)):
            walk(x, y, '%s[%d]' % (path, i), out)
    elif isinstance(a, str) and isinstance(b, str):
        if norm_str(a) != norm_str(b):
            out.append((path, '%r vs %r' % (a[:160], b[:160])))
    elif a != b and not (isinstance(a, float) and isinstance(b, float) and a != a and b != b):
        out.append((path, '%r vs %r' % (str(a)[:80], str(b)[:80])))


res = {}
for name in sorted(os.listdir(NEW)):
    pn, po = os.path.join(NEW, name, name + '.json'), os.path.join(RC8, name, name + '.json')
    if not (os.path.exists(pn) and os.path.exists(po)):
        continue
    a, b = json.load(open(po)), json.load(open(pn))
    a = {k: v for k, v in a.items() if k not in IGNORE_TOP}
    b = {k: v for k, v in b.items() if k not in IGNORE_TOP}
    expected = []
    V = b.get('validation') or {}
    if 'physics_audit' in V:
        V.pop('physics_audit')
        expected.append('validation.physics_audit (new; CHANGELOG item 3)')
    if 'design_code_reports' in V:
        V.pop('design_code_reports')
        expected.append('validation.design_code_reports (new; item 4)')
    if 'element_physics_audit' in b:
        b.pop('element_physics_audit')
        expected.append('element_physics_audit of an array (new; item 3)')
    for blk, src in ((a, OLD_SRC), (b, NEW_SRC)):     # the Pairing's new home (items 1 and 2)
        kp = (blk.get('design') or {}).get('kernel_pair')
        if isinstance(kp, dict) and kp.get('source') == src:
            kp['source'] = 'Pairing'
            if blk is b:
                expected.append('design.kernel_pair.source names the Pairing\'s new module (items 1 and 2)')
    for blk in (a, b):
        g = (blk.get('validation') or {}).get('gates') or (blk.get('refused') or {}).get('gates') or {}
        if isinstance(g.get('P8'), dict):
            g['P8'].pop('status', None)                    # the text that said the auditor was not run
        asm = blk.get('assessment') or {}
        for k in ('status', 'status_codes'):
            if isinstance(asm.get(k), dict):
                asm[k].pop('physics_auditor_P8', None)
                asm[k].pop('physics_audit', None)
        if isinstance(asm.get('checks'), list):          # the audit's rows and the design code's reports' row are new
            asm['checks'] = [c for c in asm['checks'] if not str(c.get('check', '')).startswith(
                ('physics audit', "the design code's reports"))]
        drc = (blk.get('validation') or {}).get('drc')
        if isinstance(drc, dict):
            for k in ('seconds', 'time_s', 't_tiles', 't_total'):
                drc.pop(k, None)
    out = []
    walk(a, b, name, out)
    named = []
    if name in CHANGED:
        want = CHANGED[name]['fields']
        for p_, d_ in list(out):
            if p_ in want and d_ == '%r vs %r' % want[p_]:
                named.append((p_, d_))
                out.remove((p_, d_))
        assert len(named) == len(want), (name, named)   # every named change is there, as named
    res[name] = dict(differences=out, expected=sorted(set(expected)), changed_for_named_reason=named,
                     reason=CHANGED.get(name, {}).get('reason'))
    print('%-52s %4d unexplained differences%s%s' % (name, len(out), (': ' + '; '.join('%s %s' % d for d in out[:3]))
                                                     if out else '', ' (+%d named)' % len(named) if named else ''))
summary = dict(records_compared=len(res),
               identical_apart_from_named_additions=sorted(n for n, v in res.items()
                                                           if not v['differences'] and not v['changed_for_named_reason']),
               changed_for_named_reason={n: v['reason'] for n, v in res.items() if v['changed_for_named_reason']},
               unexplained={n: len(v['differences']) for n, v in res.items() if v['differences']},
               named_additions=sorted({e for v in res.values() for e in v['expected']}))
json.dump(dict(summary=summary, records=res), open(OUT, 'w'), indent=1)
print('records compared: %d; identical apart from the named additions: %d; changed for a named reason: %d; '
      'unexplained: %d' % (len(res), len(summary['identical_apart_from_named_additions']),
                           len(summary['changed_for_named_reason']), len(summary['unexplained'])))
