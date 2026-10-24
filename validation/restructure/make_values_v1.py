"""make_values_v1.py: write the values that fill the placeholders of the release's documents (tools_v1/fill_docs.py),
each computed from the evidence, so that no number in the documents is typed by hand.

Usage: python -I make_values_v1.py TREE STANDALONE_AUDITS_DIR SUITES_DIR OUT.json [--keep KEY,KEY,...]

TREE                   the final tree (validation/RELEASE.json, validation/audit_table.csv, runs/)
STANDALONE_AUDITS_DIR  the standalone audits (<name>.json), as make_validation.py read them
SUITES_DIR             the suites' logs (<suite>.log) and summary.txt (one line per run: time, suite, rc, seconds, result)
--keep                 placeholders to keep as they are for a later pass (their value is the placeholder itself)

The code stamp, the tests stamp, the Python version and the line counts are left to fill_docs.py, which computes them
from the tree. Each value below names its source in a comment.
"""
import csv
import glob
import json
import os
import re
import sys

TREE, AUD, SUITES, OUT = sys.argv[1:5]
KEEP = set(sys.argv[sys.argv.index('--keep') + 1].split(',')) if '--keep' in sys.argv else set()
R = json.load(open(os.path.join(TREE, 'validation', 'RELEASE.json')))
C = R['counts']
PM = R['photomasks']
V = {}


def fmt(x, p=2):
    return ('%.' + str(p) + 'g') % x


def sci(x):
    """3.7e-08 -> 3.7 × 10⁻⁸"""
    m, e = ('%.1e' % x).split('e')
    sup = str.maketrans('-0123456789', '⁻⁰¹²³⁴⁵⁶⁷⁸⁹')
    return '%s × 10%s' % (m, str(int(e)).translate(sup))


def recs():
    out = {}
    for name in PM:
        p = os.path.join(TREE, 'runs', name, name + '.json')
        if os.path.exists(p):
            out[name] = json.load(open(p))
    return out


REC = recs()


def checks(rec):
    pa = (rec.get('validation') or {}).get('physics_audit') or {}
    return {c['id']: c for c in pa.get('checks', [])}


# ---- the plans, the photomasks, the records, the exit codes (validation/RELEASE.json)
V['N_SPECS'] = str(C['plans_compared'])
V['PLAN_RESULT'] = ('%d of %d identical in every compared field; the two fields compared apart are the planning time and '
                    'the text of gate P8\'s status, which 1.0.0 rewrites on all %d because it now audits every photomask'
                    % (C['plans_identical'], C['plans_compared'], C['plans_p8_status_changed']))
assert not C['plans_different'] and not C['plans_missing']
drawn = [n for n, v in PM.items() if v['dxf']]
V['N_BUILDS'] = str(C['builds'])
V['N_MASKS'] = str(C['drawn'])
V['N_ARRAYS'] = str(C['arrays'])
same = sum(1 for n in drawn if PM[n]['identical_to_rc8'])
refused = [n for n, v in PM.items() if v['exit_code'] == 5]
V['MASK_RESULT'] = ('%d of %d photomasks and %d of %d arrays identical to RC8\'s, byte for byte (sha256 in '
                    '`validation/photomasks_sha256.txt`); the refused design (%s) draws no photomask in either release'
                    % (same, len(drawn), C['arrays_identical_to_rc8'], C['arrays'], ', '.join('`%s`' % n for n in refused)))
V['MASK_RESULT_SHORT'] = '%d of %d photomasks and %d of %d arrays identical' % (same, len(drawn),
                                                                               C['arrays_identical_to_rc8'], C['arrays'])
rc = C['records']
changed = rc['changed_for_named_reason']
V['RECORD_RESULT'] = ('%d of %d identical apart from the named additions; %d changed for a named reason (%s: its validity '
                      'envelope goes from inside to outside, because the audit\'s M14 on the drawn mask finds the floor '
                      'start\'s undeclared rows above the Stokes gate, CHANGELOG item 3; its exit code, 2, is unchanged); '
                      '%d unexplained (`validation/RELEASE.json`, counts.records)'
                      % (len(rc['identical_apart_from_named_additions']), rc['records_compared'], len(changed),
                         ', '.join('`%s`' % n for n in changed), len(rc['unexplained'])))
assert not rc['unexplained']
ex = {}
for n, v in PM.items():
    ex.setdefault(v['exit_code'], []).append(n)
V['EXIT_RESULT'] = ('%d of %d as in RC8 (%s)' % (C['exit_codes_as_rc8'], C['builds'],
                                                 ', '.join('exit %s: %d' % (k, len(ex[k])) for k in sorted(ex))))

# ---- the audits (validation/RELEASE.json, the records' physics_audit, the standalone audits)
st = {}
for n in drawn:
    st.setdefault(PM[n]['physics_audit'], []).append(n)
same_st = sum(1 for n in drawn if PM[n]['physics_audit'] == PM[n]['physics_audit_standalone'])
outside = sorted(st.get('PASS_OUTSIDE_ENVELOPE', []))
failed = sorted(st.get('FAIL', []))
env_why = {}
for n in outside + failed:
    pa = REC[n]['validation']['physics_audit']
    env_why[n] = (pa.get('hard_failures') or [], pa.get('envelope_failures') or [])
V['AUDIT_BUILD'] = ('%d PASS, %d PASS outside the envelope (%s: M14 on the drawn mask), and %d FAIL (%s: %s, the design '
                    'that already fails Q2); the refused design has no photomask to audit'
                    % (len(st.get('PASS', [])), len(outside), ', '.join('`%s`' % n for n in outside), len(failed),
                       ', '.join('`%s`' % n for n in failed),
                       '; '.join('hard %s, envelope %s' % (','.join(h) or 'none', ','.join(e) or 'none')
                                 for h, e in (env_why[n] for n in failed))))
assert all(e == ['M14'] or set(e) == {'M14'} for h, e in (env_why[n] for n in outside)), env_why
V['AUDIT_STANDALONE'] = ('the same status as in the build for %d of %d photomasks (%d PASS, %d PASS outside the envelope, '
                         '%d FAIL)' % (same_st, len(drawn), len(st.get('PASS', [])), len(outside), len(failed)))
assert same_st == len(drawn)
# R0: the two rulers on the straight rows; P2: the two solvers on every device
r0_rel, r0_pos, r0_seg, r0_gap, p2 = 0.0, 0.0, 0, 0, []
for n in drawn:
    ck = checks(REC[n])
    v = ck['R0'].get('value') or {}
    r0_rel = max(r0_rel, float(v.get('max_rel_diff') or 0.0))
    r0_pos = max(r0_pos, float(v.get('gap_position_max_dev_um') or 0.0))
    r0_seg += int(v.get('segments') or 0)
    r0_gap += int(v.get('gaps') or 0)
    p2.append((float((ck['P2'].get('value') or {}).get('solver_agreement') or 0.0), n))
V['R0_SUMMARY'] = ('a relative difference of at most %s in the segment resistances and %s µm in the gap positions, over '
                   '%s straight segments and %s gaps' % (sci(r0_rel) if r0_rel else '0', sci(r0_pos) if r0_pos else '0',
                                                        '{:,}'.format(r0_seg), '{:,}'.format(r0_gap)))
pmax = max(p2)
V['PAR_DEVICES'] = '%s (the largest, on `%s`)' % (sci(pmax[0]), pmax[1])
# the strict designs: their largest audit values
strict = [n for n in drawn if ((REC[n].get('input') or {}).get('gates') or {}).get('mode', 'strict') == 'strict']
mx = dict(p4=(0.0, ''), p5hi=(0.0, ''), p5lo=(9.0, ''), p3=(0.0, ''), p6=(0.0, ''), r6=(0.0, ''), re_d=(0.0, ''))
r6um = 0.0
for n in strict:
    ck = checks(REC[n])
    v4, v5, v3, v6 = (ck[k].get('value') or {} for k in ('P4', 'P5', 'P3', 'P6'))
    vr, vm = ck['R6b'].get('value') or {}, ck['M14'].get('value') or {}
    r6um = max(r6um, float(vr.get('max_dev_along_um') or 0.0), float(vr.get('max_dev_across_um') or 0.0))
    for key, val, better in (('p4', v4.get('undeclared_max'), max), ('p5hi', v5.get('max_delta_over_delta_star'), max),
                             ('p5lo', v5.get('min_delta_over_delta_star'), min),
                             ('p3', v3.get('exact_closure_max_rel_dev'), max),
                             ('p6', abs((v6.get('F_drawn') or 0) - (v6.get('F_target') or 0)), max),
                             ('r6', vr.get('max_dev_ratio'), max), ('re_d', vm.get('Re_delta_max'), max)):
        if val is None:
            continue
        if better(val, mx[key][0]) == val and val != mx[key][0]:
            mx[key] = (float(val), n)
V['AUDIT_STRICT'] = ('over the %d strict designs: the drawn fraction of an undeclared gap within %.2f%% of the intended '
                     'one (P4, `%s`); the lamina δ/δ* from %.4f to %.4f (P5); the exact closure within %s (P3); the '
                     'drawn F within %s of its target (P6); %s; and the lamina '
                     'Reynolds number at most %s (M14, against 10⁻²)'
                     % (len(strict), 100 * mx['p4'][0], mx['p4'][1], mx['p5lo'][0], mx['p5hi'][0], sci(mx['p3'][0]),
                        sci(mx['p6'][0]), ('every post at its planned size, to the vertex (R6b)' if r6um == 0 else
                                           'every post within %s µm of its planned size (R6b)' % sci(r6um)),
                        sci(mx['re_d'][0])))
VARIANTS = ('_fold20', '_pathorder', '_start51')                  # the Fig. 4 variants (provenance/rc7_traces/variants)
presets = sorted(n for n in drawn if n.startswith(('gifford2014', 'xia2016', 'dinh2024')))
n_var = sum(1 for n in presets if n.endswith(VARIANTS))
pst = {}
for n in presets:
    pst.setdefault(PM[n]['physics_audit'], []).append(n)
LABEL = {'PASS': 'PASS', 'PASS_OUTSIDE_ENVELOPE': 'PASS outside the envelope', 'FAIL': 'FAIL'}
V['AUDIT_REPRO'] = ('%d presets and the %d Fig. 4 variants: %s; their declared departures are reported row by row, and a '
                    'reproduction outside the envelope is the paper\'s design, not opencif\'s choice '
                    '(`docs/REPRODUCTIONS.md`)'
                    % (len(presets) - n_var, n_var, ', '.join('%d %s%s' % (len(v), LABEL.get(k, k),
                                                             ' (%s)' % ', '.join('`%s`' % x for x in v)
                                                             if k != 'PASS' else '') for k, v in sorted(pst.items()))))
V['AUDIT_SUMMARY'] = ('The %d photomasks give %d PASS, %d PASS outside the envelope and %d FAIL, in the build and again '
                      'from scratch, with the same status for every photomask, and the refused design has no photomask '
                      '(`validation/audit_table.csv`).' % (len(drawn), len(st.get('PASS', [])), len(outside), len(failed)))
# ---- the design code's reports on the shipped designs (the records' validation.design_code_reports)
op5 = st5 = fab = n_d = 0
yields, tb = [], []
for n in drawn:
    d = (REC[n].get('validation') or {}).get('design_code_reports') or {}
    if not d:
        continue
    n_d += 1
    op5 += bool(((d.get('V_operating_regime') or {}).get('at_5_uL_min') or {}).get('ok'))
    st5 += bool((d.get('V_steric') or {}).get('ok'))
    fab += bool((d.get('V_fab_tolerance') or {}).get('ok'))
    y = (d.get('V_fab_tolerance') or {}).get('yield_F')
    if y is not None:
        yields.append(y)
    t = (d.get('B7_transfer_bias') or {}).get('D_true_over_D_predicted')
    if t:
        tb += list(t)
V['DCR_LIM6'] = ('the operating-regime contract holds at 5 µL/min on %d of %d designs, the steric screen passes on %d, and '
                 'the width tolerance\'s yield runs from %.2f to %.2f (ok on %d); the transfer bias, never applied, runs '
                 'from %.2f to %.2f (the records\' `validation.design_code_reports`).'
                 % (op5, n_d, st5, min(yields), max(yields), fab, min(tb), max(tb)))

# ---- the suites (SUITES_DIR): one row per suite, the composition suite's two halves together
SCOPE = {
    'test_kernels': 'every kernel equals the design code\'s, bit for bit, over grids of inputs',
    'test_teaching': 'the teaching module against the kernels',
    'test_units': 'every unit conversion and dimensionless group of the record, and the curved-duct limits',
    'test_bend_fd': 'the curved-duct closed form against a finite-difference solve',
    'test_design_code': 'the design code\'s reports and validators, each against CIF_v2\'s own function',
    'test_auditor': 'the ported auditor against CIF_v2\'s, its solver, fourteen seeded faults, the ruler\'s corners',
    'test_negative': 'refusals and spec errors: every bad input caught with its exit code',
    'test_lattice': 'the plan over a lattice of targets and processes',
    'test_presets_match_replicas': 'each preset plans its replica\'s design',
    'test_composition': 'the parts composed over 22 cases and extras',
    'test_rc7': 'the rc7 acceptance checks on the photomasks and trace results',
    'test_rc8': 'the RC8 acceptance checks, with the arrays recomposed',
    'test_figure_anchors': 'the reproductions against the papers\' printed figures (needs the papers\' page images)',
    'test_trace_anchors': 'the traces against the micrographs (needs the papers\' page images and the variant builds)',
    'test_options_sheet': 'docs/OPTIONS.md against the specification\'s keys',
    'test_appendix_trace': 'every constant and identity the code uses against the shipped derivation documents',
}
runs = {}
for line in open(os.path.join(SUITES, 'summary.txt'), encoding='utf-8'):
    m = re.match(r'^(\S+) (test_\w+) rc=(\d+) (\d+)s (.*)$', line.strip())
    if m:
        runs[m.group(2)] = dict(rc=int(m.group(3)), s=int(m.group(4)), result=m.group(5).strip())


def n_checks(name):
    p = os.path.join(SUITES, name + '.log')
    if not os.path.exists(p):
        return None
    t = open(p, encoding='utf-8', errors='replace').read().splitlines()
    last = t[-1] if t else ''
    m = re.search(r'(\d+) checks', last)
    if m:
        return int(m.group(1))
    return sum(1 for ln in t if re.match(r'^\s*ok\b', ln) or re.match(r'^\S+\s+OK\s', ln)
               or ' ok ' in ln or ln.rstrip().endswith(' ok'))


rows, n_ok = [], 0
for name in SCOPE:
    parts = [k for k in runs if k == name or k.startswith(name + '_')]
    if not parts:
        continue
    res = [runs[k] for k in sorted(parts)]
    ok = all(r['rc'] == 0 for r in res)
    n_ok += ok
    ch = [n_checks(k) for k in sorted(parts)]
    rows.append('| `%s`%s | %s | %s | %s | %s |' % (
        name, ' (two halves)' if len(parts) > 1 else '', SCOPE[name],
        ' + '.join(str(c) for c in ch if c is not None) or 'see the log',
        'ALL OK' if ok else '; '.join(r['result'] for r in res),
        ' + '.join('%d s' % r['s'] for r in res)))
V['SUITES_TABLE'] = '\n'.join(rows)
V['N_SUITES'] = str(len(rows))
V['VALIDATION_SHORT'] = ('%d of %d test suites ALL OK; %s plans identical in their compared fields; %s; %d records '
                         'identical apart from the named additions and %d changed for a named reason; every photomask '
                         'audited in the build and again from scratch with the same status (`docs/VALIDATION.md`).'
                         % (n_ok, len(rows), V['N_SPECS'], V['MASK_RESULT_SHORT'],
                            len(rc['identical_apart_from_named_additions']), len(changed)))
for k in ('DEVICES_RESULT', 'TUTORIAL_RESULT', 'REVIEW_RESULT'):
    V.setdefault(k, '{{%s}}' % k)
for k in KEEP:
    V[k] = '{{%s}}' % k
json.dump(V, open(OUT, 'w', encoding='utf-8'), indent=1, ensure_ascii=False)
for k, v in sorted(V.items()):
    print('%-18s %s' % (k, v if len(v) < 400 else v[:400] + ' ...'))
