"""test_composition.py: the composition sweep, which tests the clearance guarantee.

The generator composes a device from parts, which are a layout (linear, serpentine, or paperclip, nested or in path
order, with a fixed or fitted leg count), a head (none, pack, hairpins, arms, or narrow), an inlet (plain, filter, or
Dinh's filter), an outlet (separate, teardrop, ring, or flared), a taper, and a start rule. The sweep composes them on
four processes, and no combination may produce a photomask whose parts cross.

Every case is listed in the printed table, and the suite asserts four things.
  1. The plan stage (generate(check_only=True), which draws the outline in memory exactly as a build would and applies
     the build's collision checks) answers every case with an exit code (0, 3, 5) or a spec error (6), and never with
     a crash. A layout that the generator places itself (no spacing given) is answered 0 or 3 with its outline drawn
     and is never refused for colliding parts, and its outlet stretch is lengthened only where the outlet would meet
     the inlet end (two nested legs). A spacing given too tight is refused, and the refusal names the two parts and the
     remedies (the generator's own spacing).
  2. In the full builds of a subset, every hard geometric check of the drawn mask passes. These checks are part
     clearance, outlet separation, port spacing, one polygon with the planned holes, continuous side channels, the DRC
     with no lost wall or channel, and every other hard check apart from Q1 and Q2 (corner radius, topology, features
     found, envelope, and so on). The build's exit code is the check's, and so is its outlet stretch. A case may
     declare a FAIL on Q1 or Q2, which judge the drawn mask's physics that the check does not measure (case 00, a floor
     start in evaluate mode, like m14). The physics audit asks the same questions again from the photomask file (its P4
     is Q1's, and its P5, P6 and P6b are Q2's), so such a case may also fail the audit on those checks, and on no other
     (P0, R0, P1 to P3, P9, P10 and R6b judge the file, the code and the geometry). Its photomask is then
     <name>_FAILED.dxf, and otherwise <name>.dxf, and no <name>_UNASSESSED.dxf remains. A colliding spec is refused by the build too, with no photomask written and with
     the refusal figure drawn.
  3. A failing build's photomask is never under the deliverable name (pipeline.name_mask). The build's
     <name>_UNASSESSED.dxf takes its verdict's name, and a re-assessment that clears a FAIL renames it back.
  4. A parallel array of a clean element is composed and its photomask written, and its collision checks pass.

Run it with `python3 tests/test_composition.py`, which takes 2 to 3 hours on one core for 22 plan-stage outlines (each
drawn in memory), 8 builds, and an array. OPENCIF_TEST_QUICK=1 runs the plan stage only. OPENCIF_TEST_CASES=0-10 (or
11-21,x) runs a part of the cases (x selects the photomask-name and array checks), so that the suite can run in parts
or in parallel."""
import os
import sys
import tempfile
import time
import traceback

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
T_OUT = os.environ.get('OPENCIF_TEST_OUT') or tempfile.mkdtemp(prefix='opencif_test_')
if not os.environ.get('OPENCIF_TEST_OUT'):               # a temporary directory of our own, removed at exit
    import atexit as _atexit, shutil as _shutil
    _atexit.register(_shutil.rmtree, T_OUT, True)
sys.path.insert(0, ROOT)
from opencif import api, spec as SP                      # noqa: E402
from opencif import pipeline as PL                       # noqa: E402
print('opencif %s (code %s)' % (PL.CODE['version'], PL.CODE['sha256']))

fails = []


def check(name, cond, info=''):
    print('%-96s %s %s' % (name, 'ok ' if cond else 'FAIL', info))
    if not cond:
        fails.append(name)


QUICK = bool(os.environ.get('OPENCIF_TEST_QUICK'))
EV = dict(mode='evaluate', seal='allow')                 # heads that cannot seal still draw (gap 1 declared)
# Each case is (label, D in um, F, generate() arguments, expectation, whether to build in full[, the build's expected
# FAILs, which are Q1 or Q2 checks only]). The expectation is one of
#   'placed'  the generator places every part (no spacing given), so the case is never refused for colliding parts
#   'tight'   a leg spacing is given too tight for the head, so the case is refused for colliding parts, with remedies
CASES = [
    # A floor start in evaluate mode leaves rows below the floor, and those rows FAIL Q2.
    ('generic19 linear, no head, floor start, plain, separate', 4.0, 0.3,
     dict(process='generic19', layout='linear', meander=None, start='floor', mode='evaluate'), 'placed', True,
     ('Q2',)),
    ('generic19 linear, pack, filter, teardrop', 4.0, 0.5,
     dict(process='generic19', layout='linear', meander='pack', inlet='filter', outlet='teardrop'), 'placed', False),
    ('generic19 linear, narrow, plain, ring', 4.0, 0.5,
     dict(process='generic19', layout='linear', meander='narrow', outlet='ring'), 'placed', False),
    ('generic19 linear, pack, plain, separate, flared outlet 10 deg', 3.0, 0.5,
     dict(process='generic19', layout='linear', meander='pack', outlet=dict(style='separate', flare_deg=10.0)), 'placed', True),
    ('generic19 serpentine 3, pack, plain, separate', 4.0, 0.6,
     dict(process='generic19', layout='serpentine', legs=3, meander='pack'), 'placed', True),
    ('generic19 serpentine 5, hairpins, filter, teardrop', 4.0, 0.6,
     dict(process='generic19', layout='serpentine', legs=5, meander='hairpins', inlet='filter', outlet='teardrop', **EV), 'placed', False),
    ('generic19 serpentine auto (25 x 75 chip), narrow, filter, teardrop', 6.0, 0.4,
     dict(process='generic19', chip=(25, 75), layout='serpentine', legs='auto', meander='narrow', inlet='filter',
          outlet='teardrop', mode='evaluate'), 'placed', False),
    ('generic19 paperclip 2 nested, pack, plain, separate', 4.0, 0.6,
     dict(process='generic19', layout='paperclip', legs=2, meander='pack'), 'placed', True),
    ('generic19 paperclip 4 nested, narrow, filter, teardrop', 4.0, 0.6,
     dict(process='generic19', layout='paperclip', legs=4, meander='narrow', inlet='filter', outlet='teardrop'), 'placed', False),
    ('generic19 paperclip 4 in path order, pack, plain, separate', 4.0, 0.6,
     dict(process='generic19', layout=dict(style='paperclip', order='path'), legs=4, meander='pack'), 'placed', True),
    ('generic19 paperclip auto (25 x 75 chip), pack, plain, separate', 6.0, 0.4,
     dict(process='generic19', chip=(25, 75), layout='paperclip', legs='auto', meander='pack', mode='evaluate'), 'placed', False),
    ('generic19 paperclip 3 nested, constant-shear taper, pack, filter, ring', 4.0, 0.6,
     dict(process='generic19', layout='paperclip', legs=3, taper=('constant_shear', 100.0, 70.0), meander='pack',
          inlet='filter', outlet='ring'), 'placed', False),
    ('gifford2014 paperclip 4 in path order, pack (the Fig. 4 variant\'s parts), plain, separate', 5.2, 0.5,
     dict(process='gifford2014', layout=dict(style='paperclip', order='path'), legs=4, meander='pack'), 'placed', True),
    ('gifford2014 paperclip 3 nested, no head, floor start, plain, separate', 5.2, 0.9,
     dict(process='gifford2014', layout='paperclip', legs=3, meander=None, start='floor', mode='evaluate'), 'placed', False),
    ('gifford2014 paperclip 6 in path order, narrow, filter, teardrop', 5.2, 0.6,
     dict(process='gifford2014', layout=dict(style='paperclip', order='path'), legs=6, meander='narrow', inlet='filter',
          outlet='teardrop'), 'placed', False),
    ('dinh2024 linear, arms 3, Dinh filter, separate, flared outlet 20 deg', 2.5, 0.6,
     dict(process='dinh2024', layout='linear', meander='arms', inlet='dinh_filter',
          outlet=dict(style='separate', flare_deg=20.0), **EV), 'placed', False),
    ('dinh2024 serpentine 3, arms 3, Dinh filter, teardrop', 2.5, 0.6,
     dict(process='dinh2024', layout='serpentine', legs=3, meander='arms', inlet='dinh_filter', outlet='teardrop', **EV),
     'placed', False),
    ('dinh2024 paperclip 3 nested, pack, filter, separate', 3.0, 0.7,
     dict(process='dinh2024', layout='paperclip', legs=3, meander='pack', inlet='filter'), 'placed', True),
    ('xia2016 paperclip 4 nested, hairpins, plain, ring', 7.0, 0.6,
     dict(process='xia2016', layout='paperclip', legs=4, meander='hairpins', outlet='ring', **EV), 'placed', False),
    ('xia2016 linear, constant-shear taper, hairpins, plain, separate', 7.0, 0.6,
     dict(process='xia2016', layout='linear', taper=('constant_shear', 300.0, 150.0), meander='hairpins', **EV),
     'placed', False),
    # This is the case that A. H. Routt found (2026-10-19 16:08), Gifford Fig. 4's design at a 20 um floor with a
    # packing head and the legs in path order at the inset's 818 um. Before the guarantee, rc7 drew the head across
    # leg 2 and failed it afterwards.
    ('gifford Fig. 4 at a 20 um floor, pack head, path order, legs 818 um apart (too tight for the head)', 5.2032, 0.90644,
     dict(process=dict(depth_um=80.0, gap_um=20.5, post_length_um=25.3, post_width_um=59.5, min_feature_um=20.0,
                       fillet_um=5.0), width_in_um=100.0,
          layout=dict(style='paperclip', order='path', leg_spacing_um=818.0), legs=4, meander='pack'), 'tight', True),
    ('the same with the legs at the generator\'s own spacing', 5.2032, 0.90644,
     dict(process=dict(depth_um=80.0, gap_um=20.5, post_length_um=25.3, post_width_um=59.5, min_feature_um=20.0,
                       fillet_um=5.0), width_in_um=100.0,
          layout=dict(style='paperclip', order='path'), legs=4, meander='pack'), 'placed', False),
]

SEL = os.environ.get('OPENCIF_TEST_CASES')
PICK, EXTRAS = set(range(len(CASES))), True
if SEL:
    PICK, EXTRAS = set(), False
    for tok in SEL.split(','):
        tok = tok.strip()
        if tok == 'x':
            EXTRAS = True
        elif '-' in tok:
            a_, b_ = tok.split('-')
            PICK |= set(range(int(a_), int(b_) + 1))
        elif tok:
            PICK.add(int(tok))
    print('cases %s%s (OPENCIF_TEST_CASES)' % (sorted(PICK), ' and the extras' if EXTRAS else ''))
rows = []
for i, case_ in enumerate(CASES):
    label, D, F, kw, expect, full = case_[:6]
    q_expected = tuple(case_[6]) if len(case_) > 6 else ()
    if i not in PICK:
        continue
    name = '_t_comp%02d' % i
    out = os.path.join(T_OUT, name)
    t0 = time.time()
    crash, rec, code, verdict = '', None, None, ''
    try:
        r = api.generate(D, F, out=out, check_only=True, name=name, **kw)
        rec, code, verdict = r.record, r.exit_code, r.verdict
    except SP.SpecError as e:
        code, verdict = 6, 'SPEC ERROR: %s' % e
    except Exception as e:                               # never, because every answer is a record or a spec error
        crash = '%s: %s' % (type(e).__name__, e)
        traceback.print_exc()
    lr = (rec or {}).get('layout_refusal') or ''
    t_chk = time.time() - t0
    lay = (rec or {}).get('layout')
    ext = float((lay or {}).get('outlet_stretch_extended_um') or 0.0)
    rows.append([i, label, code, 'layout' if lr else ('refused' if code == 5 else ''), '%.0f' % ext if ext else '',
                 '%.0f' % t_chk])
    check('plan stage %02d: answered (exit %s), no crash: %s' % (i, code, label[:60]), not crash and code in (0, 3, 5, 6),
          crash[:100] or verdict[:100])
    if expect == 'placed':
        # The case must be answered 0 or 3 with its outline drawn, because a case refused before the outline, or a
        # spec error, would pass "not refused for colliding parts" without being outlined at all (code review 4, F3).
        check('plan stage %02d: the generator\'s own placement is outlined (exit 0 or 3) and never refused for colliding '
              'parts' % i, code in (0, 3) and lay is not None and not lr, (lr or verdict)[:120])
        # Two nested legs end with the last leg beside the inlet, with no loop there to clear.
        if 'paperclip 2 nested' in label:
            check('plan stage %02d: two nested legs: the outlet would meet the inlet; its stretch is lengthened instead' % i,
                  ext > 0 and not lr, '%.0f um' % ext)
        else:
            check('plan stage %02d: the outlet stretch is not lengthened (no outlet meets the inlet end)' % i, ext == 0.0,
                  '%.0f um' % ext if ext else '')
    else:
        check('plan stage %02d: a spacing too tight for the head is refused, naming the parts and the generator\'s spacing' % i,
              code == 5 and 'collide' in lr and "the generator's own rule" in lr and 'Remedies' in lr, (lr or verdict)[:120])
    if not full or QUICK or crash or code == 6:
        continue
    # ---- The full build of the same spec
    d = SP.load(os.path.join(out, name + '_spec.yaml'))
    t1 = time.time()
    try:
        rb = PL.run(d, out, verbose=False, render=(expect == 'tight'))
        crash_b = ''
    except Exception as e:
        rb, crash_b = {}, '%s: %s' % (type(e).__name__, e)
        traceback.print_exc()
    a = rb.get('assessment') or {}
    rows[-1] += ['%s' % a.get('exit_code'), '%.0f' % (time.time() - t1)]
    if expect == 'tight' or code == 5:
        check('build %02d: refused like the check (exit 5), no photomask written%s' % (
            i, ', the refusal figure drawn' if expect == 'tight' else ''),
            not crash_b and a.get('exit_code') == 5
            and not any(os.path.exists(os.path.join(out, name + s_)) for s_ in ('.dxf', '_FAILED.dxf', '_UNASSESSED.dxf'))
            and (expect != 'tight' or os.path.exists(os.path.join(out, name + '_refused_layout.png'))),
            crash_b[:100] or a.get('verdict', '')[:100])
        continue
    V = rb.get('validation') or {}
    geo = dict(part_clearance=(V.get('part_clearance') or {}).get('ok'),
               outlet_separation=(V.get('outlet_separation') or {}).get('separated'),
               port_spacing=(V.get('port_spacing') or {}).get('ok'),
               one_polygon=(V.get('topology') or {}).get('single_component'),
               holes_as_planned=(V.get('topology') or {}).get('holes_ok'),
               side_continuity=(V.get('side_continuity') or {}).get('ok'),
               drc_no_lost=((V.get('drc') or {}).get('n_lost_open') == 0 and (V.get('drc') or {}).get('n_lost_solid') == 0))
    ex_b = a.get('exit_code')
    failing = [c['check'] for c in a.get('checks', []) if not c.get('ok') and c.get('kind') == 'hard']
    q_fail = [c_ for c_ in failing if c_.startswith(('Q1', 'Q2'))]
    # The physics audit (1.0.0) asks Q1's and Q2's questions again from the file: its P4 is Q1's, and its P5, P6 and P6b
    # are Q2's. A case that declares a FAIL on Q1 or Q2 may therefore fail the audit on those checks, and on no other.
    aud = [c_ for c_ in a.get('checks', []) if str(c_.get('check', '')).startswith('physics audit of the photomask file')]
    aud_failed = sorted(((aud[0].get('value') or {}).get('failed') or [])) if aud and not aud[0].get('ok') else []
    allowed = {x for q_ in q_expected for x in {'Q1': ('P4',), 'Q2': ('P5', 'P6', 'P6b')}.get(q_, ())}
    aud_physics = bool(aud_failed) and set(aud_failed) <= allowed
    # Every other failing hard check is a geometric failure, whatever its name (corner radius, topology, features
    # found, envelope, and so on), and fails the case (code review 4, F2).
    bad = [k for k, v in geo.items() if not v] + [c_ for c_ in failing if c_ not in q_fail and
                                                   not (aud and c_ == aud[0]['check'] and aud_physics)]
    rows[-1][3] = ', '.join(failing) if failing else rows[-1][3]
    q_ok = all(c_.startswith(q_expected) for c_ in q_fail) if q_expected else not q_fail
    want = 2 if q_expected else code
    check('build %02d: exit %s: the check\'s (%s)%s' % (
        i, ex_b, code, ', with the FAIL the case declares (%s)' % ', '.join(q_expected) if q_expected else ''),
          not crash_b and ex_b == want and q_ok and (bool(q_fail) == bool(q_expected)),
          crash_b[:100] or a.get('verdict', '')[:100])
    check('build %02d: every hard geometric check of the drawn mask passes (all but Q1, Q2 and the audit\'s same '
          'questions)' % i, not bad, 'failed: %s' % bad if bad else ('the audit fails %s, as the case declares'
                                                                     % aud_failed if aud_failed else ''))
    ext_b = float((rb.get('design') or {}).get('outlet_stretch_extended_um') or 0.0)
    check('build %02d: the outlet stretch is the check\'s (+%.0f um)' % (i, ext), abs(ext_b - ext) < 1e-6,
          '%.0f um' % ext_b if ext_b != ext else '')
    has_ok, has_fail, has_un = (os.path.exists(os.path.join(out, name + s_))
                                for s_ in ('.dxf', '_FAILED.dxf', '_UNASSESSED.dxf'))
    check('build %02d: the photomask is %s, and no _UNASSESSED.dxf remains' % (
        i, '<name>_FAILED.dxf (a FAIL)' if ex_b == 2 else '<name>.dxf'),
          not has_un and ((has_fail and not has_ok) if ex_b == 2 else (has_ok and not has_fail)))

# ---- 3. A failing build's photomask is never under the deliverable name
if EXTRAS:
    nm_dir = os.path.join(T_OUT, '_t_name_mask')
    os.makedirs(nm_dir, exist_ok=True)
    def only(*names_):
        """Return True when the folder holds exactly the files names_."""
        return sorted(os.listdir(nm_dir)) == sorted(names_)
    seq = []
    for start_, code_, want_ in (('dev_UNASSESSED.dxf', 0, 'dev.dxf'), (None, 2, 'dev_FAILED.dxf'), (None, 0, 'dev.dxf'),
                                 ('dev_UNASSESSED.dxf', 2, 'dev_FAILED.dxf'), (None, 3, 'dev.dxf')):
        if start_:
            for f_ in os.listdir(nm_dir):
                os.remove(os.path.join(nm_dir, f_))
            open(os.path.join(nm_dir, start_), 'w').write('0\nEOF\n')
        got_ = PL.name_mask(nm_dir, 'dev', code_)
        seq.append(got_ == want_ and only(want_))
    check('the build\'s <name>_UNASSESSED.dxf takes its verdict\'s name: a FAIL (exit 2) <name>_FAILED.dxf, else '
          '<name>.dxf; a re-assessment that clears a FAIL renames it back', all(seq), str(seq))

# ---- 4. A parallel array of a clean element
if EXTRAS and not QUICK:
    out_a = os.path.join(T_OUT, '_t_comp_array')
    try:
        ra = api.generate(4.0, 0.5, out=out_a, name='_t_comp_array', process='generic19', layout='linear', meander='pack',
                          array=3, render=False)
        arr_ = (ra.record or {}).get('array') or {}
        ch_ = arr_.get('checks') or {}
        check('array: 3 copies of a clean element composed, the photomask written, its collision checks pass',
              str(arr_.get('verdict', '')).startswith('PASS') and os.path.exists(os.path.join(out_a, '_t_comp_array_x3.dxf'))
              and not os.path.exists(os.path.join(out_a, '_t_comp_array_x3_UNASSESSED.dxf'))
              and all((ch_.get(k) or {}).get('ok', True) for k in ('components', 'product_port_separated_in_every_element',
                                                                   'manifold_touches_elements_only_at_their_ports')),
              str(arr_.get('verdict', ''))[:100])
    except Exception as e:
        check('array: 3 copies of a clean element composed', False, '%s: %s' % (type(e).__name__, e))

print()
print('| # | case | check exit | refusal, or the build\'s failing checks | outlet stretch + um | check s | build exit | build s |')
print('|---|---|---|---|---|---|---|---|')
for r_ in rows:
    r_ = r_ + [''] * (8 - len(r_))
    print('| %s |' % ' | '.join(str(x) for x in r_))
print('ALL OK' if not fails else 'FAILED: %s' % fails)
sys.exit(1 if fails else 0)
