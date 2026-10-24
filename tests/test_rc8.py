"""test_rc8.py: the fixes made in rc8 from the third-party reviews of rc7 (2026-10-20) and their adjudication (review8).

Each fix comes with a control showing that rc7's behavior for valid input is unchanged. The letters below begin the
names of the checks.
  V   Input validation. Every length, depth, flow, radius, tolerance, and count that the design divides by or draws is
      a finite number > 0 (a fillet, the gap-entrance coefficient, corner radii, margins, and a tubing length may be
      0), D is finite, and every bad input is a SpecError naming its field. In rc7, the inputs that escaped this rule
      were a negative or zero flow, a zero gap, and a zero tubing bore (ZeroDivisionError), a negative depth (refused
      only inside the planner), a zero or negative chip, an array count of 0 or -2, D = inf, and an unknown process or
      non-numeric values (raised as ValueError). A true/false value where a number belongs (array=True,
      `layout.legs: true`; review 3, E9) is refused, and so is one passed to the numeric arguments of generate()
      (flow_uL_min, the widths, and the chip and taper tuples; the RC8 code review, B1). The command line's --array 0
      is refused rather than dropped (n1).
  C   The missed pairs of the conflict rule (rc7 LIMITATIONS 29). These are meander.enabled beside a reproduction.head
      that draws the opposite, and an explicit width that sets the head's arms (narrow_um for the narrow style,
      channel_um for pack and arms) beside a different reproduction.w_arm. Agreeing pairs, and w_arm with no explicit
      width, pass, and a style beside reproduction.head none is reported once.
  P   The kappa-plateau pre-flight. A headed strict design far below the plateau is refused before the head's sealing
      search (which took more than 20 minutes on a 25 x 75 mm chip in rc7), with the plateau's advice and a gates note
      in its record and no photomask. Head-less designs, the other modes, and a headed design on the plateau are
      unchanged.
  S   A head-less floor start outside reproduce mode says before drawing that its entrance rows are not declared and
      that, where they over-extract (the usual case), they FAIL Q2.
  Q   A plan whose F overshoots the target beyond checks.F_tolerance says before drawing that Q2 F is expected to FAIL.
  N   pipeline.name_mask names the file just assessed, whereas rc7 kept an earlier <stem>.dxf beside a new pass and
      moved it over a new FAIL. A re-assessment with no fresh file still renames.
  I   compose.inter_copy_port_spacing applies the 2 mm inlet rule between copies, passing and failing, shared and
      separate.
  L   (not QUICK) The leg fitter, unchanged from rc7, fits 5 legs to the design of review8's patch review, whose 3-leg
      outline collides (inlet against leg 3) and whose 5-leg outline is clear at a given spacing.
  D   (not QUICK) A full build records validation.drc.file_final (the verdict's name) beside file (the name audited).
  A   (not QUICK; needs runs/ with both shipped arrays' elements) compose.array clears the array's own earlier outputs
      (a planted stale photomask does not survive), and both shipped arrays compose to their recorded bytes.
The restart of the extension loop (build.plan_and_draw) is exercised by case 07 of tests/test_composition.py. The exit
codes are
  0   every check passed
  1   a check failed
  4   a check could not run (inputs missing)"""
import json
import os
import shutil
import sys
import tempfile
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import opencif                                                     # noqa: E402
from opencif import spec as SP                                     # noqa: E402
from opencif import pipeline as PL                                 # noqa: E402
from opencif import compose as CP                                  # noqa: E402

QUICK = bool(os.environ.get('OPENCIF_TEST_QUICK'))
RUNS = os.environ.get('OPENCIF_RUNS', os.path.join(ROOT, 'runs'))
TMP = tempfile.mkdtemp(prefix='opencif_rc8_')
RESULTS, MISSING = [], []


def check(name, ok, detail=''):
    RESULTS.append((name, bool(ok)))
    print('%-104s %s  %s' % (name[:104], 'ok ' if ok else 'FAIL', str(detail)[:150]), flush=True)


def outcome(D=5.0, F=0.5, **kw):
    """Plan the design and return (kind, text), where kind is 'spec' for a SpecError, 'other' for any other exception,
    and 'record' otherwise, with the exit code in the text."""
    try:
        r = opencif.generate(D, F, check_only='plan', verbose=False, out=os.path.join(TMP, 'probe'), **kw)
        return 'record', '%s %s' % (r.exit_code, r.verdict[:120])
    except SP.SpecError as e:
        # The whole message is returned, because a remedy's words may sit past 240 characters.
        return 'spec', str(e)
    except Exception as e:
        return 'other', '%s: %s' % (type(e).__name__, str(e)[:150])


GEN = dict(depth_um=100.0, gap_um=19.0, post_length_um=27.0, post_width_um=50.0, min_feature_um=19.0)


def notes_of(r):
    """Return the notes of the result r as one string, from its assessment or else from its record."""
    return ' '.join((r.record.get('assessment') or {}).get('notes') or r.record.get('notes') or [])


def main():
    """Run every check and exit with 0, 1, or 4, as the module docstring states."""
    print('code', PL.code_stamp())
    # ---- V
    for label, args, field in (
            ('negative flow', dict(flow_uL_min=-5.0), 'operation.flow_uL_min'),
            ('zero flow', dict(flow_uL_min=0.0), 'operation.flow_uL_min'),
            ('infinite flow', dict(flow_uL_min=float('inf')), 'operation.flow_uL_min'),
            ('zero gap', dict(process=dict(GEN, gap_um=0.0)), 'process.gap_um'),
            ('negative depth', dict(process=dict(GEN, depth_um=-80.0)), 'process.depth_um'),
            ('chip (0, 0)', dict(chip=(0, 0), layout='serpentine'), 'chip.width_mm'),
            ('negative chip', dict(chip=(-25, 75), layout='serpentine'), 'chip.width_mm'),
            ('array count 0', dict(array=0), 'array.count'),
            ('array count -2', dict(array=-2), 'array.count'),
            ('negative leg spacing', dict(layout='paperclip', overrides={'layout': {'leg_spacing_um': -818.0}}), 'layout.leg_spacing_um'),
            ('negative port radius', dict(inlet={'style': 'plain', 'port_radius_um': -1.0}), 'inlet.port_radius_um'),
            ('negative fillet', dict(process=dict(GEN, fillet_um=-1.0)), 'process.fillet_um'),
            ('zero tubing bore', dict(outlet={'style': 'separate', 'tubing': {'centre_id_mm': 0.0, 'side_id_mm': 0.5, 'length_m': 0.3}}), 'outlet.tubing.centre_id_mm'),
            ('zero arms', dict(meander={'style': 'arms', 'arms': 0}), 'meander.arms'),
            ('negative margin', dict(chip=(25, 75, -100), layout='serpentine'), 'chip.margin_um'),
            ('non-numeric flow', dict(flow_uL_min='fast'), 'operation.flow_uL_min'),
            ('non-numeric taper width', dict(taper=('constant', 'wide')), 'centre_channel.width_in_um')):
        kind, txt = outcome(**args)
        check('V %s is a SpecError naming %s' % (label, field), kind == 'spec' and field in txt, txt)
    for label, D, F, word in (('D = inf', float('inf'), 0.5, 'cutoff_um'), ('D = "five"', 'five', 0.5, 'cutoff_um'),
                              ('F = "half"', 5.0, 'half', 'filtrate_fraction')):
        kind, txt = outcome(D, F)
        check('V %s is a SpecError naming targets.%s' % (label, word), kind == 'spec' and word in txt, txt)
    kind, txt = outcome(process='nosuchprocess')
    check('V an unknown process is a SpecError', kind == 'spec' and 'nosuchprocess' in txt, txt)
    kind, txt = outcome(meander=None, seal='require')
    check('V a seal without a head is a SpecError', kind == 'spec' and 'seal' in txt, txt)
    # A flag is not a count, whereas rc7 and the candidate read array=True as one copy and `legs: true` as one leg
    # (review 3, E9).
    kind, txt = outcome(layout='linear', array=True)
    check('V array=True is a SpecError (not one copy)', kind == 'spec' and 'array' in txt, txt)
    kind, txt = outcome(layout='linear', overrides={'layout': {'legs': True}})
    check('V layout.legs: true is a SpecError naming layout.legs', kind == 'spec' and 'layout.legs' in txt, txt)
    # The numeric arguments of generate() refuse a flag before float() makes it 1.0 or 0.0 (the RC8 code review, B1).
    for label, kw, field in (('flow_uL_min=True', dict(flow_uL_min=True), 'operation.flow_uL_min'),
                             ('width_in_um=True', dict(width_in_um=True), 'centre_channel.width_in_um'),
                             ('chip=(25, 75, False)', dict(chip=(25, 75, False)), 'chip.margin_um'),
                             ("taper=('constant_shear', 120, True)", dict(taper=('constant_shear', 120, True)),
                              'centre_channel.width_min_um'),
                             ('flow_uL_min=numpy.True_', dict(flow_uL_min=__import__('numpy').True_),
                              'operation.flow_uL_min')):
        kind, txt = outcome(layout='linear', **kw)
        check('V generate(%s) is a SpecError naming %s' % (label, field), kind == 'spec' and field in txt, txt)
    # The command line's --array 0 reaches the validation, whereas rc7 dropped it as no array (the RC8 code review,
    # n1).
    import contextlib
    import io
    from opencif import __main__ as CLI
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        rc_ = CLI.main(['generate', '5', '0.5', '--layout', 'linear', '--array', '0', '--check', '--out',
                        os.path.join(TMP, 'cli_array0')])
    check('V the command line\'s --array 0 is a spec error (exit 6) naming array.count', rc_ == 6 and 'array.count' in buf.getvalue(),
          'exit %s: %s' % (rc_, buf.getvalue().strip()[-150:]))
    for label, proc in (('fillet 0 (square corners)', dict(GEN, fillet_um=0.0)),
                        ('gap entrance coefficient 0', dict(GEN, gap_entrance_c=0.0))):
        kind, txt = outcome(process=proc, meander=False, start='floor', flow_uL_min=10.0)
        check('V %s is accepted (only a negative value is refused)' % label, kind == 'record' and txt.startswith('0'), txt)
    d = SP.load(dict(name='t', targets=dict(cutoff_um=5.0, filtrate_fraction=0.5), process=dict(GEN),
                     centre_channel=dict(width_in_um=100), meander=dict(enabled=True, style='pack', head_width_um='auto'),
                     layout=dict(style='serpentine', legs='auto')))
    check('V "auto" values and absent optional keys pass', d is not None)
    # ---- C
    base = dict(mode='reproduce')
    for label, kw, want, needle in (
            ('meander off beside reproduction.head pack', dict(meander=False, overrides={'reproduction': {'head': 'pack'}}), True, 'meander.enabled'),
            ('meander on (no style) beside reproduction.head none', dict(meander={'enabled': True}, overrides={'reproduction': {'head': 'none'}}), True, 'meander.enabled'),
            ('meander on beside reproduction.head pack (agree)', dict(meander='pack', overrides={'reproduction': {'head': 'pack'}}), False, ''),
            ('narrow_um 19 beside reproduction.w_arm 28', dict(meander={'style': 'narrow', 'narrow_um': 19.0}, overrides={'reproduction': {'w_arm': 28.0}}), True, 'meander.narrow_um'),
            ('narrow_um 19 beside reproduction.w_arm 19 (agree)', dict(meander={'style': 'narrow', 'narrow_um': 19.0}, overrides={'reproduction': {'w_arm': 19.0}}), False, ''),
            ('channel_um 28 beside reproduction.w_arm 40 (pack)', dict(meander={'style': 'pack', 'channel_um': 28.0}, overrides={'reproduction': {'w_arm': 40.0}}), True, 'meander.channel_um'),
            ('reproduction.w_arm 40 with no explicit head width', dict(meander={'style': 'pack'}, overrides={'reproduction': {'w_arm': 40.0}}), False, '')):
        kind, txt = outcome(**dict(base, **kw))
        if want:
            check('C %s: refused, naming it' % label, kind == 'spec' and needle in txt and 'reproduction' in txt, txt)
        else:
            check('C %s: not refused by the conflict rule' % label, not (kind == 'spec' and 'replaces' in txt), txt)
    kind, txt = outcome(**dict(base, meander='pack', overrides={'reproduction': {'head': 'none'}}))
    check('C a style beside reproduction.head none is refused once (the style rule; no second message for the same head)',
          kind == 'spec' and 'meander.style' in txt and 'meander.enabled' not in txt, txt)
    # ---- P (P1 is a full build call, so that a refusal writes its record and no photomask)
    t0 = time.time()
    r = opencif.generate(0.5, 0.5, process='gifford2014', layout='linear', taper=('constant', 100), meander='pack',
                         verbose=False, out=os.path.join(TMP, 'p1'))
    dt = time.time() - t0
    check('P D 0.5 um with a pack head (u 0.005) is refused before the head\'s sealing search',
          r.exit_code == 5 and 'sealing search' in r.verdict and dt < 600, '%s (%.0f s)' % (r.verdict[:90], dt))
    ref_ = r.record.get('refused') or {}
    check('P its refusal record gives the plateau advice (the widest centre the plateau allows) and says what its gates are',
          'kappa plateau' in (ref_.get('advice') or '') and '14.3 um wide' in (ref_.get('advice') or '')
          and 'sealing search' in (ref_.get('gates_note') or ''), str(ref_.get('advice'))[:150])
    check('P no photomask is written for it', not any(f.endswith('.dxf') for f in os.listdir(os.path.join(TMP, 'p1'))),
          sorted(os.listdir(os.path.join(TMP, 'p1'))))
    r = opencif.generate(2.0, 0.3, process='generic19', layout='linear', taper=('constant', 80), meander=False, start='floor',
                         flow_uL_min=20.0, check_only='plan', verbose=False, out=os.path.join(TMP, 'p2'))
    check('P a head-less design below the plateau is refused by the post-plan gate, once', r.exit_code == 5 and
          r.verdict.count('kappa plateau') == 1 and 'sealing search' not in r.verdict, r.verdict[:100])
    r = opencif.generate(2.0, 0.3, process='generic19', layout='linear', taper=('constant', 80), meander=False, start='floor',
                         flow_uL_min=20.0, mode='evaluate', check_only='plan', verbose=False, out=os.path.join(TMP, 'p3'))
    check('P evaluate mode is unchanged (outside the envelope, exit 3)', r.exit_code == 3, r.verdict[:100])
    r = opencif.generate(4.0, 0.5, process='generic19', layout='linear', taper=('constant', 80), meander='pack',
                         check_only='plan', verbose=False, out=os.path.join(TMP, 'p4'))
    check('P a headed strict design on the plateau is untouched (exit 0)', r.exit_code == 0, r.verdict[:100])
    # ---- S
    r = opencif.generate(3.6, 0.7, process='gifford2014', layout='linear', taper=('constant', 72.0), meander=False,
                         start='floor', flow_uL_min=10.0, check_only='plan', verbose=False, out=os.path.join(TMP, 's1'))
    nts = notes_of(r)
    check('S a strict floor start notes that its entrance rows are not declared and the build FAILS Q2',
          'not declared' in nts and 'FAILS Q2' in nts, nts[-160:])
    r = opencif.generate(3.6, 0.7, process='gifford2014', layout='linear', taper=('constant', 72.0), meander=False,
                         start='floor', flow_uL_min=10.0, mode='reproduce', check_only='plan', verbose=False,
                         out=os.path.join(TMP, 's2'))
    check('S a reproduce-mode floor start does not carry that note (its rows are declared)', 'not declared' not in notes_of(r))
    kind, txt = outcome(3.6, 0.7, process='gifford2014', layout='linear', taper=('constant', 72.0), meander=False,
                        flow_uL_min=10.0)
    check('S the sealed-start refusal names the floor start\'s mode condition', kind == 'spec' and 'reproduce' in txt, txt)
    # ---- Q: a plan whose F overshoots the target by more than checks.F_tolerance says so before drawing. In review8's
    # stress design S05, rc7's check said exit 0 and its build FAILED Q2 F, with the plan's last row carrying F 0.00101
    # past the target.
    DEEP = dict(depth_um=120.0, gap_um=25.0, post_length_um=30.0, post_width_um=60.0, min_feature_um=20.0)
    r = opencif.generate(7.01, 0.334, process=DEEP, chip=(25, 75), layout='linear', taper=('constant_shear', 87.6, 61.3),
                         meander='hairpins', inlet='filter', outlet='ring', flow_uL_min=10.0, check_only='plan',
                         verbose=False, out=os.path.join(TMP, 'q1'))
    check('Q a plan that overshoots F beyond the tolerance notes the expected Q2 F FAIL and the remedy',
          'expected to FAIL Q2 F' in notes_of(r) and 'trim_F' in notes_of(r), notes_of(r)[-150:])
    r = opencif.generate(4.0, 0.5, process='generic19', layout='linear', taper=('constant', 80), meander='pack',
                         check_only='plan', verbose=False, out=os.path.join(TMP, 'q2'))
    check('Q a plan within the tolerance carries no such note', 'expected to FAIL Q2 F' not in notes_of(r))
    # ---- N
    nd = os.path.join(TMP, 'names')
    os.makedirs(nd)

    def put(fn, txt):
        """Write txt to the file fn in the names folder."""
        with open(os.path.join(nd, fn), 'w') as fh:
            fh.write(txt)

    def get(fn):
        """Return the text of the file fn in the names folder, or None if there is no such file."""
        p = os.path.join(nd, fn)
        return open(p).read() if os.path.exists(p) else None
    put('a.dxf', 'OLD'); put('a_UNASSESSED.dxf', 'NEW')
    fn = PL.name_mask(nd, 'a', 0)
    check('N pass over an earlier pass: the file just assessed is named', fn == 'a.dxf' and get('a.dxf') == 'NEW'
          and get('a_UNASSESSED.dxf') is None, (fn, get('a.dxf')))
    put('b.dxf', 'OLD'); put('b_UNASSESSED.dxf', 'NEWFAIL')
    fn = PL.name_mask(nd, 'b', 2)
    check('N FAIL over an earlier pass: the new mask is _FAILED, no deliverable name left', fn == 'b_FAILED.dxf' and
          get('b_FAILED.dxf') == 'NEWFAIL' and get('b.dxf') is None, (fn, get('b_FAILED.dxf'), get('b.dxf')))
    put('c_FAILED.dxf', 'OLDFAIL'); put('c_UNASSESSED.dxf', 'NEWPASS')
    fn = PL.name_mask(nd, 'c', 0)
    check('N pass over an earlier FAIL: no stale _FAILED left beside it', fn == 'c.dxf' and get('c.dxf') == 'NEWPASS'
          and get('c_FAILED.dxf') is None, fn)
    put('e_FAILED.dxf', 'X')
    fn = PL.name_mask(nd, 'e', 0)
    check('N a re-assessment (no fresh file) that clears a FAIL renames it, as rc7', fn == 'e.dxf' and get('e.dxf') == 'X', fn)
    put('f.dxf', 'Y')
    fn = PL.name_mask(nd, 'f', 2)
    check('N a re-assessment (no fresh file) that FAILS renames it _FAILED, as rc7', fn == 'f_FAILED.dxf' and get('f_FAILED.dxf') == 'Y'
          and get('f.dxf') is None, fn)
    # ---- I
    ports = {'A': (20000.0, 3000.0), 'B': (20000.0, -300.0), 'C': (21500.0, 1400.0)}
    radii = {'A': 750.0, 'B': 750.0, 'C': 750.0}
    far = CP.inter_copy_port_spacing([((-4000.0, j * 5600.0), 750.0, j) for j in range(3)], ports, radii, 3, 5600.0, 2000.0)
    check('I separate inlets far from every other copy\'s outlets pass', far['ok'] and far['inlet_to_other_copies_outlet_ports_min_um'] > 20000, far)
    near = CP.inter_copy_port_spacing([((19000.0, j * 3000.0), 1020.0, j) for j in range(2)], ports, radii, 2, 3000.0, 2000.0)
    check('I an inlet within 2 mm of another copy\'s outlet port fails, with the pair named', not near['ok']
          and near['closest']['inlet_copy'] != near['closest']['outlet_copy'], near['closest'])
    sh = CP.inter_copy_port_spacing([((-6000.0, 2800.0), 750.0, None)], ports, radii, 2, 5600.0, 2000.0)
    check('I a shared port is checked against every copy\'s outlets', sh['ok'] and sh['closest']['inlet_copy'] == 'shared', sh['closest'])
    one = CP.inter_copy_port_spacing([((-4000.0, 0.0), 750.0, 0)], ports, radii, 1, 5600.0, 2000.0)
    check('I one copy: nothing to check, passes', one['ok'] and one['inlet_to_other_copies_outlet_ports_min_um'] is None, one)
    if not QUICK:
        # ---- L
        t0 = time.time()
        r = opencif.generate(4.0, 0.5, process='generic19', layout='paperclip', taper=('constant', 80), meander=None,
                             start='floor', flow_uL_min=10.0, chip=(10, 19), overrides={'layout': {'leg_spacing_um': 905.0}},
                             check_only=True, verbose=False, out=os.path.join(TMP, 'l1'))
        lay = r.record.get('layout') or {}
        check('L a given spacing whose 3-leg outline collides and whose 5-leg outline is clear fits 5 legs (as rc7)',
              r.exit_code == 0 and lay.get('legs') == 5, '%s legs %s (%.0f s)' % (r.verdict[:70], lay.get('legs'), time.time() - t0))
        # ---- D
        r = opencif.generate(4.0, 0.3, process='generic19', layout='linear', taper=('constant', 80), meander='pack',
                             flow_uL_min=20.0, render=False, verbose=False, out=os.path.join(TMP, 'd1'), name='rc8_drc_name')
        drc = (r.record.get('validation') or {}).get('drc') or {}
        check('D a build records the DRC\'s file and its final name', drc.get('file') == 'rc8_drc_name_UNASSESSED.dxf' and
              drc.get('file_final') == r.record['files']['dxf'] == os.path.basename(r.dxf), (drc.get('file'), drc.get('file_final')))
        # ---- A
        arrays = (('m07_serpentine_const_pack_filter_separate', 3, 'shared'),
                  ('m13_paperclip_const_nomeander_sealed_fine', 4, 'separate'))
        for el, cnt, inlet in arrays:
            an = '%s_x%d' % (el, cnt)
            src_el, src_ar = os.path.join(RUNS, el), os.path.join(RUNS, an, an + '.json')
            if not (os.path.exists(os.path.join(src_el, el + '.pkl')) and os.path.exists(src_ar)):
                MISSING.append(an)
                print('A %s: not run (needs runs/%s with its pickle and runs/%s)' % (an, el, an), flush=True)
                continue
            ed = os.path.join(TMP, 'arr_' + el)
            os.makedirs(ed)
            for ext in ('.json', '.pkl', '.dxf'):
                shutil.copy(os.path.join(src_el, el + ext), ed)
            rec = json.load(open(src_ar))
            stale = os.path.join(ed, an + '.dxf')
            with open(stale, 'w') as fh:
                fh.write('STALE')                              # an earlier array's photomask under the deliverable name
            out = CP.array(ed, el, cnt, ed, rec.get('pitch_um'), inlet, chip=None, do_drc=True)
            from opencif.geometry import dxf as io_dxf
            sha = io_dxf.sha256(os.path.join(ed, out['dxf']))
            icp = out['checks'].get('inter_copy_port_spacing', {})
            check('A %s into a used folder: the stale file is gone, the bytes are the recorded array\'s' % an,
                  out['dxf'] == an + '.dxf' and sha == rec['dxf_sha256'], sha[:12])
            check('A %s: the 2 mm rule between copies holds (%s um; neighbouring outlet ports %s um, reported)'
                  % (an, icp.get('inlet_to_other_copies_outlet_ports_min_um'), icp.get('outlet_ports_of_neighbouring_copies_min_um')),
                  icp.get('ok') is True, icp.get('closest'))
    n_fail = sum(1 for r in RESULTS if not r[1])
    print('\n%d checks, %d failed%s' % (len(RESULTS), n_fail, (', %d not run (inputs missing)' % len(MISSING)) if MISSING else ''))
    shutil.rmtree(TMP, ignore_errors=True)
    if n_fail:
        print('FAILED: %s' % [r[0] for r in RESULTS if not r[1]])
        sys.exit(1)
    if MISSING:
        print('INCOMPLETE')
        sys.exit(4)
    print('ALL OK')
    sys.exit(0)


if __name__ == '__main__':
    main()
