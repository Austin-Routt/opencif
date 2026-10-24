"""test_negative.py: seeded defects, which show that the validators listed here catch the defects they exist for.

A check that cannot fail proves nothing. The suite therefore seeds defects and synthetic clearance violations, each of
which must be caught, with clean controls and follow-up checks (the count is printed at the end). The controls show
what must not be enforced, which is gamma_GS without a head, a design over its chip area (reported and never failed,
App. B C8 CQ5), and, outside strict mode, a channel deeper than the process aspect ratio allows (reported,
docs/LIMITATIONS.md 7). The DRC, topology, continuity, corner-radius, and head-envelope checks are not seeded here (see
docs/EVALUATION.md). Run it with `python3 tests/test_negative.py` (15 to 25 minutes, in which three devices are drawn
and a few are planned and refused)."""
import copy, os, sys
import numpy as np
from shapely import affinity
from shapely.geometry import box, Polygon
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import tempfile as _tempfile
# Scratch output goes outside the package tree, so that a read-only install works, and OPENCIF_TEST_OUT overrides
# the location (review R3 of rc4).
T_OUT = os.environ.get('OPENCIF_TEST_OUT') or _tempfile.mkdtemp(prefix='opencif_test_')
if not os.environ.get('OPENCIF_TEST_OUT'):               # a temporary directory of our own, removed at exit
    import atexit as _atexit, shutil as _shutil
    _atexit.register(_shutil.rmtree, T_OUT, True)

sys.path.insert(0, ROOT)
from opencif import spec as SP, pipeline as PL, build as BU
from opencif import pipeline as _PL_CODE                # noqa: E402  (every suite's log names the code it ran on)
print('opencif %s (code %s)' % (_PL_CODE.CODE['version'], _PL_CODE.CODE['sha256']))
from opencif.physics import schedule as S
from opencif.geometry import layout as L

fails = []
n_checks = [0]
def check(name, cond, info=''):
    n_checks[0] += 1
    print('%-74s %s %s' % (name, 'ok ' if cond else 'FAIL', info))
    if not cond:
        fails.append(name)

# 1. Part clearance on synthetic parts (the flow order a - b - c, in which a joins b, b joins c, and a and c are
#    unrelated)
a, b = box(0, 0, 100, 10), box(90, 0, 200, 10)
# The tolerance is 0.1 um, as in the DRC.
for gap, want in ((-5.0, False), (5.0, False), (18.85, False), (18.95, True), (40.0, True)):
    c = box(0, 10 + gap, 100, 20 + gap)
    r = L.part_clearance([('a', a), ('b', b), ('c', c)], 19.0)
    check('part clearance: unrelated parts %+.2f um apart (floor 19, tol 0.1) -> %s' % (gap, 'ok' if want else 'caught'), r['ok'] == want,
          str(r['closest'][:1]))
r = L.part_clearance([('a', a), ('b', box(50, 0, 150, 10)), ('c', box(0, 60, 10, 70))], 19.0)
check('part clearance: parts joined in flow order may overlap', r['ok'])

# 2. A real drawing (m15: linear, no head, filter inlet, teardrop outlet)
d = SP.load(os.path.join(ROOT, 'examples', 'matrix', 'm15_linear_shear_nomeander_filter_teardrop.yaml'))
spec, opts, sc, g0, env_bad, notes = PL.prepare(d)
sc_l, legs = L.plan_layout(sc, 'linear', 1, None)
fluid, sc_d, info = L.assemble(sc_l, L.make_legs(sc_l))
pc = L.part_clearance(info['parts'], spec.F_min)
check('drawn m15: parts keep the floor apart', pc['ok'], str(pc['closest'][:1]))
moved = [(n, affinity.translate(p, info['inlet']['x_port_c'] - p.centroid.x, -p.centroid.y) if n == 'outlet' else p)
         for n, p in info['parts']]
pc2 = L.part_clearance(moved, spec.F_min)
check('seeded: the outlet moved onto the inlet port is caught by name', (not pc2['ok']) and
      set(pc2['closest'][0]['parts']) == {'inlet and head', 'outlet'}, str(pc2['closest'][:1]))
ps = L.port_spacing(info, spec)
check('drawn m15: outlet ports clear of the inlet port', ps['ok'], '%.0f um' % ps['inlet_to_nearest_outlet_port_um'])
info2 = dict(info); info2['outlet'] = dict(info['outlet'])
ip = (float(info['inlet']['x_port_c']), 0.0)
k0 = next(iter(info['outlet']['ports']))
info2['outlet']['ports'] = dict(info['outlet']['ports'])
info2['outlet']['ports'][k0] = (ip[0], ip[1] + float(info['inlet'].get('port_r', spec.port_r_in)) + spec.port_r + 300.0)
ps2 = L.port_spacing(info2, spec)
check('seeded: an outlet port 0.3 mm from the inlet port is caught', not ps2['ok'], '%.0f um' % ps2['inlet_to_nearest_outlet_port_um'])
cd = BU.components_drawn(fluid, sc_d, info)
bad = PL.requested_vs_drawn(spec, opts, cd)
check('drawn m15: components as requested (filter posts counted on the drawing)', not bad and cd['filter_posts_counted'] > 0,
      'posts %d' % cd['filter_posts_counted'])
inl = info['inlet']
region = box(inl['x_f0'] - 1, -inl['Hf'], inl['x_f1'] + 1, inl['Hf'])
empty = Polygon(fluid.exterior, [r for r in fluid.interiors if not region.contains(Polygon(r).centroid)])
cd2 = BU.components_drawn(empty, sc_d, info)
bad2 = PL.requested_vs_drawn(spec, opts, cd2)
check('seeded: an empty filter bulb is caught', any('filter' in x for x in bad2), str(bad2))

# 3. The request-vs-drawing comparison on seeded records
opts_plain = dict(opts, inlet='plain')
check('seeded: a plain inlet asked, the drawing has filter posts -> caught (counted on the drawing)',
      any('plain inlet' in x for x in PL.requested_vs_drawn(spec, opts_plain, cd)))
o_ = info['outlet']; k_ = next(iter(o_['ports'])); (px, py) = o_['ports'][k_]
from shapely.geometry import Point
no_port = fluid.difference(Point(float(px), float(py)).buffer(float((o_.get('port_r') or {}).get(k_, spec.port_r)) * 0.95))
cd3 = BU.components_drawn(no_port, sc_d, info)
check('seeded: an outlet port missing from the drawing -> caught', any('outlet ports' in x for x in PL.requested_vs_drawn(spec, opts, cd3)),
      'ports on the drawing %s' % cd3['outlet_ports'])
ok_cd = dict(head_units=72, head_arms=0, arm_lane_pairs=[], head_unit_widths_um=[28.0], gap1_f_over_fstar=1.0, gap1_sealed=True,
             filter_posts_counted=0, outlet_style='separate', outlet_ports=3, legs=3, u_turns=2)
sp28 = copy.copy(spec); sp28.GS = 28.0; sp28.w_narrow = None; sp28.n_arms_fixed = 3
o = dict(style='hairpins', seal='require', inlet='plain', outlet='separate', layout='paperclip', legs='auto')
check('clean hairpins record passes', PL.requested_vs_drawn(sp28, o, ok_cd) == [])
# These are the plan-against-request rules. The plan's own values are injected here, and the rules catch a generator
# that ignores the request.
seeds = [('hairpins asked, narrow units drawn (the old silent fallback)', dict(head_unit_widths_um=[19.0, 28.0], head_arms=2)),
         ('seal: require, gap 1 unsealed', dict(gap1_sealed=False, gap1_f_over_fstar=2.5)),
         ('ring asked, separate drawn', None),
         ('3 legs asked, 2 drawn', None),
         ('legs and U-turns inconsistent', dict(u_turns=1))]
for what, mod in seeds:
    cdx = dict(ok_cd); ox = dict(o)
    if mod:
        cdx.update(mod)
    if what.startswith('ring'):
        ox['outlet'] = 'ring'
    if what.startswith('3 legs'):
        ox['legs'] = 3; cdx['legs'] = 2; cdx['u_turns'] = 1
    r = PL.requested_vs_drawn(sp28, ox, cdx)
    check('seeded: %s -> caught' % what, bool(r), str(r[:1]))
check('seeded: the packing head asked, narrow units drawn -> caught (whatever spec.w_narrow says)',
      bool(PL.requested_vs_drawn(sp28, dict(o, style='pack'), dict(ok_cd, head_unit_widths_um=[19.0, 28.0], head_arms=2))))
oa = dict(o, style='arms')
check("seeded: Dinh's arms asked, folded (packed) arms drawn -> caught",
      bool(PL.requested_vs_drawn(sp28, oa, dict(ok_cd, head_arms=3, arm_lane_pairs=[4, 2, 1]))))

# 4. No head is ever substituted, so hairpins that cannot seal gap 1 in strict mode are refused with the numbers
dm = SP.load(os.path.join(ROOT, 'examples', 'matrix', 'm10_paperclip_shear_hairpins_filter_separate.yaml'))
dm['meander']['seal'] = 'require'; dm['gates'] = dict(mode='strict')
sp_, op_, sc_, g_, e_, n_ = PL.prepare(dm)
check('seeded: hairpins + seal: require that cannot seal -> refused, nothing substituted',
      bool(op_.get('seal_refusal')) and sp_.head == 'dinh3' and not sp_.w_narrow, (op_.get('seal_refusal') or '')[:90])
dm2 = copy.deepcopy(dm); dm2['meander']['seal'] = 'allow'
try:
    SP.load(dm2); SP.to_engine(dm2); ok = False
except SP.SpecError:
    ok = True
check('seeded: seal: allow in strict mode -> refused', ok)
dn = SP.load(os.path.join(ROOT, 'examples', 'matrix', 'm09_paperclip_const_narrow_plain_teardrop.yaml'))
dn['meander']['head_width_um'] = 600.0
recn = PL.run(dn, outdir=T_OUT, verbose=False, render=False)
check('seeded: a head that cannot even be planned -> a refusal record (exit 5), not a crash',
      recn['assessment']['exit_code'] == 5 and os.path.exists(os.path.join(T_OUT, dn['name'] + '.json')),
      recn['assessment']['verdict'][:80])
import subprocess, tempfile, yaml
with tempfile.NamedTemporaryFile('w', suffix='.yaml', delete=False) as fh:
    yaml.safe_dump(dn, fh)
pr = subprocess.run([sys.executable, '-m', 'opencif', 'check', fh.name], cwd=ROOT, capture_output=True, text=True)
check('seeded: `opencif check` on that head reports the refusal (exit 5, no traceback)',
      pr.returncode == 5 and 'Traceback' not in pr.stderr and 'seal_refusal' in pr.stdout, pr.stderr[-120:])
os.remove(fh.name)

# 5. Strict mode refuses a design outside the validity envelope before drawing (App. B thesis; CIF_v2 M13 and M14;
#    App. B C8 CQ7)
gl = SP.load(os.path.join(ROOT, 'examples', 'gifford_physics_linear.yaml'))
def strict_refusal(mod, what, key):
    """Build the example changed by mod, check that strict mode refuses it (exit 5) with key in its verdict and no mask
    drawn, and return its record."""
    dd = copy.deepcopy(gl)
    mod(dd)
    rec = PL.run(dd, outdir=T_OUT, verbose=False, render=False)
    a_ = rec['assessment']
    check('seeded: %s -> REFUSED in strict mode (exit 5), no mask drawn' % what,
          a_['exit_code'] == 5 and key in a_['verdict'] and not os.path.exists(os.path.join(T_OUT, dd['name'] + '.dxf')),
          a_['verdict'][:100])
    return rec
# The example's explicit 2 um fillet is not below half a 3 um floor, which is refused as a spec error before the
# plateau gate is reached (that refusal is section 13's). The fillet is therefore lowered with the floor, so that this
# case still reaches the gate it seeds.
strict_refusal(lambda dd: (dd['targets'].update(cutoff_um=3.2), dd['process'].update(min_feature_um=3.0, min_radius_um=0.75,
                                                                                     fillet_um=1.0)),
               'u_min 0.032 below the kappa plateau 0.035', 'kappa plateau')
strict_refusal(lambda dd: dd.update(operation=dict(flow_uL_min=200.0)), 'channel Re 37 > 10 (200 uL/min)', 'Reynolds')
r_ = strict_refusal(lambda dd: dd.update(operation=dict(flow_uL_min=30.0)), 'lamina Re 1.05e-2 > 1e-2 (30 uL/min)', 'Stokes')
g_ = r_['refused']['gates']
dl_ = gl['targets']['cutoff_um'] / S.KAPPA_G
check('  ... the refusal reports both Stokes numbers and the plateau band (kappa_g interval x delta*)',
      g_['Re_delta_max'] > 1e-2 and g_['Re_gap_max'] < 1e-2 and g_['cutoff_band_um']['passage_edge_below_knee'] is None
      and abs(g_['cutoff_band_um']['band_um'][0] - S.KAPPA_CI[0] * dl_) < 1e-9
      and 'cutoff band' in r_['refused']['advice'], 'Re_delta %.3g, Re_gap %.3g' % (g_['Re_delta_max'], g_['Re_gap_max']))
dd = copy.deepcopy(gl); dd['gates'] = dict(mode='evaluate'); dd['operation'] = dict(flow_uL_min=30.0)
sp_, op_, sc_, g5, e5, n5 = PL.prepare(dd)
check('seeded: the same in evaluate mode is not refused, and the gate failure is listed', any('Stokes' in x for x in e5), str(e5)[:90])
for D_, what, key in ((250.0, 'delta*/w_c >= 1/2 (D 250 um on a 100 um channel) -> SpecError (A implementers step 3)', 'half'),
                      (100.0, 'a closure with no solution (D 100 um on 100 um, f_gap 0.11) -> SpecError, not a crash',
                       'cannot be planned'),
                      (0.25, 'a march of more than 400,000 rows (D 0.25 um on 100 um, f_gap 7.8e-7) -> SpecError, not a '
                       'crash', 'march did not reach')):
    dd = copy.deepcopy(gl); dd['targets']['cutoff_um'] = D_
    try:
        PL.prepare(SP.load(dd)); ok = False; msg = ''
    except SP.SpecError as e:
        ok = True; msg = str(e)
    check('seeded: %s' % what, ok and key in msg, msg[:80])
dd = copy.deepcopy(gl); dd['targets']['cutoff_um'] = 100.0
with tempfile.NamedTemporaryFile('w', suffix='.yaml', delete=False) as fh:
    yaml.safe_dump(dd, fh)
pr = subprocess.run([sys.executable, '-m', 'opencif', 'build', fh.name, '--out', T_OUT], cwd=ROOT,
                    capture_output=True, text=True)
check('seeded: the CLI reports a spec error (exit 6, no traceback)', pr.returncode == 6 and 'Traceback' not in pr.stderr
      and 'SPEC ERROR' in pr.stdout, (pr.stdout + pr.stderr)[-100:])
os.remove(fh.name)

# 6. The pre-flight gamma_GS of spec v5 section 18 (and section 14.4; CIF_v2 design()) on a meander head whose fold
#    width over gap lies outside [1.4, 1.5]
dg = SP.load(os.path.join(ROOT, 'examples', 'matrix', 'm01_linear_const_pack_plain_separate.yaml'))
dg['meander']['channel_um'] = 30.0                       # GS/G = 30/19 = 1.58
sp_, op_, sc_, g6, e6, n6 = PL.prepare(copy.deepcopy(dg))
refused_, reasons_, tag_ = PL.refusal(op_, e6)
check('seeded: meander GS/G 1.58 outside gamma_GS -> refused in strict mode, with the reason',
      refused_ and any('gamma_GS' in r for r in reasons_) and op_['mode'] == 'strict', (reasons_ or [''])[0][:80])
dg2 = copy.deepcopy(dg); dg2['gates'] = dict(mode='evaluate')
sp_, op_, sc_, g7, e7, n7 = PL.prepare(dg2)
check('  ... in evaluate mode it is reported in the notes, not refused',
      not PL.refusal(op_, e7)[0] and any('gamma_GS' in n for n in n7))
dg3 = copy.deepcopy(gl); dg3.setdefault('meander', {})['channel_um'] = 60.0     # no head drawn, so GS is not drawn
sp_, op_, sc_, g8, e8, n8 = PL.prepare(dg3)
check('control: without a meander head GS is not drawn and gamma_GS is not applied', not op_.get('preflight_refusal'))

# 7. The chip area is a packaging budget, reported and never a pass condition (App. B thesis; C8 CQ5 "without failing")
dc = SP.load(os.path.join(ROOT, 'examples', 'decoupling', 'dc_D5_F0p3.yaml'))
dc['name'] = '_t_chip_budget'; dc['layout'] = dict(style='linear'); dc['chip'] = dict(width_mm=5.0, length_mm=10.0)
recc = PL.run(dc, outdir=T_OUT, verbose=False, render=False)
ac = recc['assessment']
cc = [c for c in ac['checks'] if c['check'].startswith('fits the chip area')]
check('control: a design over its chip area is drawn and reported (exit 0, verdict says so, not FAIL)',
      ac['exit_code'] == 0 and ac.get('over_chip_area') and 'OVER THE CHIP AREA' in ac['verdict'] and cc and
      not cc[0]['ok'] and cc[0]['kind'] == 'budget', ac['verdict'][-90:])
from opencif import api
rc_ = api.Result(record=recc, spec=dc, out='', verdict=ac['verdict'], exit_code=ac['exit_code'])
check('  ... and the API says so apart from ok: Result.ok and Result.over_chip_area both true; status chip_area OVER',
      rc_.ok and rc_.over_chip_area and not rc_.over_aspect_ratio and ac['status']['chip_area'] == 'OVER')

# 8. CIF_v2's E11 and M16 (march()) and its pre-flight w_min = max(F_min, h/AR_max) (design()). Outside strict mode, the
#    aspect ratio of the drawn CIF channels is reported rather than a pass condition (docs/LIMITATIONS.md 7).
de = copy.deepcopy(gl)                                  # 11 um deep, so delta*/h = 0.111, with every other gate inside
de['process']['depth_um'] = 11.0; de['centre_channel']['width_in_um'] = 150.0
de['operation'] = dict(flow_uL_min=1.0); de['start'] = dict(style='floor')
sp_, op_, sc_, g9, e9, n9 = PL.prepare(SP.load(de))
rf_, rs_, _ = PL.refusal(op_, e9)
check('seeded: delta*/h 0.111 > eps_lam 0.10 (E11) -> refused in strict mode, E11 the only reason',
      rf_ and len(rs_) == 1 and 'thin lamina E11' in rs_[0], (rs_ or [''])[0][:80])
dt = copy.deepcopy(gl)                                  # a steep constant-shear taper (about 0.59 um per row)
dt['targets'].update(cutoff_um=16.0, filtrate_fraction=0.5)
dt['centre_channel'].update(taper='constant_shear', width_min_um=60.0)
sp_, op_, sc_, g10, e10, n10 = PL.prepare(SP.load(dt))
rf_, rs_, _ = PL.refusal(op_, e10)
check('seeded: |dw_c|/(G + g_p) %.2e > eps_taper 1e-2 (M16) -> refused in strict mode (with M13 and Stokes)'
      % g10['taper_rate_max'], rf_ and any('taper rate M16' in r for r in rs_) and g10['taper_rate_max'] > 1e-2,
      '; '.join(r[:28] for r in rs_))
dt['gates'] = dict(mode='evaluate')
with tempfile.NamedTemporaryFile('w', suffix='.yaml', delete=False) as fh:
    yaml.safe_dump(dt, fh)
pr = subprocess.run([sys.executable, '-m', 'opencif', 'check', fh.name], cwd=ROOT, capture_output=True, text=True)
check('seeded: the same in evaluate mode -> `opencif check` exits 3 (outside the envelope), M16 listed',
      pr.returncode == 3 and 'taper rate M16' in pr.stdout and 'Traceback' not in pr.stderr, 'exit %d' % pr.returncode)
os.remove(fh.name)
dw = copy.deepcopy(gl); dw['process']['depth_um'] = 300.0          # h/AR_max = 30 um > G = 20.5 um
dw['process']['aspect_ratio_max'] = 10.0                            # CIF_v2's default (the example declares 25)
sp_, op_, sc_, g11, e11, n11 = PL.prepare(SP.load(dw))
check('seeded: G 20.5 um below w_min = max(F_min, h/AR_max) = 30 um -> refused in strict mode (CIF_v2 pre-flight)',
      PL.refusal(op_, e11)[0] and 'gap G = 20.50 um below w_min' in (op_.get('preflight_refusal') or ''),
      (op_.get('preflight_refusal') or '')[:60])
dw['gates'] = dict(mode='evaluate')
sp_, op_, sc_, g12, e12, n12 = PL.prepare(SP.load(dw))
check('  ... in evaluate mode it is reported in the notes, not refused',
      not op_.get('preflight_refusal') and any('gap G = 20.50 um below w_min' in x for x in n12))
# Spec v5 G7 also applies to the straight side channels. G11 would fold them, whereas opencif does not, so a straight
# side channel below w_min is refused in strict mode, and drawn and reported in the other modes (docs/LIMITATIONS.md 7).
g10_ = copy.deepcopy(gl); g10_['process']['aspect_ratio_max'] = 10.0  # the 5.35 um side channel 80 um deep is 15 : 1
sp_, op_, sc_, g13, e13, n13 = PL.prepare(SP.load(copy.deepcopy(g10_)))
check('seeded: a 5.35 um straight side channel 80 um deep (15 : 1 > 10) -> refused in strict mode before drawing (G7)',
      PL.refusal(op_, e13)[0] and 'straight side channel' in (op_.get('preflight_refusal') or ''),
      (op_.get('preflight_refusal') or '')[:70])
g10_['gates'] = dict(mode='evaluate')
ra = PL.run(copy.deepcopy(g10_), outdir=T_OUT, verbose=False, render=False)
aa = ra['assessment']
ca = [c for c in aa['checks'] if c['check'].startswith('aspect ratio of the narrowest CIF channel')]
check('control: the same in evaluate mode is drawn and reported (exit 0, verdict and API say so, kind fabrication)',
      aa['exit_code'] == 0 and 'ABOVE THE PROCESS ASPECT-RATIO LIMIT' in aa['verdict'] and ca and not ca[0]['ok']
      and ca[0]['kind'] == 'fabrication' and api.Result(record=ra, spec=g10_, out='', exit_code=0).over_aspect_ratio
      and aa['status']['aspect_ratio_of_CIF_channels'].startswith('ABOVE THE LIMIT'), aa['verdict'][-80:])

# 9. An outlet that cannot close is refused by name before drawing (by the check and the build alike), the fitter does
#    not stop at the first leg count at which an outlet fails, and no library text reaches the user
out_t = T_OUT
r9 = api.generate(4.0, 0.9, chip=(25, 75), layout='paperclip', outlet='ring', check_only=True, out=out_t)
check('seeded: a ring outlet that cannot close at F 0.9 -> check_only refuses it by name (exit 5), no library text',
      r9.exit_code == 5 and 'ring outlet' in r9.verdict and 'f(a)' not in r9.verdict
      and 'ring outlet' in (r9.record.get('outlet_refusal') or ''), r9.verdict[:90])
r9b = api.generate(4.0, 0.9, chip=(25, 75), layout='paperclip', outlet='ring', name='_t_ring', out=out_t, render=False)
check('  ... and its build writes the refusal record (exit 5, REFUSED (outlet), advice names teardrop), no mask',
      r9b.exit_code == 5 and r9b.verdict.startswith('REFUSED (outlet)') and 'teardrop' in r9b.record['refused']['advice']
      and not os.path.exists(os.path.join(out_t, '_t_ring.dxf')), r9b.verdict[:90])
sp9, op9, sc9, g9_, e9_, n9_ = PL.prepare(SP.load(api.build_spec(4.0, 0.9, chip=(25, 75), layout='paperclip', outlet='ring')))
tried9 = []
L_ok = L.plan_layout
def _flaky(sc_, kind_, n_, *a_, **k_):
    """Stand in for the layout planner so that the outlet fails at the first count only."""
    if not tried9:
        tried9.append(n_)
        raise L.OutletInfeasible('ring outlet: seeded failure at %d legs' % n_)
    tried9.append(n_)
    return L_ok(sc_, kind_, n_, *a_, **k_)
L.plan_layout = _flaky
try:
    sp9b, op9b, sc9b, g9b, e9b, n9b = PL.prepare(SP.load(api.build_spec(4.0, 0.5, chip=(25, 75), layout='paperclip',
                                                                        outlet='teardrop')))
finally:
    L.plan_layout = L_ok
check('seeded: the outlet fails at the first leg count only -> the next count is tried, nothing refused',
      len(tried9) >= 2 and not op9b.get('outlet_refusal'), 'counts tried %s' % tried9)
calls9 = []
def _stub(sc_, kind_, n_, *a_, **k_):
    """Stand in for the layout planner so that the fitter meets an outlet failure, then a count too long to fold."""
    calls9.append(n_)
    if len(calls9) == 1:
        raise L.OutletInfeasible('ring outlet: seeded failure at %d legs' % n_)
    raise ValueError('device too short for %d paperclip legs (seeded)' % n_)
L.plan_layout = _stub
try:
    n9, best9, tried9b = PL._fit_once(sp9b, sc9b, 'paperclip', op9b, 73000.0, 23000.0, lambda m: None)
finally:
    L.plan_layout = L_ok
check('seeded: the leg fitter meets an outlet failure -> it goes on to the next count (and stops at "too short")',
      calls9[:2] == [calls9[0], calls9[0] + 1] and len(tried9b) == 2 and 'outlet' in tried9b[0]
      and 'too short' in tried9b[1] and n9 is None, '; '.join(t[:40] for t in tried9b))
try:
    api.generate(2.0, 0.5, width_in_um=100, meander=None, check_only=True, out=out_t,
                 process=dict(depth_um=80, gap_um=20.5, post_length_um=25.3, post_width_um=59.5, min_feature_um=4.5,
                              min_radius_um=1.8))
    ok9, msg9 = False, ''
except SP.SpecError as e:
    ok9, msg9 = True, str(e)
check('seeded: a sealed start below the floor -> the spec error also names every other finding (w_min, the plateau)',
      ok9 and 'Also found before drawing' in msg9 and 'kappa plateau' in msg9 and 'straight side channel' in msg9,
      msg9[-90:])

# 10. E11 on every planned row (spec v5 "for all i"), and M16 alone, on the gate function itself. On a constructible
#     design M16 cannot fail alone here, because a taper steep enough for it makes M13 fail first (test 8).
sp_, op_, sc_, g14, e14, n14 = PL.prepare(SP.load(copy.deepcopy(gl)))
sc_e = copy.copy(sc_); sc_e.delta = np.array(sc_.delta, float).copy()
# One planned row's lamina is set to 0.11 h, and the delta* of the cutoff rows is left untouched.
sc_e.delta[sc_.N // 2] = 0.11 * sp_.h
ge = BU.gates(sp_, sc_e, S.KAPPA_G * sp_.delta_star)
check('seeded: one planned row with delta/h 0.11 -> E11 fails on every-row reading, the cutoff-row reading passes',
      ge['lamina_ok'] and not ge['lamina_all_rows_ok'] and abs(ge['lamina_over_depth_planned_max'] - 0.11) < 1e-9,
      'planned max %.3f' % ge['lamina_over_depth_planned_max'])
sc_t = copy.copy(sc_); sc_t.wc = np.array(sc_.wc, float).copy()
k_ = sc_.N // 2
sc_t.wc[k_:] = sc_t.wc[k_:] - 0.5                        # one step of 0.5 um on the G + g_p pitch (> 1e-2)
gt = BU.gates(sp_, sc_t, S.KAPPA_G * sp_.delta_star)
others = gt['M13_ok'] and gt['Re_ok'] and gt['stokes_ok'] and gt['plateau_ok'] and gt['lamina_ok']
check('seeded: one 0.5 um step of the centre width (0.5 / %.1f = %.2e > 1e-2) -> M16 fails alone on the gate function'
      % (sp_.G + sp_.g_p, 0.5 / (sp_.G + sp_.g_p)),
      not gt['taper_ok'] and others and gt['taper_rate_row'] == k_ + 1, 'rate %.3e at row %d' % (gt['taper_rate_max'], gt['taper_rate_row']))

# 11. The command line documents itself
pr = subprocess.run([sys.executable, '-m', 'opencif', 'generate', '--help'], cwd=ROOT, capture_output=True, text=True)
pv = subprocess.run([sys.executable, '-m', 'opencif', '--version'], cwd=ROOT, capture_output=True, text=True)
check('control: `generate --help` names units and choices; --version gives the code hash',
      pr.returncode == 0 and 'target cutoff diameter D (um)' in pr.stdout and '{separate,teardrop,ring}' in pr.stdout
      and PL.CODE['sha256'] in pv.stdout, pv.stdout.strip())
# 12. The spec rules of the flared outlet lead (the flare is drawn for the teardrop and separate outlets only)
def _spec_err(mut):
    """Apply mut to a copy of the example and return whether loading it raises a SpecError, with the message."""
    d_ = copy.deepcopy(gl)
    mut(d_)
    try:
        SP.to_engine(SP.load(d_))
        return False, ''
    except SP.SpecError as e:
        return True, str(e)
ok12a, m12a = _spec_err(lambda d_: d_.setdefault('outlet', {}).update(style='ring', flare_deg=14.5))
ok12b, m12b = _spec_err(lambda d_: d_.setdefault('outlet', {}).update(style='teardrop', flare_start_um=67))
ok12c, m12c = _spec_err(lambda d_: d_.setdefault('outlet', {}).update(style='teardrop', flare_deg=60))
check('seeded: a flare on the ring outlet, a flare_start without a flare, a 60 deg flare -> each a spec error',
      ok12a and 'teardrop and separate' in m12a and ok12b and 'flare_deg' in m12b and ok12c and '[1, 45]' in m12c, m12a[:60])
ok12d, m12d = _spec_err(lambda d_: d_.setdefault('outlet', {}).update(style='teardrop', flare_deg=14.5, stretch_um=300))
ok12e, m12e = _spec_err(lambda d_: (d_.setdefault('outlet', {}).update(style='teardrop', flare_deg=14.5),
                                    d_.setdefault('physics', {}).update(fail_safe_margin=True)))
ok12f, m12f = _spec_err(lambda d_: d_.setdefault('outlet', {}).update(style='teardrop', flare_deg=14.5))
check('seeded: a flare with outlet.stretch_um, a flare with the fail-safe margin -> each a spec error; control: the flare '
      'alone loads', ok12d and 'flare_start_um' in m12d and ok12e and 'fail_safe_margin' in m12e and not ok12f, m12d[:60])

# 12b. Two items deferred from the rc5 verification. First, reproduction.L_out with a flare is refused by name, because
#      the flare sets the stretch. Second, the domain of the outlet stretch-end tables is tested with the deep step
#      tables cut at R 1.5 (they end at R 8, and no shipped design comes near, since the largest ratio is 3.87), so that
#      m15's teardrop lead (R 2.87) lies beyond them. Strict mode then refuses it by name before drawing, whereas
#      evaluate mode plans it.
dl_ = SP.load(os.path.join(ROOT, 'presets', 'dinh2024_v1_depicted.yaml'))
dl_['reproduction']['L_out'] = 300.0
try:
    SP.to_engine(SP.load(dl_))
    ok12g, m12g = False, ''
except SP.SpecError as e:
    ok12g, m12g = True, str(e)
del dl_['reproduction']['L_out']
SP.to_engine(SP.load(dl_))                               # control: the preset itself loads
check('seeded: reproduction.L_out with outlet.flare_deg (Dinh V1) -> a spec error naming both; control: the preset loads',
      ok12g and 'reproduction.L_out' in m12g and 'flare_deg' in m12g, m12g[:70])
d12 = SP.load(os.path.join(ROOT, 'examples', 'matrix', 'm15_linear_shear_nomeander_filter_teardrop.yaml'))
one_, sym_ = S.STEP_ONE_DEEP, S.STEP_SYM_DEEP
S.STEP_ONE_DEEP = tuple(r_ for r_ in one_ if r_[0] <= 1.5)
S.STEP_SYM_DEEP = tuple(r_ for r_ in sym_ if r_[0] <= 1.5)
try:
    sp12, op12, sc12, g12_, e12_, n12_ = PL.prepare(copy.deepcopy(d12))
    d12e = copy.deepcopy(d12)
    d12e['gates'] = dict(mode='evaluate')
    sp12e, op12e, sc12e, g12e, e12e, n12e = PL.prepare(d12e)
    dom12 = (L.plan_layout(sc12e, op12e['layout'], 1, op12e.get('leg_spacing'))[0].info or {}).get('exp_domain')
finally:
    S.STEP_ONE_DEEP, S.STEP_SYM_DEEP = one_, sym_
check('seeded: the step tables cut at R 1.5 -> m15 refused before drawing in strict mode, by name (the ratio, the table '
      'end, the other modes)', PL.refusal(op12, e12_)[0] and 'outside the solved 2-D tables' in (op12.get('outlet_refusal') or '')
      and 'beyond the deep table\'s 1.5' in (op12.get('outlet_refusal') or '') and 'evaluate or reproduce mode'
      in (op12.get('outlet_refusal') or ''), (op12.get('outlet_refusal') or '')[:90])
check('  ... in evaluate mode it plans (nothing refused) and the accepted plan carries the status (the ratio, the table end)',
      not op12e.get('outlet_refusal') and not PL.refusal(op12e, e12e)[0] and dom12 is not None and '1.5' in dom12
      and '2.8' in dom12, str(dom12)[:70])

# 13. A user process of 12 um gaps, a 10 um floor, and 60 um depth with the narrow head (A's fresh device u05) crashed
#     rc5's drawn-mask measurement with a ZeroDivisionError in measure (Reviewer A's E5). The cause was the drawing,
#     because the head and outlet smoothing (closing then opening, with radius process.fillet_um) erases every channel
#     no wider than twice the radius, and the 5 um default did so at a 10 um floor. The default is therefore min(5,
#     0.4 x the narrowest drawn feature, 0.3 x the narrowest head unit), an explicit fillet at or above half the
#     narrowest feature is refused, and a feature that the measurement cannot find on the drawing is a named FAIL
#     rather than a crash.
u05 = api.build_spec(1.5, 0.50, process=dict(depth_um=60.0, gap_um=12.0, post_length_um=20.0, post_width_um=40.0,
                                              min_feature_um=10.0), chip=(25, 75), layout='linear', taper='constant_shear',
                     width_in_um=40.0, width_min_um=30.0, meander='narrow', inlet='plain', outlet='teardrop',
                     flow_uL_min=5.0, name='_t_e5')
u05f = copy.deepcopy(u05)
u05f['process']['fillet_um'] = 5.0
try:
    SP.to_engine(SP.load(u05f))
    ok13a, m13a = False, ''
except SP.SpecError as e:
    ok13a, m13a = True, str(e)
sp13, _ = SP.to_engine(SP.load(copy.deepcopy(u05)))
check('seeded: an explicit 5 um fillet at a 10 um floor -> a spec error naming the field and the feature; the default '
      'is 3 um (0.3 x the 10 um head units)', ok13a and 'process.fillet_um' in m13a and abs(sp13.r_fillet - 3.0) < 1e-9,
      m13a[:80])
prep_ok_ = PL.prepare
def _rc5_fillet(d_, say=lambda m: None):                 # rc5's 5 um fillet, forced past the spec rule
    out_ = prep_ok_(d_, say)
    out_[0].r_fillet = 5.0
    return out_
PL.prepare = _rc5_fillet
try:
    r13 = PL.run(copy.deepcopy(u05), outdir=T_OUT, verbose=False, render=False)
    crash13 = ''
except Exception as e:                                   # rc5 raised ZeroDivisionError here
    r13, crash13 = {}, '%s: %s' % (type(e).__name__, e)
finally:
    PL.prepare = prep_ok_
ff13 = (r13.get('validation') or {}).get('features_found') or {}
k13 = {q.split(' ')[0] + (' unit' if q.startswith('head') else '')  # every kind lost (the list keeps the first 20)
       for q in (ff13.get('by_kind') or {q_['kind']: 1 for q_ in ff13.get('not_found', [])})}
c13 = [c for c in (r13.get('assessment') or {}).get('checks', []) if c['check'] == 'every planned feature found on the drawing']
check('seeded: rc5\'s 5 um fillet forced past the rule -> drawn and measured to the end (no crash), the vanished head '
      'units and gaps listed by name, the build FAILs on that check', not crash13 and not ff13.get('ok', True)
      and {'head unit', 'gap'} <= k13 and c13 and not c13[0]['ok'] and (r13.get('assessment') or {}).get('exit_code') == 2,
      crash13[:80] or '%d not found, first %s' % (ff13.get('n_not_found', 0), (ff13.get('not_found') or [{}])[0]))
r13d = PL.run(copy.deepcopy(u05), outdir=T_OUT, verbose=False, render=False)
V13 = r13d.get('validation') or {}
t13 = V13.get('topology') or {}
check('control: the same at the default 3 um fillet -> every planned feature found, topology as planned (single '
      'component, holes as planned, every row kind drawn)', (V13.get('features_found') or {}).get('ok')
      and t13.get('single_component') and t13.get('holes_ok') and t13.get('row_kinds_survive'),
      (r13d.get('assessment') or {}).get('verdict', '')[:80])
leg_, gx_ = info['legs'][0], float(info['gx'][0])          # m15's drawing (test 2), with one gap filled on side A
k_g = sc_d.N // 2
Xg = gx_ + (float(sc_d.x[k_g]) - float(leg_['p0'] or 0.0)) * leg_['dirn']
Yg = float(leg_['Y']) + leg_['dirn'] * (float(sc_d.wc[k_g]) / 2.0 + spec.l_g / 2.0)
hole = box(Xg - spec.G / 2 - 2.0, Yg - spec.l_g / 2 + 1.0, Xg + spec.G / 2 + 2.0, Yg + spec.l_g / 2 - 1.0)
from opencif.checks import measure as M2
m_clean = M2.measure_device(fluid, sc_d, info)
m_fill = M2.measure_device(fluid.difference(hole), sc_d, info)
nf = m_fill['not_found']
check('seeded: one gap filled on m15\'s drawing -> the measurement names it (row, gap, side A) and only it; control: '
      'the clean drawing finds every feature', not m_clean['not_found'] and len(nf) == 1 and nf[0]['kind'] == 'gap'
      and nf[0]['row'] == k_g + 1 and nf[0]['side'] == 'A', str(nf[:1])[:90])

print('%d checks' % n_checks[0])
print('ALL OK' if not fails else 'FAILED: %s' % fails)
sys.exit(1 if fails else 0)
