"""test_units.py: fast unit tests, which draw nothing.

Run it with `python3 tests/test_units.py`."""
import copy, os, sys
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import tempfile as _tempfile
# Scratch output goes outside the package tree, so that a read-only install works, and OPENCIF_TEST_OUT overrides
# the location (review R3 of rc4).
T_OUT = os.environ.get('OPENCIF_TEST_OUT') or _tempfile.mkdtemp(prefix='opencif_test_')
if not os.environ.get('OPENCIF_TEST_OUT'):               # a temporary directory of our own, removed at exit
    import atexit as _atexit, shutil as _shutil
    _atexit.register(_shutil.rmtree, T_OUT, True)

sys.path.insert(0, ROOT)
from opencif import spec as SP, pipeline as PL
from opencif import pipeline as _PL_CODE                # noqa: E402  (every suite's log names the code it ran on)
print('opencif %s (code %s)' % (_PL_CODE.CODE['version'], _PL_CODE.CODE['sha256']))
from opencif.physics import schedule as S
from opencif.geometry import layout as L

fails = []
def check(name, cond, info=''):
    print('%-62s %s %s' % (name, 'ok ' if cond else 'FAIL', info))
    if not cond:
        fails.append(name)

base = SP.load(os.path.join(ROOT, 'examples', 'gifford_physics_linear.yaml'))
# 1. The spec rejects unknown keys and values, and missing required inputs
for bad, what in ((dict(process=dict(depth_um=80, gap_um=20, min_feature_um=4.5, colour=1)), 'unknown key'),
                  (dict(layout=dict(style='spiral')), 'unknown layout'),
                  (dict(centre_channel=dict(taper='linear', width_in_um=100)), 'unknown taper')):
    d = copy.deepcopy(base); d.update(bad)
    try:
        SP.load(d); check('spec rejects %s' % what, False)
    except SP.SpecError:
        check('spec rejects %s' % what, True)
d = copy.deepcopy(base); del d['process']['min_feature_um']
try:
    SP.load(d); check('spec requires the process floor', False)
except SP.SpecError:
    check('spec requires the process floor', True)
# 2. D enters only through delta* = D / kappa_g and F only through the stop, so changing D at fixed F keeps F, and
#    vice versa
sp, _ = SP.to_engine(base)
sc1 = S.plan(sp)
check('cutoff of every designed gap = D (constant w_c)', np.allclose(S.KAPPA_G * np.array([S.x_of(f) for f in sc1.f]) * sc1.wc, 6.0, rtol=1e-12))
d2 = copy.deepcopy(base); d2['targets']['cutoff_um'] = 7.0
sc2 = S.plan(SP.to_engine(d2)[0])
check('D changes N, not F', abs(sc2.F - 0.8) < 2 * np.max(sc2.f) and sc2.N != sc1.N, 'N %d -> %d' % (sc1.N, sc2.N))
d3 = copy.deepcopy(base); d3['targets']['filtrate_fraction'] = 0.6
sc3 = S.plan(SP.to_engine(d3)[0])
check('F changes N, not the per-gap cutoff', np.allclose(sc3.f, sc1.f[:sc3.N]) and sc3.N < sc1.N, 'N %d -> %d' % (sc1.N, sc3.N))
# 3. The resistance-exact wall, in which a segment's drawn profile has exactly the closure resistance
for n0, n1, w, ell, h in ((4.857, 5.5, 4.857, 46.1, 80), (20.0, 23.0, 21.0, 46.1, 80), (300.0, 290.0, 296.0, 56.25, 140)):
    m = L.solve_mid_width(n0, n1, w, ell, h)
    rel = L.seg_profile_R(n0, m, n1, ell, h) / (S.rr(w, h) * ell) - 1
    check('wall profile resistance exact (w %.1f)' % w, abs(rel) < 1e-10, 'rel %.1e, mid %.3f' % (rel, m))
# 4. The floor start begins the design at the requested side width
gf = SP.load(os.path.join(ROOT, 'presets', 'gifford2014_fig4_as_fabricated_20um.yaml'))
spf, of = SP.to_engine(gf)
spf2 = PL.floor_start(spf, of, [], lambda m: None)
scf = S.plan(spf2)
check('floor start draws w_s(1) at the floor', abs(scf.ws[0] - 20.0) < 1e-6, 'w_s(1) %.6f, S0 %.3e' % (scf.ws[0], spf2.S0))
check('floor start stops on the DELIVERED F (Gifford stops on the final width)', abs(scf.F - 0.90644) < 2 * scf.f[-1] and
      abs(scf.ws[scf.N - 2] - 300.0) < 0.5, 'N %d, F %.6f, final w_s %.2f um (printed 300)' % (scf.N, scf.F, scf.ws[scf.N - 2]))
# 5. A sealed start without a head below the floor is refused with advice
d4 = copy.deepcopy(base); d4['process']['min_feature_um'] = 19.0
try:
    PL.run(d4, outdir=T_OUT, verbose=False, render=False); check('sub-floor sealed start refused', False)
except SP.SpecError as e:
    check('sub-floor sealed start refused', 'meander' in str(e))
# 6. The curved-duct correction: the exact Stokes solution against its two closed-form limits and the straight limit
import math
from opencif import build as BU, compose as CO


def planar(r1, r2):
    """Return the deep (planar-annulus) limit of the bend ratio between the radii r1 and r2, that is, the closed-form
    flux of the annulus over that of the straight channel with the same centerline gradient."""
    c = -0.5
    Mx = np.array([[r1, 1 / r1], [r2, 1 / r2]])
    A, B = np.linalg.solve(Mx, -np.array([c * r1 * math.log(r1), c * r2 * math.log(r2)]))
    F_ = lambda r: c * (r * r / 2 * math.log(r) - r * r / 4) + A * r * r / 2 + B * math.log(r)
    return (F_(r2) - F_(r1)) / ((2.0 / (r1 + r2)) * (r2 - r1) ** 3 / 12.0)


for r1, r2 in ((431.7, 531.7), (20.5, 43.3)):
    w = r2 - r1
    deep, shallow = planar(r1, r2), 0.5 * (r1 + r2) * math.log(r2 / r1) / w
    e_deep, e_sh = S.bend_ratio(r1, r2, 1e4 * w, n_modes=3000), S.bend_ratio(r1, r2, 1e-5 * w)
    check('bend ratio -> planar annulus (deep), Hele-Shaw (shallow), r1 %.1f' % r1,
          abs(e_deep / deep - 1) < 2e-5 and abs(e_sh / shallow - 1) < 2e-5,
          'deep %.6f/%.6f shallow %.6f/%.6f' % (e_deep, deep, e_sh, shallow))
check('bend ratio -> 1 for a straight duct', abs(S.bend_ratio(1e6, 1e6 + 100, 80) - 1) < 1e-8)
# 7. The realized-cutoff metric, which passes f* everywhere and is tested on one row off, a reverse-flow row, and a
#    declared row
N = 50
wc = np.full(N, 100.0)
fs = S.f_of(1.0, 100.0)
dd = S.KAPPA_G * S.x_of(fs) * 100.0
f = np.full(N, fs)
mk = lambda fA, fB: dict(fA=fA, fB=fB, dA=np.array([S.KAPPA_G * S.x_of(v) * 100 if 0 < v < .5 else np.nan for v in fA]),
                         dB=np.array([S.KAPPA_G * S.x_of(v) * 100 if 0 < v < .5 else np.nan for v in fB]), F=0.5)
r0 = BU.realised_cutoff(None, mk(f, f), dd, 0.01)
fA = f.copy(); fA[10] *= 1.05; fA[20] = -1e-5
r1 = BU.realised_cutoff(None, mk(fA, f), dd, 0.01, declared={'seed': [20]})
check('realised cutoff: clean rows pass, a 5% f error and a reverse row are caught, declared rows reported apart',
      r0['undeclared_rows_outside'] == 0 and r1['undeclared_rows_outside'] == 1 and r1['undeclared_rows_list'] == [11]
      and r1['declared_by_option']['seed']['undefined_rows'] == 1 and r1['reverse_flow_rows'] == [21],
      '%d / %d, %s' % (r0['undeclared_rows_outside'], r1['undeclared_rows_outside'], r1['undeclared_rows_list']))
fid = BU.drawing_fidelity(mk(fA, f), mk(fA, f), np.full(N, fs), N)
fid2 = BU.drawing_fidelity(mk(fA, f), mk(f, f), np.full(N, fs), N)
check('drawing fidelity: identical networks 0; a 5% f error at row 11 is found', fid['cutoff_dev_max'] == 0 and
      fid2['cutoff_dev_row'] == 11 and fid2['cutoff_dev_max'] > 0.02, '%.4f at %d' % (fid2['cutoff_dev_max'], fid2['cutoff_dev_row']))
# 8. The compose module takes only a clean pass whose photomask passed the physics auditor
for v, c, pa, want in (('PASS', 0, 'PASS', True), ('REPRODUCED with declared departures (seed)', 0, 'PASS', False),
                       ('PASS -- OUTSIDE THE VALIDITY ENVELOPE (u_min 0.03)', 3, 'PASS', False),
                       ('PASS with an allowed departure (seed)', 0, 'PASS', False), ('FAIL (Q1)', 2, 'FAIL', False),
                       ('PASS', 0, None, False), ('PASS', 0, 'PASS_OUTSIDE_ENVELOPE', False)):
    rec_ = dict(assessment=dict(verdict=v, exit_code=c), validation=dict(physics_audit=dict(status=pa)) if pa else {})
    check('compose element gate: %s, audit %s' % (v[:40], pa), CO.element_ok(rec_)[0] == want)
# 9. Declared rows are attributed rather than zoned, so an option exempts only the rows it moves, and a fault that no
#    option explains is not exempt
from types import SimpleNamespace
spx = SimpleNamespace(f_head=None, transition='paper', head='pack', loop_pull_in=0.0, S0=0.0)
fx = np.full(N, fs); fx[0] = 3 * fs
full_f = fx.copy(); full_f[[5, 6, 7, 30]] *= 1.3
nopaper_f = fx.copy(); nopaper_f[30] *= 1.3
nets = dict(full=mk(full_f, full_f), no_paper=mk(nopaper_f, nopaper_f), none=mk(nopaper_f, nopaper_f))
dec = BU.declared_rows(spx, fx, np.full(N, fs), [], dd, 0.01, nets, 'reproduce', 'require')
check('declared rows: seed = gap 1 only, transition = its own rows, a fault at row 31 undeclared',
      dec.get(BU.OPT_SEED) == [0] and dec.get(BU.OPT_PAPER) == [5, 6, 7] and all(30 not in v for v in dec.values()), str(dec))
dec_e = BU.declared_rows(spx, fx, np.full(N, fs), [], dd, 0.01, nets, 'evaluate', 'require')
check('declared rows: nothing in evaluate mode without seal: allow', dec_e == {}, str(dec_e))
# 10. The filtrate-to-retentate flow ratio r = F/(1-F) (which is not the glossary's FR) gives F = r / (1 + r)
dfr = copy.deepcopy(base); del dfr['targets']['filtrate_fraction']; dfr['targets']['filtrate_to_retentate_ratio'] = 4.0
check('targets.filtrate_to_retentate_ratio 4 -> F 0.8', abs(SP.load(dfr)['targets']['filtrate_fraction'] - 0.8) < 1e-15)
dfr2 = copy.deepcopy(base); dfr2['targets']['filtrate_to_retentate_ratio'] = 3.0
try:
    SP.load(dfr2); check('F and FR that disagree are refused', False)
except SP.SpecError:
    check('F and FR that disagree are refused', True)
# 11. A meander head needs at least one rectilinear row under its transition unit's drop (review row 40)
drr = copy.deepcopy(base); drr['meander'] = dict(enabled=True, style='narrow', rectilinear_rows=0)
try:
    SP.to_engine(SP.load(drr)); check('meander.rectilinear_rows 0 is refused', False)
except SP.SpecError:
    check('meander.rectilinear_rows 0 is refused', True)
# 12. The length of the flared stretch (layout.stretch_total). The straight run is flare_start where no U-turn precedes
#     the outlet, the layout's clearance stretch where that is longer, and the spec's stretch when the lead does not
#     widen.
from types import SimpleNamespace as _NS
spf = _NS(flare_deg=14.5, flare_start=67.0)
if True:
    Lf_ = L.flare_len(spf, 44.0, 105.75)
    tail_ = Lf_ + L.flare_kink(spf, Lf_)[1]
    check('flared stretch: linear = flare_start + flare + end arc; clearance = the longer of the two; no widening = the '
          'spec\'s stretch', abs(L.stretch_total(spf, 'linear', 200.0, 44.0, 105.75) - (67.0 + tail_)) < 1e-9 and
          abs(L.stretch_total(spf, 'clearance', 1205.0, 44.0, 105.75) - (1205.0 + tail_)) < 1e-9 and
          abs(L.stretch_total(spf, 'clearance', 30.0, 44.0, 105.75) - (67.0 + tail_)) < 1e-9 and
          L.stretch_total(spf, 'linear', 200.0, 44.0, 44.0) == 200.0 and L.stretch_total(spf, 'linear', 200.0, 50.0, 44.0) == 200.0,
          'L_f %.2f um, tail %.2f um' % (Lf_, tail_))
# 12b. The layout's choice of rule (layout._flare_fixed_point). One leg of a serpentine is 'linear' (flare_start) and
#      three legs are 'clearance' (the loop's stretch), both planned without drawing on m05 with a flare added.
dfl = SP.load(os.path.join(ROOT, 'examples', 'matrix', 'm05_serpentine_const_arms_plain_teardrop.yaml'))
dfl.setdefault('outlet', {})['flare_deg'] = 14.5
spl, _ = SP.to_engine(SP.load(dfl))
scl = S.plan(spl)
res_ = {}
for legs_ in (1, 3):
    sc_l, _ = L.plan_layout(scl, 'serpentine', legs_, None)
    i_ = sc_l.info
    tl_ = i_['flare_len'] + L.flare_kink(sc_l.spec, i_['flare_len'])[1]
    res_[legs_] = (i_.get('flare_kind'), sc_l.spec.L_out, i_['L_out_base'], tl_)
check('flared stretch on a serpentine: one leg flare_start + flare + arc, three legs the loop\'s stretch + flare + arc',
      res_[1][0] == 'linear' and abs(res_[1][1] - (spl.flare_start + res_[1][3])) < 1e-6 and res_[3][0] == 'clearance' and
      res_[3][2] > spl.flare_start and abs(res_[3][1] - (res_[3][2] + res_[3][3])) < 1e-6,
      '1 leg %.2f um, 3 legs %.2f um (base %.2f)' % (res_[1][1], res_[3][1], res_[3][2]))
# 13. The domain of the stretch-end tables (deferred from the rc5 verification). The deep step table ends at a width
#     ratio of 8 and the flare tables at 4 (the spec limits flare_deg to [1, 45]). Beyond them np.interp holds the end
#     value, so the status names it, and strict mode refuses (spec.exp_domain_refuse) whereas the other modes report.
d_one = S.STEP_ONE_DEEP[-1][0]
d_fl = S.FLARE_RATIOS[-1]
sn = _NS(junction_band=None, flare_deg=None)
sf = _NS(junction_band=None, flare_deg=14.5)
dom_ = [S.outlet_exp_domain(sn, 10.0, 10.0 * d_one), S.outlet_exp_domain(sn, 10.0, 10.0 * (d_one + 0.01)),
        S.outlet_exp_domain(sf, 10.0, 10.0 * d_fl), S.outlet_exp_domain(sf, 10.0, 10.0 * (d_fl + 0.01)),
        S.outlet_exp_domain(_NS(junction_band='none', flare_deg=None), 10.0, 200.0)]
check('stretch-end domain: step R %.1f inside, %.2f beyond; flare R %.1f inside, %.2f beyond; no junction band: none'
      % (d_one, d_one + 0.01, d_fl, d_fl + 0.01), dom_[0] is None and dom_[1] is not None and 'deep table' in dom_[1]
      and dom_[2] is None and dom_[3] is not None and 'flare width ratio' in dom_[3] and dom_[4] is None, str(dom_[1])[:60])
sr = _NS(junction_band=None, flare_deg=None, exp_domain_refuse=True)
se = _NS(junction_band=None, flare_deg=None, exp_domain_refuse=False)
try:
    L.check_exp_domain(sr, 10.0, 10.0 * (d_one + 0.01))
    raised_ = False
except L.OutletOutsideTables:
    raised_ = True
check('check_exp_domain: refuses beyond the table in strict mode, returns the status otherwise, None inside',
      raised_ and 'deep table' in (L.check_exp_domain(se, 10.0, 10.0 * (d_one + 0.01)) or '')
      and L.check_exp_domain(sr, 10.0, 10.0 * d_one) is None)
print('ALL OK' if not fails else 'FAILED: %s' % fails)
sys.exit(1 if fails else 0)
