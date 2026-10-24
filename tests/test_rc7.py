"""test_rc7.py: fast tests of the additions made in rc7, none of which draws a whole device.

Run it with `python3 tests/test_rc7.py`. The suite checks these items:
  BUILD-1     The memo caches key on rounded inputs (a nearby input returns the cached value until they are cleared),
              and pipeline.run clears them before anything else.
  XIA-1       centre_channel.shear_formula gives the parallel-plate shear, its table's inverse, and a taper w_c
              proportional to the center flow, and 'exact' (the default) is unchanged.
  DINH-1      meander.turns set to square gives the corner term c_corner (table nodes, depth weight, and band), a
              square Path's length through its corner points, the drawn corner's roundings, and the spec's checks,
              and the arc Path is unchanged.
  defaults    The options added in rc7 are listed among a record's defaults only where they apply.
  ADJ-C4..C9  The adjudication items of the rc6 reviews (not App. B's claims C1 to C9) are the post-dimension aliases,
              the fluid section, the status codes and the source trace, the environment stamp, the defaults' origins
              and effective values, and the passage report."""
import copy
import json
import math
import os
import subprocess
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from opencif import spec as SP, pipeline as PL           # noqa: E402
print('opencif %s (code %s)' % (PL.CODE['version'], PL.CODE['sha256']))
from opencif.physics import schedule as S  # noqa: E402
from opencif.geometry import shapes as Gm  # noqa: E402

fails = []


def check(name, cond, info=''):
    print('%-66s %s %s' % (name, 'ok ' if cond else 'FAIL', info))
    if not cond:
        fails.append(name)


PRESETS = os.path.join(ROOT, 'presets')

# ---------------------------------------------------------------------------------------------------- BUILD-1
S.bend_ratio(100.0, 128.0, 140.0)
S.depth_lambda(140.0, 28.0)
check('BUILD-1: the caches hold entries before clearing', len(S._BEND_CACHE) > 0 and len(S._LAM_CACHE) > 0)
S.clear_caches()
check('BUILD-1: clear_caches empties both', len(S._BEND_CACHE) == 0 and len(S._LAM_CACHE) == 0)
src = open(os.path.join(ROOT, 'opencif', 'pipeline.py')).read()
i = src.index('def run(user, outdir=None, verbose=True, render=True):')
check('BUILD-1: pipeline.run clears the caches first (source)', 'S.clear_caches()' in src[i:i + 300])
# The caches key on rounded inputs, so a nearby input returns the value first computed (code review R12, behavior).
S.clear_caches()
v1 = S.bend_ratio(100.0, 128.0, 140.0)
v2 = S.bend_ratio(100.00004, 128.0, 140.0)
S.clear_caches()
v3 = S.bend_ratio(100.00004, 128.0, 140.0)
check('BUILD-1: a nearby input returns the cached value until the caches are cleared', v2 == v1 and v3 != v1,
      'cached %r, fresh %r' % (v1, v3))


class _Stop(Exception):
    """The exception that the stand-in for clear_caches raises."""
    pass


def _stop():
    """Stand in for clear_caches and stop the build at once."""
    raise _Stop()


_orig = S.clear_caches
S.clear_caches = _stop
try:
    import tempfile as _tf                                # a failure must not build into runs/ (codereview2 N7)
    PL.run(SP.load(os.path.join(PRESETS, 'gifford2014_fig4_design.yaml')), outdir=_tf.mkdtemp(prefix='build1_'),
           verbose=False, render=False)
    stopped = False
except _Stop:
    stopped = True
finally:
    S.clear_caches = _orig
check('BUILD-1: pipeline.run calls clear_caches before building (behaviour)', stopped)
# The earlier check of rc7 that "a plan after another equals a fresh process" is removed, because it cleared the
# caches itself and the planner does not depend on them, so it could not fail (codereview2). The cache behavior and
# run()'s call are tested above, and the order independence of the full builds is tested by the rebuild against rc6.
sp_a, _ = SP.to_engine(SP.load(os.path.join(PRESETS, 'gifford2014_fig4_design.yaml')))
sp_b, _ = SP.to_engine(SP.load(os.path.join(PRESETS, 'gifford2014_fig6_design.yaml')))

# ---------------------------------------------------------------------------------------------------- XIA-1
check('XIA-1: parallel-plate shear 6 Q / (a b^2)', S.shear(1.0, 200.0, 140.0, 'parallel_plates') == 6.0 / (200.0 * 140.0 ** 2))
check('XIA-1: narrow channel uses a = h, b = w', S.shear(1.0, 50.0, 140.0, 'parallel_plates') == 6.0 / (140.0 * 50.0 ** 2))
check('XIA-1: exact is the default and unchanged', S.shear(1.0, 200.0, 140.0) == S.C.C_shear(200.0 / 140.0) / (200.0 * 140.0 ** 2))
st = S.ShearTable(140.0, 150.0, 300.0, formula='parallel_plates')
w_ = st.wc_for(S.shear(1.0, 222.0, 140.0, 'parallel_plates'))
check('XIA-1: the table inverts the parallel-plate shear', abs(w_ - 222.0) < 1e-6, '%.9f' % w_)
d = SP.load(os.path.join(PRESETS, 'xia2016_rs7_aligned.yaml'))
sp, _ = SP.to_engine(d)
check('XIA-1: the aligned Xia preset takes parallel_plates', sp.shear_formula == 'parallel_plates')
sc = S.plan(sp)
wc, P = np.asarray(sc.wc), np.asarray(sc.P)
mid = [i for i in range(sc.NM + 1, sc.N) if sp.wc_min + 1e-6 < wc[i] < sp.wc0 - 1e-6]
ratio = wc[mid] / P[np.array(mid) - 1]
check('XIA-1: w_c / (centre flow) constant along the taper', len(mid) > 100 and float(np.ptp(ratio)) < 1e-6 * float(ratio.mean()),
      'rows %d, spread %.2e' % (len(mid), float(np.ptp(ratio))))
try:
    bad = copy.deepcopy(SP.load(os.path.join(PRESETS, 'gifford2014_fig4_design.yaml')))
    bad.setdefault('centre_channel', {})['shear_formula'] = 'parallel_plates'
    SP.to_engine(SP.load(bad))
    check('XIA-1: shear_formula refused without the constant_shear taper', False)
except SP.SpecError:
    check('XIA-1: shear_formula refused without the constant_shear taper', True)
bad2 = copy.deepcopy(SP.load(os.path.join(PRESETS, 'xia2016_rs7_aligned.yaml')))
bad2['gates'] = dict(mode='evaluate')
try:
    SP.to_engine(SP.load(bad2)); check('XIA-1: parallel_plates refused outside reproduce mode (spec v5 S16-b)', False)
except SP.SpecError:
    check('XIA-1: parallel_plates refused outside reproduce mode (spec v5 S16-b)', True)
check('XIA-1: parallel_plates recorded as a reproduction-only departure',
      'shear_formula: parallel_plates' in SP.to_engine(SP.load(os.path.join(PRESETS, 'xia2016_rs7_aligned.yaml')))[1]['repro_only'])

# ---------------------------------------------------------------------------------------------------- DINH-1
dd = SP.load(os.path.join(PRESETS, 'dinh2024_v1_depicted.yaml'))
sp, _ = SP.to_engine(dd)
check('DINH-1: the Dinh preset takes square turns', sp.turns == 'square' and sp.r_corner_out == 7.0 and sp.r_corner_in == 10.0)
sp0 = copy.copy(sp)
for ri, deep, shal in zip(S.CORNER_RIN, S.CORNER_DEEP, S.CORNER_SHALLOW):
    sp0.r_corner_in = ri * 28.0
    sp0.junction_band = 'deep'
    c_d = S.c_corner(sp0, 28.0)
    sp0.junction_band = 'shallow'
    c_s = S.c_corner(sp0, 28.0)
    check('DINH-1: c_corner at r_in/w %.2f: the table (deep %.4f, shallow %.4f)' % (ri, deep, shal),
          abs(c_d - deep) < 1e-12 and abs(c_s - shal) < 1e-12)
sp0.junction_band = None
sp0.r_corner_in = 10.0
c0, dp, sh, lam, held = S.c_corner(sp0, 28.0, parts=True)
check('DINH-1: c_corner between its limits by the depth weight', abs(c0 - (dp + lam * (sh - dp))) < 1e-12 and not held)
sp0.junction_band = 'lo'; c_lo = S.c_corner(sp0, 28.0)
sp0.junction_band = 'hi'; c_hi = S.c_corner(sp0, 28.0)
check('DINH-1: the band brackets the value', c_hi < c0 < c_lo, '%.4f < %.4f < %.4f' % (c_hi, c0, c_lo))
sp0.junction_band = 'none'
check('DINH-1: no corner term without junction terms', S.c_corner(sp0, 28.0) == 0.0)
sp0.junction_band = None
sp0.bend_model = 'none'
check('DINH-1: no corner term under the plain 1-D rule (bend_model none)', S.c_corner(sp0, 28.0) == 0.0)
sp0.bend_model = 'stokes'
p = Gm.Path(0.0, 0.0, 90.0, square=True)
p.line(100.0); p.turn(27.5, 90.0); p.line(50.0)
check('DINH-1: a square turn runs r to the corner and r on', abs(p.length() - (100.0 + 55.0 + 50.0)) < 1e-9 and p.corners == [2])
q = Gm.Path(0.0, 0.0, 90.0)
q.line(100.0); q.turn(27.5, 90.0); q.line(50.0)
r_ = Gm.Path(0.0, 0.0, 90.0)
r_.line(100.0); r_.arc(27.5, 90.0); r_.line(50.0)
check('DINH-1: turn() on an arc path draws rc6\'s arc exactly', q.p == r_.p and q.arcs == r_.arcs and not q.corners)
check('DINH-1: the bend excess counts each corner once', abs(Gm.bend_excess(sp, p, 28.0) - S.c_corner(sp, 28.0) * 28.0) < 1e-12)
check('DINH-1: quarter-turn lengths', Gm.qturn_len(sp, 10.0) == 20.0 and Gm.qturn_len(sp_b, 10.0) == math.pi * 10.0 / 2.0)
g = Gm.ribbon(p, 28.0, sp)
g_sharp = Gm.ribbon(p, 28.0, None)
# The corner point is C = (0, 127.5), where the path heading +y turns to -x (a left turn). The outer corner is at
# C + (w/2)(u - v) = (14, 141.5) and the inner corner at (-14, 113.5). The rounding removes area at the outer corner,
# and the fillet adds area at the inner one.
out_loss = g_sharp.difference(g).area
in_gain = g.difference(g_sharp).area
check('DINH-1: the outer corner loses (1 - pi/4) r_out^2', abs(out_loss - (1 - math.pi / 4) * 49.0) < 0.2, '%.2f um2' % out_loss)
check('DINH-1: the inner corner gains (1 - pi/4) r_in^2', abs(in_gain - (1 - math.pi / 4) * 100.0) < 0.3, '%.2f um2' % in_gain)
for bad_, what in ((dict(meander=dict(enabled=False, turns='square')), 'square without a head'),):
    d2 = copy.deepcopy(SP.load(os.path.join(PRESETS, 'gifford2014_fig4_design.yaml')))
    d2.update(bad_)
    try:
        SP.to_engine(SP.load(d2)); check('DINH-1: spec refuses %s' % what, False)
    except SP.SpecError:
        check('DINH-1: spec refuses %s' % what, True)
d3 = copy.deepcopy(dd); d3['meander'] = dict(d3['meander']); d3['meander']['turns'] = 'arc'
try:
    SP.to_engine(SP.load(d3)); check('DINH-1: spec refuses corner radii with arc turns', False)
except SP.SpecError:
    check('DINH-1: spec refuses corner radii with arc turns', True)

# ---------------------------------------------------------------------------------------------------- defaults
fig4 = SP.load(os.path.join(PRESETS, 'gifford2014_fig4_design.yaml'))
check('defaults: shear_formula not listed for a constant taper', not any('shear_formula' in x for x in SP.defaults_used(fig4)))
xia = SP.load(os.path.join(PRESETS, 'xia2016_rs7_depicted.yaml'))
check('defaults: shear_formula listed for the constant_shear taper', any('shear_formula' in x for x in SP.defaults_used(xia)))

# --------------------------------------------------------------------------------------------------- DINH-1 measurement
check('DINH-1: corner_round_area = the drawn ribbon\'s net rounding area', abs((g.area - g_sharp.area) - Gm.corner_round_area(sp)) < 1e-6,
      '%.4f um2' % Gm.corner_round_area(sp))
check('DINH-1: no rounding area without square-turn radii', Gm.corner_round_area(sp_b) == 0.0)
m2src = open(os.path.join(ROOT, 'opencif', 'checks', 'measure.py')).read()
check('DINH-1: measure takes the rounding area off before area -> length', 'Lm = (area - a_round) / wm + excl' in m2src)

# ---------------------------------------------------------------------------------------------------- C4 aliases
base = SP.load(os.path.join(PRESETS, 'gifford2014_fig4_design.yaml'))
d4 = copy.deepcopy(base); pr = d4['process']
gp_, lg_ = pr.get('post_length_um', 25.0), pr.get('post_width_um', 50.0)
pr.pop('post_length_um', None); pr.pop('post_width_um', None)
pr['wall_along_flow_um'] = gp_; pr['wall_cross_flow_um'] = lg_
l4 = SP.load(d4)
check('ADJ-C4: the aliases map onto post_length_um / post_width_um', l4['process'].get('post_length_um') == gp_ and
      l4['process'].get('post_width_um') == lg_ and 'wall_along_flow_um' not in l4['process'])
sp4, _ = SP.to_engine(l4); sp4b, _ = SP.to_engine(base)
check('ADJ-C4: the same engine spec as the YAML keys', sp4.g_p == sp4b.g_p and sp4.l_g == sp4b.l_g)
d4b = copy.deepcopy(base); d4b['process']['wall_along_flow_um'] = float(base['process'].get('post_length_um', 25.0)) + 1.0
try:
    SP.load(d4b); check('ADJ-C4: an alias that disagrees with its key is refused', False)
except SP.SpecError:
    check('ADJ-C4: an alias that disagrees with its key is refused', True)
d4c = copy.deepcopy(base); d4c['process']['wall_along_flow_um'] = float(base['process'].get('post_length_um', 25.0))
check('ADJ-C4: an alias equal to its key is accepted', SP.load(d4c)['process']['post_length_um'] == base['process'].get('post_length_um', 25.0))

# ---------------------------------------------------------------------------------------------------- C5 fluid
from opencif import build as BU                           # noqa: E402
check('ADJ-C5: water by default', sp_a.rho_f == 1000.0 and sp_a.mu_f == 1.0e-3)
d5 = copy.deepcopy(base); d5['fluid'] = dict(density_kg_m3=1025.0, viscosity_Pa_s=2.0e-3)
sp5, _ = SP.to_engine(SP.load(d5))
check('ADJ-C5: the fluid section reaches the engine spec', sp5.rho_f == 1025.0 and sp5.mu_f == 2.0e-3)
S.clear_caches(); sc_w = S.plan(sp4b)
S.clear_caches(); sc_5 = S.plan(sp5)
check('ADJ-C5: the schedule does not depend on the fluid', list(sc_w.ws[:-1]) == list(sc_5.ws[:-1]) and list(sc_w.wc) == list(sc_5.wc))
g_w = BU.gates(sp4b, sc_w, S.KAPPA_G * sp4b.delta_star)
g_5 = BU.gates(sp5, sc_5, S.KAPPA_G * sp5.delta_star)
ratio_ = g_5['Re_max'] / g_w['Re_max']
check('ADJ-C5: the Reynolds gates scale with rho/mu', abs(ratio_ - 1025.0 / 1000.0 / 2.0) < 1e-12, 'Re ratio %.6f' % ratio_)
check('ADJ-C5: the record\'s fluid names it', g_5['fluid']['mu_Pa_s'] == 2.0e-3 and 'spec' in g_5['fluid']['note'] and
      'water' in g_w['fluid']['note'])
for bad5, what in ((dict(density_kg_m3=0.0), 'a zero density'), (dict(rheology='carreau'), 'a non-Newtonian rheology')):
    d5b = copy.deepcopy(base); d5b['fluid'] = bad5
    try:
        SP.load(d5b); check('ADJ-C5: refuses %s' % what, False)
    except SP.SpecError:
        check('ADJ-C5: refuses %s' % what, True)

# ---------------------------------------------------------------------------------------------------- C6 status codes
check('ADJ-C6: every status field has its enumerated codes', set(PL.STATUS_CODES) >= {
    'hard_checks', 'validity_envelope', 'chip_area', 'aspect_ratio_of_CIF_channels', 'length_budget_G9',
    'drawn_mask_assessment', 'physics_audit', 'fabricated_or_measured', 'reproduction', 'source_trace'})
d6 = copy.deepcopy(base); d6.setdefault('paper', {})['source_trace'] = 'provenance/rc7_traces/x.json'
check('ADJ-C6: paper.source_trace reaches the options', SP.to_engine(SP.load(d6))[1]['source_trace'] == 'provenance/rc7_traces/x.json'
      and SP.to_engine(SP.load(os.path.join(ROOT, 'examples', 'gifford_physics_linear.yaml')))[1]['source_trace'] is None)
check('ADJ-C6 (R10): the traced presets point to their rc7 trace records', all(
    str(SP.to_engine(SP.load(os.path.join(PRESETS, n_ + '.yaml')))[1]['source_trace']).startswith('provenance/rc7_traces/results/')
    for n_ in ('dinh2024_v1_depicted', 'dinh2024_v1_figure_aligned', 'gifford2014_fig4_design', 'gifford2014_fig4_aligned',
               'gifford2014_fig6_design', 'gifford2014_fig6_aligned', 'xia2016_rs7_depicted', 'xia2016_rs7_aligned')))
# Every pointer resolves to trace records of the shipped mask, a check that runs once the rc7 traces exist
# (codereview2 N1).
TRACE_FILES = {'dinh2024_v1_depicted': ['dinh/dinh_results_rc7.json'], 'dinh2024_v1_figure_aligned': ['dinh/dinh_results_fa.json'],
               'gifford2014_fig4_design': ['gifford/gifford_results_rc7.json'],
               'gifford2014_fig4_aligned': ['gifford/gifford_results_aligned.json'],
               'gifford2014_fig6_design': ['gifford/gifford6_results_design.json', 'gifford/gifford6_results_design_onemap.json'],
               'gifford2014_fig6_aligned': ['gifford/gifford6_results_aligned.json', 'gifford/gifford6_results_aligned_onemap.json'],
               'xia2016_rs7_depicted': ['xia/xia_results_rc7.json'] + ['xia/xia_results_fig3%s_rc7.json' % L_ for L_ in 'cdef']
               + ['xia/xia_fig3a_device_rc7.json'],
               'xia2016_rs7_aligned': ['xia/xia_results_pp.json'] + ['xia/xia_results_fig3%s_pp.json' % L_ for L_ in 'cdef']
               + ['xia/xia_fig3a_device_pp.json']}


def _braces(t):
    """Expand brace groups ('a{b,c}d' -> 'abd acd'), so that a pointer may name several files in one (codereview3
    M1)."""
    import re as re_
    m_ = re_.search(r'\{([^{}]*)\}', t)
    if not m_:
        return t
    return ' '.join(_braces(t[:m_.start()] + x_ + t[m_.end():]) for x_ in m_.group(1).split(','))
RES_DIR = os.path.join(ROOT, 'provenance', 'rc7_traces', 'results')


def _mask_shas(o, out):
    """Add to the set out every mask_sha256 and dxf_sha256 string found anywhere in the JSON object o, and return
    out."""
    if isinstance(o, dict):
        for k_, v_ in o.items():
            # Fig. 3a's whole-device record names its mask by dxf_sha256.
            if k_ in ('mask_sha256', 'dxf_sha256') and isinstance(v_, str):
                out.add(v_)
            else:
                _mask_shas(v_, out)
    elif isinstance(o, list):
        for v_ in o:
            _mask_shas(v_, out)
    return out


if os.path.isdir(RES_DIR):
    import hashlib
    bad_ = []
    for n_, fl_ in TRACE_FILES.items():
        st_ = _braces(SP.to_engine(SP.load(os.path.join(PRESETS, n_ + '.yaml')))[1]['source_trace'] or '')
        dxf_ = os.path.join(ROOT, 'runs', n_, n_ + '.dxf')
        sha_ = hashlib.sha256(open(dxf_, 'rb').read()).hexdigest() if os.path.exists(dxf_) else None
        for f_ in fl_:
            p_ = os.path.join(RES_DIR, f_)
            if os.path.basename(f_) not in st_ or not os.path.exists(p_):
                bad_.append('%s: %s %s' % (n_, f_, 'not named' if os.path.basename(f_) not in st_ else 'missing'))
                continue
            with open(p_) as fh_:
                shas_ = _mask_shas(json.load(fh_), set())
            if sha_ is None or shas_ != {sha_}:
                bad_.append('%s: %s traces %s, the shipped DXF is %s' % (n_, f_, sorted(x[:12] for x in shas_), (sha_ or 'unbuilt')[:12]))
    check('ADJ-C6 (codereview2 N1): every source_trace pointer resolves to traces of the shipped mask', not bad_, '; '.join(bad_[:3]))
else:
    print('%-66s %s' % ('ADJ-C6 (codereview2 N1): every source_trace pointer resolves', 'not run (no provenance/rc7_traces/results yet)'))
d6b = copy.deepcopy(base); d6b.setdefault('paper', {})['source_trace'] = 3
o6b = SP.to_engine(SP.load(d6b))[1]
check('ADJ-C6 (R16): a value that is not a pointer is noted, not refused (paper is free-form)', o6b['source_trace'] is None and
      any('paper.source_trace 3 is not a string' in n for n in o6b['notes_spec']))
d6c = copy.deepcopy(base); d6c.setdefault('paper', {})['source_trace'] = ['a.json', 'b.json']
check('ADJ-C6 (R16): a list of pointers is joined', SP.to_engine(SP.load(d6c))[1]['source_trace'] == 'a.json; b.json')
src_pl = open(os.path.join(ROOT, 'opencif', 'pipeline.py')).read()
check('ADJ-C6: assess writes status_codes and source_trace', 'status_codes=codes' in src_pl and "status['source_trace']" in src_pl)
check('ADJ-C6: the CLI prints them', '_print_status(rec)' in open(os.path.join(ROOT, 'opencif', '__main__.py')).read())

# ---------------------------------------------------------------------------------------------------- C7 environment
E = PL.environment_stamp()
check('ADJ-C7: tests and lock hashed, versions recorded', E['tests_sha256'] and E['lock_sha256'] and E['installed'].get('numpy')
      and isinstance(E['matches_lock'], bool) and E['python'], 'tests %s lock %s matches %s' % (E['tests_sha256'], E['lock_sha256'], E['matches_lock']))
check('ADJ-C7: the code stamp is unchanged in definition (no tests in it)', 'tests' not in str(PL.code_stamp.__code__.co_consts))

# ---------------------------------------------------------------------------------------------------- C8 origins
du = SP.defaults_used(base, sp4b)
# A meander-less preset lists no meander value and no inlet stem, because nothing draws them, although its
# reproduction block sets wall_M, W, and stem_half (code review R5(a)).
check('ADJ-C8: no meander or inlet-stem values listed without a head or a filter',
      not any(x.startswith('meander.') or x.startswith('inlet.stem_half') for x in du), str([x[:40] for x in du if x.startswith('meander')]))
d_mh = SP.load(os.path.join(ROOT, 'examples', 'gifford_targets_meander_head_19um.yaml'))
sp_mh, _ = SP.to_engine(d_mh)
du_mh = SP.defaults_used(d_mh, sp_mh)
check('ADJ-C8: a derived value carries its effective value', any(x.startswith('meander.wall_um = max(min_feature, channel/2) (derived) -> ') for x in du_mh)
      or any(' (derived) -> ' in x and x.startswith('meander.') for x in du_mh), str([x[:60] for x in du_mh if 'derived' in x][:3]))
check('ADJ-C8: without the engine spec the rc6 wording', any(x.startswith('meander.') and x.endswith('(derived)') for x in SP.defaults_used(d_mh)))
# Under the clearance guarantee, a value that the spec gives but a reproduction block replaces with another is
# refused, whereas rc7 before the guarantee only named it in defaults_used. The depicted head's own rows are named as
# drawn (code review R5(b)).
du_x = SP.defaults_used(xia, SP.to_engine(xia)[0])
check('ADJ-C8: the depicted head\'s rectilinear rows are named as the drawn value (Xia: 76)',
      any(x.startswith('meander.rectilinear_rows: not the default 2') and 'n_rect 76' in x for x in du_x), str([x[:60] for x in du_x if 'rectilinear' in x]))
def _refused(d_, *needles):
    """Return whether loading d_ raises a SpecError whose message holds every needle, with the message (or
    'not refused')."""
    try:
        SP.load(d_)
        return False, 'not refused'
    except SP.SpecError as e_:
        return all(n_ in str(e_) for n_ in needles), str(e_)[:120]
d_x9 = copy.deepcopy(xia); d_x9['meander'] = dict(d_x9['meander'], rectilinear_rows=9)
check('clearance guarantee: rectilinear rows 9 beside a depicted head (76) is refused, naming both',
      *_refused(d_x9, 'meander.rectilinear_rows = 9', 'depicted_head'))
d_rf = copy.deepcopy(xia); d_rf['process'] = dict(d_rf['process']); d_rf['process']['fillet_um'] = 3.0
d_rf['reproduction'] = dict(d_rf.get('reproduction') or {}); d_rf['reproduction']['r_fillet'] = 2.0
check('clearance guarantee: a given fillet a reproduction block replaces is refused', *_refused(d_rf, 'process.fillet_um = 3.0', 'reproduction.r_fillet'))
d_rs = copy.deepcopy(d_rf); d_rs['reproduction']['r_fillet'] = 3.0
try:
    SP.load(d_rs); ok_rs = True
except SP.SpecError:
    ok_rs = False
check('clearance guarantee: the same value given in both places is accepted', ok_rs)
d_hw = copy.deepcopy(xia); d_hw['reproduction'] = dict(d_hw.get('reproduction') or {}, W=2400.0)
check('clearance guarantee: a head width beside a different reproduction.W is refused', *_refused(d_hw, 'meander.head_width_um', 'reproduction.W'))

# ---------------------------------------------------------------------------------------------------- C9 passage
nt = []
rp = PL.passage_report({'gap opening G': 20.5, 'straight side channel': 4.86, 'centre channel': 100.0, 'outlet lead': 9.0},
                       5.2, nt)
check('ADJ-C9: the classes against D and 2D', rp['below_D'] == ['straight side channel'] and
      set(rp['below_2D']) == {'straight side channel', 'outlet lead'} and rp['classes']['centre channel']['over_D'] == round(100 / 5.2, 3))
check('ADJ-C9: a warning note, the centre channel not warned', len(nt) == 1 and 'straight side channel 4.86' in nt[0] and
      'outlet lead' in nt[0] and 'centre channel' not in nt[0])
d9 = copy.deepcopy(base); d9.setdefault('checks', {})['passage_report'] = False
check('ADJ-C9: the report can be turned off', SP.to_engine(SP.load(d9))[1]['passage_report'] is False and
      SP.to_engine(base)[1]['passage_report'] is True)

# ---------------------------------------------------------------------------------------------------- DINH-7
check('DINH-7: the Dinh presets round the head junctions 11 um', sp.r_junction == 11.0)
dj = copy.deepcopy(dd); dj['meander'] = dict(dj['meander']); dj['meander']['junction_radius_um'] = 30.0
try:
    SP.to_engine(SP.load(dj)); check('DINH-7: a junction radius above half the wall refused', False)
except SP.SpecError:
    check('DINH-7: a junction radius above half the wall refused', True)
dj2 = copy.deepcopy(SP.load(os.path.join(PRESETS, 'gifford2014_fig4_design.yaml'))); dj2.setdefault('meander', {})['junction_radius_um'] = 5.0
try:
    SP.to_engine(SP.load(dj2)); check('DINH-7: a junction radius refused without a head', False)
except SP.SpecError:
    check('DINH-7: a junction radius refused without a head', True)
from opencif.geometry import layout as L  # noqa: E402
_units = [dict(X=1000.0, w=28.0), dict(X=1060.0, w=28.0)]
_jr = L.junction_roundings(sp, _units)
_a1 = (1.0 - math.pi / 4.0) * 11.0 ** 2
check('DINH-7: two roundings per unit, each (1 - pi/4) r^2, at the stem edges on the trunk wall',
      len(_jr) == 4 and all(abs(q.area - _a1) < 0.05 for q in _jr) and
      sorted(round(q.bounds[0], 6) for q in _jr) == [1000.0 - 14.0 - 11.0, 1000.0 + 14.0, 1060.0 - 14.0 - 11.0, 1060.0 + 14.0]
      and all(abs(q.bounds[1] - sp.wc0 / 2.0) < 1e-9 and abs(q.bounds[3] - (sp.wc0 / 2.0 + 11.0)) < 1e-9 for q in _jr),
      '%.3f um2 each' % (_jr[0].area if _jr else 0.0))
_sp0j = copy.copy(sp); _sp0j.r_junction = 0.0
check('DINH-7: none without meander.junction_radius_um', L.junction_roundings(_sp0j, _units) == [])

# ---------------------------------------------------------------------------------------------------- GIF4-1 B
# Reading B (legs in path order) is the variant gifford2014_fig4_pathorder, and the reproduction
# gifford2014_fig4_aligned is the nested reading (D1, D6).
dB = SP.load(os.path.join(ROOT, 'provenance', 'rc7_traces', 'variants', 'gifford2014_fig4_pathorder.yaml'))
spB, oB = SP.to_engine(dB)
check('GIF4-1 B: the path-order variant (reading B) takes legs in path order', spB.leg_order == 'path' and oB['layout'] == 'paperclip')
spA8, oA8 = SP.to_engine(SP.load(os.path.join(PRESETS, 'gifford2014_fig4_aligned.yaml')))
check('GIF4-1 (RC8): the aligned Fig. 4 reproduction is nested, leg 4 ending 2 mm short, the 8 um start, the inset gap',
      spA8.leg_order == 'nested' and spA8.pc_last_short == 2000.0 and spA8.G == 17.2 and oA8['layout'] == 'paperclip'
      and oA8.get('start') == 'floor' and oA8.get('start_width') == 8.0,
      '%s %s %s %s' % (spA8.leg_order, spA8.pc_last_short, spA8.G, oA8.get('start_width')))
S.clear_caches(); scB = S.plan(spB)
scLB, legsB = L.plan_layout(scB, 'paperclip', 4, 818.0)
check('GIF4-1 B: path order puts the legs at 0, d, 2d, 3d', [round(l['Y'], 6) for l in legsB] == [0.0, 818.0, 1636.0, 2454.0])
spN = copy.copy(spB); spN.leg_order = 'nested'
scN = copy.copy(scB); scN.spec = spN
_, legsN = L.plan_layout(scN, 'paperclip', 4, 818.0)
check('GIF4-1 B: nested order unchanged (0, +d, -d, +2d)', [round(l['Y'], 6) for l in legsN] == [0.0, 818.0, -818.0, 1636.0])
dBad = copy.deepcopy(dB); dBad['layout'] = dict(style='linear', order='path')
try:
    SP.to_engine(SP.load(dBad)); check('GIF4-1 B: layout.order refused outside the paperclip', False)
except SP.SpecError:
    check('GIF4-1 B: layout.order refused outside the paperclip', True)
check('GIF4-1 B: layout.order listed among the defaults only for a paperclip',
      any('layout.order' in x for x in SP.defaults_used(SP.load(os.path.join(PRESETS, 'gifford2014_fig4_design.yaml')))) and
      not any('layout.order' in x for x in SP.defaults_used(SP.load(os.path.join(PRESETS, 'dinh2024_v1_depicted.yaml')))))

# ---------------------------------------------------------------------------------------------------- code review (rc7)
# Square corners below the head fillet are drawn, measured, and modeled at the fillet (R1).
d1 = copy.deepcopy(dd); d1['meander'] = dict(d1['meander'])
d1['meander'].pop('corner_radius_um'); d1['meander'].pop('inner_corner_radius_um')
sp1, o1 = SP.to_engine(SP.load(d1))
check('R1: default radii drawn and modelled at the head fillet', sp1.r_corner_out == sp1.r_fillet == 5.0 and sp1.r_corner_in == 5.0)
check('R1: the derived radii are listed with their effective value',
      any(x.startswith('meander.corner_radius_um = the head fillet') and x.endswith('-> 5') for x in o1['defaults_used']) and
      any(x.startswith('meander.inner_corner_radius_um = the head fillet') and x.endswith('-> 5') for x in o1['defaults_used']))
d1b = copy.deepcopy(dd); d1b['meander'] = dict(d1b['meander']); d1b['meander']['corner_radius_um'] = 2.0
sp1b, o1b = SP.to_engine(SP.load(d1b))
check('R1: a given radius below the fillet is raised to it, with a note',
      sp1b.r_corner_out == 5.0 and any('corner_radius_um 2.00 um is below the head fillet' in n for n in o1b['notes_spec']))
sp1c = copy.copy(sp1); sp1c.r_corner_in = 0.0
check('R1: c_corner is taken at the drawn (fillet) radius, not 0', S.c_corner(sp1, 28.0) != S.c_corner(sp1c, 28.0) and
      abs(S.c_corner(sp1, 28.0) - (-0.3916)) < 2e-3, '%.4f' % S.c_corner(sp1, 28.0))
pR = Gm.Path(0.0, 0.0, 90.0, square=True); pR.line(100.0); pR.turn(27.5, 90.0); pR.line(100.0)
sp_lo = copy.copy(sp1); sp_lo.r_corner_out, sp_lo.r_corner_in = 2.0, 2.0
a_eff = Gm.smooth(Gm.ribbon(pR, 28.0, sp1), 5.0, 1500.0).area
a_lo = Gm.smooth(Gm.ribbon(pR, 28.0, sp_lo), 5.0, 1500.0).area
check('R1: the smoothing draws a 2 um corner as the fillet-radius corner', abs(a_eff - a_lo) < 0.3, '%.3f um2' % (a_eff - a_lo))
# There is no spike at a square corner in any orientation, and the smoothing keeps the rounding (R13).
worst, worst_s = 0.0, 0.0
end_loss = 4.0 * (1.0 - math.pi / 4.0) * 25.0                 # the opening rounds the ribbon's four flat-end corners
for hd in range(0, 360, 7):
    for tn in (90.0, -90.0):
        pq = Gm.Path(10.0, -20.0, float(hd), square=True); pq.line(100.0); pq.turn(27.5, tn); pq.line(100.0)
        gq = Gm.ribbon(pq, 28.0, sp); gs = Gm.ribbon(pq, 28.0, None)
        worst = max(worst, abs((gq.area - gs.area) - Gm.corner_round_area(sp)))
        worst_s = max(worst_s, abs(Gm.smooth(gq, 5.0, 1500.0).area - (gq.area - end_loss)))
check('R13: the drawn roundings equal corner_round_area at every heading (both turns)', worst < 0.02, 'worst %.4f um2' % worst)
check('R13: smoothing refills no rounding (no spike) at any heading', worst_s < 0.3, 'worst %.3f um2' % worst_s)
# These checks test the domain of the radii on the drawn unit widths and the record's parts under bend_model none
# (R14).
for k_, v_, what in (('corner_radius_um', 9.0, 'an outer rounding above 0.307 w'),
                     ('inner_corner_radius_um', 14.0, 'an inner rounding above half the wall')):
    d14 = copy.deepcopy(dd); d14['meander'] = dict(d14['meander']); d14['meander'][k_] = v_
    try:
        SP.to_engine(SP.load(d14)); check('R14: refuses %s' % what, False)
    except SP.SpecError:
        check('R14: refuses %s' % what, True)
sp14 = copy.copy(sp); sp14.bend_model = 'none'
c14, dp14, sh14, lam14, held14 = S.c_corner(sp14, 28.0, parts=True)
check('R14: under bend_model none no term, the table values reported', c14 == 0.0 and dp14 < -0.3 and sh14 < -0.3,
      'deep %.4f shallow %.4f' % (dp14, sh14))
# The depicted Dinh head is derived for square turns, and the figure-aligned head is the depicted plan's (R2).
with open(os.path.join(ROOT, 'opencif', 'data', 'depicted.json')) as fh_:
    dj_ = json.load(fh_)
check('R2: data/depicted.json[dinh] is derived for square turns (7 / 10 um)',
      dj_['dinh'].get('turns', {}).get('turns') == 'square' and dj_['dinh']['turns']['r_corner_out_um'] == 7.0 and
      dj_['dinh']['turns']['r_corner_in_um'] == 10.0)
from opencif.geometry import depicted as DP  # noqa: E402
from opencif import published as DV  # noqa: E402
S.clear_caches()
sp2 = DV.dinh_fig2('V1'); sp2 = sp2.__class__(**{**sp2.__dict__, 'turns': 'square', 'r_corner_out': 7.0, 'r_corner_in': 10.0})
f2_ = DP.report('x', sp2, DP.dinh_depicted_geometry(sp2))['f_head']
check('R2: the stored fractions are the square-turn derivation', max(abs(a / b - 1.0) for a, b in zip(f2_, dj_['dinh']['f_head'])) < 1e-9)
r2c = subprocess.run([sys.executable, '-I', os.path.join(ROOT, 'provenance', 'rc7_scripts', 'make_dinh_figure_aligned_head.py'),
                      '--check'], capture_output=True, text=True)
check('R2: the figure-aligned preset carries the depicted preset\'s planned head', r2c.returncode == 0 and 'unchanged' in r2c.stdout,
      r2c.stdout.strip().split('\n')[-1][:80])
# The passage report is made against the largest realized cutoff, and repeated calls do not repeat the notes (R3 and
# R4).
nt3 = []
w3 = {'gap opening G': 17.2, 'straight side channel': 8.0, 'centre channel': 100.0}
rp3 = PL.passage_report(w3, 5.2, nt3, D_max=10.23, G=17.2)
check('R3: classes against min(D_realised_max, G), with a warning', rp3['D_pass_max_um'] == 10.23 and
      rp3['below_D_pass_max'] == ['straight side channel'] and any(n.startswith(PL.PASSAGE_NOTE_MAX) for n in nt3))
PL.passage_report(w3, 5.2, nt3, D_max=10.23, G=17.2)
check('R4: a second report does not repeat its notes', len(nt3) == 2 and
      sum(n.startswith(PL.PASSAGE_NOTE) for n in nt3) == 1 and sum(n.startswith(PL.PASSAGE_NOTE_MAX) for n in nt3) == 1)
rp3b = PL.passage_report(w3, 5.2, [], D_max=4.0, G=17.2)
check('R3: no comparison when no row is realised above D', rp3b['D_pass_max_um'] is None and rp3b['below_D_pass_max'] == [])
rec_dir = os.path.join(ROOT, 'runs', 'gifford2014_fig4_aligned')
if os.path.exists(os.path.join(rec_dir, 'gifford2014_fig4_aligned.pkl')):
    import shutil
    import tempfile
    tmp_ = tempfile.mkdtemp()
    for ext in ('.json', '.pkl'):
        shutil.copy(os.path.join(rec_dir, 'gifford2014_fig4_aligned' + ext), tmp_)
    n_ = []
    for _ in range(2):
        r_ = PL.reassess(tmp_, 'gifford2014_fig4_aligned', render=False)
        n_.append(len(r_['assessment']['notes']))
    shutil.rmtree(tmp_)
    check('R4: reassess twice gives the same notes', n_[0] == n_[1], str(n_))
else:
    print('%-66s %s' % ('R4: reassess twice gives the same notes', 'not run (no rebuilt record in runs/)'))
# The plateau minimum is taken over the rows that carry flow (R6).
xa = SP.load(os.path.join(PRESETS, 'xia2016_rs7_aligned.yaml'))
spx, _ = SP.to_engine(xa); S.clear_caches(); scx = S.plan(spx)
gx = BU.gates(spx, scx, S.KAPPA_G * spx.delta_star)
check('R6: Xia aligned: no zero from rows that carry no flow (its rectilinear rows at f/f* ~ 2e-7 left out)',
      gx['plateau_u_min_all_planned_rows'] is not None and gx['plateau_u_min_all_planned_rows'] > 0.005,
      '%.4f' % (gx['plateau_u_min_all_planned_rows'] or 0.0))
fa = SP.load(os.path.join(PRESETS, 'dinh2024_v1_figure_aligned.yaml'))
spf, _ = SP.to_engine(fa); S.clear_caches(); scf = S.plan(spf)
gf = BU.gates(spf, scf, S.KAPPA_G * spf.delta_star)
check('R6: Dinh figure-aligned: the head rows at 2.17 um below the plateau, with their share',
      gf['plateau_ok'] and 0.02 < gf['plateau_u_min_all_planned_rows'] < S.U_LO and gf['plateau_F_share_below_all_planned_rows'] > 0.0,
      'u %.4f share %.3f' % (gf['plateau_u_min_all_planned_rows'], gf['plateau_F_share_below_all_planned_rows']))
# The capacity counts a square turn's corner, and arc heads are unchanged (R7).
check('R7: square-turn capacity = geometric + c_corner w', abs(sp.L_cap - sp.cap_frac * (sp.L_taper(sp.R_max) + S.c_corner(sp, sp.GS) * sp.GS)) < 1e-9
      and sp.corner_term(28.0) < 0.0)
spA, _ = SP.to_engine(SP.load(os.path.join(PRESETS, 'dinh2024_v1_specified.yaml')))
check('R7: arc-turn capacity unchanged (geometric)', spA.L_cap == spA.cap_frac * spA.L_taper(spA.R_max) and spA.corner_term(28.0) == 0.0)
# Arc-turn heads plan and draw as in rc6, by the fingerprints of the planner and the head geometry computed on rc6
# (4bbdaabef1225a24) (R12).
RC6_FP = {"dinh2024_v1_specified.yaml": "2729c8144e919b4a", "gifford_targets_meander_head_19um.yaml": "edecd55bfc2cabeb",
          "m01_linear_const_pack_plain_separate.yaml": "a529f9f42009a5a3", "m09_paperclip_const_narrow_plain_teardrop.yaml": "1d91857d4345b90f",
          "xia2016_rs7_depicted.yaml": "7f6fb1b22e6161a8"}
fp_ = subprocess.run([sys.executable, '-I', os.path.join(ROOT, 'provenance', 'rc7_scripts', 'fingerprint_head.py'), ROOT,
                      'presets/dinh2024_v1_specified.yaml', 'presets/xia2016_rs7_depicted.yaml',
                      'examples/matrix/m01_linear_const_pack_plain_separate.yaml',
                      'examples/matrix/m09_paperclip_const_narrow_plain_teardrop.yaml',
                      'examples/gifford_targets_meander_head_19um.yaml'], capture_output=True, text=True, cwd=ROOT)
got_ = json.loads(fp_.stdout.strip().split('\n')[-1]) if fp_.returncode == 0 else {}
check('R12: arc-turn heads (arms, hairpins, pack, narrow) plan and draw as in rc6', got_ == RC6_FP,
      ', '.join(k for k in RC6_FP if got_.get(k) != RC6_FP[k]) or '5 specs')

# ------------------------------------------------------------------------------------------- codereview3 (T1-T10)
sp_f = gf['shear_peak']
check('T1 (B1): an exact constant-shear taper holds its peak shear (head rows outside the taper)',
      abs(sp_f['over_the_taper']['max_over_min'] - 1.0) < 1e-9 and sp_f['taper_shear_formula'] == 'exact',
      'max/min - 1 = %.2e over %d rows' % (sp_f['over_the_taper']['max_over_min'] - 1.0, sp_f['over_the_taper']['rows']))
kind_x = np.asarray(scx.kind)
floor_x = max(spx.wc_min, 0.3 * spx.wc0)
rows_x = int(np.sum((np.asarray(scx.wc, float) > floor_x + 1e-6) & (kind_x == 'L')))
c0_ = S.C.C_shear(max(spx.wc0, spx.h) / min(spx.wc0, spx.h)) / 6.0 - 1.0
cf_ = S.C.C_shear(max(floor_x, spx.h) / min(floor_x, spx.h)) / 6.0 - 1.0
tp_x = gx['shear_peak']['true_over_parallel_plates_minus_1']
check('T2 (B1): parallel-plate taper: the taper is the L rows above the floor; true/pp - 1 spans C(alpha)/6 - 1 from w_c0 to the floor',
      gx['shear_peak']['over_the_taper']['rows'] == rows_x and gx['shear_peak']['over_the_taper']['max_over_min'] > 1.1
      and abs(tp_x[1] - max(c0_, cf_)) < 1e-9 and abs(tp_x[0] - min(c0_, cf_)) < 0.02,
      'rows %d (L rows above %.1f um: %d); range %s; C/6 - 1 at w_c0 %.4f, at the floor %.4f' % (
          gx['shear_peak']['over_the_taper']['rows'], floor_x, rows_x, [round(v, 4) for v in tp_x], c0_, cf_))
S.clear_caches(); sc_a = S.plan(sp_a)
ga = BU.gates(sp_a, sc_a, S.KAPPA_G * sp_a.delta_star)
a_, b_ = max(sp_a.wc0, sp_a.h) * 1e-6, min(sp_a.wc0, sp_a.h) * 1e-6
want_ = S.C.C_shear(a_ / b_) * sp_a.Q_in_uL_min * 1e-9 / 60.0 / (a_ * b_ ** 2)
check('T3: shear_peak row 1 = C(alpha) Q / (a b^2) (constant centre channel, s^-1)',
      abs(ga['shear_peak']['row1_per_s'] / want_ - 1.0) < 1e-12 and 'over_the_taper' not in ga['shear_peak'],
      '%.3f s^-1 (by hand %.3f)' % (ga['shear_peak']['row1_per_s'], want_))
nan_ = float('nan')
nw4 = dict(fA=[0.0, -0.01, 0.6, 0.1], fB=[0.1, 0.1, 0.1, 0.0], dA=[nan_, nan_, nan_, 5.0], dB=[5.0, 5.0, 5.0, nan_])
rc4 = BU.realised_cutoff(None, nw4, 5.0)
check('T4 (B2): rows_f_half_or_more counts only rows with f >= 1/2 on a side (not f = 0, not reverse)',
      rc4['rows_f_half_or_more'] == 1, str(rc4['rows_f_half_or_more']))
nt5 = []
rp5 = PL.passage_report({'gap opening G': 17.2, 'straight side channel': 8.0, 'centre channel': 100.0}, 5.2, nt5,
                        D_max=None, G=17.2, rows_no_cutoff=2)
check('T5 (m1, m2): f >= 1/2 rows: d_pass = G, the note names them; D_max absent does not fail',
      rp5['D_pass_max_um'] == 17.2 and rp5['rows_f_half_or_more'] == 2 and
      any('2 rows extract f >= 1/2' in n for n in nt5 if n.startswith(PL.PASSAGE_NOTE_MAX)), ' | '.join(nt5)[:120])
check('T6: the source_trace status codes', PL.STATUS_CODES['source_trace'] == ('NAMED', 'NOT_NAMED'))
check('T7 (M1): brace groups in a pointer expand', set(_braces('a{b,c}d').split()) == {'abd', 'acd'} and _braces('xy') == 'xy')
for src_, d8 in (('reproduction.r_fillet', dict(reproduction=dict(r_fillet=9.0))), ('process.fillet_um', dict(process=dict(fillet_um=9.0)))):
    d14 = copy.deepcopy(dd)
    for sec_, kv_ in d8.items():
        d14[sec_] = dict(d14.get(sec_) or {}, **kv_)
    try:
        SP.to_engine(SP.load(d14)); check('T8: R14 names %s when the fillet comes from it' % src_, False, 'not refused')
    except SP.SpecError as e_:
        check('T8: R14 names %s when the fillet comes from it' % src_, src_ in str(e_), str(e_)[:100])
check('T9: under bend_model none the record reports the depth weight lambda(h, w)', abs(lam14 - S.depth_lambda(sp14.h, 28.0)) < 1e-15,
      '%.6f' % lam14)
ox_ = SP.to_engine(xa)[1]
check('T10: paper.targets_fitted reaches the reproduction status (the aligned Xia preset fits D and F)',
      'F 0.94072' in (ox_.get('targets_fitted') or '') and SP.to_engine(SP.load(os.path.join(PRESETS, 'gifford2014_fig4_aligned.yaml')))[1].get('targets_fitted') is None
      and 'paper.targets_fitted: %s' in open(os.path.join(ROOT, 'opencif', 'pipeline.py')).read())

print('ALL OK' if not fails else 'FAILED: %s' % fails)
sys.exit(1 if fails else 0)
