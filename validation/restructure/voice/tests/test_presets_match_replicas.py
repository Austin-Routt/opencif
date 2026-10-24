"""test_presets_match_replicas.py: the paper presets against the replica generator's designs.

The paper presets (the decoupled D and F, plus the papers' other parameters) must produce exactly the replica
generator's designs. They must give the same engine Spec (apart from opencif's trim_F and the defaults it now derives),
the same N, the same f at every row except the trimmed last gap (to 1e-15, the round trip D -> delta* = D/kappa_g),
and the same F (identical, because there is no trim). Run it with `python3 tests/test_presets_match_replicas.py`."""
import json, os, sys
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from opencif import spec as SP
from opencif import pipeline as _PL_CODE                # noqa: E402  (every suite's log names the code it ran on)
print('opencif %s (code %s)' % (_PL_CODE.CODE['version'], _PL_CODE.CODE['sha256']))
from opencif.physics import schedule as S
from opencif import published as D

HERE = os.path.join(ROOT, 'opencif', 'data')
dep = json.load(open(os.path.join(HERE, 'depicted.json')))
CASES = {
    'gifford2014_fig4_design': lambda: D.gifford_fig4(),
    'dinh2024_v1_specified': lambda: D.dinh_fig2('V1'),
    # The depicted head is drawn with square quarter-turns (outer 7, inner 10 um), and its fractions are derived for
    # them (provenance/rc7_scripts/regen_depicted_dinh_square.py), so the replica takes the same turns (DINH-1,
    # codereview2 N2).
    'dinh2024_v1_depicted': lambda: (lambda s_: s_.__class__(**{**s_.__dict__, 'turns': 'square', 'r_corner_out': 7.0,
                                                                'r_corner_in': 10.0}))(
        D.dinh_fig2('V1', f_head=tuple(dep['dinh']['f_head']))),
    'xia2016_rs7_specified': lambda: D.xia_fig(),
    'xia2016_rs7_depicted': lambda: D.xia_fig(f_head=tuple(dep['xia']['f_head']), NM_force=dep['xia']['NM'], n_rect=dep['xia']['n_rect']),
}
# IGNORE lists operating and bookkeeping values rather than device geometry, and the replica factories predate them.
# They are the flow rate for the Reynolds gate and the tubing I.D. kept for the tubing scenarios, whose resistance
# R_ext_c is still compared. The flared outlet lead (flare_deg, flare_start) is a drawing option of the outlet stretch,
# set in the two Dinh presets after the replica factories were written. It does not enter the march compared here (N,
# every f, and F), and its geometry is checked against the paper's figure by tests/test_figure_anchors.py. The
# outlet-refusal switch of the stretch-end table domain (exp_domain_refuse) is not compared here either, because it is
# the gates mode rather than device geometry. The switch is set by spec.to_engine from gates.mode, whereas the factories
# set no mode. Finally, r_junction (meander.junction_radius_um, DINH-7) is a drawing option of the entrance corners of
# the head gaps, and it does not enter the march compared here.
IGNORE = {'trim_F', 'name', 'note', 'Q_in_uL_min', 'tubing_centre_id_mm', 'flare_deg', 'flare_start', 'exp_domain_refuse',
          'r_junction'}
ok_all = True
for name, fac in CASES.items():
    d = SP.load(os.path.join(ROOT, 'presets', name + '.yaml'))
    sp, opts = SP.to_engine(d)
    ref = fac()
    diffs = []
    for k in set(ref.__dict__) | set(sp.__dict__):
        if k in IGNORE:
            continue
        a, b = ref.__dict__.get(k), sp.__dict__.get(k)
        if isinstance(a, float) or isinstance(b, float):
            if a is None or b is None or abs(float(a) - float(b)) > 1e-9 * max(1, abs(float(a))):
                diffs.append((k, a, b))
        elif (tuple(a) if isinstance(a, (list, tuple)) else a) != (tuple(b) if isinstance(b, (list, tuple)) else b):
            diffs.append((k, a, b))
    s1, s2 = S.plan(ref), S.plan(sp)
    same_N = s1.N == s2.N and s1.NM == s2.NM
    df = float(np.max(np.abs(s1.f[:-1] - s2.f[:-1]))) if same_N else float('nan')
    ok = not diffs and same_N and df <= 1e-15 and abs(s2.F - s1.F) < 1e-12
    ok_all &= ok
    print('%-24s %s  spec diffs %s; N %d/%d NM %d/%d; max |df| rows 1..N-1 %.1e; F %.8f (ref %.8f, target %.5f)' % (
        name, 'OK ' if ok else 'BAD', diffs[:6], s1.N, s2.N, s1.NM, s2.NM, df, s2.F, s1.F, sp.F_target))
print('ALL OK' if ok_all else 'MISMATCH')
sys.exit(0 if ok_all else 1)
