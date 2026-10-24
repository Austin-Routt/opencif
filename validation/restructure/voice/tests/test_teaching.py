"""test_teaching.py: the teaching core (opencif/teaching.py) against the full generator and against its own photomask.

The suite checks two things.
1. Its design of a straight ladder is the generator's plan of the same device, number for number (the march, the
   closure, and the side widths), on the example examples/gifford_physics_linear.yaml.
2. Its photomask, read back with no knowledge of the design, realizes the design, that is, the cutoff at every gap and
   F.
Run it with `python3 tests/test_teaching.py`."""
import os
import sys
import tempfile

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from opencif import pipeline as PL, spec as SP            # noqa: E402
print('opencif %s (code %s)' % (PL.CODE['version'], PL.CODE['sha256']))
from opencif import teaching as T                         # noqa: E402
from opencif.physics import schedule as S                 # noqa: E402

fails = []


def check(name, cond, info=''):
    print('%-80s %s %s' % (name, 'ok ' if cond else 'FAIL', info))
    if not cond:
        fails.append(name)


# 1. The generator's plan of examples/gifford_physics_linear.yaml (no head, constant center, sealed start)
d = SP.load(os.path.join(ROOT, 'examples', 'gifford_physics_linear.yaml'))
spec_u, opts = SP.to_engine(d)
S.clear_caches()
sc = S.plan(spec_u)
p = d['process']
lad = T.design(D_um=d['targets']['cutoff_um'], F=d['targets']['filtrate_fraction'], depth_um=p['depth_um'],
               gap_um=p['gap_um'], post_um=p['post_length_um'], wall_um=p['post_width_um'],
               centre_um=d['centre_channel']['width_in_um'], outlet_um=spec_u.L_out)
check('the same lamina delta* and fraction f as the generator (bit for bit)', lad.delta == spec_u.delta_star and
      np.all(sc.f == lad.f), 'delta* %.6f um, f %.6e' % (lad.delta, lad.f))
check('the same number of gaps N = ceil(ln(1-F)/ln(1-2f)) (App. A eq. Nconst)',
      lad.N == sc.N == int(np.ceil(np.log(1 - lad.F_target) / np.log(1 - 2 * lad.f))), 'N %d' % lad.N)
check('the same march P(i), S(i) and F = 1 - P(N) (bit for bit)', np.array_equal(lad.P, sc.P) and
      np.array_equal(lad.S, sc.S) and lad.F == sc.F, 'F %.6f' % lad.F)
check('the same centre and gap resistances (bit for bit)', np.all(sc.Rc[:-1] == lad.Rc) and np.all(sc.Rg == lad.Rg))
check('the same exact closure R_s(i), every segment (bit for bit)', np.array_equal(lad.Rs, sc.Rs))
check('the same side widths w_s(i), every segment (bit for bit)', np.array_equal(lad.ws, sc.ws[:sc.N - 1]),
      'w_s %.3f -> %.3f um' % (lad.ws[0], lad.ws[-1]))
m13 = lad.Rg * lad.f / lad.Rc                             # gap dominance, R_g f / R_c (spec v5 M13)
pred = -2.0 * lad.f * m13 / (1.0 - 2.0 * lad.f)
dev = float(np.max(np.abs((lad.Rs / lad.Rs_pe - 1.0) / pred - 1.0)))
check('exact closure = PE closure x (1 - 2 f M13/(1 - 2f)), every segment (the gap drop\'s change)', dev < 1e-6,
      'R_s/R_s,PE - 1 = %.3e; %.1e off the identity' % (pred, dev))
check('the mass identity P(N) + 2 S(N) = 1 (App. A)', abs(lad.P[-1] + 2 * lad.S[-1] - 1) < 1e-12)
# 2. A smaller ladder drawn, written, and audited blind
small = T.design(D_um=6.0, F=0.3, depth_um=80, gap_um=20.5, post_um=25.3, wall_um=59.5, centre_um=100)
fluid = T.draw(small)
out = tempfile.mkdtemp(prefix='opencif_teaching_')
path = os.path.join(out, 'ladder.dxf')
T.write_photomask(fluid, path)
a = T.audit_photomask(path, depth_um=80, D_um=6.0)
check('the drawn ladder has N gaps and N - 1 posts per side',
      a['gaps_per_side'] == small.N and a['posts_per_side'] == small.N - 1,
      'N %d' % small.N)
check('read back: the wall, the gaps and the centre as designed', abs(a['wall_um'] - 59.5) < 1e-9 and
      max(abs(g - 20.5) for g in a['gap_um']) < 1e-9 and max(abs(w - 100) for w in a['centre_um']) < 1e-9)
check('read back: the side widths as designed',
      abs(a['side_um'][0] - small.ws[0]) < 1e-9 and abs(a['side_um'][1] - small.ws[-1]) < 1e-9)
check('the drawn network realises D at every gap, both sides (< 1e-6)', a['D_rel_dev_max'] < 1e-6,
      '%.1e' % a['D_rel_dev_max'])
check('the drawn network realises F (< 1e-9)', abs(a['F'] - small.F) < 1e-9, 'F %.6f' % a['F'])
print('ALL OK' if not fails else 'FAILED: %s' % fails)
sys.exit(1 if fails else 0)
