"""test_kernels.py: the kernels that opencif uses (opencif/physics/kernels.py) against the design code they restate.

The design code is reference/CIF_v2.py v2.1.0, kept byte-frozen. The kernels must return the same floating-point
numbers as the design code, bit for bit, over grids of inputs, and the frozen file must be the one that Appendix A
names. Run it with `python3 tests/test_kernels.py`."""
import hashlib
import os
import sys
from dataclasses import fields

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, 'reference'))
from opencif import pipeline as _PL                      # noqa: E402  (every suite's log names the code it ran on)
print('opencif %s (code %s)' % (_PL.CODE['version'], _PL.CODE['sha256']))
from opencif.physics import kernels as K                 # noqa: E402
import CIF_v2 as C                                       # noqa: E402  the frozen design code

fails = []


def check(name, cond, info=''):
    print('%-74s %s %s' % (name, 'ok ' if cond else 'FAIL', info))
    if not cond:
        fails.append(name)


def same(a, b):
    """Return True when two floats, or two sequences of floats, are equal bit for bit, so that nan equals nan."""
    a, b = np.atleast_1d(np.asarray(a, float)), np.atleast_1d(np.asarray(b, float))
    return a.shape == b.shape and a.tobytes() == b.tobytes()


# 0. The design code is the frozen v2.1.0 file
md5 = hashlib.md5(open(os.path.join(ROOT, 'reference', 'CIF_v2.py'), 'rb').read()).hexdigest()
check('reference/CIF_v2.py is the frozen v2.1.0 file (md5 0e0961a5...)', md5 == '0e0961a50441afb94671b4470d41bace', md5)
check('the series length: 50 odd terms in both', K.N_ODD_TERMS == C._N_ODD == 50)

# 1. The duct law (App. A Result 14): B, r, and r_inverse
alphas = np.concatenate([[1.0], np.geomspace(1.0 + 1e-9, 400.0, 3000)])
check('duct_bracket = CIF_v2.B, bit for bit, %d aspect ratios in [1, 400]' % len(alphas),
      same([K.duct_bracket(a) for a in alphas], [C.B(a) for a in alphas]))
UM = 1e-6
ws = np.geomspace(0.5, 3000.0, 400) * UM
hs = np.array([20.0, 50.0, 80.0, 100.0, 140.0, 400.0]) * UM
rk = [K.r(w, h, 1e-3) for h in hs for w in ws]
check('r = CIF_v2.r, bit for bit, %d widths x %d depths (both orders of the sides)' % (len(ws), len(hs)),
      same(rk, [C.r(w, h, 1e-3) for h in hs for w in ws]))
ys = [K.r(w, h, 1e-3) * fac for h in hs[::2] for w in ws[::8] for fac in (1.0, 1.0137)]
hh = [h for h in hs[::2] for w in ws[::8] for fac in (1.0, 1.0137)]
check('r_inverse = CIF_v2.r_inverse, bit for bit, %d resistances' % len(ys),
      same([K.r_inverse(y, h, 1e-3) for y, h in zip(ys, hh)], [C.r_inverse(y, h, 1e-3) for y, h in zip(ys, hh)]))
# 2. The skim cubic and its inverse (App. A Results 1 and 2)
xs = np.linspace(0.0, 0.5, 20001)
check('g_cubic = CIF_v2.g_cubic, bit for bit, 20001 points of [0, 1/2]',
      same([K.g_cubic(x) for x in xs], [C.g_cubic(x) for x in xs]))
fs = np.concatenate([np.geomspace(1e-9, 0.49, 4000), [1.04e-4, 5.76e-4, 0.25, 0.4999999]])
check('x_of_f = CIF_v2.x_of_f(f, 1.0, "cubic"), bit for bit, %d fractions' % len(fs),
      same([K.x_of_f(f) for f in fs], [C.x_of_f(f, 1.0, 'cubic') for f in fs]))
bad = []
for f in (0.0, 0.5, -1e-6, 0.7):
    try:
        K.x_of_f(f)
        bad.append(f)
    except ValueError:
        pass
check('x_of_f refuses f outside (0, 1/2) as CIF_v2 does', not bad, str(bad))
# 3. The cutoff pairing (App. A Result 9; App. B C3b)
Pk, Pc = K.Pairing(), C.Pairing()
common = [f.name for f in fields(Pk) if f.name in {g.name for g in fields(Pc)}]
check('Pairing: every field equal to CIF_v2\'s (%s)' % ', '.join(common),
      all(getattr(Pk, n) == getattr(Pc, n) for n in common) and Pk.kappa == 4.9138 and Pk.kernel == 'cubic')
check('Pairing.D(delta) = kappa delta = CIF_v2',
      same([Pk.D(d) for d in (0.1, 0.44, 1.0)], [Pc.D(d) for d in (0.1, 0.44, 1.0)])
      if hasattr(Pc, 'D') else same([Pk.D(d) for d in (0.1, 0.44)], [Pc.kappa * d for d in (0.1, 0.44)]))
# 4. Spec v5: the peak shear, the entrance length, and the post outline
al = np.concatenate([[1.0], np.geomspace(1.01, 50.0, 200)])
check('C_shear = CIF_v2.C_shear, bit for bit, %d aspect ratios' % len(al),
      same([K.C_shear(a) for a in al], [C.C_shear(a) for a in al]))
res = np.geomspace(1e-4, 100.0, 300)
check('entrance_length = CIF_v2.entrance_length, bit for bit', same([K.entrance_length(r_, 50e-6) for r_ in res],
                                                                    [C.entrance_length(r_, 50e-6) for r_ in res]))
check('stadium_post = CIF_v2.stadium_post, bit for bit (3 posts)',
      all(same(K.stadium_post(1.0, 2.0, Lp, Wp), C.stadium_post(1.0, 2.0, Lp, Wp))
          for Lp, Wp in ((27.0, 50.0), (26.1, 50.0), (25.6, 50.6))))
try:
    K.stadium_post(0.0, 0.0, 50.0, 27.0)
    check('stadium_post refuses crossing caps', False)
except K.GateError:
    check('stadium_post refuses crossing caps', True)
# 5. The tolerances and the process defaults
Tk, Tc = K.Tolerances(), C.Tolerances()
check('Tolerances: every field equal to CIF_v2\'s', [(f.name, getattr(Tk, f.name)) for f in fields(Tk)] ==
      [(f.name, getattr(Tc, f.name)) for f in fields(Tc)])
Pk_, Pc_ = K.Process(), C.Process()
pf = [f.name for f in fields(Pk_)]
check('Process: the fields opencif keeps equal CIF_v2\'s (%s), w_min and R_gap too' % ', '.join(pf),
      all(getattr(Pk_, n) == getattr(Pc_, n) for n in pf) and Pk_.w_min == Pc_.w_min and same(Pk_.R_gap, Pc_.R_gap))
# 6. The schedule uses these kernels (and nothing else of the design code)
from opencif.physics import schedule as S                # noqa: E402
check('the schedule\'s kernels are opencif.physics.kernels', S.C is K and S.KAPPA_G == K.Pairing().kappa)
check('no module of the package imports the design code', not any(
    'import CIF_v2' in open(os.path.join(dp, f)).read() for dp, _, fs_ in os.walk(os.path.join(ROOT, 'opencif'))
    for f in fs_ if f.endswith('.py')))
print('ALL OK' if not fails else 'FAILED: %s' % fails)
sys.exit(1 if fails else 0)
