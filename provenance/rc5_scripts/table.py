"""c_step(ratio) one-sided, deep (FEM, two finest levels + Richardson-free: report finest) and shallow (FEM and the
conformal closed form); c_flare(deg) at ratio 2.4 and 1.5 / 3 for the ratio dependence."""
import sys, math
sys.path.insert(0, '.')
import expfem as E
def conformal_one_sided(R):
    r = 1.0 / R
    return ((r + 1 / r) * math.log((1 + r) / (1 - r)) - 2 * math.log(4 * r / (1 - r * r))) / math.pi
what = sys.argv[1]
if what == 'step':
    for W in [float(v) for v in sys.argv[2:]]:
        hs = max(1.0, W / 3.0)
        d, ne = E.excess('step', W, hmin=0.005, hmax=0.08 * hs, limit='deep')
        s, _ = E.excess('step', W, hmin=0.005, hmax=0.08 * hs, limit='shallow')
        print('step ratio %.2f  deep %.4f  shallow %.4f  conformal %.4f  (elems %d)' % (W, d, s, conformal_one_sided(W), ne), flush=True)
else:
    W = float(sys.argv[2])
    for deg in [float(v) for v in sys.argv[3:]]:
        d, ne = E.excess('flare', W, deg, hmin=0.005, hmax=0.08, limit='deep')
        s, _ = E.excess('flare', W, deg, hmin=0.005, hmax=0.08, limit='shallow')
        print('flare ratio %.2f deg %.1f  deep %.4f  shallow %.4f  (elems %d)' % (W, deg, d, s, ne), flush=True)
