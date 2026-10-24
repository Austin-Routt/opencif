"""mesh convergence of the deep-limit excess (rc5): the one-sided step at ratio 2.71 and the 14.5 deg flare at ratio
2.4, at the table's mesh (hmin 0.005, hmax 0.08 x scale) and at two refinements; and the shallow step against the
conformal closed form at each."""
import sys, math, time
sys.path.insert(0, '.')
import expfem as E
def conformal_one_sided(R):
    r = 1.0 / R
    return ((r + 1 / r) * math.log((1 + r) / (1 - r)) - 2 * math.log(4 * r / (1 - r * r))) / math.pi
for fac in [float(v) for v in (sys.argv[1:] or (1.0, 0.7, 0.5))]:
    t = time.time()
    W = 2.71; hs = max(1.0, W / 3.0)
    d, ne = E.excess('step', W, hmin=0.005 * fac, hmax=0.08 * hs * fac, limit='deep')
    s, _ = E.excess('step', W, hmin=0.005 * fac, hmax=0.08 * hs * fac, limit='shallow')
    print('step 2.71 mesh x%.2f: deep %.5f shallow %.5f (conformal %.5f) elems %d (%.0f s)' % (fac, d, s, conformal_one_sided(W), ne, time.time() - t), flush=True)
    t = time.time()
    d, ne = E.excess('flare', 2.4, 14.5, hmin=0.005 * fac, hmax=0.08 * fac, limit='deep')
    s, _ = E.excess('flare', 2.4, 14.5, hmin=0.005 * fac, hmax=0.08 * fac, limit='shallow')
    print('flare 2.4 14.5 mesh x%.2f: deep %.5f shallow %.5f elems %d (%.0f s)' % (fac, d, s, ne, time.time() - t), flush=True)
