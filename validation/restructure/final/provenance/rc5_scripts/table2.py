"""rc5 coefficient tables (deep limit by FEM; shallow by FEM, checked against the conformal closed form for steps)."""
import sys, math, json
sys.path.insert(0, '.')
import expfem as E
what = sys.argv[1]
res = []
if what in ('sym', 'one'):
    for W in [float(v) for v in sys.argv[2:]]:
        hs = max(1.0, W / 3.0)
        d, ne = E.excess('step', W, hmin=0.005, hmax=0.08 * hs, limit='deep', sym=(what == 'sym'))
        s, _ = E.excess('step', W, hmin=0.005, hmax=0.08 * hs, limit='shallow', sym=(what == 'sym'))
        print('%s step ratio %.2f deep %.4f shallow %.4f elems %d' % (what, W, d, s, ne), flush=True)
        res.append([W, d, s])
else:                                                     # flare grid: ratio, then angles
    W = float(sys.argv[2])
    for deg in [float(v) for v in sys.argv[3:]]:
        d, ne = E.excess('flare', W, deg, hmin=0.005, hmax=0.08, limit='deep')
        s, _ = E.excess('flare', W, deg, hmin=0.005, hmax=0.08, limit='shallow')
        print('flare ratio %.2f deg %.1f deep %.4f shallow %.4f elems %d' % (W, deg, d, s, ne), flush=True)
        res.append([W, deg, d, s])
json.dump(res, open('tab_%s_%s.json' % (what, '_'.join(sys.argv[2:3])), 'w'))
