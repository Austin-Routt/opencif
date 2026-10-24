"""rc5: opencif's stretch-end tables (cifsched) against the solver outputs saved beside this script.
Usage: python3 check_tables.py <opencif tree>   (run from this folder)"""
import json, sys, os
T = sys.argv[1]
sys.path[:0] = [T, os.path.join(T, 'opencif'), os.path.join(T, 'opencif', 'engine')]
import cifsched as S
C = json.load(open('coefficients_all.json'))
sym = {float(k): v for k, v in C['sym'].items()}
for line in open('sym_24_8.out'):                       # the two centred entries re-solved in rc5
    w = line.split()
    sym[float(w[3])] = [float(w[5]), float(w[7])]
one = {float(k): v for k, v in C['one'].items()}
worst = 0.0
worst_step = 0.0
for tab, src, sh_one in ((S.STEP_ONE_DEEP, one, True), (S.STEP_SYM_DEEP, sym, False)):
    for R, v in tab:
        if R == 1.0:
            continue
        d, s = src[R]
        e1, e2 = abs(d - v), abs(S.step_shallow(R, sh_one) - s)
        worst = max(worst, e1)
        worst_step = max(worst_step, e1)
        print('%s R %.2f  table %.4f  solver %.5f  |d| %.1e   shallow closed form %.5f  solver %.5f  |d| %.1e'
              % ('one' if sh_one else 'sym', R, v, d, e1, S.step_shallow(R, sh_one), s, e2))
print('step tables (deep entries): largest |table - solver| %.1e' % worst_step)
worst_flare = 0.0
for i, R in enumerate(S.FLARE_RATIOS):
    if R == 1.0:
        continue
    t = json.load(open('tab_flare_%s.json' % ('%g' % R if R not in (2.0, 3.0, 4.0) else '%.1f' % R)))
    for RR, deg, d, s in t:
        j = S.FLARE_DEGS.index(deg)
        e = max(abs(S.FLARE_DEEP[i][j] - d), abs(S.FLARE_SHALLOW[i][j] - s))
        worst = max(worst, e)
        worst_flare = max(worst_flare, e)
print('flare tables (deep and shallow entries): every entry checked; largest |table - solver| %.1e' % worst_flare)
print('all tables: largest |table - solver| %.1e' % worst)
