"""R3: convergence and self-check of r3_junction2d: (a) th = 2pi minus th = pi isolates the developed bend (junctions
cancel) and must equal the planar-annulus / Hele-Shaw formula; (b) resolution 40/80/160 for the Dinh hairpin."""
import sys, math
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from r3_junction2d import deep, shallow
import numpy as np
w, rm = 28.0, 27.5
r1, r2 = rm - w / 2, rm + w / 2
c = -0.5
Mx = np.array([[r1, 1 / r1], [r2, 1 / r2]])
A_, B_ = np.linalg.solve(Mx, -np.array([c * r1 * math.log(r1), c * r2 * math.log(r2)]))
F_ = lambda r: c * (r * r / 2 * math.log(r) - r * r / 4) + A_ * r * r / 2 + B_ * math.log(r)
ratio = {'deep': (F_(r2) - F_(r1)) / ((1 / rm) * w ** 3 / 12.0), 'shallow': rm * math.log(r2 / r1) / w}
rs = {'deep': 12.0 / w ** 3, 'shallow': 1.0 / w}
Ls = 2.0 * w
for kind, fn in (('deep', deep), ('shallow', shallow)):
    ex = {}
    for nres in (40, 80, 160):
        for th in (math.pi / 2, math.pi, 2 * math.pi):
            dP, Lt, La = fn(w, rm, th, Ls, w / nres)
            ex[(nres, th)] = dP / rs[kind] - Lt
        print(kind, 'res', nres, 'excess/w  90: %.5f  180: %.5f  360: %.5f   (360-180)/w %.5f vs developed per pi %.5f' % (
            ex[(nres, math.pi/2)] / w, ex[(nres, math.pi)] / w, ex[(nres, 2*math.pi)] / w,
            (ex[(nres, 2*math.pi)] - ex[(nres, math.pi)]) / w, math.pi * rm * (1 / ratio[kind] - 1) / w), flush=True)
    for th, lab in ((math.pi / 2, '90'), (math.pi, '180')):
        e1, e2, e3 = ex[(40, th)], ex[(80, th)], ex[(160, th)]
        p = math.log(abs((e1 - e2) / (e2 - e3))) / math.log(2) if e2 != e3 else float('nan')
        rich = e3 + (e3 - e2) / (2 ** p - 1)
        dev = th * rm * (1 / ratio[kind] - 1)
        junc_iso = 2 * e3 - ex[(160, 2 * math.pi)] if th == math.pi else float('nan')
        print('  %s %s deg: observed order %.2f, extrapolated excess %.5f w, developed %.5f w -> junction %.0f%% of developed' % (
            kind, lab, p, rich / w, dev / w, 100 * (rich - dev) / dev))
