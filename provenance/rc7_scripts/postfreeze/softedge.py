"""Dinh V1's bead data against Gifford 2014 Fig. 5's measured soft edge.
Gifford Fig. 5 (provenance/gifford_fig5_d50.out, the duct-law ceiling): retention R of beads of diameter d in four long
devices of known design cutoff D (R = 1: fully retained).  Dinh ESI Fig. S1B (App. B P-C7.2): 4.7 um beads reach the
filtrate at 1.8 / 3.5 / 6.8 % against the fluid's 50 / 52 / 48 % (106 / 799 / 1714 beads per uL), i.e. a passage of
3.6 / 6.7 / 14.2 % of a tracer's.  Question 1: at which d/D do Gifford's devices retain that much (R = 1 - passage)?
That gives the D Dinh's 4.7 um beads imply.  Question 2: what passage do Gifford's curves predict for 4.7 um beads at
the figure-aligned reading's D (3.67 um since the printed-length refit; 4.29 um before) and at the printed 2.17 um?"""
import math
G5 = {'5A': (4.220, {4: 0.366, 6: 0.753, 8: 0.944, 10: 1.025}),
      '5B': (5.327, {4: 0.098, 6: 0.719, 8: 0.961, 10: 1.029}),
      '5C': (6.982, {4: 0.043, 6: 0.098, 8: 0.821, 10: 0.965}),
      '5D': (9.027, {4: 0.043, 6: 0.043, 8: 0.043, 10: 0.51})}
def curve(D, R):
    pts = sorted((d / D, min(r, 1.0)) for d, r in R.items())
    return pts
def x_at(pts, r):                       # d/D where R first reaches r (linear in log d/D)
    for (x0, r0), (x1, r1) in zip(pts, pts[1:]):
        if r0 < r <= r1:
            t = (r - r0) / (r1 - r0)
            return math.exp(math.log(x0) + t * (math.log(x1) - math.log(x0)))
    return None
def r_at(pts, x):
    if x <= pts[0][0]: return pts[0][1]
    for (x0, r0), (x1, r1) in zip(pts, pts[1:]):
        if x0 <= x <= x1:
            t = (math.log(x) - math.log(x0)) / (math.log(x1) - math.log(x0))
            return r0 + t * (r1 - r0)
    return pts[-1][1]
for conc, bead, fluid in ((106, 1.8, 50), (799, 3.5, 52), (1714, 6.8, 48)):
    passage = bead / fluid; R = 1 - passage
    xs = {k: x_at(curve(*v), R) for k, v in G5.items()}
    Ds = {k: (4.7 / x if x else None) for k, x in xs.items()}
    print('Dinh 4.7 um beads at %4d /uL: passage %.1f %% of the fluid -> R %.3f; Gifford d/D at that R: %s; implied D: %s' % (
        conc, 100 * passage, R, {k: (round(x, 2) if x else 'not reached') for k, x in xs.items()},
        {k: (round(d, 2) if d else '-') for k, d in Ds.items()}))
for D in (2.17, 3.0, 3.67, 4.29, 4.52):
    x = 4.7 / D
    rs = {k: round(r_at(curve(*v), x), 2) for k, v in G5.items()}
    print('D %.2f um: 4.7 um = %.2f D; Gifford Fig. 5 retention there %s -> passage %s %% of the fluid (Dinh: 3.6-14 %%)' % (
        D, x, rs, {k: round(100 * (1 - r), 0) for k, r in rs.items()}))
