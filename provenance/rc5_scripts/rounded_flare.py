"""rc5 (after the rc5 reviews): the stretch-end excess of the flare AS DRAWN (Dinh V1: 25 um tangent arcs at both kinks,
the polyline layout.flare_profile draws and stretch_R integrates) against the sharp-cornered flare the tables were
solved on.  Same Taylor-Hood solver and limits as expfem.py, on a mapped tensor mesh y = eta * w(x) whose wall is the
profile w(x); the 1-D reference integrates the same w(x).  Units of the narrow width (Dinh V1: w 43.994 um, lead
105.75 um, ratio 2.4037, 14.5 deg, arc radius 25 um = 0.568 w).
Usage: python3 rounded_flare.py OPENCIF_ROOT"""
import math, os, sys, time
from types import SimpleNamespace
import numpy as np
R = sys.argv[1]
sys.path[:0] = [os.path.dirname(os.path.abspath(__file__)), os.path.join(R, 'opencif', 'engine')]
import expfem as E
import layout as L
from skfem import MeshTri

w_um, lead_um, deg = 43.99442, 105.75009, 14.5
W = lead_um / w_um
spec = SimpleNamespace(flare_deg=deg, flare_start=67.0)
Lf = L.flare_len(spec, w_um, lead_um)
_, tt = L.flare_kink(spec, Lf)
L_out = 100.0 + Lf + tt                                  # a straight run of 100 um before the flare's first arc
pts = np.array(L.flare_profile(spec, w_um, lead_um, L_out))
s1 = L_out - tt - Lf                                     # the first virtual corner
xs_r, ws_r = (pts[:, 0] - s1) / w_um, pts[:, 1] / w_um   # x = 0 at the first virtual corner, units of w
t = math.tan(math.radians(deg))
Lfw = (W - 1.0) / t


def w_round(x):
    return np.interp(x, xs_r, ws_r, left=1.0, right=W)


def w_sharp(x):
    return np.where(x < 0, 1.0, np.where(x > Lfw, W, 1.0 + x * t))


def run(wfun, hmin, hmax, L1=5.0, L2f=5.0):
    L2 = L2f * W
    x = E.graded(-L1, Lfw + L2, [0.0, Lfw] + [float(v) for v in xs_r[1:-1:4]], hmin, hmax)
    eta = E.graded(0.0, 1.0, [1.0], hmin, hmax)
    m = MeshTri.init_tensor(x, eta)
    X, Y = m.p
    m = MeshTri(np.vstack([X, Y * wfun(X)]), m.t)
    xin, xout = -L1, Lfw + L2
    q = np.linspace(xin, xout, 400001)                   # the 1-D reference on the mesh's own wall (piecewise linear)
    wq = np.interp(q, x, wfun(x))
    R1D_deep = 12.0 * np.trapezoid(1.0 / wq ** 3, q)
    R1D_sh = np.trapezoid(1.0 / wq, q)
    prof_in = lambda y: 6.0 * y * (1.0 - y)
    prof_out = lambda y: 6.0 / W * (y / W) * (1.0 - y / W)
    Phi, _, _ = E.stokes(m, xin, xout, prof_in, prof_out)
    Rsh = E.hele_shaw(m, xin, xout)
    return (Phi - R1D_deep) / 12.0, Rsh - R1D_sh, m.t.shape[1]


print('Dinh V1 flare: ratio %.4f, %.1f deg, arc radius %.3f w, tangent length %.4f w, flare %.3f w' % (W, deg, 25.0 / w_um, tt / w_um, Lfw))
for hmin, hmax in ((0.005, 0.08), (0.0025, 0.06)):
    for lab, f in (('sharp (the tables)', w_sharp), ('as drawn (25 um arcs)', w_round)):
        t0 = time.time()
        d, s, ne = run(f, hmin, hmax)
        print('%-22s mesh hmin %.4f hmax %.3f: deep %.5f shallow %.5f elems %d (%.0f s)' % (lab, hmin, hmax, d, s, ne, time.time() - t0), flush=True)
