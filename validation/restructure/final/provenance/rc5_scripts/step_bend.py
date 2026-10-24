"""rc5 (after the rc5 reviews; R1's E2): does the lead's first quarter arc, which starts right at the stretch's end,
change the step's excess?  The tables solve the step with a STRAIGHT wide channel after it; opencif adds the step's
excess and the arc's own bend term separately.  Here the same Taylor-Hood solver as expfem.py (deep limit, planar
Stokes; shallow limit, Hele-Shaw) solves four channels and isolates the interaction:
  (a) narrow inlet -> one-sided step -> outward quarter arc of the wide lead -> straight exit
  (b) wide inlet (no step)            -> the same arc                          -> straight exit
  (c) narrow inlet -> step -> straight wide piece as long as the arc's centreline -> straight exit
  (d) wide inlet  -> straight wide piece                                          -> straight exit
interaction = [(a) - (b)] - [(c) - (d)]: the step's excess with the arc after it minus its excess with a straight
channel after it, in narrow widths (the unit of the tables' c).  The arc turns toward the step's outer wall (opencif's
teardrop and separate leads turn outward, away from the axis), with opencif's centreline radius 2 lead widths.
Usage: [HMIN=0.005 HMAX=0.08] python3 step_bend.py W [W ...]   (W = lead width / side-channel width)"""
import math, os, sys, time
import numpy as np
import scipy.sparse as sp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import expfem as E
from skfem import MeshTri, Basis, ElementTriP2, ElementTriP1, ElementVector, BilinearForm, asm, condense, solve
from skfem.helpers import ddot, grad
from skfem.models.general import divergence
from skfem.models.poisson import laplace


def join(parts):
    P = np.hstack([p for p, t in parts]); T = []; off = 0
    for p, t in parts:
        T.append(t + off); off += p.shape[1]
    T = np.hstack(T)
    key = np.round(P.T, 9)
    _, idx, inv = np.unique(key, axis=0, return_index=True, return_inverse=True)
    P2 = P[:, idx]; T2 = inv.ravel()[T]
    a = P2[:, T2]
    area = (a[0, 1] - a[0, 0]) * (a[1, 2] - a[1, 0]) - (a[0, 2] - a[0, 0]) * (a[1, 1] - a[1, 0])
    neg = area < 0
    T2[[1, 2]][:, neg] = T2[[2, 1]][:, neg]
    T2[1, neg], T2[2, neg] = T2[2, neg].copy(), T2[1, neg].copy()
    return MeshTri(P2, T2)


def tensor(xs, ys):
    m = MeshTri.init_tensor(xs, ys)
    return m.p, m.t


def build(W, bend, wide_inlet, hmin=0.005, hmax=0.08, L1=5.0):
    s = max(1.0, W / 3.0)
    yW = E.graded(0.0, W, [1.0], hmin, hmax * s)                  # lateral nodes of the wide lead, node at y = 1
    yin = yW if wide_inlet else yW[yW <= 1.0 + 1e-12]
    xin = E.graded(-L1, 0.0, [0.0], hmin, hmax * s)
    parts = [tensor(xin, yin)]
    rq = 2.0 * W                                                  # opencif: centreline radius 2 lead widths
    Lb = math.pi * rq / 2.0
    L2 = 5.0 * W
    if bend:
        Yc = W / 2.0 + rq                                          # arc centre on the step's outer-wall side
        phis = E.graded(-math.pi / 2.0, 0.0, [-math.pi / 2.0], hmin / rq, hmax * s / rq)
        p, t = tensor(phis, yW)                                   # (phi, y): radius r = Yc - y
        ph, yy = p
        r = Yc - yy
        parts.append((np.vstack([r * np.cos(ph), Yc + r * np.sin(ph)]), t))
        ex = E.graded(0.0, L2, [], hmin, hmax * s)
        p, t = tensor(ex, yW)                                     # exit: x = r, y = Yc + s
        e_, yy = p
        parts.append((np.vstack([Yc - yy, Yc + e_]), t))
        out = ('y', Yc + L2, (Yc - W, Yc))
    else:
        xs = E.graded(0.0, Lb + L2, [0.0], hmin, hmax * s)
        parts.append(tensor(xs, yW))
        out = ('x', Lb + L2, (0.0, W))
    return join(parts), -L1, (W if wide_inlet else 1.0), out


def stokes_dissipation(m, xin, win, out):
    e = ElementVector(ElementTriP2()); ub = Basis(m, e, intorder=4); pb = ub.with_element(ElementTriP1())
    A = asm(BilinearForm(lambda u, v, w: ddot(grad(u), grad(v))), ub)
    B = -asm(divergence, ub, pb)
    K = sp.bmat([[A, B.T], [B, None]], 'csr'); nu = A.shape[0]
    D = ub.get_dofs(); ux = D.all('u^1'); uy = D.all('u^2'); loc = ub.doflocs
    xd = np.zeros(nu)
    X, Y = loc[0, ux], loc[1, ux]
    ii = np.abs(X - xin) < 1e-9
    xd[ux[ii]] = 6.0 / win * (Y[ii] / win) * (1.0 - Y[ii] / win)
    ax, val, (a, b) = out
    wo = b - a
    if ax == 'x':
        oo = np.abs(X - val) < 1e-9
        xd[ux[oo]] = 6.0 / wo * ((Y[oo] - a) / wo) * (1.0 - (Y[oo] - a) / wo)
    else:
        X2, Y2 = loc[0, uy], loc[1, uy]
        oo = np.abs(Y2 - val) < 1e-9
        xd[uy[oo]] = 6.0 / wo * ((X2[oo] - a) / wo) * (1.0 - (X2[oo] - a) / wo)
    xf = np.zeros(K.shape[0]); xf[:nu] = xd
    I = np.setdiff1d(np.arange(K.shape[0]), np.concatenate([D.all(), [nu]]))
    sol = solve(*condense(K, np.zeros(K.shape[0]), x=xf, I=I))
    u = sol[:nu]
    return float(u @ (A @ u))


def hele_shaw_R(m, xin, out):
    b = Basis(m, ElementTriP2(), intorder=4); A = asm(laplace, b)
    Din = b.get_dofs(lambda x: np.abs(x[0] - xin) < 1e-9).all()
    ax, val, _ = out
    Dout = b.get_dofs(lambda x: np.abs(x[0 if ax == 'x' else 1] - val) < 1e-9).all()
    p = np.zeros(A.shape[0]); p[Din] = 1.0
    I = np.setdiff1d(np.arange(A.shape[0]), np.concatenate([Din, Dout]))
    p = solve(*condense(A, np.zeros(A.shape[0]), x=p, I=I))
    return 1.0 / float(p @ (A @ p))


if __name__ == '__main__':
    for W in [float(v) for v in sys.argv[1:]]:
        t0 = time.time(); res = {}
        hmin, hmax = float(os.environ.get('HMIN', '0.005')), float(os.environ.get('HMAX', '0.08'))
        for lab, bend, wide in (('a', True, False), ('b', True, True), ('c', False, False), ('d', False, True)):
            m, xin, win, out = build(W, bend, wide, hmin=hmin, hmax=hmax)
            res[lab] = (stokes_dissipation(m, xin, win, out) / 12.0, hele_shaw_R(m, xin, out), m.t.shape[1])
        deep_step_bend = res['a'][0] - res['b'][0]; deep_step = res['c'][0] - res['d'][0]
        sh_step_bend = res['a'][1] - res['b'][1]; sh_step = res['c'][1] - res['d'][1]
        # (c) - (d) also holds the inlet's narrow-vs-wide length: L1 (1 - 1/W^3) deep, L1 (1 - 1/W) shallow
        L1 = 5.0
        c_deep = deep_step - L1 * (1.0 - 1.0 / W ** 3); c_sh = sh_step - L1 * (1.0 - 1.0 / W)
        print('mesh hmin %.4f hmax %.3f  W %.3f (arc centreline radius %.2f w): step excess with a straight lead: deep %.4f shallow %.4f;'
              ' interaction with the arc: deep %+.4f (%+.1f %% of c) shallow %+.4f (%+.1f %% of c); elems %d-%d (%.0f s)'
              % (hmin, hmax, W, 2 * W, c_deep, c_sh, deep_step_bend - deep_step, 100 * (deep_step_bend - deep_step) / c_deep,
                 sh_step_bend - sh_step, 100 * (sh_step_bend - sh_step) / c_sh, min(v[2] for v in res.values()),
                 max(v[2] for v in res.values()), time.time() - t0), flush=True)
