"""Independent 2-D check of an outlet expansion's excess resistance (rc5 verification of Reviewer 1's FLARE note).
Taylor-Hood P2/P1 finite elements (scikit-fem) on graded tensor meshes; NOT the reviewer's MAC finite-difference
solver and not the package's streamfunction-vorticity solver.  Two limits of a rectangular duct of uniform depth:
  deep    (planar Stokes, h >> w): dissipation Phi = int |grad u|^2 = dP Q with developed profiles at both ends
  shallow (Hele-Shaw, h << w): Laplace for p, no flux on the walls, conductance per length proportional to w
Excess = (R_2D - R_1D) / r(w_narrow), in units of the narrow width w (w = 1 here).
Geometries (one-sided: the inner wall y = 0 straight, only the outer wall moves, as in opencif's outlet stretch):
  step(W)        : y < 1 for x < 0, y < W for x > 0 (sharp corners)
  flare(W, deg)  : y < w(x), w = 1 up to x = 0, then w = 1 + x tan(deg) up to W, then W
  sym_step(W)    : the symmetric expansion (both walls step), for comparison with the literature limits"""
import math, sys, time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl
from skfem import MeshTri, Basis, ElementTriP2, ElementTriP1, ElementVector, BilinearForm, asm, condense, solve
from skfem.helpers import ddot, grad, dot
from skfem.models.general import divergence
from skfem.models.poisson import laplace


def graded(a, b, foci, hmin, hmax, g=0.08):
    """coordinates a..b, spacing hmin at each focus growing linearly (rate g) to hmax"""
    xs = [a]
    while xs[-1] < b - 1e-12:
        x = xs[-1]
        d = min(abs(x - f) for f in foci) if foci else 1e9
        h = min(hmax, hmin + g * d)
        # do not step over a focus
        for f in foci:
            if x < f < x + h - 1e-12:
                h = f - x
        xs.append(min(x + h, b))
    return np.array(xs)


def mesh_step(W, L1, L2, hmin, hmax, sym=False):
    """one-sided (or symmetric) step: narrow channel width 1 (x<0), wide W (x>0); tensor mesh minus the solid"""
    if sym:
        y0, y1 = -(W - 1.0) / 2.0, (W + 1.0) / 2.0          # narrow channel 0..1, wide y0..y1
        fy = [0.0, 1.0]
    else:
        y0, y1 = 0.0, W
        fy = [1.0]
    x = graded(-L1, L2, [0.0], hmin, hmax)
    y = graded(y0, y1, fy, hmin, hmax)
    m = MeshTri.init_tensor(x, y)
    c = m.p[:, m.t].mean(axis=1)
    solid = (c[0] < 0) & ((c[1] > 1.0) | (c[1] < 0.0))
    m = m.remove_elements(np.nonzero(solid)[0])
    return m


def mesh_flare(W, deg, L1, L2, hmin, hmax, nyx=None):
    """one-sided linear flare: mapped tensor mesh, y = eta * w(x)"""
    t = math.tan(math.radians(deg))
    Lf = (W - 1.0) / t
    x = graded(-L1, Lf + L2, [0.0, Lf], hmin, hmax)
    eta = graded(0.0, 1.0, [1.0], hmin, hmax)
    m = MeshTri.init_tensor(x, eta)
    X, E = m.p
    wx = np.where(X < 0, 1.0, np.where(X > Lf, W, 1.0 + X * t))
    p = np.vstack([X, E * wx])
    return MeshTri(p, m.t), Lf


def wall_of(x, W, kind, deg=None):
    if kind == 'flare':
        t = math.tan(math.radians(deg))
        Lf = (W - 1.0) / t
        return np.where(x < 0, 1.0, np.where(x > Lf, W, 1.0 + x * t))
    return np.where(x < 0, 1.0, W)


def stokes(m, xin, xout, prof_in, prof_out):
    """Taylor-Hood Stokes, Dirichlet velocity on every boundary (developed profiles at the ends, no slip elsewhere);
    returns the dissipation int |grad u|^2 (mu = 1, Q = 1 per unit depth)"""
    e = ElementVector(ElementTriP2())
    ub = Basis(m, e, intorder=4)
    pb = ub.with_element(ElementTriP1())

    @BilinearForm
    def vlap(u, v, w):
        return ddot(grad(u), grad(v))

    A = asm(vlap, ub)
    B = -asm(divergence, ub, pb)
    K = sp.bmat([[A, B.T], [B, None]], 'csr')
    nu = A.shape[0]
    # Dirichlet: all boundary velocity dofs
    D = ub.get_dofs()                       # all boundary facets
    dofs = D.all()
    xd = np.zeros(nu)
    # positions of the velocity dofs: use the basis' doflocs
    loc = ub.doflocs
    for comp, key in ((0, 'u^1'), (1, 'u^2')):
        pass
    # component dofs: ElementVector interleaves; use nodal dofs by name
    ux = D.all('u^1'); uy = D.all('u^2')
    X, Y = loc[0, ux], loc[1, ux]
    vals = np.zeros(len(ux))
    tol = 1e-9
    ii = np.abs(X - xin) < tol
    vals[ii] = prof_in(Y[ii])
    oo = np.abs(X - xout) < tol
    vals[oo] = prof_out(Y[oo])
    xd[ux] = vals
    xd[uy] = 0.0
    # pin one pressure dof
    pin = nu + 0
    xfull = np.zeros(K.shape[0])
    xfull[:nu] = xd
    I = np.setdiff1d(np.arange(K.shape[0]), np.concatenate([dofs, [pin]]))
    sol = solve(*condense(K, np.zeros(K.shape[0]), x=xfull, I=I))
    u = sol[:nu]
    Phi = float(u @ (A @ u))
    return Phi, u, ub


def hele_shaw(m, xin, xout):
    """Laplace for p (mobility 1), p = 1 at x = xin, 0 at xout, no flux elsewhere; returns the resistance 1/Q"""
    b = Basis(m, ElementTriP2(), intorder=4)
    A = asm(laplace, b)
    loc = b.doflocs
    D = b.get_dofs(lambda x: np.abs(x[0] - xin) < 1e-9).all()
    E = b.get_dofs(lambda x: np.abs(x[0] - xout) < 1e-9).all()
    p = np.zeros(A.shape[0])
    p[D] = 1.0
    I = np.setdiff1d(np.arange(A.shape[0]), np.concatenate([D, E]))
    p = solve(*condense(A, np.zeros(A.shape[0]), x=p, I=I))
    Q = float(p @ (A @ p))                  # = flux * (1 - 0)
    return 1.0 / Q


def excess(kind, W, deg=None, hmin=0.004, hmax=0.08, L1=5.0, L2f=5.0, limit='deep', sym=False):
    L2 = L2f * W
    if kind == 'step':
        m = mesh_step(W, L1, L2, hmin, hmax, sym=sym)
        xin, xout = -L1, L2
        if sym:
            y0 = -(W - 1.0) / 2.0
            prof_in = lambda y: 6.0 * y * (1.0 - y)
            prof_out = lambda y: 6.0 / W * ((y - y0) / W) * (1.0 - (y - y0) / W)
        else:
            prof_in = lambda y: 6.0 * y * (1.0 - y)
            prof_out = lambda y: 6.0 / W * (y / W) * (1.0 - y / W)
        R1D_deep = 12.0 * (L1 + L2 / W ** 3)
        R1D_sh = L1 + L2 / W
    elif kind == 'straight':
        m = mesh_step(1.0, L1, L2, hmin, hmax)
        xin, xout = -L1, L2
        prof_in = prof_out = lambda y: 6.0 * y * (1.0 - y)
        R1D_deep = 12.0 * (L1 + L2)
        R1D_sh = L1 + L2
    else:
        m, Lf = mesh_flare(W, deg, L1, L2, hmin, hmax)
        xin, xout = -L1, Lf + L2
        prof_in = lambda y: 6.0 * y * (1.0 - y)
        prof_out = lambda y: 6.0 / W * (y / W) * (1.0 - y / W)
        t = math.tan(math.radians(deg))
        # 1-D integrals over the flare: int 12/w^3 dx = 12 * (1/(2t)) (1 - 1/W^2) ; int dx/w = ln(W)/t
        R1D_deep = 12.0 * (L1 + L2 / W ** 3 + (1.0 - 1.0 / W ** 2) / (2.0 * t))
        R1D_sh = L1 + L2 / W + math.log(W) / t
    if limit == 'deep':
        Phi, _, _ = stokes(m, xin, xout, prof_in, prof_out)
        return (Phi - R1D_deep) / 12.0, m.t.shape[1]
    R = hele_shaw(m, xin, xout)
    return R - R1D_sh, m.t.shape[1]


if __name__ == '__main__':
    args = sys.argv[1:]
    kind = args[0]
    W = float(args[1])
    deg = float(args[2]) if kind == 'flare' else None
    limit = args[3] if len(args) > 3 else 'deep'
    sym = len(args) > 4 and args[4] == 'sym'
    import os
    hs = float(os.environ.get('HMAXSCALE', '1'))
    levels = ((0.01, 0.12), (0.005, 0.08), (0.0025, 0.06))[:int(os.environ.get('NLEV', '3'))]
    for hmin, hmax in levels:
        hmax *= hs
        t0 = time.time()
        c, ne = excess(kind, W, deg, hmin=hmin, hmax=hmax, limit=limit, sym=sym)
        print('%-8s W %.3f deg %s %-7s %s hmin %.4f hmax %.3f elems %7d  excess %+.5f w   (%.0f s)'
              % (kind, W, deg, limit, 'sym' if sym else 'one-sided', hmin, hmax, ne, c, time.time() - t0), flush=True)
