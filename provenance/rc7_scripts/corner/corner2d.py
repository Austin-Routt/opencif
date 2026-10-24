"""rc7 (DINH-1): the head's square quarter-turn against the 1-D convention (centreline through the corner point J), in
the two limits of a rectangular duct, as rc6's junction fits (provenance/junction_fits/r4_tjunc2d.py, r4b_fillet.py):
  shallow (Hele-Shaw): cell-centred finite volumes for p;
  deep (planar Stokes): psi-omega on a node grid, Thom wall vorticity, the pressure drop integrated along the centreline.
Geometry (w = 1): in along +x (y in [-1/2, 1/2]) for Ld, a 90-degree turn, out along +y (x in [-1/2, 1/2]) for Lu.
  r_out: rounding of the turn's OUTER (convex) corner at (1/2, -1/2), the fluid's corner (Dinh Fig. 2B: 6.0 / 8.6 um at
         w 28 um, r/w 0.21 / 0.31);
  r_in : fillet of the INNER (re-entrant) corner at (-1/2, 1/2) (the drawing's process fillet);
  arc  : instead, the constant-width arc of centreline radius rm about (-1/2 - (rm - 1/2), 1/2 + ...) -- the rc6 drawing,
         used to validate the solver against rc6's exact curved-duct correction (junction2d.log 'Dinh hairpin quarter
         arc': deep +0.0431 w, shallow -0.1112 w relative to the arc's centreline length).
Excess c = (2-D resistance - 1-D resistance of the centreline length) / r, in units of w; the centreline length is
Ld + Lu for the square turn (through J) and Ld' + Lu' + rm pi/2 for the arc.
Usage: python3 -I corner2d.py [n1 n2]"""
import math
import sys
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl


def make_inside(kind, r_out=0.0, r_in=0.0, rm=None, Ld=3.0, Lu=3.0, method=None):
    if kind == 'square':
        R = [(-Ld, 0.5, -0.5, 0.5), (-0.5, 0.5, -0.5, Lu)]

        def ins(x, y, eps=1e-9):
            if not any(a - eps <= x <= b + eps and c - eps <= y <= d + eps for a, b, c, d in R):
                # inner fillet: solid quadrant x < -1/2, y > 1/2 gains fluid? no: a fillet ADDS solid at a re-entrant
                # corner of the fluid; here the re-entrant fluid corner is at (-1/2, 1/2) seen from the fluid: the fluid
                # occupies x > -1/2 above y = 1/2 and y < 1/2 left of x = -1/2, so the solid quadrant's corner is the
                # convex solid corner; the fillet rounds it by adding FLUID? -> see below
                if r_in > 0 and -0.5 - r_in - eps <= x <= -0.5 + eps and 0.5 - eps <= y <= 0.5 + r_in + eps \
                        and math.hypot(x + 0.5 + r_in, y - 0.5 - r_in) >= r_in - eps:
                    return True          # fillet: fluid in the solid quadrant's corner (rounds the solid's convex corner)
                return False
            if r_out > 0 and x >= 0.5 - r_out - eps and y <= -0.5 + r_out + eps \
                    and math.hypot(x - 0.5 + r_out, y + 0.5 - r_out) > r_out + eps:
                return False             # outer corner rounding: the fluid's convex corner cut to radius r_out
            return True
        return ins, Ld + Lu
    # arc: centre at (cx, cy) = (-Ld0, Lu0)?  Build in a frame where the straight in-line ends at x = 0 and the arc
    # (centre (0, rm)) turns from +x to +y, ending at (rm, rm); the out line continues up from y = rm.
    a1, a2 = rm - 0.5, rm + 0.5
    def ins(x, y, eps=1e-9):
        if x <= 0 and -0.5 - eps <= y <= 0.5 + eps and x >= -Ld - eps:
            return True                                    # in line, y = 0 is its centreline
        if y >= rm and rm - 0.5 - eps <= x <= rm + 0.5 + eps and y <= rm + Lu + eps:
            return True                                    # out line, x = rm its centreline
        if x >= 0 and y <= rm:
            r = math.hypot(x, y - rm)
            return a1 - eps <= r <= a2 + eps
        return False
    return ins, Ld + Lu + rm * math.pi / 2


def hele_shaw(n, kind, **kw):
    h = 1.0 / n
    ins, L1 = make_inside(kind, **kw)
    Ld, Lu = kw.get('Ld', 3.0), kw.get('Lu', 3.0)
    if kind == 'square':
        x0, x1, y0, y1 = -Ld, 0.5, -0.5, Lu
        out_x = (-0.5, 0.5)
    else:
        rm = kw['rm']
        x0, x1, y0, y1 = -Ld, rm + 0.5, -0.5, rm + Lu
        out_x = (rm - 0.5, rm + 0.5)
    xs = np.arange(x0 + h / 2, x1, h); ys = np.arange(y0 + h / 2, y1, h)
    fl = np.array([[ins(x, y, 0.0) for y in ys] for x in xs])
    idx = -np.ones(fl.shape, int); idx[fl] = np.arange(fl.sum()); nf = int(fl.sum())
    I, J, V = [], [], []; b = np.zeros(nf); nx, ny = fl.shape
    for i in range(nx):
        for j in range(ny):
            if not fl[i, j]:
                continue
            r = idx[i, j]; diag = 0.0
            for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                ii, jj = i + di, j + dj
                if 0 <= ii < nx and 0 <= jj < ny and fl[ii, jj]:
                    I.append(r); J.append(idx[ii, jj]); V.append(1.0); diag -= 1.0
                elif di == -1 and ii < 0:                      # inlet face x = x0: p = 1
                    diag -= 2.0; b[r] -= 2.0
                elif dj == 1 and jj >= ny and out_x[0] < xs[i] < out_x[1]:   # outlet face y = y1: p = 0
                    diag -= 2.0
            I.append(r); J.append(r); V.append(diag)
    p = spl.spsolve(sp.csr_matrix((V, (I, J)), shape=(nf, nf)).tocsc(), b)
    P = np.full(fl.shape, np.nan); P[fl] = p
    q = sum((1.0 - P[0, j]) / (h / 2) * h for j in range(ny) if fl[0, j])
    return 1.0 / q - L1


def deep(n, kind, **kw):
    """planar Stokes, psi-omega (Thom), the pressure drop from the vorticity's normal derivative along the centreline"""
    h = 1.0 / n
    ins, L1 = make_inside(kind, **kw)
    Ld, Lu = kw.get('Ld', 3.0), kw.get('Lu', 3.0)
    if kind == 'square':
        x0, x1, y0, y1 = -Ld, 0.5, -0.5, Lu
        xc_out = 0.0
    else:
        rm = kw['rm']
        x0, x1, y0, y1 = -Ld, rm + 0.5, -0.5, rm + Lu
        xc_out = rm
    xs = np.round(np.arange(x0, x1 + h / 2, h), 12); ys = np.round(np.arange(y0, y1 + h / 2, h), 12)
    nx, ny = len(xs), len(ys)
    fl = np.array([[ins(x, y) for y in ys] for x in xs])
    isfl = lambda i, j: 0 <= i < nx and 0 <= j < ny and fl[i, j]
    interior = np.array([[fl[i, j] and all(isfl(i + a, j + c) for a in (-1, 0, 1) for c in (-1, 0, 1))
                          for j in range(ny)] for i in range(nx)])
    inlet = np.zeros_like(fl); inlet[0, :] = fl[0, :]
    outlet = np.zeros_like(fl); outlet[:, -1] = fl[:, -1]
    prof = lambda s: 0.5 + 1.5 * s - 2.0 * s ** 3
    omg = lambda s: 12.0 * s
    # wall psi: the inner (left/top) wall chain psi = 1, the outer (bottom/right) chain psi = 0.  A wall node belongs to
    # the inner chain when it lies on the inner side of the centreline path (nearer the turn's centre).
    if kind == 'square':
        def inner_side(x, y):
            return (y > 0 and x < 0) if True else False
    else:
        def inner_side(x, y):
            if x <= 0:
                return y > 0
            if y >= rm:
                return x < rm
            return math.hypot(x, y - rm) < rm
    def psi_known(i, j):
        if inlet[i, j]:
            return prof(ys[j])
        if outlet[i, j]:
            return prof(-(xs[i] - xc_out))
        return 1.0 if inner_side(xs[i], ys[j]) else 0.0
    uid = {}
    for i in range(nx):
        for j in range(ny):
            if interior[i, j]:
                uid[('p', i, j)] = len(uid)
    for i in range(nx):
        for j in range(ny):
            if fl[i, j]:
                uid[('w', i, j)] = len(uid)
    NU = len(uid); I, J, V = [], [], []; b = np.zeros(NU)
    for i in range(nx):
        for j in range(ny):
            if interior[i, j]:
                r = uid[('p', i, j)]
                for (a, c) in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                    ii, jj = i + a, j + c
                    if interior[ii, jj]:
                        I.append(r); J.append(uid[('p', ii, jj)]); V.append(1 / h ** 2)
                    else:
                        b[r] -= psi_known(ii, jj) / h ** 2
                I.append(r); J.append(r); V.append(-4 / h ** 2)
                I.append(r); J.append(uid[('w', i, j)]); V.append(1.0)
                r2 = uid[('w', i, j)]
                for (a, c) in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                    I.append(r2); J.append(uid[('w', i + a, j + c)]); V.append(1 / h ** 2)
                I.append(r2); J.append(r2); V.append(-4 / h ** 2)
            elif fl[i, j]:
                r2 = uid[('w', i, j)]
                if inlet[i, j] or outlet[i, j]:
                    s = ys[j] if inlet[i, j] else -(xs[i] - xc_out)
                    I.append(r2); J.append(r2); V.append(1.0); b[r2] = omg(s); continue
                nb = [(a, c) for (a, c) in ((1, 0), (-1, 0), (0, 1), (0, -1))
                      if 0 <= i + a < nx and 0 <= j + c < ny and interior[i + a, j + c]]
                I.append(r2); J.append(r2); V.append(1.0)
                if not nb:
                    continue
                pw = psi_known(i, j)
                for (a, c) in nb:
                    I.append(r2); J.append(uid[('p', i + a, j + c)]); V.append(2.0 / h ** 2 / len(nb))
                    b[r2] += 2.0 * pw / h ** 2 / len(nb)
    x = spl.spsolve(sp.csr_matrix((V, (I, J)), shape=(NU, NU)).tocsc(), b)
    W = np.full((nx, ny), np.nan)
    for (t, i, j), u in uid.items():
        if t == 'w':
            W[i, j] = x[u]
    if kw.get('method') == 'dissipation':
        # Phi = mu int omega^2 dA = dP Q (no-slip walls, fully developed ends); node quadrature, a node's weight the
        # share of its four surrounding cells that lie in the fluid
        wgt = np.zeros((nx, ny))
        for i in range(nx - 1):
            for j in range(ny - 1):
                if fl[i, j] and fl[i + 1, j] and fl[i, j + 1] and fl[i + 1, j + 1]:
                    for (a, c) in ((0, 0), (1, 0), (0, 1), (1, 1)):
                        wgt[i + a, j + c] += 0.25 * h * h
        Phi = float(np.nansum(np.where(fl, W, 0.0) ** 2 * wgt))
        return Phi / 12.0 - L1
    # dp/ds = -d(omega)/dn (mu = 1) along the centreline: the in line (y = 0, x0 -> 0), [the arc], the out line
    j0 = int(np.argmin(np.abs(ys)))
    if kind == 'square':
        i0 = int(np.argmin(np.abs(xs)))
        seg1 = -np.trapezoid(((W[:, j0 + 1] - W[:, j0 - 1]) / (2 * h))[: i0 + 1], xs[: i0 + 1])
        seg2 = np.trapezoid(((W[i0 + 1, :] - W[i0 - 1, :]) / (2 * h))[j0:], ys[j0:])
        return -(seg1 + seg2) / 12.0 - L1
    # arc: integrate along the in line to x = 0, along the arc (polar about (0, rm)), and the out line from y = rm
    i_0 = int(np.argmin(np.abs(xs - 0.0)))
    seg1 = -np.trapezoid(((W[:, j0 + 1] - W[:, j0 - 1]) / (2 * h))[: i_0 + 1], xs[: i_0 + 1])
    from scipy.interpolate import RegularGridInterpolator
    Wf = np.where(np.isnan(W), 0.0, W)
    G = RegularGridInterpolator((xs, ys), Wf)
    th = np.linspace(-math.pi / 2, 0.0, 400)
    dn = 1e-3
    vals = []
    for t in th:
        ux, uy = math.cos(t), math.sin(t)                  # outward radial unit vector from the centre (0, rm)
        pa = (rm + dn) * ux, rm + (rm + dn) * uy
        pb = (rm - dn) * ux, rm + (rm - dn) * uy
        vals.append(float(G([pa])[0] - G([pb])[0]) / (2 * dn))
    seg_arc = -np.trapezoid(np.array(vals) * rm, th) * 1.0
    i_r = int(np.argmin(np.abs(xs - rm)))
    jr = int(np.argmin(np.abs(ys - rm)))
    seg2 = np.trapezoid(((W[i_r + 1, :] - W[i_r - 1, :]) / (2 * h))[jr:], ys[jr:])
    return -(seg1 + seg_arc + seg2) / 12.0 - L1


if __name__ == '__main__':
    ns = [int(v) for v in sys.argv[1:3]] or [20, 40]
    print('sharp L (rc6 r4_tjunc2d stem 0: deep -0.2419 extrapolated, shallow -0.4403):')
    for n in ns:
        print('  n %d: shallow %.4f  deep %.4f' % (n, hele_shaw(n, 'square'), deep(n, 'square')), flush=True)
