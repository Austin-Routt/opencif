"""R4: the head's T-junction (drop of unit k -> out line of unit k+1, with the dead-ended mouth below) against the 1-D
convention (centreline through the junction centre J), in the two limits of a rectangular duct:
  shallow (Hele-Shaw): cell-centred finite volumes for p, exact for the rectilinear geometry;
  deep (planar Stokes): psi-omega on a node grid, Thom wall vorticity, pressure drop integrated along the centreline.
Excess = (2-D resistance - 1-D resistance of the centreline length) / r, in units of w.  Geometries:
  'L'  : a plain sharp 90-degree corner (no stem)        -- a check against known values
  'T2' : the head junction, dead stem of 2 w below J (Dinh l_g/w = 1.9, Gifford-head 2-3)"""
import sys, math
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl

def rects(kind, Ld, Lu, stem):
    # w = 1; drop: x in [-Ld, 0.5], y in [-0.5, 0.5]; out line: x in [-0.5, 0.5], y in [-0.5 - stem, Lu]
    R = [(-Ld, 0.5, -0.5, 0.5), (-0.5, 0.5, -0.5 - stem, Lu)]
    return R

def inside(R, x, y, eps=1e-9):
    return any(a - eps <= x <= b + eps and c - eps <= y <= d + eps for a, b, c, d in R)

def hele_shaw(n, Ld=4.0, Lu=4.0, stem=2.0):
    h = 1.0 / n
    R = rects(None, Ld, Lu, stem)
    xs = np.arange(-Ld + h / 2, 0.5, h); ys = np.arange(-0.5 - stem + h / 2, Lu, h)
    fl = np.array([[inside(R, x, y) for y in ys] for x in xs])
    idx = -np.ones(fl.shape, int); idx[fl] = np.arange(fl.sum())
    nf = fl.sum(); I, J, V = [], [], []; b = np.zeros(nf)
    nx, ny = fl.shape
    for i in range(nx):
        for j in range(ny):
            if not fl[i, j]:
                continue
            r = idx[i, j]; diag = 0.0
            for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                ii, jj = i + di, j + dj
                if 0 <= ii < nx and 0 <= jj < ny and fl[ii, jj]:
                    I.append(r); J.append(idx[ii, jj]); V.append(1.0); diag -= 1.0
                elif di == -1 and ii < 0:                     # inlet face (x = -Ld): p = 1, half cell
                    diag -= 2.0; b[r] -= 2.0 * 1.0
                elif dj == 1 and jj >= ny and abs(xs[i]) < 0.5:  # outlet face (y = Lu): p = 0
                    diag -= 2.0
            I.append(r); J.append(r); V.append(diag)
    A = sp.csr_matrix((V, (I, J)), shape=(nf, nf))
    p = spl.spsolve(A.tocsc(), b)
    P = np.full(fl.shape, np.nan); P[fl] = p
    # flux through the inlet face: sum over inlet cells of (1 - p)/(h/2) * h
    q = sum((1.0 - P[0, j]) / (h / 2) * h for j in range(ny) if fl[0, j])
    # straight resistance per length = 1/w = 1 (w = 1); 1-D length through J (0,0): Ld + Lu
    return 1.0 / q - (Ld + Lu)

def deep(n, Ld=4.0, Lu=4.0, stem=2.0):
    h = 1.0 / n
    R = rects(None, Ld, Lu, stem)
    xs = np.round(np.arange(-Ld, 0.5 + h / 2, h), 12); ys = np.round(np.arange(-0.5 - stem, Lu + h / 2, h), 12)
    nx, ny = len(xs), len(ys)
    fl = np.array([[inside(R, x, y) for y in ys] for x in xs])
    def isfl(i, j): return 0 <= i < nx and 0 <= j < ny and fl[i, j]
    interior = np.array([[fl[i, j] and all(isfl(i + a, j + c) for a in (-1, 0, 1) for c in (-1, 0, 1)) for j in range(ny)] for i in range(nx)])
    inlet = np.array([[fl[i, j] and i == 0 for j in range(ny)] for i in range(nx)])
    outlet = np.array([[fl[i, j] and j == ny - 1 for j in range(ny)] for i in range(nx)])
    # Poiseuille, Q = 1, psi from 0 (lower/right wall chain) to 1 (upper/left wall chain)
    prof = lambda s: 0.5 + 1.5 * s - 2.0 * s ** 3           # s in [-0.5, 0.5]
    omg = lambda s: 12.0 * s                               # omega = -psi'' = 12 s (w = 1)
    # wall psi: lower/right chain = 0 (drop bottom y=-0.5, stem walls and bottom, out line right x=0.5);
    # upper/left chain = 1 (drop top y=0.5 for x<=-0.5, out line left x=-0.5 for y>=0.5)
    def wall_psi(x, y):
        if (abs(y - 0.5) < 1e-9 and x <= -0.5 + 1e-9) or (abs(x + 0.5) < 1e-9 and y >= 0.5 - 1e-9):
            return 1.0
        return 0.0
    bnd = fl & ~interior
    uid = {}
    for i in range(nx):
        for j in range(ny):
            if interior[i, j]:
                uid[('p', i, j)] = len(uid)
    for i in range(nx):
        for j in range(ny):
            if fl[i, j]:
                uid[('w', i, j)] = len(uid)                 # omega everywhere in the fluid
    NU = len(uid); I, J, V = [], [], []; b = np.zeros(NU)
    def psi_known(i, j):
        if inlet[i, j]: return prof(ys[j])
        if outlet[i, j]: return prof(-xs[i])                 # outlet: s across = -x (psi=1 at x=-0.5)
        return wall_psi(xs[i], ys[j])
    row = 0
    eqs = []
    for i in range(nx):
        for j in range(ny):
            if interior[i, j]:
                # Laplace psi + omega = 0
                r = uid[('p', i, j)]
                for (a, c) in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                    ii, jj = i + a, j + c
                    if interior[ii, jj]:
                        I.append(r); J.append(uid[('p', ii, jj)]); V.append(1 / h ** 2)
                    else:
                        b[r] -= psi_known(ii, jj) / h ** 2
                I.append(r); J.append(r); V.append(-4 / h ** 2)
                I.append(r); J.append(uid[('w', i, j)]); V.append(1.0)
                # Laplace omega = 0
                r2 = uid[('w', i, j)]
                for (a, c) in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                    I.append(r2); J.append(uid[('w', i + a, j + c)]); V.append(1 / h ** 2)
                I.append(r2); J.append(r2); V.append(-4 / h ** 2)
            elif fl[i, j]:
                r2 = uid[('w', i, j)]
                if inlet[i, j] or outlet[i, j]:
                    s = ys[j] if inlet[i, j] else -xs[i]
                    I.append(r2); J.append(r2); V.append(1.0); b[r2] = omg(s) if inlet[i, j] else omg(s)
                    continue
                # wall node: Thom, averaged over the interior neighbours along the axes
                nb = [(a, c) for (a, c) in ((1, 0), (-1, 0), (0, 1), (0, -1)) if 0 <= i + a < nx and 0 <= j + c < ny and interior[i + a, j + c]]
                I.append(r2); J.append(r2); V.append(1.0)
                if not nb:
                    continue                                 # convex corner: unused
                pw = psi_known(i, j)
                for (a, c) in nb:
                    I.append(r2); J.append(uid[('p', i + a, j + c)]); V.append(2.0 / h ** 2 / len(nb))
                    b[r2] += 2.0 * pw / h ** 2 / len(nb)
    A = sp.csr_matrix((V, (I, J)), shape=(NU, NU))
    x = spl.spsolve(A.tocsc(), b)
    W = np.full((nx, ny), np.nan)
    for (t, i, j), u in uid.items():
        if t == 'w': W[i, j] = x[u]
    # pressure drop along the centreline: dp/dx = -domega/dy (y = 0, x from -Ld to 0); dp/dy = domega/dx (x = 0, y 0..Lu)
    j0 = int(np.argmin(np.abs(ys))); i0 = int(np.argmin(np.abs(xs)))
    dwdy = (W[:, j0 + 1] - W[:, j0 - 1]) / (2 * h)
    seg1 = -np.trapezoid(dwdy[: i0 + 1], xs[: i0 + 1])
    dwdx = (W[i0 + 1, :] - W[i0 - 1, :]) / (2 * h)
    seg2 = np.trapezoid(dwdx[j0:], ys[j0:])
    dP = -(seg1 + seg2)            # p_in - p_out
    return dP / 12.0 - (Ld + Lu)   # straight: dp/ds = 12 Q / w^3 = 12

if __name__ == '__main__':
    for stem in (0.0, 2.0):
        lab = 'L corner (no stem)' if stem == 0 else 'head T (dead stem 2w)'
        hs = [hele_shaw(n, stem=stem) for n in (10, 20, 40)]
        print('%-24s shallow: excess/w %s -> Richardson %.4f' % (lab, ['%.4f' % v for v in hs], hs[2] + (hs[2] - hs[1]) / 3), flush=True)
        dp = [deep(n, stem=stem) for n in (10, 20, 40)]
        r1 = dp[2] + (dp[2] - dp[1]) / 3; r0 = dp[2] + (dp[2] - dp[1])
        print('%-24s deep   : excess/w %s -> Richardson %.4f (1st order %.4f)' % (lab, ['%.4f' % v for v in dp], r1, r0), flush=True)
