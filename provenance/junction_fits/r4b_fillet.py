"""R4b: effect of the drawing's fillets (r_f, Gm.smooth closing+opening) at the head T-junction's two re-entrant corners
(inner corner of the turn at (-1/2, 1/2) and drop/stem corner at (-1/2, -1/2); w = 1) on c_T, Hele-Shaw (cell FV) and
deep planar Stokes (psi-omega, staircase walls).  Same geometry as r4_tjunc2d 'T2' (dead stem 2 w)."""
import sys, math
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl

def make_inside(rf, Ld=3.0, Lu=3.0, stem=2.0):
    R = [(-Ld, 0.5, -0.5, 0.5), (-0.5, 0.5, -0.5 - stem, Lu)]
    def ins(x, y, eps=1e-9):
        if any(a - eps <= x <= b + eps and c - eps <= y <= d + eps for a, b, c, d in R):
            return True
        if rf > 0:
            # top-left re-entrant corner: solid quadrant x < -0.5, y > 0.5; fillet centre (-0.5 - rf, 0.5 + rf)
            if -0.5 - rf - eps <= x <= -0.5 + eps and 0.5 - eps <= y <= 0.5 + rf + eps and math.hypot(x + 0.5 + rf, y - 0.5 - rf) >= rf - eps:
                return True
            # bottom-left: solid quadrant x < -0.5, y < -0.5; centre (-0.5 - rf, -0.5 - rf)
            if -0.5 - rf - eps <= x <= -0.5 + eps and -0.5 - rf - eps <= y <= -0.5 + eps and math.hypot(x + 0.5 + rf, y + 0.5 + rf) >= rf - eps:
                return True
        return False
    return ins

def hele_shaw(n, rf, Ld=3.0, Lu=3.0, stem=2.0):
    h = 1.0 / n; ins = make_inside(rf, Ld, Lu, stem)
    xs = np.arange(-Ld + h / 2, 0.5, h); ys = np.arange(-0.5 - stem + h / 2, Lu, h)
    fl = np.array([[ins(x, y, 0.0) for y in ys] for x in xs])
    idx = -np.ones(fl.shape, int); idx[fl] = np.arange(fl.sum()); nf = fl.sum()
    I, J, V = [], [], []; b = np.zeros(nf); nx, ny = fl.shape
    for i in range(nx):
        for j in range(ny):
            if not fl[i, j]: continue
            r = idx[i, j]; diag = 0.0
            for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                ii, jj = i + di, j + dj
                if 0 <= ii < nx and 0 <= jj < ny and fl[ii, jj]:
                    I.append(r); J.append(idx[ii, jj]); V.append(1.0); diag -= 1.0
                elif di == -1 and ii < 0:
                    diag -= 2.0; b[r] -= 2.0
                elif dj == 1 and jj >= ny and abs(xs[i]) < 0.5:
                    diag -= 2.0
            I.append(r); J.append(r); V.append(diag)
    p = spl.spsolve(sp.csr_matrix((V, (I, J)), shape=(nf, nf)).tocsc(), b)
    P = np.full(fl.shape, np.nan); P[fl] = p
    q = sum((1.0 - P[0, j]) / (h / 2) * h for j in range(ny) if fl[0, j])
    return 1.0 / q - (Ld + Lu)

def deep(n, rf, Ld=3.0, Lu=3.0, stem=2.0):
    h = 1.0 / n; ins = make_inside(rf, Ld, Lu, stem)
    xs = np.round(np.arange(-Ld, 0.5 + h / 2, h), 12); ys = np.round(np.arange(-0.5 - stem, Lu + h / 2, h), 12)
    nx, ny = len(xs), len(ys)
    fl = np.array([[ins(x, y) for y in ys] for x in xs])
    isfl = lambda i, j: 0 <= i < nx and 0 <= j < ny and fl[i, j]
    interior = np.array([[fl[i, j] and all(isfl(i + a, j + c) for a in (-1, 0, 1) for c in (-1, 0, 1)) for j in range(ny)] for i in range(nx)])
    inlet = np.zeros_like(fl); inlet[0, :] = fl[0, :]
    outlet = np.zeros_like(fl); outlet[:, -1] = fl[:, -1]
    prof = lambda s: 0.5 + 1.5 * s - 2.0 * s ** 3; omg = lambda s: 12.0 * s
    # wall psi: 1 on the chain drop-top -> fillet -> out-line-left; 0 elsewhere (lower/right chain incl. the bottom fillet)
    def wall_psi(x, y):
        if y >= 0.5 - 1e-9 and x <= -0.5 + 1e-9:
            return 1.0
        if abs(x + 0.5) < 1e-9 and y >= 0.5 - 1e-9:
            return 1.0
        if rf > 0 and x <= -0.5 + 1e-9 and y > 0.0:      # top-left fillet arc
            return 1.0
        return 0.0
    def psi_known(i, j):
        if inlet[i, j]: return prof(ys[j])
        if outlet[i, j]: return prof(-xs[i])
        return wall_psi(xs[i], ys[j])
    uid = {}
    for i in range(nx):
        for j in range(ny):
            if interior[i, j]: uid[('p', i, j)] = len(uid)
    for i in range(nx):
        for j in range(ny):
            if fl[i, j]: uid[('w', i, j)] = len(uid)
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
                I.append(r); J.append(r); V.append(-4 / h ** 2); I.append(r); J.append(uid[('w', i, j)]); V.append(1.0)
                r2 = uid[('w', i, j)]
                for (a, c) in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                    I.append(r2); J.append(uid[('w', i + a, j + c)]); V.append(1 / h ** 2)
                I.append(r2); J.append(r2); V.append(-4 / h ** 2)
            elif fl[i, j]:
                r2 = uid[('w', i, j)]
                if inlet[i, j] or outlet[i, j]:
                    s = ys[j] if inlet[i, j] else -xs[i]
                    I.append(r2); J.append(r2); V.append(1.0); b[r2] = omg(s); continue
                nb = [(a, c) for (a, c) in ((1, 0), (-1, 0), (0, 1), (0, -1)) if 0 <= i + a < nx and 0 <= j + c < ny and interior[i + a, j + c]]
                I.append(r2); J.append(r2); V.append(1.0)
                if not nb: continue
                pw = psi_known(i, j)
                for (a, c) in nb:
                    I.append(r2); J.append(uid[('p', i + a, j + c)]); V.append(2.0 / h ** 2 / len(nb)); b[r2] += 2.0 * pw / h ** 2 / len(nb)
    x = spl.spsolve(sp.csr_matrix((V, (I, J)), shape=(NU, NU)).tocsc(), b)
    W = np.full((nx, ny), np.nan)
    for (t, i, j), u in uid.items():
        if t == 'w': W[i, j] = x[u]
    j0 = int(np.argmin(np.abs(ys))); i0 = int(np.argmin(np.abs(xs)))
    seg1 = -np.trapezoid(((W[:, j0 + 1] - W[:, j0 - 1]) / (2 * h))[: i0 + 1], xs[: i0 + 1])
    seg2 = np.trapezoid(((W[i0 + 1, :] - W[i0 - 1, :]) / (2 * h))[j0:], ys[j0:])
    return -(seg1 + seg2) / 12.0 - (Ld + Lu)

if __name__ == '__main__':
    for rf in (0.0, 0.18, 0.26):
        hs = [hele_shaw(n, rf) for n in (40, 80)]
        dp = [deep(n, rf) for n in (40, 80)]
        print('r_f/w %.2f: shallow c_T %.4f (n40 %.4f), deep c_T %.4f (n40 %.4f; first-order extrapolation %.4f)'
              % (rf, hs[1], hs[0], dp[1], dp[0], 2 * dp[1] - dp[0]), flush=True)
