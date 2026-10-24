"""test_bend_fd.py: the curved-duct correction against an independent finite-difference solve.

The suite compares the curved-duct correction (schedule.bend_ratio, a Fourier-Bessel series) with an independent
second-order finite-difference solve of the same azimuthal Stokes problem, u_rr + u_r/r - u/r^2 + u_zz = -1/r on the
duct section with u = 0 on its walls. The finite-difference flux is divided by that of the straight duct with the same
centerline gradient on the same grid, and the ratio is extrapolated from two grids with a Richardson step. The
governing equation itself is checked through its planar-annulus, Hele-Shaw, and straight limits in tests/test_units.py.
Run it with `python3 tests/test_bend_fd.py` (about a minute)."""
import os, sys
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from opencif import pipeline as _PL_CODE                # noqa: E402  (every suite's log names the code it ran on)
print('opencif %s (code %s)' % (_PL_CODE.CODE['version'], _PL_CODE.CODE['sha256']))
from opencif.physics import schedule as S

fails = []
def check(name, cond, info=''):
    print('%-74s %s %s' % (name, 'ok ' if cond else 'FAIL', info))
    if not cond:
        fails.append(name)


def fd_flux(r1, r2, h, n, curved=True):
    """Return the flux through the section r1 < r < r2, 0 < z < h, solved by second-order finite differences on n
    radial points, for the curved (azimuthal) problem or for the straight duct with the same centerline gradient."""
    nr, nz = n, max(8, int(round(n * h / (r2 - r1))))
    dr, dz = (r2 - r1) / (nr + 1), h / (nz + 1)
    r = r1 + dr * np.arange(1, nr + 1)
    rm = 0.5 * (r1 + r2)
    D2r = sp.diags([1.0, -2.0, 1.0], [-1, 0, 1], shape=(nr, nr)) / dr ** 2
    D1r = sp.diags([-1.0, 1.0], [-1, 1], shape=(nr, nr)) / (2 * dr)
    D2z = sp.diags([1.0, -2.0, 1.0], [-1, 0, 1], shape=(nz, nz)) / dz ** 2
    if curved:
        Lr = D2r + sp.diags(1.0 / r) @ D1r - sp.diags(1.0 / r ** 2)
        rhs = np.kron(1.0 / r, np.ones(nz))
    else:
        Lr = D2r
        rhs = np.full(nr * nz, 1.0 / rm)
    A = sp.kron(Lr, sp.identity(nz)) + sp.kron(sp.identity(nr), D2z)
    return spl.spsolve(A.tocsc(), -rhs).sum() * dr * dz


# The three bends are a deep head hairpin (28 um wide, 140 um deep, inner radius 30 um), a shallow wide loop, and a
# tight, nearly square bend.
for r1, r2, h in ((30.0, 58.0, 140.0), (100.0, 128.0, 28.0), (20.0, 60.0, 40.0)):
    q = [fd_flux(r1, r2, h, n, True) / fd_flux(r1, r2, h, n, False) for n in (80, 160)]
    fd = q[1] + (q[1] - q[0]) / 3.0
    ser = S.bend_ratio(r1, r2, h)
    check('bend_ratio = finite differences (r %g-%g um, h %g um): %.7f vs %.7f' % (r1, r2, h, ser, fd),
          abs(fd / ser - 1) < 2e-6, '%.1e' % (fd / ser - 1))
print('ALL OK' if not fails else 'FAILED: %s' % fails)
sys.exit(1 if fails else 0)
