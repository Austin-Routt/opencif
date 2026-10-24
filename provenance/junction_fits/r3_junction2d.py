"""R3: how much does a FULLY DEVELOPED curved-duct correction leave out at the straight-arc junctions?
2-D solves in curvilinear (s, n) coordinates of a channel  straight Ls -> arc (radius rm, angle th) -> straight Ls,
in the two limits that bracket a rectangular duct:
  deep  (planar Stokes, h >> w): streamfunction-vorticity, dissipation  Phi = mu int omega^2 dA = dP Q
  shallow (Hele-Shaw, h << w): Laplace for p, no-flux walls
Excess resistance expressed as an equivalent straight length; compared with the developed correction
pi*rm*(1/ratio - 1) (deep: planar annulus; shallow: rm ln(r2/r1)/w)."""
import math, sys
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl

def geom(w, rm, th, Ls, d):
    La = rm * th
    Lt = 2 * Ls + La
    ns = int(round(Lt / d)) + 1
    nn = int(round(w / d)) + 1
    s = np.linspace(0, Lt, ns); n = np.linspace(-w / 2, w / 2, nn)
    ds, dn = s[1] - s[0], n[1] - n[0]
    kap = lambda ss: np.where((ss > Ls) & (ss < Ls + La), 1.0 / rm, 0.0)
    return s, n, ds, dn, kap, Lt, La

def lap_matrix(s, n, ds, dn, kap):
    ns, nn = len(s), len(n)
    idx = lambda i, j: i * nn + j
    I, J, V = [], [], []
    rows_int = []
    for i in range(1, ns - 1):
        kp, km, k0 = kap(s[i] + ds / 2), kap(s[i] - ds / 2), kap(s[i])
        for j in range(1, nn - 1):
            h0 = 1 - k0 * n[j]
            hp_s, hm_s = 1 - kp * n[j], 1 - km * n[j]
            hp_n, hm_n = 1 - k0 * (n[j] + dn / 2), 1 - k0 * (n[j] - dn / 2)
            c_ip, c_im = 1 / (hp_s * ds * ds) / h0, 1 / (hm_s * ds * ds) / h0
            c_jp, c_jm = hp_n / (dn * dn) / h0, hm_n / (dn * dn) / h0
            r = idx(i, j)
            I += [r] * 5; J += [r, idx(i + 1, j), idx(i - 1, j), idx(i, j + 1), idx(i, j - 1)]
            V += [-(c_ip + c_im + c_jp + c_jm), c_ip, c_im, c_jp, c_jm]
    return I, J, V

def deep(w, rm, th, Ls, d):
    s, n, ds, dn, kap, Lt, La = geom(w, rm, th, Ls, d)
    ns, nn = len(s), len(n)
    N = ns * nn
    idx = lambda i, j: i * nn + j
    Q = 1.0
    # Poiseuille in 2-D: u = 6Q/w (1/4 - (n/w)^2)... psi(n) = Q (1/2 + 3/2 eta - 2 eta^3), eta = n/w
    eta = n / w
    psi_p = Q * (0.5 + 1.5 * eta - 2 * eta ** 3)
    omg_p = -(Q * (-12 * eta) / w ** 2)          # omega = -psi_nn ; psi_nn = Q(-12 eta)/w^2
    I, J, V = lap_matrix(s, n, ds, dn, kap)
    # unknown vector [psi (N), omega (N)]
    II, JJ, VV, b = [], [], [], np.zeros(2 * N)
    # interior: L psi + omega = 0  (row r), L omega = 0 (row N + r)
    II += I; JJ += J; VV += V
    II += [N + r for r in I]; JJ += [N + c for c in J]; VV += V
    for i in range(1, ns - 1):
        for j in range(1, nn - 1):
            r = idx(i, j); II.append(r); JJ.append(N + r); VV.append(1.0)
    def dirichlet(r, val):
        II.append(r); JJ.append(r); VV.append(1.0); b[r] = val
    for i in range(ns):
        for j in (0, nn - 1):
            r = idx(i, j)
            dirichlet(r, 0.0 if j == 0 else Q)
            if i in (0, ns - 1):
                dirichlet(N + r, omg_p[j])
            else:
                # Jensen wall vorticity: omega_w = -(8 psi_1 - psi_2 - 7 psi_w)/(2 dn^2) (wall at j, interior j+-1, j+-2)
                jj = 1 if j == 0 else nn - 2; j2 = 2 if j == 0 else nn - 3
                II += [N + r, N + r, N + r, N + r]; JJ += [N + r, idx(i, jj), idx(i, j2), r]
                VV += [1.0, 8 / (2 * dn * dn), -1 / (2 * dn * dn), -7 / (2 * dn * dn)]
                # (wall curvature term: at the wall psi_s = 0, psi_n = 0 -> omega = -psi_nn exactly)
    for i in (0, ns - 1):
        for j in range(1, nn - 1):
            r = idx(i, j)
            dirichlet(r, psi_p[j]); dirichlet(N + r, omg_p[j])
    A = sp.csr_matrix((VV, (II, JJ)), shape=(2 * N, 2 * N))
    x = spl.spsolve(A.tocsc(), b)
    om = x[N:].reshape(ns, nn)
    S_, Nn = np.meshgrid(s, n, indexing='ij')
    hs = 1 - kap(S_) * Nn
    integrand = om ** 2 * hs
    Phi = np.trapezoid(np.trapezoid(integrand, n, axis=1), s)   # mu = 1
    dP = Phi / Q ** 2
    return dP, Lt, La

def shallow(w, rm, th, Ls, d):
    s, n, ds, dn, kap, Lt, La = geom(w, rm, th, Ls, d)
    ns, nn = len(s), len(n)
    N = ns * nn
    idx = lambda i, j: i * nn + j
    I, J, V = lap_matrix(s, n, ds, dn, kap)
    II, JJ, VV, b = list(I), list(J), list(V), np.zeros(N)
    # walls: no flux, second-order one-sided: -3p0 + 4p1 - p2 = 0
    for i in range(1, ns - 1):
        for j, j1, j2 in ((0, 1, 2), (nn - 1, nn - 2, nn - 3)):
            r = idx(i, j); II += [r, r, r]; JJ += [r, idx(i, j1), idx(i, j2)]; VV += [-3.0, 4.0, -1.0]
    for i, val in ((0, 1.0), (ns - 1, 0.0)):
        for j in range(nn):
            r = idx(i, j); II.append(r); JJ.append(r); VV.append(1.0); b[r] = val
    A = sp.csr_matrix((VV, (II, JJ)), shape=(N, N))
    p = spl.spsolve(A.tocsc(), b).reshape(ns, nn)
    # flux through the mid-section of the first straight (h_s = 1 there): q = -int dp/ds dn  (conductance h^3/12 = 1)
    im = int(round((Ls / 2) / ds))
    q = -np.trapezoid((p[im + 1] - p[im - 1]) / (2 * ds), n)
    return 1.0 / q, Lt, La

def report(lab, w, rm, th, Ls_w=4.0, ds_list=(40, 80)):
    Ls = Ls_w * w
    r1, r2 = rm - w / 2, rm + w / 2
    # developed predictions
    c = -0.5
    Mx = np.array([[r1, 1 / r1], [r2, 1 / r2]])
    A_, B_ = np.linalg.solve(Mx, -np.array([c * r1 * math.log(r1), c * r2 * math.log(r2)]))
    F_ = lambda r: c * (r * r / 2 * math.log(r) - r * r / 4) + A_ * r * r / 2 + B_ * math.log(r)
    ratio_deep = (F_(r2) - F_(r1)) / ((1 / rm) * w ** 3 / 12.0)
    ratio_sh = rm * math.log(r2 / r1) / w
    for kind, fn, ratio, rstraight in (('deep', deep, ratio_deep, 12.0 / w ** 3), ('shallow', shallow, ratio_sh, 1.0 / w)):
        vals = []
        for nres in ds_list:
            dP, Lt, La = fn(w, rm, th, Ls, w / nres)
            vals.append(dP / rstraight - (2 * Ls + La))        # excess equivalent length over the plain 1-D rule
        ex = vals[-1] + (vals[-1] - vals[-2]) / 3.0              # Richardson (2nd order)
        dev = La * (1 / ratio - 1)                               # developed (fully developed curved-duct) excess
        print('%-34s %-7s w %.1f rm %.1f th %.0f deg: 1-D excess, 2-D solve %+.4f w (res %s: %s), developed %+.4f w, '
              'junction part %+.4f w = %+.0f%% of the developed correction'
              % (lab, kind, w, rm, math.degrees(th), ex / w, ds_list, ['%+.4f' % (v / w) for v in vals], dev / w,
                 (ex - dev) / w, 100 * (ex - dev) / dev))
        sys.stdout.flush()

if __name__ == '__main__':
    report('straight check (rm huge)', 28.0, 1e6, 1e-4, ds_list=(20, 40))
    report('Dinh hairpin U-turn', 28.0, 27.5, math.pi)
    report('Dinh hairpin quarter arc', 28.0, 27.5, math.pi / 2)
    report('Xia hairpin U-turn', 22.8, 31.9, math.pi)
    report('Gifford loop 1 centre', 100.0, 481.7, math.pi, Ls_w=3.0, ds_list=(20, 40))
