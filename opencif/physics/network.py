"""opencif.physics.network: the symmetric 3-rail Kirchhoff network of a CIF device.

The network has a center rail and one side rail, the other side following by symmetry, with gap resistances that are
set row by row. A unit inlet flow enters center node 1, and both outlets are held at pressure 0. The unknowns are
P_c(1..N) and P_s(1..N), and the solve returns the gap flow of each row, the emergent fraction
f(i) = Q_gap(i)/Q_c,in(i), and F."""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl


def solve(Rc, Rs, Rg, Rc_out, Rs_out):
    """Solve the network for the resistances below and return its pressures and flows as a dict (Pc, Ps, Qg, Qc_in,
    f, F, Qs_out).

        Rc[i]            center segment i -> i+1 (i = 0..N-2)
        Rs[i]            side segment i -> i+1
        Rg[i]            gap i (i = 0..N-1)
        Rc_out, Rs_out   from node N to the outlets (per side)"""
    N = len(Rg)
    Rc = np.asarray(Rc, float); Rs = np.asarray(Rs, float); Rg = np.asarray(Rg, float)
    n = 2 * N
    I, J, V = [], [], []
    b = np.zeros(n)

    def add(i, j, v):
        I.append(i); J.append(j); V.append(v)
    for i in range(N):
        c, s = i, N + i
        gc = 1.0 / Rg[i]
        # The center node balances the inflow from upstream (the inlet for i = 0) against the outflow downstream and
        # the two gap flows.
        if i == 0:
            b[c] += 1.0
        else:
            g = 1.0 / Rc[i - 1]; add(c, c, g); add(c, c - 1, -g)
        g = 1.0 / (Rc[i] if i < N - 1 else Rc_out); add(c, c, g)
        if i < N - 1:
            add(c, c + 1, -g)
        add(c, c, 2 * gc); add(c, s, -2 * gc)
        # The side node balances the gap inflow and the upstream side inflow against the downstream side outflow.
        add(s, s, gc); add(s, c, -gc)
        if i > 0:
            g = 1.0 / Rs[i - 1]; add(s, s, g); add(s, s - 1, -g)
        g = 1.0 / (Rs[i] if i < N - 1 else Rs_out); add(s, s, g)
        if i < N - 1:
            add(s, s + 1, -g)
    A = sp.csr_matrix((V, (I, J)), shape=(n, n))
    x = spl.spsolve(A.tocsc(), b)
    Pc, Ps = x[:N], x[N:]
    Qg = (Pc - Ps) / Rg
    Qc_in = np.concatenate([[1.0], (Pc[:-1] - Pc[1:]) / Rc[:N - 1]])
    Qc_out_last = Pc[-1] / Rc_out
    f = Qg / Qc_in
    F = 1.0 - Qc_out_last
    return dict(Pc=Pc, Ps=Ps, Qg=Qg, Qc_in=Qc_in, f=f, F=F, Qs_out=Ps[-1] / Rs_out)
