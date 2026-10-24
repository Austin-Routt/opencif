"""test_auditor.py: the physics auditor (opencif/checks/auditor.py).

The suite has three parts.
1. The port against the original. CIF_v2.py's own physics auditor (its embedded ladder_physics module, `audit()`) is
   run on its three demonstration devices, and the port's functions, given the same measured geometry and the same
   kernel, must return the same numbers (the network solution, the terminal window, the bulk, P3, P4, P5, P6, P7, and
   P8) and the same verdicts. The original's tolerances are read from its frozen source text.
2. The port's network solver against the build's (measure.solve3) on random ladders.
3. A drawn opencif device. The photomask file read back equals the drawn polygon, and the audit passes. Then faults
   are seeded one at a time (a vertex of the file moved, a kernel corrupted, a gap widened, a side channel narrowed,
   a post thickened, a post moved along the flow, a side segment 2% narrower, another device's photomask, the flow
   raised, the DRC or the corners failing, the target F raised, an aspect-ratio limit below the drawing, and an
   intended network away from the drawn one for P4's reference), and each must be caught by the check meant to catch
   it.
Run it with `python3 tests/test_auditor.py` (about two minutes)."""
import copy
import os
import re
import sys
import tempfile

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, 'reference'))
from opencif import pipeline as PL                        # noqa: E402
print('opencif %s (code %s)' % (PL.CODE['version'], PL.CODE['sha256']))
from opencif import spec as SP, build as BU               # noqa: E402
from opencif.checks import auditor as AU, measure as M    # noqa: E402
from opencif.physics import kernels as K, schedule as S   # noqa: E402
from opencif.geometry import dxf as io_dxf                # noqa: E402
import CIF_v2 as C2                                       # noqa: E402  the frozen design code

fails = []


def check(name, cond, info=''):
    print('%-86s %s %s' % (name, 'ok ' if cond else 'FAIL', info))
    if not cond:
        fails.append(name)


def num(txt, pat):
    """Return the first number after the pattern pat in one of CIF_v2's detail strings, or None if there is none."""
    m_ = re.search(pat + r'\s*=?\s*([-+0-9.eE]+)', txt)
    return float(m_.group(1)) if m_ else None


# ------------------------------------------------------------------------------------- 1. the port against the original
pc = C2.physics_core
src = C2._PHYSICS_SRC
import dataclasses                                         # noqa: E402
spec0 = type('Defaults', (), {f.name: f.default for f in dataclasses.fields(pc.PhysicsSpec)})   # its default values
check('tolerances = CIF_v2 PhysicsSpec: tol_f_bulk %g, tol_delta %g, eps_net %g, Re_max %g, Q_in %g uL/min'
      % (spec0.tol_f_bulk, spec0.tol_delta, spec0.eps_net, spec0.Re_max, spec0.Q_in_uLmin),
      (AU.TOL_F_BULK, AU.TOL_DELTA, AU.EPS_NET, AU.STRICT_RE_MAX, AU.STRICT_Q_UL_MIN) ==
      (spec0.tol_f_bulk, spec0.tol_delta, spec0.eps_net, spec0.Re_max, spec0.Q_in_uLmin))
check('tolerances = CIF_v2 audit() source: P3 0.15, P6 1e-3 and 0.02, P6b max(8e-3, 4 f), P7 x3, P8 L_e <= Lambda/2',
      'worst <= %s' % AU.TOL_CLOSURE in src and 'spec.F_target - 1e-3' in src and 'abs(F_em - F_march) < 0.02' in src
      and 'max(8e-3, 4.0 * fbar)' in src and 'mismatch=3.0' in src and 'Le = 0.05 * Re * Dh' in src
      and 'Le <= 0.5 * geom["Lambda"]' in src and AU.TOL_F_TARGET == 1e-3 and AU.TOL_F_MARCH == 0.02
      and AU.TOL_F6_MIN == 8e-3 and AU.PORT_MISMATCH == 3.0)
check('the design envelope = CIF_v2 Tolerances (M14): Re_max 10, Re_stokes 1e-2',
      (AU.RE_MAX, AU.RE_STOKES) == (C2.Tolerances().Re_max, C2.Tolerances().Re_stokes) == (10.0, 1e-2))


def half_ladder(geom, sol, mismatch=1.0):
    """Return CIF_v2's own network (ladder_physics.solve_network) as an edge list. The list holds rung k, the trunk
    doubled (the one-rail half model), the excursion, the run-outs, the ground tie (1.0), and the common outlet
    (1e-9)."""
    U = geom['U']
    Rc, Rr, Re = sol['Rc'], sol['R_rung'], list(sol['R_exc'])
    Re[U - 1] *= mismatch
    E = []
    for k in range(U):
        E.append((2 * k, 2 * k + 1, Rr[k]))
        if k < U - 1:
            E.append((2 * k, 2 * k + 2, 2.0 * Rc[k]))
            E.append((2 * k + 1, 2 * k + 3, Re[k]))
    Po, Jo = 2 * U, 2 * U + 1
    E += [(2 * (U - 1), Po, 2.0 * Rc[U - 1]), (2 * (U - 1) + 1, Jo, Re[U - 1]), (Po, -1, 1.0), (Jo, Po, 1e-9)]
    inj = np.zeros(2 * U + 2)
    inj[0] = 0.5
    return 2 * U + 2, E, inj


def half_ladder_f(x, sol, U):
    """Return CIF_v2's f_em and F_em from the node pressures x."""
    Rc, Rr = sol['Rc'], sol['R_rung']
    Qh = [0.5] + [(x[2 * (k - 1)] - x[2 * k]) / (2.0 * Rc[k - 1]) for k in range(1, U)]
    f = np.array([float((x[2 * k] - x[2 * k + 1]) / Rr[k] / (2.0 * Qh[k])) for k in range(U)])
    return f, float(1.0 - 2.0 * (x[2 * (U - 1)] - x[2 * U]) / (2.0 * Rc[U - 1]))


def solve_long_double(n, E, inj):
    """Solve the same network in 80-bit extended precision (numpy longdouble, banded Gaussian elimination). The
    solution is the reference that both double-precision solvers are measured against."""
    bw = max(abs(i - j) for i, j, _ in E if j >= 0)
    A = np.zeros((n, 2 * bw + 1), dtype=np.longdouble)
    for i, j, R in E:
        g = np.longdouble(1) / np.longdouble(R)
        A[i, bw] += g
        if j >= 0:
            A[j, bw] += g
            A[i, bw + (j - i)] -= g
            A[j, bw + (i - j)] -= g
    b = np.array(inj, dtype=np.longdouble)
    for i in range(n):
        for jj in range(1, bw + 1):
            k_ = i + jj
            if k_ >= n:
                break
            fct = A[k_, bw - jj] / A[i, bw]
            A[k_, 0:2 * bw + 1 - jj] -= fct * A[i, jj:2 * bw + 1]
            b[k_] -= fct * b[i]
    x = np.zeros(n, dtype=np.longdouble)
    for i in range(n - 1, -1, -1):
        s_ = b[i]
        for jj in range(1, bw + 1):
            if i + jj < n:
                s_ -= A[i, bw + jj] * x[i + jj]
        x[i] = s_ / A[i, bw]
    return x


for mk in ('demo_spec_A', 'demo_spec_B', 'demo_spec_C'):
    spec = getattr(pc, mk)()
    pat, plan, aud = pc.compile_and_audit(spec)
    det = {c.name: c for c in aud.checks}
    geom = pc._measure_geometry(pat)
    sol = pc.solve_network(geom, spec.h)
    U = geom['U']
    n_, E_, inj_ = half_ladder(geom, sol)
    x_p, kres = AU.kirchhoff(n_, E_, inj_)
    f_p, F_p = half_ladder_f(x_p, sol, U)
    f_r, F_r = half_ladder_f(solve_long_double(n_, E_, inj_), sol, U)
    err_p = float(np.max(np.abs(f_p / f_r - 1.0)))
    err_c = float(np.max(np.abs(np.array(sol['f_em']) / f_r - 1.0)))
    # CIF_v2's network floats at the pressure that its 1-ohm ground tie gives the outlet (0.5, against drops of about
    # 1e-4 in its micrometer units), so double precision keeps about 1e-8 of a gap's flow in either solver.
    check('%s (%d units): the port\'s solver = CIF_v2\'s network to the precision either reaches' % (mk, U),
          err_p <= 3.0 * err_c + 1e-12 and err_p < 1e-6 and abs(F_p - F_r) < 1e-9 and abs(sol['F_em'] - F_r) < 1e-9,
          'vs extended precision: port %.1e, CIF_v2 %.1e' % (err_p, err_c))
    n3, E3, inj3 = half_ladder(geom, sol, mismatch=3.0)
    F3 = half_ladder_f(AU.kirchhoff(n3, E3, inj3)[0], sol, U)[1]
    check('%s: P7 (run-out x3) F as CIF_v2' % mk, abs(F3 - pc.solve_network(geom, spec.h, mismatch=3.0)['F_em']) < 1e-9)
    qmin, lam, n_term = AU.terminal_window(sol['q'], spec.eps_net)
    check('%s: terminal window lambda %.4f, n_term %d as CIF_v2' % (mk, lam, n_term),
          abs(lam - aud.lam_terminal) < 1e-12 and n_term == aud.n_term)
    bulk = AU.bulk_rows(U, n_term)
    check('%s: the bulk is CIF_v2\'s range(1, U - n_term)' % mk, list(np.where(bulk)[0]) == list(range(1, U - n_term)))
    e4 = float(np.max(np.abs(np.array(aud.f_em)[bulk] / np.array(aud.f_kernel)[bulk] - 1.0)))
    ref4 = num(det['P4_network_bulk'].detail, r'max bulk \|f_em/f - 1\|')
    check('%s: P4 %.3e as CIF_v2 (%s), same verdict' % (mk, e4, det['P4_network_bulk'].detail[:40]),
          abs(e4 - ref4) <= 0.006 * ref4 and (e4 <= AU.TOL_F_BULK) == det['P4_network_bulk'].ok)
    dstar = spec.D_target / spec.kappa
    dr = [AU.lamina_ratio(min(max(fe, 1e-9), 0.4999), geom['wc_rows'][2 * k], dstar,
                          x_of=lambda f, b=spec.h / geom['wc_rows'][2 * k]: pc.Flam_inv(f, b))
          for k, fe in enumerate(aud.f_em)]
    check('%s: P5 lamina ratios = CIF_v2\'s delta_ratio (all %d), same verdict' % (mk, U),
          np.array_equal(np.array(dr), np.array(aud.delta_ratio))
          and (max(dr) <= 1 + AU.TOL_DELTA) == det['P5_lamina_uniform'].ok)
    tgt = AU.pe_closure(aud.f_kernel, sol['Rc'])
    w3 = float(np.max(np.abs(np.array(sol['R_exc']) / tgt - 1.0)))
    ref3 = num(det['P3_closure_drawn'].detail, r'max \|R_exc/\(rho\*Rc\) - 1\|')
    check('%s: P3 PE-closure deviation %.3e as CIF_v2, same verdict' % (mk, w3),
          abs(w3 - ref3) <= 0.006 * ref3 and (w3 <= AU.TOL_CLOSURE) == det['P3_closure_drawn'].ok)
    Pp = float(np.prod([1 - 2 * f for f in aud.f_kernel]))
    ok6 = aud.F_em >= spec.F_target - AU.TOL_F_TARGET and abs(aud.F_em - (1 - Pp)) < AU.TOL_F_MARCH
    check('%s: P6 verdict as CIF_v2 (F_em %.4f)' % (mk, aud.F_em),
          ok6 == det['P6_F_target'].ok and abs(F_p - aud.F_em) < 1e-9)
    p8 = AU.strict_inlet_check(spec.Q_in_uLmin, geom['wc_rows'][0], spec.h, geom['Lambda'], spec.rho_SI, spec.mu_SI)
    ref8 = num(det['P8_stokes_gates'].detail, r'Re')
    check('%s: P8 inlet Re %.3f as CIF_v2 (%s), same verdict'
          % (mk, p8['Re_exact'], det['P8_stokes_gates'].detail[:22]),
          abs(p8['Re_exact'] - ref8) < 6e-4 and p8['passes'] == det['P8_stokes_gates'].ok)

# --------------------------------------------------------------------------------- 2. the solver against measure.solve3
rng = np.random.default_rng(20261022)
worst = 0.0
for N in (2, 7, 50, 400):
    for _ in range(3):
        Rc = rng.uniform(1.0, 2.0, N)
        RsA, RsB = rng.uniform(0.5, 3e3, N), rng.uniform(0.5, 3e3, N)
        RgA, RgB = rng.uniform(50.0, 500.0, N), rng.uniform(50.0, 500.0, N)
        a = AU.solve_ladder(Rc, RsA, RsB, RgA, RgB, 3.0, 9.0, 11.0)
        b = M.solve3(Rc[:-1], RsA[:-1], RsB[:-1], RgA, RgB, 3.0, 9.0, 11.0)
        worst = max(worst, float(np.max(np.abs(a['fA'] - b['fA']) / np.abs(b['fA']))),
                    float(np.max(np.abs(a['fB'] - b['fB']) / np.abs(b['fB']))), abs(a['F'] - b['F']),
                    a['kcl_residual_max'])
check('the port\'s solver (banded Cholesky) = measure.solve3 (sparse LU) on 12 random ladders', worst < 1e-9,
      '%.1e' % worst)

# ---------------------------------------------------------------------------------- 3. a drawn device and seeded faults
T = tempfile.mkdtemp(prefix='opencif_audit_test_')
d = SP.load(os.path.join(ROOT, 'examples', 'matrix', 'm15_linear_shear_nomeander_filter_teardrop.yaml'))
S.clear_caches()
spec, opts, sc, g, env_bad, notes = PL.prepare(d)
kind, n_legs, spec, sc = PL.fit_layout(spec, sc, opts, notes, lambda m_: None)
fluid, sc_f, info, _, _, _, _ = BU.plan_and_draw(spec, kind, n_legs, opts['leg_spacing'], opts['failsafe'],
                                                 lambda m_: None, 'm15', 0.0)
dxf = os.path.join(T, 'm15.dxf')
io_dxf.write_dxf(fluid, dxf)
pm, facts = AU.read_photomask(dxf)
check('the photomask file read back is the drawn polygon, vertex for vertex', AU.same_polygon(pm, fluid),
      '%d outline, %d islands, %d vertices' % (facts['fluid_outlines'], facts['islands'], facts['vertices']))


def blocks(poly, m):
    """Return the validator blocks that the build computes for the polygon poly and its measurement m (rounding, side
    continuity, components, outlet and loops, drawn parameters, and the DRC)."""
    V = dict(rounding=BU.rounding_check(poly, spec), side_continuity=BU.side_continuity(poly, sc_f, info),
             components_drawn=BU.components_drawn(poly, sc_f, info))
    BU.outlet_and_loops_measured(V, m, spec)
    kL = int(next(k for k in range(sc_f.NM, sc_f.N) if sc_f.kind[k] == 'L'))
    V['drawn_parameters'] = BU.drawn_parameters(m, sc_f, poly, info, kL)
    V['drc'] = BU.run_drc(dxf, T, 'm15', spec)
    return V


m0 = M.measure_device(pm, sc_f, info)
V0 = blocks(pm, m0)
inten, _, decl = BU.intended_and_declared(sc_f, info, m0, spec, S.KAPPA_G * spec.delta_star, opts['cut_tol'],
                                          opts['mode'], opts['seal'])
ARM = opts.get('AR_max', BU.AR_MAX_DEFAULT)
a0 = AU.audit(dxf, sc_f, info, m0, V0, decl, opts['mode'], ARM, fluid=fluid, nw=M.network(sc_f, m0))
n_posts = next(c for c in a0['checks'] if c['id'] == 'R6b')['value']['posts']
r0_0 = next(c for c in a0['checks'] if c['id'] == 'R0')['value']
check('the drawn device passes every check (P0-P10, R0, M14, R6b: %d posts; the two rulers within %s on %d segments)'
      % (n_posts, r0_0['max_rel_diff'], r0_0['segments']),
      a0['status'] == 'PASS' and all(c['ok'] for c in a0['checks']), AU.summary_line(a0))


def failed(a):
    """Return the sorted IDs of the checks that failed in the audit a."""
    return sorted(c['id'] for c in a['checks'] if not c['ok'])


# This block tests the reference of P4. Given the intended network, P4 compares the drawn fractions with it and reports
# the deviation from the kernel beside it, and on a design without declared departures the two references agree. An
# intended network 1% away from the drawn one fails P4 although the kernel deviation is unchanged, which shows that P4
# uses its reference.
a_i = AU.audit(dxf, sc_f, info, m0, V0, decl, opts['mode'], ARM, fluid=fluid, intended=inten)
p4i = next(c for c in a_i['checks'] if c['id'] == 'P4')['value']
p4k = next(c for c in a0['checks'] if c['id'] == 'P4')['value']
check('P4 with the intended network: the reference named, the kernel deviation reported and equal to the kernel run\'s',
      a_i['status'] == 'PASS' and p4i['reference'] == 'the intended network' and p4k['reference'] == 'the kernel'
      and p4i['kernel_bulk_max'] == p4k['bulk_max'],
      'bulk %s against the intended network, %s against the kernel' % (p4i['bulk_max'], p4k['bulk_max']))
off_net = dict(inten, fA=np.asarray(inten['fA']) * 1.01, fB=np.asarray(inten['fB']) * 1.01)
a = AU.audit(dxf, sc_f, info, m0, V0, decl, opts['mode'], ARM, fluid=fluid, intended=off_net)
check('fault: an intended network 1% away from the drawn one -> P4 (and nothing else)', failed(a) == ['P4'],
      str(failed(a)))


# 3a. One vertex of the file moved by 1 nm
import ezdxf                                               # noqa: E402
doc = ezdxf.readfile(dxf)
e = next(iter(doc.modelspace().query('LWPOLYLINE[layer=="ISLANDS"]')))
pts = list(e.get_points('xyseb'))
pts[3] = (pts[3][0] + 1e-3,) + tuple(pts[3][1:])
e.set_points(pts, format='xyseb')
dxf_b = os.path.join(T, 'm15_moved.dxf')
doc.saveas(dxf_b)
a = AU.audit(dxf_b, sc_f, info, m0, V0, decl, opts['mode'], ARM, fluid=fluid)
check('fault: one vertex of the file moved by 1 nm -> P0', failed(a) == ['P0'], str(failed(a)))
# 3b. A corrupted kernel
g_ok = K.g_cubic
K.g_cubic = lambda x: 3.0 * x ** 2 - 2.0 * x ** 3 + 1e-9
try:
    a = AU.audit(dxf, sc_f, info, m0, V0, decl, opts['mode'], ARM, fluid=fluid)
finally:
    K.g_cubic = g_ok
check('fault: the cubic kernel off by 1e-9 -> P1', 'P1' in failed(a), str(failed(a)))
# 3c. One gap widened by 0.6 um (row 400, side A), with the drawn polygon edited and measured again
from shapely.geometry import Point, Polygon, box           # noqa: E402
k = 400
leg = M.leg_of(info['legs'], k)
gx0 = info['gx'][leg['j']]
off = leg['p0'] if leg['p0'] is not None else 0.0
hw = float(sc_f.wc[k]) / 2.0
xl = float(sc_f.x[k]) - off
G = float(sc_f.wL[k])
p0, p1 = M.l2g(leg, gx0, xl + G / 2.0 - 0.01, hw + 0.5), M.l2g(leg, gx0, xl + G / 2.0 + 0.6, hw + spec.l_g - 0.5)
poly_g = fluid.union(box(min(p0[0], p1[0]), min(p0[1], p1[1]), max(p0[0], p1[0]), max(p0[1], p1[1])))
dxf_g = os.path.join(T, 'm15_gap.dxf')
io_dxf.write_dxf(poly_g, dxf_g)
pm_g, _ = AU.read_photomask(dxf_g)
m_g = M.measure_device(pm_g, sc_f, info)
a = AU.audit(dxf_g, sc_f, info, m_g, V0, decl, opts['mode'], ARM, fluid=poly_g)
p4 = next(c for c in a['checks'] if c['id'] == 'P4')['value']
check('fault: gap %d widened by 0.6 um -> P4 and P5 (the drawn f rises %.1f %%)' % (k + 1, 100 * p4['bulk_max']),
      {'P4', 'P5'} <= set(failed(a)) and p4['bulk_max_row'] == k + 1, str(failed(a)))
# 3d. One side segment narrowed to 60% (row 300, side B)
k = 300
leg = M.leg_of(info['legs'], k)
gx0 = info['gx'][leg['j']]
off = leg['p0'] if leg['p0'] is not None else 0.0
xa_, xb_ = float(sc_f.x[k]) + G / 2.0 + 2.0 - off, float(sc_f.x[k + 1]) - G / 2.0 - 2.0 - off
y0_ = float(max(sc_f.wc[k], sc_f.wc[k + 1])) / 2.0 + spec.l_g + 0.6 * float(sc_f.ws[k])
q0, q1 = M.l2g(leg, gx0, xa_, -y0_), M.l2g(leg, gx0, xb_, -(y0_ + float(sc_f.ws[k]) + 5.0))
poly_s = fluid.difference(box(min(q0[0], q1[0]), min(q0[1], q1[1]), max(q0[0], q1[0]), max(q0[1], q1[1])))
dxf_s = os.path.join(T, 'm15_side.dxf')
io_dxf.write_dxf(poly_s, dxf_s)
pm_s, _ = AU.read_photomask(dxf_s)
m_s = M.measure_device(pm_s, sc_f, info)
a = AU.audit(dxf_s, sc_f, info, m_s, V0, decl, opts['mode'], ARM, fluid=poly_s)
p3 = next(c for c in a['checks'] if c['id'] == 'P3')['value']
check('fault: side segment %d narrowed to 60 %% -> P3 (and P4) at that segment' % (k + 1),
      {'P3', 'P4'} <= set(failed(a)) and tuple(p3['at']) == (k + 1, 'B'), '%s %s' % (failed(a), p3))
# 3d'. One post 0.6 um thicker across the flow (the post after gap 501 on side A, with its outer face moved into the
# side channel), with the drawn polygon edited and measured again
k = 500
leg = M.leg_of(info['legs'], k)
gx0 = info['gx'][leg['j']]
off = leg['p0'] if leg['p0'] is not None else 0.0
xa_ = float(sc_f.x[k]) + float(sc_f.wL[k]) / 2.0 + 2.0 - off                 # 2 um inside the post's two faces
xb_ = float(sc_f.x[k + 1]) - float(sc_f.wL[k + 1]) / 2.0 - 2.0 - off
y0_ = float(np.interp(0.5 * (sc_f.x[k] + sc_f.x[k + 1]), sc_f.x, sc_f.wc)) / 2.0 + spec.l_g
q0, q1 = M.l2g(leg, gx0, xa_, y0_ - 1.0), M.l2g(leg, gx0, xb_, y0_ + 0.6)
poly_p = fluid.difference(box(min(q0[0], q1[0]), min(q0[1], q1[1]), max(q0[0], q1[0]), max(q0[1], q1[1])))
dxf_p = os.path.join(T, 'm15_post.dxf')
io_dxf.write_dxf(poly_p, dxf_p)
pm_p, _ = AU.read_photomask(dxf_p)
m_p = M.measure_device(pm_p, sc_f, info)
a = AU.audit(dxf_p, sc_f, info, m_p, V0, decl, opts['mode'], ARM, fluid=poly_p)
p9p = next(c for c in a['checks'] if c['id'] == 'R6b')['value']
check('fault: the post after gap %d 0.6 um thicker across the flow -> R6b at that post' % (k + 1),
      'R6b' in failed(a) and p9p['worst']['row'] == k + 1 and p9p['worst']['side'] == 'A'
      and abs(p9p['max_dev_across_um'] - 0.6) < 0.01 and p9p['max_dev_along_um'] < 1e-6,
      '%s %s' % (failed(a), p9p['worst']))
# 3d''. One post moved half a gap along the flow (the post after gap 601, side A). The generator's measurement reads
# the widths where the plan puts the parts and the lengths from the plan, whereas the auditor's ruler reads them from
# the posts.
from shapely import affinity                                # noqa: E402
k = 600
post = next(Polygon(r) for r in fluid.interiors
            if Polygon(r).contains(Point(0.5 * (sc_f.x[k] + sc_f.x[k + 1]), float(sc_f.wc[k]) / 2.0 + spec.l_g / 2.0)))
poly_s = fluid.union(post).difference(affinity.translate(post, xoff=0.5 * float(sc_f.wL[k])))
dxf_s = os.path.join(T, 'm15_shifted.dxf')
io_dxf.write_dxf(poly_s, dxf_s)
pm_s, _ = AU.read_photomask(dxf_s)
m_s = M.measure_device(pm_s, sc_f, info)
a = AU.audit(dxf_s, sc_f, info, m_s, V0, decl, opts['mode'], ARM, fluid=poly_s)
r0 = next(c for c in a['checks'] if c['id'] == 'R0')['value']
check('fault: the post after gap %d moved half a gap along the flow -> R0 (the rulers disagree) and P4/P5' % (k + 1),
      {'R0', 'P4', 'P5'} <= set(failed(a)) and r0['max_rel_diff'] > 0.05,
      '%s; R0 max %s at %s' % (failed(a), r0['max_rel_diff'], r0['by_kind']))
# 3d'''. One side segment 2% narrower over its middle half (row 301, side B), with steps at both ends. This is a real
# defect rather than a drawing that the generator makes. The network checks catch what it does to the cutoff, and R0
# reports how far the generator's measurement (which samples a wall at fixed points) and the ruler's (which integrates
# between the polygon's vertices) differ on a stepped wall.
k = 300
y_out = float(max(sc_f.wc[k], sc_f.wc[k + 1])) / 2.0 + spec.l_g + float(sc_f.ws[k])
poly_n = fluid.difference(box(float(sc_f.x[k]) + 12.0, -(y_out + 5.0), float(sc_f.x[k + 1]) - 12.0,
                              -(y_out - 0.02 * float(sc_f.ws[k]))))
dxf_n = os.path.join(T, 'm15_narrow2.dxf')
io_dxf.write_dxf(poly_n, dxf_n)
pm_n, _ = AU.read_photomask(dxf_n)
m_n = M.measure_device(pm_n, sc_f, info)
a = AU.audit(dxf_n, sc_f, info, m_n, V0, decl, opts['mode'], ARM, fluid=poly_n)
r0 = next(c for c in a['checks'] if c['id'] == 'R0')['value']
p4 = next(c for c in a['checks'] if c['id'] == 'P4')['value']
p5 = next(c for c in a['checks'] if c['id'] == 'P5')['value']
check('fault: side segment %d 2 %% narrower over half its length -> P4 and P5 at that row' % (k + 1),
      {'P4', 'P5'} <= set(failed(a)) and p4['bulk_max_row'] in (k + 1, k + 2),
      '%s; P4 %s at row %s; P5 max %s; R0 on RsB %s at row %s'
      % (failed(a), p4['bulk_max'], p4['bulk_max_row'], p5['max_delta_over_delta_star'],
         r0['by_kind']['RsB']['max_rel_diff'], r0['by_kind']['RsB']['row']))
# 3d''''. Another device's photomask under this build's name (the teaching module's straight ladder)
from opencif import teaching as TCH                         # noqa: E402
dxf_o = os.path.join(T, 'other_device.dxf')
TCH.write_photomask(TCH.draw(TCH.design(D_um=6.0, F=0.5, depth_um=80, gap_um=20.5, post_um=25.3, wall_um=59.5,
                                        centre_um=100)), dxf_o)
a = AU.audit(dxf_o, sc_f, info, m0, V0, decl, opts['mode'], ARM, fluid=fluid)
check('fault: another device\'s photomask under this build\'s name -> P0 (and the ruler finds none of its posts)',
      {'P0', 'R0', 'R6b'} <= set(failed(a)), str(failed(a)))
# 3e. The flow raised three-fold, which leaves the envelope (Re > 10) without a FAIL
sp3 = copy.copy(sc_f.spec)
sp3.Q_in_uL_min = 3.0 * float(sp3.Q_in_uL_min)
sc3 = copy.copy(sc_f)
sc3.spec = sp3
a = AU.audit(dxf, sc3, info, m0, V0, decl, opts['mode'], ARM, fluid=fluid)
check('fault: the flow x3 -> M14 (envelope): PASS_OUTSIDE_ENVELOPE, no hard failure',
      failed(a) == ['M14'] and a['status'] == 'PASS_OUTSIDE_ENVELOPE', '%s %s' % (failed(a), a['status']))
# 3f to 3h. The DRC, the corner radius, and the target F
V1 = copy.deepcopy(V0)
V1['drc']['n_lost_open'] = 1
a = AU.audit(dxf, sc_f, info, m0, V1, decl, opts['mode'], ARM, fluid=fluid)
check('fault: the DRC finds a lost feature -> P9', failed(a) == ['P9'], str(failed(a)))
V1 = copy.deepcopy(V0)
V1['rounding']['ok'] = False
a = AU.audit(dxf, sc_f, info, m0, V1, decl, opts['mode'], ARM, fluid=fluid)
check('fault: a corner sharper than the criterion -> P10', failed(a) == ['P10'], str(failed(a)))
spF = copy.copy(sc_f.spec)
spF.F_target = float(spF.F_target) + 0.01
scF = copy.copy(sc_f)
scF.spec = spF
a = AU.audit(dxf, scF, info, m0, V0, decl, opts['mode'], ARM, fluid=fluid)
check('fault: the target F 0.01 above the drawn F -> P6', failed(a) == ['P6'], str(failed(a)))
w_ar = next(c for c in a0['checks'] if c['id'] == 'P9-AR')['value']['width_um']
a = AU.audit(dxf, sc_f, info, m0, V0, decl, opts['mode'], 0.5 * float(spec.h) / w_ar, fluid=fluid)
check('fault: an aspect-ratio limit below the drawing -> P9-AR (a pass condition in strict mode)',
      failed(a) == ['P9-AR'] and a['status'] == 'FAIL', str(failed(a)))
print('ALL OK' if not fails else 'FAILED: %s' % fails)
sys.exit(1 if fails else 0)
