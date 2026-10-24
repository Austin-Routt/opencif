"""opencif.checks.auditor: the physics auditor, which reads the photomask back from its DXF file and checks the
physics of Appendices A and B on the device that the file contains.

The auditor comes from the design code (reference/CIF_v2.py v2.1.0), which carries a physics auditor, the `audit()` of
its embedded ladder_physics module, with ten checks, P1 to P10. Appendix B relies on it twice. It states that "the
physics auditor separately checks the drawn device's inlet against its own, stricter channel limit" (App. B
:4349-4350), and that the targets-to-photomask tool and the validity envelope are "warranted by their acceptance suite
and the physics auditor" (:4375-4377). CIF_v2's auditor reads CIF_v2's own meander patterns (their declared stroke
widths and trunk outline), and it cannot read a DXF file or opencif's parts. This module is that auditor ported to
opencif and evaluated on the photomask file rather than on a pattern in memory. Every check keeps its purpose and its
threshold. Where what a DXF can show forces another formula (the cubic pair's anchors in P1, App. A's outlet identity
in P6b, and measurement on the polygon in P9 and P10), or App. B another role (P8 is reported), docs/PHYSICS_AUDIT.md
says how and why. The port adds three checks, P0 (the file is the drawing), R0 (an independent ruler), and R6b (the
design code's acceptance check of the posts), and it evaluates the design code's gate M14 on the drawn mask beside
P8. The suite tests/test_auditor.py runs the port and CIF_v2's own auditor side by side on CIF_v2's demonstration
devices, requires the same numbers for the network, the terminal window, P3 (its PE part), P4, P5, P7, and P8, and
seeds a fault for each check.

The auditor reads the DXF file back (layer CHANNELS minus layer ISLANDS), and the result must be the drawn polygon,
vertex for vertex (P0). The straight rows (most of a device's gaps) are measured twice on that polygon, once by the
generator's own measurement (opencif.checks.measure, which reads every part where the plan puts it) and once by the
auditor's own ruler (opencif.checks.ruler, which finds the posts on the file and reads every position, width, and
length from them). The two must agree (R0), and the auditor's network uses the ruler's resistances on those rows.
The head units, U-turns, and outlet are taken from the generator's measurement. The auditor assembles and solves its
own Kirchhoff network, with its own node numbering and solver, independent of measure.solve3, and applies the
following checks:

  P0   the file is the drawing; the DXF read back equals the drawn polygon exactly (new, because CIF_v2 has no file)
  R0   the auditor's own ruler, which finds every post of every ordinary straight row on the file and measures every
       straight row's center segment, side segment, and gap from the file alone within 1e-3 of the generator's
       measurement
  P1   the kernel anchors, which are the duct bracket B(1) = 0.4217310, the square duct's R w h^3/(mu L) = 28.454,
       the cubic's g(0.3) = 0.216 and g(1/2) = 1/2, its inverse, and kappa_g = 4.9138 (the cubic pair)
       (App. A R1, R2; App. B C3b)
  P2   mass, through the plan's identity P(N) + 2[S(N) - sigma_0] = 1 (App. A, mass identity) and, on the drawn
       network, Kirchhoff's current law at every node and 1 - F = prod (1 - f_A - f_B) (App. A eq. Frec)
  P3   the closure realized by the drawn geometry, which holds every drawn side resistance against the design's exact
       closure recomputed from the drawn center and gap resistances, and against the PE closure rho R_c (App. A eq.
       closure), within CIF_v2's sanity bound 0.15
  P4   the network bulk, |f_drawn / f_int - 1| <= 5e-3 on the bulk rows (CIF_v2 tol_f_bulk), where f_int is the
       intended network's fraction (build.intended_network, solved on the design's own resistances). Without declared
       departures f_int is the plan's f(i) = g(delta*/w_c(i)) at the design width; with them it also carries their
       effect on the rows downstream. The deviation from the kernel at each gap's drawn center width, CIF_v2's own
       reference, is reported beside it (kernel_bulk_max), and it is the reference when no intended network is given.
       The bulk leaves out (i) row 1, (ii) CIF_v2's terminal window n_term = max(2, ceil(lambda ln(1/eps_net)) + 1),
       where lambda = 1/ln z+, z+ = [(2 + q) + sqrt((2 + q)^2 - 4)]/2, q = min (2 R_c + R_s)/R_g, and
       eps_net = 1e-3, and (iii) a window on both sides of every declared departure, lambda ln(A/eps_net) rows long,
       where A is the largest deviation of its rows from the P4 reference, and never shorter than n_term
       (departure_window). Because the intended network carries the departure, A measures the drawing's error there,
       and the window is n_term whenever A <= 1
  P5   lamina uniformity, the cutoff guarantee, delta_i/delta* = x(f_i) w_c(i)/delta* <= 1.02 at every undeclared
       gap on both sides (CIF_v2 tol_delta)
  P6   the filtrate fraction, F_drawn >= F_target - 1e-3 and |F_drawn - F_march| < 0.02 (CIF_v2)
  P6b  a closed form for F, App. A's march-outlet identity F = 2 R_c/(R_s + 2 R_c) at the drawn outlet (with both
       sides, F = (s_A + s_B)/(1 + s_A + s_B) and s = R_c,out/R_s,out), against the drawn network's F within CIF_v2's
       max(8e-3, 4 f_mean). The identity's own premise (PE at the outlet row) is not exact (App. A), and its
       residual is reported
  P7   port matching, with both side outlets' resistance x3 and the F it gives (reported, as in CIF_v2)
  P8   the auditor's own, stricter channel limit (App. B :4349-4350; CIF_v2 P8, PhysicsSpec.Re_max = 1.0), which is
       the drawn inlet's Re <= 1 and the naive entrance length 0.05 Re D_h <= Lambda/2, at the design flow and at
       CIF_v2's default 5 uL/min. It is reported rather than a pass condition, because the design envelope is the
       design code's M14 (docs/LIMITATIONS.md 6)
  M14  the design code's Reynolds and Stokes gates (CIF_v2 M14, which App. B gives the design code) on the drawn mask
       at the design flow, which are the center channel's Re <= 10 at every row and Re_gap and Re_delta <= 1e-2 at
       every cutoff gap, from the drawn widths and the drawn network's flows (role 'envelope')
  P9   the fabrication floors on the drawn geometry, which are a clean DRC at the process floor, the narrowest drawn
       CIF channel >= F_min, and depth/width <= AR_max (CIF_v2 w_min = max(F_min, h/AR_max)). The aspect ratio,
       P9-AR, is a pass condition in strict mode and is reported in evaluate and reproduce modes, as pipeline.assess
       does
  R6b  the posts read back from the file (CIF_v2's acceptance check R6b, at every straight-row post rather than at
       one), whose extents across and along the flow are the wall thickness l_g and the post length g_p within
       0.5 um, with their ratio l_g/g_p within 0.05
  P10  rounded corners, with the smallest local radius of the drawn boundary >= the spec's criterion (r_min)

Each check has one of three roles. A failure of a 'hard' check fails the build (FAIL, exit 2) in every mode. A failure
of an 'envelope' check puts the drawn device outside the validity envelope, which gives exit 3 in evaluate and reproduce
modes, whereas strict mode refuses such designs before drawing. A 'report' check does not change the status. Declared
departures (in reproduce mode a depicted head, the papers' transition, a loop pull-in, or a floor start, and an allowed
seed) are the rows that build.declared_rows names before the drawing is measured. P3 to P5 and M14's gap terms are
evaluated on the other rows and report the declared ones apart.
"""
import hashlib
import math
import os
import time

import numpy as np

from opencif.checks import ruler as RU
from opencif.physics import kernels as K
from opencif.physics import schedule as S

# The tolerances of CIF_v2's ladder_physics PhysicsSpec and audit(), unchanged. The suite tests/test_auditor.py reads
# them from the frozen source and requires these values.
TOL_F_BULK = 5.0e-3          # P4  |f_em/f - 1| in the bulk
TOL_DELTA = 0.02             # P5  max delta_em/delta* - 1
EPS_NET = 1.0e-3             # P4  terminal-layer accounting tolerance
TOL_CLOSURE = 0.15           # P3  sanity bound on the drawn closure
TOL_F_TARGET = 1.0e-3        # P6  F_em >= F_target - 1e-3
TOL_F_MARCH = 0.02           # P6  |F_em - F_march| < 0.02
TOL_F6_MIN = 8.0e-3          # P6b max(8e-3, 4 f_mean)
STRICT_RE_MAX = 1.0          # P8  the auditor's own channel limit
STRICT_Q_UL_MIN = 5.0        # P8  its default operating point
PORT_MISMATCH = 3.0          # P7  the run-out mismatch factor
# the extents and ratio of CIF_v2's acceptance validator R6b (the DXF's posts are the process's pill)
TOL_POST_UM = 0.5            # |along - g_p| and |across - l_g| < 0.5 um
TOL_POST_RATIO = 0.05        # |across/along - l_g/g_p| < 0.05
# R0, the generator's measurement of a straight row against the auditor's own ruler (opencif.checks.ruler)
TOL_RULER = 1.0e-3           # |R_measure/R_ruler - 1| on every center segment, side segment, and gap
# M14, the channel and Stokes gates (CIF_v2 Tolerances, the design envelope that App. B gives the design code)
RE_MAX = float(K.Tolerances().Re_max)          # 10
RE_STOKES = float(K.Tolerances().Re_stokes)    # 1e-2
# P1's anchors, the published values that each kernel must reproduce
ANCHORS = dict(duct_bracket_1=(0.4217310, 5e-7), square_duct_R=(28.454154, 5e-6), g_0p3=(0.216, 1e-15),
               g_half=(0.5, 1e-15), kappa_g=(4.9138, 0.0))


# ---------------------------------------------------------------------------------------------------------------------
# reading the photomask
# ---------------------------------------------------------------------------------------------------------------------

def read_photomask(path):
    """Read the fluid region of an opencif photomask and return (geometry, facts about the file).

    Every closed polyline on layer CHANNELS is the outline of a fluid polygon, and every closed polyline on layer
    ISLANDS that follows it is a solid island (a post) inside it, because opencif.geometry.dxf writes them in that
    order."""
    import ezdxf
    from shapely.geometry import Polygon, MultiPolygon
    doc = ezdxf.readfile(path)
    polys, other_layers, n_bulge, n_open, n_vert = [], {}, 0, 0, 0
    for e in doc.modelspace():
        if e.dxftype() != 'LWPOLYLINE':
            other_layers[e.dxftype()] = other_layers.get(e.dxftype(), 0) + 1
            continue
        if not e.closed:
            n_open += 1
            continue
        if any(abs(b[0]) > 0 for b in e.get_points('b')):
            n_bulge += 1                                 # an arc segment, which opencif's writer never draws
        pts = [(float(x), float(y)) for x, y in e.get_points('xy')]
        n_vert += len(pts)
        lay = e.dxf.layer.upper()
        if lay == 'CHANNELS':
            polys.append([pts, []])
        elif lay == 'ISLANDS':
            if not polys:
                raise ValueError('%s: an ISLANDS outline before any CHANNELS outline' % path)
            polys[-1][1].append(pts)
        else:
            other_layers['LWPOLYLINE on ' + lay] = other_layers.get('LWPOLYLINE on ' + lay, 0) + 1
    gs = [Polygon(ext + [ext[0]], [h + [h[0]] for h in holes]) for ext, holes in polys]
    geom = gs[0] if len(gs) == 1 else MultiPolygon(gs)
    facts = dict(file=os.path.basename(path), sha256=hashlib.sha256(open(path, 'rb').read()).hexdigest(),
                 fluid_outlines=len(polys), islands=int(sum(len(h) for _, h in polys)), vertices=int(n_vert),
                 arc_segments=n_bulge, open_polylines=n_open, other_entities=other_layers,
                 units=int(doc.header.get('$INSUNITS', 0)), valid=bool(geom.is_valid))
    return geom, facts


def same_polygon(a, b):
    """Return True when two polygons are identical, with the same rings in the same order and every vertex equal
    (no tolerance)."""
    return a.geom_type == b.geom_type and a.equals_exact(b, 0.0)


# ---------------------------------------------------------------------------------------------------------------------
# the Kirchhoff network
# ---------------------------------------------------------------------------------------------------------------------

def kirchhoff(n, edges, inject):
    """Return the node pressures of a resistor network from Kirchhoff's current law at every node.

    n       the number of nodes
    edges   a list of (i, j, R), each a resistor R between nodes i and j, where j = -1 ties node i to ground
            (pressure 0) through R
    inject  the flow injected at each node (length n)

    The conductance matrix is assembled in banded form, so the caller numbers the nodes such that every resistor joins
    nearby nodes. The matrix is solved by banded Cholesky factorization (scipy.linalg.solveh_banded), which applies
    because the matrix is symmetric positive definite as soon as one node is tied to ground. The function returns
    (pressures, the residual of the current law at each node)."""
    from scipy.linalg import solveh_banded
    bw = max([abs(i - j) for i, j, _ in edges if j >= 0] + [0])
    ab = np.zeros((bw + 1, n))                # the upper form, ab[bw + i - j, j] = A[i, j] for i <= j
    for i, j, R in edges:
        g = 1.0 / R
        ab[bw, i] += g
        if j >= 0:
            ab[bw, j] += g
            lo, hi = min(i, j), max(i, j)
            ab[bw + lo - hi, hi] -= g
    p = solveh_banded(ab, np.asarray(inject, float))
    res = -np.asarray(inject, float).copy()
    for i, j, R in edges:                       # the flow out of each node through each resistor
        q = (p[i] - (p[j] if j >= 0 else 0.0)) / R
        res[i] += q
        if j >= 0:
            res[j] -= q
    return p, res


def solve_ladder(Rc, RsA, RsB, RgA, RgB, Rc_out, RsA_out, RsB_out):
    """Solve the drawn device as a ladder network.

    The ladder has a center rail c_0..c_{N-1} and two side rails a_i and b_i, and gap i joins c_i to a_i (R_gA) and to
    b_i (R_gB). Unit flow enters at c_0, and the three outlets end at pressure 0. The nodes are numbered c_i = 3i,
    a_i = 3i + 1, and b_i = 3i + 2 (bandwidth 3). The function returns the emergent fractions f_A and f_B (each gap's
    flow over the center flow arriving at its row), F, the center flow entering each row, the gap flows, the outlet
    flows, and the largest residual of the current law."""
    N = len(RgA)
    sc = float(np.median(Rc[:max(N - 1, 1)]))          # resistances in units of a typical center segment
    E = []
    for i in range(N):
        E.append((3 * i, 3 * i + 1, RgA[i] / sc))
        E.append((3 * i, 3 * i + 2, RgB[i] / sc))
        if i < N - 1:
            E.append((3 * i, 3 * i + 3, Rc[i] / sc))
            E.append((3 * i + 1, 3 * i + 4, RsA[i] / sc))
            E.append((3 * i + 2, 3 * i + 5, RsB[i] / sc))
    E += [(3 * (N - 1), -1, Rc_out / sc), (3 * (N - 1) + 1, -1, RsA_out / sc), (3 * (N - 1) + 2, -1, RsB_out / sc)]
    inj = np.zeros(3 * N)
    inj[0] = 1.0
    p, res = kirchhoff(3 * N, E, inj)
    Pc, PA, PB = p[0::3], p[1::3], p[2::3]
    QA = (Pc - PA) / (np.asarray(RgA) / sc)
    QB = (Pc - PB) / (np.asarray(RgB) / sc)
    Qin = np.concatenate([[1.0], (Pc[:-1] - Pc[1:]) / (np.asarray(Rc[:N - 1]) / sc)])
    out = dict(C=Pc[-1] / (Rc_out / sc), A=PA[-1] / (RsA_out / sc), B=PB[-1] / (RsB_out / sc))
    return dict(fA=QA / Qin, fB=QB / Qin, F=1.0 - out['C'], Qc_in=Qin, QgA=QA, QgB=QB, outlets=out,
                kcl_residual_max=float(np.max(np.abs(res))))


def terminal_window(q, eps_net=EPS_NET):
    """Return CIF_v2's terminal window as (q_min, lambda, n_term).

    A disturbance of the ladder decays row by row like z+^(-k), where z+ is the root above 1 of z + 1/z = 2 + q for the
    loop ratio q = (2 R_c + R_s)/R_g. Thus it falls below eps_net within lambda ln(1/eps_net) rows, with
    lambda = 1/ln z+ taken at the smallest q, which gives the slowest decay."""
    qmin = float(np.min(q))
    zp = ((2.0 + qmin) + math.sqrt((2.0 + qmin) ** 2 - 4.0)) / 2.0
    lam = 1.0 / math.log(zp)
    return qmin, lam, max(2, int(math.ceil(lam * math.log(1.0 / eps_net))) + 1)


def runs_of(rows):
    """Return the consecutive runs in a set of row indices as [(first, last), ...]."""
    rows = sorted(set(int(r) for r in rows))
    out = []
    for r in rows:
        if out and r == out[-1][1] + 1:
            out[-1][1] = r
        else:
            out.append([r, r])
    return [tuple(x) for x in out]


def exact_closure(Rc, RgX, P, Qg, Sx, Rc_out):
    """Return the side resistances that the design's closure asks for (opencif.physics.schedule.build), evaluated with
    the given (drawn) R_c and R_g and the plan's flows.

    The side and center pressure drops between two gaps differ by the change of the gap drop R_g Q_g, so
        R_s(i) = [R_c(i) P(i) - R_g(i) Q_g(i) + R_g(i+1) Q_g(i+1)] / S(i),
    and at the outlet R_s,out = [R_c,out P(N) - R_g(N) Q_g(N)] / S(N). The PE closure (App. A eq. closure) is the
    same without the gap terms, R_s = (P/S) R_c."""
    N = len(P)
    Pc = np.zeros(N)
    for k in range(N - 1):
        Pc[k + 1] = Pc[k] - Rc[k] * P[k]
    Ps = Pc - RgX * Qg
    Rs = (Ps[:-1] - Ps[1:]) / Sx[:-1]
    P_out = Pc[N - 1] - Rc_out * P[N - 1]
    return Rs, (Ps[N - 1] - P_out) / Sx[N - 1]


def pe_closure(f, Rc, sigma0=0.0):
    """Return the PE closure's side resistances (App. A eq. sigma and eq. closure) in the form of CIF_v2's P3.

    The side-flow ratio marches as sigma(i) = (sigma(i-1) + f(i))/(1 - 2 f(i)) from sigma_0, which is 0 for a sealed
    inlet, and R_s(i) = R_c(i)/sigma(i). The audit takes sigma from the plan's own P and S, which is this march of the
    plan's f."""
    sig, out = float(sigma0), []
    for fi, rc in zip(f, Rc):
        sig = (sig + fi) / (1.0 - 2.0 * fi)
        out.append(rc / sig)
    return np.array(out)


def strict_inlet_check(q_ul_min, w_um, h_um, pitch_um, rho=1000.0, mu=1e-3):
    """Evaluate CIF_v2 P8's own criteria at the inlet, which are Re = rho u D_h/mu within its limit 1.0, with
    u = Q/(w h) and D_h = 2 w h/(w + h), and the naive entrance length 0.05 Re D_h within half the pitch. Spec v5
    S16-c's entrance length (the form that keeps its Stokes floor) is reported beside them."""
    UM = S.UM
    q_ = q_ul_min * 1e-9 / 60.0
    dh = 2.0 * w_um * h_um / (w_um + h_um) * UM
    re_ = rho * (q_ / (w_um * h_um * UM * UM)) * dh / mu
    le_ = 0.05 * re_ * dh / UM
    return dict(flow_uL_min=float(q_ul_min), inlet_Re=_f(re_, 4), within_Re_limit=bool(re_ <= STRICT_RE_MAX),
                naive_Le_um=_f(le_, 4), naive_Le_over_half_pitch=_f(le_ / (0.5 * pitch_um), 4),
                entrance_length_um=_f(K.entrance_length(re_, dh) / UM, 4),
                passes=bool(re_ <= STRICT_RE_MAX and le_ <= 0.5 * pitch_um), Re_exact=float(re_))


def bulk_rows(N, n_term, declared_runs=()):
    """Return CIF_v2's bulk (P4) as a mask over the rows.

    The bulk is every row but the first and the terminal window (the last n_term rows), less a window on each side of
    every declared departure. Each departure is a run of rows (first, last[, window]), 0-based, and its window is
    n_term rows when the run gives none of its own. A departure disturbs the network as a boundary does, and the
    disturbance decays over a window of the same kind (departure_window)."""
    b = np.ones(N, bool)
    b[0] = False
    b[max(N - n_term, 0):] = False
    for run in declared_runs:
        a, c = run[0], run[1]
        w = run[2] if len(run) > 2 else n_term
        b[max(a - w, 0):min(c + w + 1, N)] = False
    return b


def departure_window(lam, amplitude, n_term, eps_net=EPS_NET):
    """Return the rows that a declared departure's disturbance needs to decay below eps_net. CIF_v2's terminal window
    is sized for a disturbance of order one, lambda ln(1/eps_net), so a departure whose fractions are off by A
    (relative) needs lambda ln(A/eps_net), and never fewer rows than the terminal window."""
    if not np.isfinite(amplitude) or amplitude <= 1.0:
        return n_term
    return max(n_term, int(math.ceil(lam * math.log(amplitude / eps_net))) + 1)


def lamina_ratio(f, w, dstar, x_of=None):
    """Return delta/delta* of a gap that extracts the fraction f from a center channel of width w. The lamina is
    x(f) w, where x is the kernel's inverse (the cubic's, Appendix A Result 2, unless another is given). The ratio is
    nan when f is outside (0, 1/2), where the gap has no cutoff."""
    x_of = x_of or S.x_of
    return x_of(f) * w / dstar if 0.0 < f < 0.5 else float('nan')


# ---------------------------------------------------------------------------------------------------------------------
# the inventory of drawn widths (P9; pipeline.assess's aspect-ratio check reads the same inventory)
# ---------------------------------------------------------------------------------------------------------------------

def drawn_widths(V, spec):
    """Return the narrowest measured width of each kind of CIF channel on the drawn mask. The kinds are the gap
    openings, the center channel (both ends of its taper), every straight side section, the U-turn sides, the head
    units, and the outlet leads."""
    dp = V.get('drawn_parameters') or {}
    widths = {'gap opening G': float(dp.get('G_meas_min', spec.G))}
    if dp.get('wc_in_meas') and dp.get('wc_out_meas'):
        widths['centre channel'] = min(float(dp['wc_in_meas']), float(dp['wc_out_meas']))
    sc_ = V.get('side_continuity')
    if sc_:
        widths['straight side channel'] = sc_['narrowest_um']
    lw = [min(lp['A']['w_min'], lp['B']['w_min']) for lp in V.get('loops_measured', [])]
    if lw:
        widths['U-turn side channel'] = min(lw)
    cd = V.get('components_drawn')
    if cd is not None and cd.get('head_unit_widths_um'):
        widths['meander head unit'] = min(cd['head_unit_widths_um'])
    leads = list(((V.get('outlet_measured') or {}).get('lead_width_min_um') or {}).values())
    if leads:
        widths['outlet lead'] = min(leads)
    return widths


def posts_check(posts, sc):
    """Apply CIF_v2's R6b on the file, at every post of every straight row.

    The check takes each post that the ruler found among the islands of the polygon read from the file
    (opencif.checks.ruler.find_posts, which leaves out the head's rows and the rows next to a U-turn). A post's
    extent across the flow must be the wall thickness l_g, which spec v5's gap resistance takes as the gap's length,
    and its extent along the flow must be the post length g_p. The function returns the counts, the largest
    deviations, and the worst post (the 1-based row of the gap before it, and its side)."""
    found, missing, _ = posts
    lg, gp = float(sc.spec.l_g), float(sc.spec.g_p)
    d_across, d_along, d_ratio, worst, worst_dev = 0.0, 0.0, 0.0, None, -1.0
    for (k, side), (xa, xb, ya, yb) in sorted(found.items()):
        along, across = xb - xa, yb - ya
        da, dg, dq = abs(across - lg), abs(along - gp), abs(across / along - lg / gp)
        d_across, d_along, d_ratio = max(d_across, da), max(d_along, dg), max(d_ratio, dq)
        dev = max(da / TOL_POST_UM, dg / TOL_POST_UM, dq / TOL_POST_RATIO)
        if dev > worst_dev:
            worst_dev, worst = dev, dict(row=k + 1, side=side, across_um=round(across, 4), along_um=round(along, 4))
    n = len(found)
    ok = n > 0 and not missing and d_across < TOL_POST_UM and d_along < TOL_POST_UM and d_ratio < TOL_POST_RATIO
    return dict(ok=ok, posts=n, missing=len(missing), missing_at=missing[:10], max_dev_across_um=_f(d_across, 3),
                max_dev_along_um=_f(d_along, 3), max_dev_ratio=_f(d_ratio, 3), worst=worst, l_g_um=lg, g_p_um=gp)


def rulers_agree(ru, m, N):
    """Hold the generator's measurement (opencif.checks.measure) of every straight row that the auditor's own ruler
    read against the ruler (R0), and return the largest relative difference of each kind of resistance and where it
    occurs."""
    worst, rows = {}, 0
    nan = []
    for key, mkey in (('Rc', 'Rc'), ('RsA', 'RsA'), ('RsB', 'RsB'), ('RgA', 'RgA'), ('RgB', 'RgB')):
        arr = np.asarray(m[mkey], float)
        best = (0.0, None)
        for k, v in ru[key].items():
            if not np.isfinite(v) or v <= 0:
                nan.append((key, k + 1))
                continue
            d = abs(arr[k] / v - 1.0)
            if d > best[0]:
                best = (d, k + 1)
        worst[key] = best
        rows = max(rows, len(ru[key]))
    dmax = max(v[0] for v in worst.values())
    return dict(ok=bool(rows > 0 and not nan and not ru['missing'] and dmax <= TOL_RULER),
                segments=len(ru['Rc']), gaps=len(ru['RgA']) + len(ru['RgB']),
                gap_share=_f(len(ru['RgA']) / max(N, 1), 3),
                max_rel_diff=_f(dmax, 3),
                by_kind={k: dict(max_rel_diff=_f(v[0], 3), row=v[1]) for k, v in worst.items()},
                unmeasurable=nan[:10], posts_missing=len(ru['missing']),
                gap_position_max_dev_um=_f(ru['x_dev_max_um'], 3))


# ---------------------------------------------------------------------------------------------------------------------
# the audit
# ---------------------------------------------------------------------------------------------------------------------

def _f(x, n=4):
    """Round x to n significant figures for the record, or return None when x is None or not finite."""
    return None if x is None or not np.isfinite(x) else float('%.*g' % (n, x))


def audit(dxf_path, sc, info, m, V, declared, mode, ar_max, fluid=None, nw=None, intended=None):
    """Audit one photomask and return the audit block of the record.

    dxf_path  the photomask file, read back here (P0)
    sc        the design as drawn (opencif.physics.schedule.Sched, build.plan_and_draw's sc_f)
    info      the drawing's layout (build.plan_and_draw's info), which says where each leg lies, so that the ruler
              starts looking for each post there
    m         the result of opencif.checks.measure.measure_device on the drawn polygon, which holds the resistances of
              every drawn part
    V         the build's validation block, with its DRC on the same file, its corner radii, and its width inventory
    declared  build.declared_rows as {option: [0-based rows]}, empty for a design without declared departures
    mode      'strict', 'evaluate', or 'reproduce'
    ar_max    the process's aspect-ratio limit
    fluid     the drawn polygon that m was measured on; when it is given, P0 requires the file to equal it exactly
    nw        the build's own network solution (measure.network), compared with this auditor's (P2)
    intended  the intended network (build.intended_network: the closure with the declared departures as intended),
              the reference of P4 when it is given; without it, P4 compares with the kernel"""
    t0 = time.time()
    spec = sc.spec
    N = int(sc.N)
    checks = []

    def ck(cid, name, ok, role, value, need, detail=None):
        checks.append(dict(id=cid, check=name, ok=bool(ok), role=role, value=value, requirement=need,
                           detail=detail or {}))

    # ---- P0 the file is the drawing
    pm, facts = read_photomask(dxf_path)
    clean = facts['arc_segments'] == 0 and facts['open_polylines'] == 0 and not facts['other_entities'] \
        and facts['units'] == 13 and facts['valid']
    same = same_polygon(pm, fluid) if fluid is not None else None
    ck('P0', 'the photomask file is the drawing', clean and same is not False, 'hard',
       dict(identical_to_drawn_polygon=same, fluid_outlines=facts['fluid_outlines'], islands=facts['islands'],
            vertices=facts['vertices']),
       'the DXF read back (CHANNELS minus ISLANDS, micrometres, closed straight-segment polylines only) is a valid '
       'polygon equal to the drawn polygon vertex for vertex', facts)
    # ---- P1 the kernel anchors, which guard the ports of Appendix A's kernels
    a1 = K.duct_bracket(1.0)
    sq = K.r(1e-4, 1e-4, 1e-3) * (1e-4 * 1e-12) / 1e-3          # R w h^3/(mu L) of a square duct
    x_rt = max(abs(K.g_cubic(K.x_of_f(f)) - f) / f for f in (1e-6, 1.04e-4, 5.76e-4, 0.01, 0.2, 0.49))
    kap = K.Pairing().kappa
    v1 = dict(duct_bracket_1=a1, square_duct_R=sq, g_0p3=K.g_cubic(0.3), g_half=K.g_cubic(0.5), kappa_g=kap)
    ok1 = all(abs(v1[k] - ANCHORS[k][0]) <= ANCHORS[k][1] for k in ANCHORS) and x_rt < 1e-12 and S.KAPPA_G == kap
    ck('P1', 'kernel anchors', ok1, 'hard', dict({k: float(v) for k, v in v1.items()}, x_of_f_round_trip=x_rt),
       'B(1) = 0.4217310; square duct R w h^3/(mu L) = 12/B(1) = 28.454154; g(0.3) = 0.216; g(1/2) = 1/2; '
       'g(x(f)) = f; kappa_g = 4.9138, the prefactor the schedule uses (App. A Results 1-2; App. B C3b)')
    # ---- R0 the auditor's own ruler (opencif.checks.ruler), which measures every straight row again from the file
    # with code that shares nothing with the generator's measurement, and holds the generator's measurement to it
    posts = RU.find_posts(pm, sc, info)
    ru = RU.measure_straight_rows(pm, sc, info, posts)
    r0 = rulers_agree(ru, m, N)
    ck('R0', 'the auditor\'s own ruler and the generator\'s measurement agree', r0.pop('ok'), 'hard', r0,
       'every post of every ordinary straight row found on the file; every centre segment, side segment and gap the '
       'auditor\'s ruler measures from the file alone (opencif.checks.ruler: posts, openings and cross-sections; no '
       'planned width, position or length) within %g of the generator\'s measurement (opencif.checks.measure); the '
       'auditor\'s network below uses the ruler\'s resistances on those rows and the generator\'s on the head units, '
       'U-turns and outlet' % TOL_RULER)
    # ---- the drawn network, from the generator's measurement with the ruler's resistances on the straight rows
    o = m['outlet']
    Rc_m = np.asarray(m['Rc'], float)
    RsA_m, RsB_m = np.asarray(m['RsA'], float), np.asarray(m['RsB'], float)
    RgA_m, RgB_m = np.asarray(m['RgA'], float), np.asarray(m['RgB'], float)
    Rc, RsA, RsB, RgA, RgB = Rc_m.copy(), RsA_m.copy(), RsB_m.copy(), RgA_m.copy(), RgB_m.copy()
    wc = np.array(m.get('wc_gap', sc.wc), float)                    # the drawn center width at each gap
    for arr, key in ((Rc, 'Rc'), (RsA, 'RsA'), (RsB, 'RsB'), (RgA, 'RgA'), (RgB, 'RgB'), (wc, 'wc_gap')):
        for k, v in ru[key].items():
            if np.isfinite(v) and v > 0:
                arr[k] = v
    net = solve_ladder(Rc, RsA, RsB, RgA, RgB, o['Rc_out'], o['RsA_out'], o['RsB_out'])
    fA, fB, F_em = net['fA'], net['fB'], float(net['F'])
    dstar = float(spec.delta_star)
    # ---- the declared rows, the terminal window, and the bulk
    dec = np.zeros(N, bool)
    for rows in (declared or {}).values():
        dec[[r for r in rows if 0 <= r < N]] = True
    Rc_full = np.concatenate([Rc[:N - 1], [o['Rc_out']]])
    q = np.minimum((2.0 * Rc_full + np.concatenate([RsA[:N - 1], [o['RsA_out']]])) / RgA,
                   (2.0 * Rc_full + np.concatenate([RsB[:N - 1], [o['RsB_out']]])) / RgB)
    qmin, lam, n_term = terminal_window(q)
    # P4 compares the drawn fractions with the intended network, which is the kernel at every row of a design without
    # declared departures. A declared departure, such as Xia's loop pull-in, can shift the fractions of the rows
    # downstream of it without decaying, and the intended network carries that shift, so that it counts as declared
    # rather than as a failure of the drawing. The deviation from the kernel itself is reported beside it. Each declared
    # departure gets its own window, sized for how far its drawn fractions are off from that reference.
    fk = np.array([S.f_of(dstar, w) for w in wc])
    e4k = np.maximum(np.abs(fA / fk - 1.0), np.abs(fB / fk - 1.0))
    if intended is not None:
        fiA, fiB = np.asarray(intended['fA'], float), np.asarray(intended['fB'], float)
        with np.errstate(divide='ignore', invalid='ignore'):
            e4 = np.maximum(np.abs(fA / fiA - 1.0), np.abs(fB / fiB - 1.0))
        e4 = np.where(np.isfinite(e4), e4, e4k)
    else:
        e4 = e4k
    dec_runs = [(a, b, departure_window(lam, float(np.nanmax(e4[a:b + 1])), n_term))
                for a, b in runs_of(np.where(dec)[0])]
    bulk = bulk_rows(N, n_term, dec_runs)
    und = ~dec
    # ---- P2 mass
    S0 = float(getattr(spec, 'S0', 0.0) or 0.0)
    P_, S_, Qg_ = np.asarray(sc.P, float), np.asarray(sc.S, float), np.asarray(sc.Qg, float)
    mass = abs(P_[-1] + 2.0 * (S_[-1] - S0) - 1.0)
    prod = abs((1.0 - F_em) - float(np.prod(1.0 - fA - fB)))
    outl = abs(net['outlets']['C'] + net['outlets']['A'] + net['outlets']['B'] - 1.0)
    # The network's conservation is required to a ten-thousandth of the smallest flow it is used to measure, which is
    # the smallest gap flow at an undeclared row (in units of the inlet flow).
    qg = np.minimum(np.abs(net['QgA']), np.abs(net['QgB']))
    q_small = float(np.min(qg[und & (qg > 0)])) if np.any(und & (qg > 0)) else float(np.median(qg))
    tol_mass = 1e-4 * q_small
    # The outlets' sum differs from the inlet flow by exactly the sum of the node residuals, so it is held to the node
    # tolerance times the number of nodes (3N), which a long device's rounding reaches before any one node does.
    tol_out = 3 * len(net['QgA']) * tol_mass
    par = None
    if nw is not None:
        # The build's solution is compared with this solver's on the same resistances (the generator's measurement),
        # in units of the larger of |f| and f*, so that a row extracting almost nothing is not magnified.
        net_m = solve_ladder(Rc_m, RsA_m, RsB_m, RgA_m, RgB_m, o['Rc_out'], o['RsA_out'], o['RsB_out'])
        fs = np.maximum(np.abs(np.asarray(sc.f, float)), float(sc.info.get('fstar', 0.0)) or 1e-300)
        par = float(max(np.max(np.abs(np.asarray(nw['fA']) - net_m['fA']) / fs),
                        np.max(np.abs(np.asarray(nw['fB']) - net_m['fB']) / fs),
                        abs(float(nw['F']) - float(net_m['F']))))
    ok2 = (mass < 1e-10 and net['kcl_residual_max'] < tol_mass and prod < 1e-9 and outl < tol_out
           and (par is None or par < 1e-6))
    ck('P2', 'mass: the march identity and the drawn network\'s current law', ok2, 'hard',
       dict(march_identity=_f(mass, 3), kcl_residual_max=_f(net['kcl_residual_max'], 3), product_form=_f(prod, 3),
            outlets_sum_minus_1=_f(outl, 3), smallest_gap_flow=_f(q_small, 3), node_tolerance=_f(tol_mass, 3),
            outlets_tolerance=_f(tol_out, 3), solver_agreement=_f(par, 3)),
       '|P(N) + 2[S(N) - sigma_0] - 1| < 1e-10 (App. A); on the drawn network the flows balance at every node to 1e-4 '
       'of the smallest gap flow at an undeclared row (in units of the inlet flow), and the three outlets carry the '
       'inlet flow to that tolerance times the 3N nodes, because their sum differs from the inlet flow by the sum of '
       'the node residuals; 1 - F = prod(1 - f_A - f_B) (App. A eq. Frec) to 1e-9; on the generator\'s measured '
       'resistances the build\'s own network solution (measure.solve3) agrees with this solver\'s (< 1e-6 of the '
       'larger of f and f* at every gap)')
    # ---- P3 the closure realized by the drawn geometry
    worst_x, worst_pe, where = 0.0, 0.0, None
    seg_ok = und[:-1] & und[1:]                                     # a segment joins two undeclared gaps
    rel = {}
    for side, RsX, RgX, RsXo in (('A', RsA, RgA, o['RsA_out']), ('B', RsB, RgB, o['RsB_out'])):
        tgt, tgt_o = exact_closure(Rc, RgX, P_, Qg_, S_, o['Rc_out'])
        dx = np.abs(RsX[:N - 1] / tgt - 1.0)
        pe = np.abs(RsX[:N - 1] / (P_[:N - 1] / S_[:N - 1] * Rc[:N - 1]) - 1.0)
        dxo = abs(RsXo / tgt_o - 1.0) if und[-1] else 0.0
        peo = abs(RsXo / (P_[-1] / S_[-1] * o['Rc_out']) - 1.0) if und[-1] else 0.0
        if np.any(seg_ok):
            i = int(np.argmax(np.where(seg_ok, dx, -1.0)))
            if dx[i] > worst_x:
                worst_x, where = float(dx[i]), (i + 1, side)
            worst_pe = max(worst_pe, float(np.max(pe[seg_ok])))
        if dxo > worst_x:
            worst_x, where = float(dxo), ('outlet', side)
        worst_pe = max(worst_pe, float(peo))
        rel[side] = dict(exact_max=_f(float(np.max(dx[seg_ok])) if np.any(seg_ok) else 0.0, 3), outlet=_f(dxo, 3),
                         pe_max=_f(float(np.max(pe[seg_ok])) if np.any(seg_ok) else 0.0, 3), pe_outlet=_f(peo, 3))
    ck('P3', 'the closure realised by the drawn geometry', worst_x <= TOL_CLOSURE and worst_pe <= TOL_CLOSURE, 'hard',
       dict(exact_closure_max_rel_dev=_f(worst_x, 3), at=where, pe_closure_max_rel_dev=_f(worst_pe, 3)),
       'every drawn side resistance (segments between undeclared gaps, and the outlets) within %g of the design\'s '
       'exact closure recomputed from the drawn centre and gap resistances, and of the PE closure (P/S) R_c (App. A '
       'eq. closure); %g is CIF_v2\'s sanity bound' % (TOL_CLOSURE, TOL_CLOSURE), rel)
    # ---- P4 the network bulk's agreement with the intended network (the kernel without declared departures)
    e4_bulk = float(np.max(e4[bulk])) if np.any(bulk) else None
    e4_und = float(np.max(e4[und])) if np.any(und) else None
    row_of = lambda mask, e: (int(np.argmax(np.where(mask, e, -1.0))) + 1) if np.any(mask) else None
    v4 = dict(reference='the intended network' if intended is not None else 'the kernel',
              bulk_max=_f(e4_bulk, 3), bulk_max_row=row_of(bulk, e4), bulk_rows=int(np.sum(bulk)),
              undeclared_max=_f(e4_und, 3), undeclared_max_row=row_of(und, e4))
    if intended is not None:
        v4.update(kernel_bulk_max=_f(float(np.max(e4k[bulk])) if np.any(bulk) else None, 3),
                  kernel_bulk_max_row=row_of(bulk, e4k),
                  kernel_undeclared_max=_f(float(np.max(e4k[und])) if np.any(und) else None, 3))
    ck('P4', 'network bulk: the drawn fractions against the intended network', e4_bulk is not None and
       e4_bulk <= TOL_F_BULK, 'hard', v4,
       'max |f_i/f_i,int - 1| <= %g over the bulk rows, where f_int is the intended network\'s fraction, which equals '
       'the kernel g(delta*/w_c(i)) at every row of a design without declared departures and carries the effect of the '
       'declared departures on the rows downstream of them otherwise (build.intended_network; the deviation from the '
       'kernel is reported as kernel_bulk_max). The bulk is every row but row 1, the terminal window (the last %d '
       'rows: CIF_v2 tol_f_bulk and its window, P4) and, on each side of every declared departure, the rows its '
       'disturbance needs to decay below eps_net (lambda ln(A/eps_net), A the departure\'s largest deviation; at least '
       '%d)' % (TOL_F_BULK, n_term, n_term),
       dict(q_min=_f(qmin, 4), lambda_rows=_f(lam, 4), n_term=n_term, eps_net=EPS_NET,
            declared_runs=[[a + 1, b + 1, w] for a, b, w in dec_runs]))
    # ---- P5 lamina uniformity, the cutoff guarantee
    dr = np.array([[lamina_ratio(a, w, dstar), lamina_ratio(b, w, dstar)] for a, b, w in zip(fA, fB, wc)])
    undef = und & ~np.all(np.isfinite(dr), axis=1)
    dmax = float(np.nanmax(dr[und])) if np.any(und) else None
    dmin = float(np.nanmin(dr[und])) if np.any(und) else None
    ddec = float(np.nanmax(dr[dec])) if np.any(dec) and np.any(np.isfinite(dr[dec])) else None
    ok5 = dmax is not None and dmax <= 1.0 + TOL_DELTA and not undef.any()
    ck('P5', 'lamina uniformity: the cutoff guarantee', ok5, 'hard',
       dict(max_delta_over_delta_star=_f(dmax, 5), min_delta_over_delta_star=_f(dmin, 5),
            rows_without_cutoff=int(np.sum(undef)), declared_rows=int(np.sum(dec)), declared_max=_f(ddec, 4)),
       'delta_i/delta* = x(f_i) w_c(i)/delta* <= 1 + %g at every undeclared gap, both sides, each with 0 < f < 1/2 '
       '(CIF_v2 tol_delta, P5: no gap passes particles larger than %.0f %% above D)' % (TOL_DELTA, 100 * TOL_DELTA))
    # ---- P6 the filtrate fraction
    F_t, F_march = float(spec.F_target), float(sc.F)
    ck('P6', 'filtrate fraction', F_em >= F_t - TOL_F_TARGET and abs(F_em - F_march) < TOL_F_MARCH, 'hard',
       dict(F_drawn=_f(F_em, 7), F_march=_f(F_march, 7), F_target=F_t),
       'F_drawn >= F_target - %g and |F_drawn - F_march| < %g (CIF_v2 P6)' % (TOL_F_TARGET, TOL_F_MARCH))
    # ---- P6b App. A's march-outlet identity on the drawn outlet
    sA, sB = o['Rc_out'] / o['RsA_out'], o['Rc_out'] / o['RsB_out']
    F_id = (sA + sB) / (1.0 + sA + sB)
    fbar = float(np.mean(np.asarray(sc.f, float)))
    tol6 = max(TOL_F6_MIN, 4.0 * fbar)
    ol = net['outlets']                                             # the identity's premise, S_X/P = R_c,out/R_sX,out
    pe_res = float(max(abs(ol['A'] / ol['C'] - sA), abs(ol['B'] / ol['C'] - sB)) / max(sA, sB))
    ck('P6b', 'a closed form for F: the march-outlet identity at the drawn outlet', abs(F_id - F_em) < tol6, 'hard',
       dict(F_identity=_f(F_id, 7), F_drawn=_f(F_em, 7), difference=_f(F_id - F_em, 3),
            pe_residual_at_outlet=_f(pe_res, 3)),
       '|F_identity - F_drawn| < max(%g, 4 f_mean) = %.3g, F_identity = 2R_c/(R_s + 2R_c) per side at the drawn '
       'outlet resistances (App. A Result: the march and the outlet identity give the same F, for every schedule; its '
       'premise, PE at the outlet row, is not exact in a real network: the residual is reported)' % (TOL_F6_MIN, tol6))
    # ---- P7 port matching (sensitivity, reported)
    net3 = solve_ladder(Rc, RsA, RsB, RgA, RgB, o['Rc_out'], PORT_MISMATCH * o['RsA_out'], PORT_MISMATCH * o['RsB_out'])
    d3 = np.array([[lamina_ratio(a, w, dstar), lamina_ratio(b, w, dstar)]
                   for a, b, w in zip(net3['fA'], net3['fB'], wc)])
    ck('P7', 'port matching (reported)', True, 'report',
       dict(F_drawn=_f(F_em, 6), F_side_outlets_x3=_f(float(net3['F']), 6),
            max_delta_over_delta_star_x3=_f(float(np.nanmax(d3[und])) if np.any(und) else None, 4)),
       'both side outlets\' resistance x%g: how far F and the cutoff move (CIF_v2 P7: the outlet\'s matching is '
       'load-bearing; the record\'s model band gives the +-1 %% case)' % PORT_MISMATCH)
    # ---- M14 the Stokes and Reynolds gates of the design envelope (the design code's M14), on the drawn mask, at the
    # design flow (a check of its own, beside the auditor's P8)
    rho = float(getattr(spec, 'rho_f', 1000.0))
    mu = float(getattr(spec, 'mu_f', S.MU))
    h = float(spec.h)
    Q = float(spec.Q_in_uL_min) * 1e-9 / 60.0
    UM = S.UM
    Re_c = 2.0 * rho * Q * net['Qc_in'] / (mu * (wc + h) * UM)           # rho u D_h/mu = 2 rho Q/(mu (w + h))
    gaps = {(int(k), s_): float(wg) for k, s_, wg, _ in m['gaps']}
    GA = np.array([gaps.get((i, 'A'), float(sc.wL[i])) for i in range(N)])
    GB = np.array([gaps.get((i, 'B'), float(sc.wL[i])) for i in range(N)])
    Re_g = np.maximum(2.0 * rho * Q * np.abs(net['QgA']) / (mu * (GA + h) * UM),
                      2.0 * rho * Q * np.abs(net['QgB']) / (mu * (GB + h) * UM))
    a_, b_ = np.maximum(wc, h) * UM, np.minimum(wc, h) * UM
    Cs = {}
    C_ = np.array([Cs.setdefault(round(float(x), 9), K.C_shear(float(x))) for x in a_ / b_])
    gam = C_ * Q * net['Qc_in'] / (a_ * b_ ** 2)                         # the peak wall shear (spec v5 S16-b)
    dl = np.fmax(dr[:, 0], dr[:, 1]) * dstar                          # the wider lamina of the two sides (nan if none)
    Re_d = rho * gam * (np.where(np.isfinite(dl), np.maximum(dl, dstar), dstar) * UM) ** 2 / mu
    cut = und & np.all(np.isfinite(dr), axis=1)                         # the cutoff rows, undeclared and with a cutoff
    # CIF_v2's literal row-1 term (its M14 at the full inlet flow, gap 1 at f*) is included, as build.gates does.
    re_g1 = 2.0 * rho * (S.f_of(dstar, float(spec.wc0)) * Q) / (mu * (float(spec.G) + h) * UM)
    g8 = dict(Re_channel_max=_f(float(np.max(Re_c)), 4), Re_channel_row=int(np.argmax(Re_c)) + 1,
              Re_gap_max=_f(max(float(np.max(Re_g[cut])) if np.any(cut) else 0.0, re_g1), 4),
              Re_delta_max=_f(float(np.max(Re_d[cut])) if np.any(cut) else 0.0, 4),
              Re_gap_declared_max=_f(float(np.max(Re_g[dec])) if np.any(dec) else None, 3),
              Re_delta_declared_max=_f(float(np.max(Re_d[dec])) if np.any(dec) else None, 3),
              flow_uL_min=float(spec.Q_in_uL_min), fluid=dict(rho_kg_m3=rho, mu_Pa_s=mu))
    ok8 = g8['Re_channel_max'] <= RE_MAX and g8['Re_gap_max'] <= RE_STOKES and g8['Re_delta_max'] <= RE_STOKES
    ck('M14', 'the design code\'s Reynolds and Stokes gates (M14) on the drawn mask', ok8, 'envelope', g8,
       'at the design flow, from the drawn widths and the drawn network\'s flows: the centre channel\'s Re = 2 rho Q_c/'
       '(mu (w_c + h)) <= %g at every row; Re_gap = 2 rho Q_g/(mu (G + h)) and Re_delta = rho gamma_max delta^2/mu '
       '<= %g at every cutoff gap (CIF_v2 M14; App. B: "the gates are Re_stokes = 10^-2 on the gap and lamina and '
       'Re_max = 10 on the channel"); declared rows reported' % (RE_MAX, RE_STOKES))
    # ---- P8 the auditor's own, stricter channel limit (CIF_v2's P8, reported)
    at_d = strict_inlet_check(float(spec.Q_in_uL_min), float(wc[0]), h, float(spec.Lam_L), rho, mu)
    at_5 = strict_inlet_check(STRICT_Q_UL_MIN, float(wc[0]), h, float(spec.Lam_L), rho, mu)
    ck('P8', 'CIF_v2\'s P8: the auditor\'s own, stricter inlet limit (reported)', True, 'report',
       dict(at_design_flow=at_d, at_auditor_default_flow=at_5,
            result=('PASSES' if at_d['passes'] else 'EXCEEDS') + ' at the design flow, ' +
                   ('PASSES' if at_5['passes'] else 'EXCEEDS') + ' at %g uL/min' % STRICT_Q_UL_MIN),
       'the drawn inlet (centre width at gap 1) at Re <= %g and the naive entrance length 0.05 Re D_h <= Lambda/2 '
       '(CIF_v2 PhysicsSpec.Re_max and P8_stokes_gates); reported, not a pass condition: the design envelope is '
       'Re <= %g (the design code\'s M14: check M14), and 0.05 Re D_h is the form spec v5, the glossary and CIF_v2\'s '
       'entrance_length call wrong in this regime (docs/LIMITATIONS.md 6); entrance_length_um is S16-c\'s form'
       % (STRICT_RE_MAX, RE_MAX))
    # ---- P9 the fabrication floors
    drc = V.get('drc') or {}
    drc_ok = bool(drc) and 'error' not in drc and drc.get('n_lost_open', 1) == 0 and drc.get('n_lost_solid', 1) == 0 \
        and drc.get('n_components', 2) == 1 and drc.get('invalid_rings', 1) == 0
    widths = drawn_widths(V, spec)
    w_kind = min(widths, key=widths.get)
    w_n = float(widths[w_kind])
    floor_ok = w_n >= float(spec.F_min) - 0.05
    ck('P9', 'fabrication floors on the drawn geometry', drc_ok and floor_ok, 'hard',
       dict(drc_clean=drc_ok, narrowest_um=_f(w_n, 5), narrowest=w_kind, floor_um=float(spec.F_min),
            drc={k: drc.get(k) for k in ('n_lost_open', 'n_lost_solid', 'n_components', 'invalid_rings', 'd_um')}),
       'the DRC at the process floor finds nothing lost and one fluid component, and the narrowest drawn CIF channel '
       'is at least F_min (within the 0.05 um of the section sampling): CIF_v2 P9\'s width and wall floors')
    pe = posts_check(posts, sc)
    ck('R6b', 'the posts read back from the file (CIF_v2 R6b)', pe.pop('ok'), 'hard', pe,
       'every post of every straight row read back from the photomask (an island of the file; the head\'s rows and the '
       'rows next to a U-turn left out): its extent across the flow within %g um of the wall thickness l_g = %g um, '
       'along the flow within %g um of g_p = %g um, and across/along within %g of l_g/g_p (CIF_v2 R6b, there on one '
       'post of its own DXF)' % (TOL_POST_UM, float(spec.l_g), TOL_POST_UM, float(spec.g_p), TOL_POST_RATIO))
    ar = float(spec.h) / w_n
    ck('P9-AR', 'aspect ratio of the narrowest CIF channel', ar <= float(ar_max) + 1e-9,
       'hard' if mode == 'strict' else 'report', dict(depth_over_width=_f(ar, 4), at=w_kind, width_um=_f(w_n, 5)),
       'depth/width <= AR_max = %g (CIF_v2 w_min = max(F_min, h/AR_max) = %.2f um): a pass condition in strict mode, '
       'reported in evaluate and reproduce modes (docs/LIMITATIONS.md 7)'
       % (float(ar_max), max(float(spec.F_min), float(spec.h) / float(ar_max))))
    # ---- P10 rounded corners
    rd = V.get('rounding') or {}
    ck('P10', 'rounded corners', bool(rd.get('ok')), 'hard',
       dict(min_local_radius_um=_f(rd.get('min_local_radius_um'), 4), criterion_um=rd.get('criterion_um'),
            at=rd.get('at')),
       'the smallest local radius of curvature of the drawn boundary >= the spec\'s r_min (CIF_v2 P10: no sharp edge '
       'on the ink boundary)')
    order = ['P0', 'R0', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P6b', 'P7', 'P8', 'M14', 'P9', 'P9-AR', 'R6b', 'P10']
    checks.sort(key=lambda c: order.index(c['id']))
    hard_fail = [c['id'] for c in checks if c['role'] == 'hard' and not c['ok']]
    env_fail = [c['id'] for c in checks if c['role'] == 'envelope' and not c['ok']]
    status = 'FAIL' if hard_fail else ('PASS_OUTSIDE_ENVELOPE' if env_fail else 'PASS')
    return dict(auditor='opencif.checks.auditor: CIF_v2.py v2.1.0\'s physics auditor (ladder_physics audit(), P1-P10, '
                        'its gate M14 and its acceptance check R6b) ported to opencif and run on the photomask file '
                        'read back, with its own ruler on the straight rows (R0)',
                photomask=dict(file=facts['file'], sha256=facts['sha256']),
                status=status, hard_failures=hard_fail, envelope_failures=env_fail,
                checks=checks, rows=dict(N=N, declared=int(np.sum(dec)), bulk=int(np.sum(bulk)), n_term=n_term),
                seconds=round(time.time() - t0, 2))


def summary_line(a):
    """Return one line for the build log."""
    bad = ', '.join(a['hard_failures'] + a['envelope_failures'])
    return 'physics audit of %s: %s%s' % (a['photomask']['file'], a['status'], (' (%s)' % bad) if bad else '')


# ---------------------------------------------------------------------------------------------------------------------
# the standalone audit of a finished build: python -m opencif audit RUN_DIR
# ---------------------------------------------------------------------------------------------------------------------

def audit_run(run_dir, name=None, say=print, drc=True):
    """Audit a finished build from its folder, using the photomask file alone.

    The record (<name>.json) gives the specification that the build was made from (its `input`). The design is
    planned and its outline drawn again in memory with the build's own planning and drawing code, and nothing is
    written. The photomask is read back from the DXF file and must equal that outline vertex for vertex (P0).
    Everything the audit uses is then measured on the polygon read from the file, which covers the resistances of
    every part, the side-channel continuity and the narrowest widths, and the corner radii, and the DRC runs on the
    file. Nothing is taken from the build's record but its input, and no saved Python object (pickle) is loaded. The
    function returns the audit block, as in the record's validation.physics_audit."""
    import json
    from opencif import pipeline as PL
    from opencif import build as BU
    from opencif.checks import measure as M
    t0 = time.time()
    if name is None:
        cands = [f[:-5] for f in os.listdir(run_dir) if f.endswith('.json') and not f.endswith('_audit.json')
                 and not f.endswith('_physics_audit.json')]
        if len(cands) != 1:
            raise ValueError('%s: name the build (more than one record, or none: %s)' % (run_dir, cands))
        name = cands[0]
    rec = json.load(open(os.path.join(run_dir, name + '.json')))
    if 'input' not in rec or not (rec.get('files') or {}).get('dxf'):
        raise ValueError('%s: the record has no input or no photomask (a refused build has nothing to audit)' % name)
    dxf = os.path.join(run_dir, rec['files']['dxf'])
    S.clear_caches()                                     # as in pipeline.run
    spec, opts, sc, g, env_bad, notes = PL.prepare(rec['input'], lambda m_: None)
    kind, n_legs, spec, sc = PL.fit_layout(spec, sc, opts, notes, lambda m_: None)
    fluid, sc_f, info, _, _, _, _ = BU.plan_and_draw(spec, kind, n_legs, opts['leg_spacing'], opts['failsafe'],
                                                     lambda m_: None, name, t0)
    say('[%s] planned and drawn again in memory (%.0fs)' % (name, time.time() - t0))
    pm, facts = read_photomask(dxf)
    say('[%s] photomask read back: %d outline(s), %d islands, %d vertices (%.0fs)' % (
        name, facts['fluid_outlines'], facts['islands'], facts['vertices'], time.time() - t0))
    m = M.measure_device(pm, sc_f, info)                 # measured on the polygon read from the file
    nw = M.network(sc_f, m)
    say('[%s] measured on the file\'s polygon (%.0fs)' % (name, time.time() - t0))
    V = dict(rounding=BU.rounding_check(pm, spec), side_continuity=BU.side_continuity(pm, sc_f, info),
             components_drawn=BU.components_drawn(pm, sc_f, info))
    BU.outlet_and_loops_measured(V, m, spec)
    kL = int(next(k for k in range(sc_f.NM, sc_f.N) if sc_f.kind[k] == 'L'))
    V['drawn_parameters'] = BU.drawn_parameters(m, sc_f, pm, info, kL)
    if drc:
        import tempfile
        with tempfile.TemporaryDirectory() as td:
            V['drc'] = BU.run_drc(dxf, td, name, spec)
    inten, _, decl = BU.intended_and_declared(sc_f, info, m, spec, S.KAPPA_G * spec.delta_star, opts['cut_tol'],
                                              opts['mode'], opts['seal'])
    a = audit(dxf, sc_f, info, m, V, decl, opts['mode'], opts.get('AR_max', BU.AR_MAX_DEFAULT), fluid=fluid, nw=nw,
              intended=inten)
    a['standalone'] = dict(record=name + '.json', seconds=round(time.time() - t0, 1),
                           method='planned and drawn again from the record\'s input; every measurement on the polygon '
                                  'read from the photomask file' + ('' if drc else '; the DRC not run (P9 incomplete)'))
    return a
