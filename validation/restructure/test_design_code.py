"""test_design_code.py: the design code's reports and acceptance validators (opencif/checks/design_code.py), and two
properties that the design code checks of itself and that opencif had shown only in its release records.

The suite checks six things.
1. Each report and validator is compared with the design code's own function (reference/CIF_v2.py) on the design
   code's own reference device. The laminar coefficient c2 and the transfer bias agree bit for bit, and the
   sensitivity map, the operating regime, and the steric screen agree to rounding.
2. The Monte-Carlo fabrication tolerance is deterministic, exact at zero spread, and spreads with sigma.
3. The reports on a drawn device carry every key, and on a linear device without a head the drawn network's pressure
   drop is the trunk's.
4. Determinism (the design code's R1) holds, so the same specification drawn in two separate processes gives the same
   photomask, byte for byte.
5. Invariance (the design code's M9) holds, so the schedule does not depend on the flow or the fluid's viscosity, and
   N, F, and the fractions do not depend on the pitch. The GS column of M9 concerns the design code's folded side
   channels, which opencif does not draw.
6. Gifford 2014, Fig. 2 (the design code's M19) follows from opencif's planner. For a constant fraction f_gap and the
   endpoint w_s = w_c (F = 2/3), N is Gifford's i_T = ln 3 / (-ln(1 - 2 f_gap)) rounded up, the flow ratio follows the
   closed form [(1 - 2 f_gap)^(-i) - 1]/2, the ratio does not depend on the depth whereas the widths do, and the six
   values of panel A, k x 9.6e-5, end where the curves digitized from the published figure end.
Run it with `python3 tests/test_design_code.py` (about two minutes)."""
import copy
import hashlib
import os
import subprocess
import sys
import tempfile

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, 'reference'))
from opencif import pipeline as PL, spec as SP, build as BU    # noqa: E402
print('opencif %s (code %s)' % (PL.CODE['version'], PL.CODE['sha256']))
from opencif.checks import design_code as DCR, measure as M     # noqa: E402
from opencif.physics import kernels as K, schedule as S          # noqa: E402
import CIF_v2 as C                                               # noqa: E402  the frozen design code

fails = []


def check(name, cond, info=''):
    print('%-84s %s %s' % (name, 'ok ' if cond else 'FAIL', info))
    if not cond:
        fails.append(name)


def rel(a, b):
    """Return the relative difference of a from b."""
    return abs(a - b) / max(abs(b), 1e-300)


UM = 1e-6
# 1. Against the design code's own functions, on its reference device
betas = np.geomspace(0.02, 200.0, 300)
check('c2(beta) = the design code\'s c2, bit for bit, 300 aspect ratios in [0.02, 200]',
      np.array([K.c2(b) for b in betas]).tobytes() == np.array([C.c2(b) for b in betas]).tobytes())
check('c2_exact = the design code\'s _c2_exact, bit for bit (beta 0.3, 1, 2.33, 10)',
      all(K.c2_exact(b) == C._c2_exact(b) for b in (0.3, 1.0, 2.3337282112032494, 10.0)))
check('c2 tends to the cubic\'s 3 for a wide slot (beta = 1e4: within 1e-3)', abs(K.c2(1e4) - 3.0) < 1e-3,
      '%.6f' % K.c2(1e4))
dev = C.design(C.reference_request())
proc, pair = dev.req.proc, dev.req.pair
rb = [rw.beta for rw in dev.rows]
check('the transfer bias = the design code\'s Pairing.transfer_bias on its reference device, bit for bit',
      K.Pairing().transfer_bias(rb) == C.Pairing().transfer_bias(rb) == dev.bias, '%s' % (dev.bias,))
sm = DCR.sensitivity_map(dev.rows[-1].w_c / UM, dev.rows[-1].w_s / UM, dev.rows[0].w_c / UM, dev.rows[0].f,
                         proc.h / UM, pair.kappa, proc.mu)
check('the sensitivity map = the design code\'s sensitivity_map (S16-e), every entry within 1e-6',
      all(rel(sm[k], v) < 1e-6 if v else abs(sm[k]) < 1e-12 for k, v in dev.sens.items()),
      ', '.join('%s %.4f' % (k.split('/')[-1], v) for k, v in sm.items()))
vo = DCR.operating_regime([rw.w_c / UM for rw in dev.rows], [rw.f for rw in dev.rows], proc.h / UM,
                          (proc.G + proc.g_p) / UM, 5.0, proc.mu, proc.rho_f)
vc = C.V_operating_regime(dev, proc)
check('V_operating_regime = the design code\'s (Re, L_e, peak shear, trunk pressure drop, the four gates)',
      rel(vo['Re'], vc['Re']) < 1e-3 and rel(vo['Le_um'], vc['Le_um']) < 1e-3
      and rel(vo['gdot_per_s'], vc['gdot']) < 1e-3
      and rel(vo['dP_kPa'], float(vc['dP_kPa'])) < 1e-3 and vo['gates'] == {k: bool(v) for k, v in vc['gates'].items()},
      'Re %s, gdot %s 1/s, dP %s kPa' % (vo['Re'], vo['gdot_per_s'], vo['dP_kPa']))
vs = DCR.steric(dev.req.D_target / UM, proc.G / UM, proc.h / UM, dev.req.D_target / pair.kappa / UM)
vsc = C.V_steric(dev, proc)
check('V_steric = the design code\'s (D/G, D/h, delta*/h and the verdict)',
      vs['ok'] == vsc['ok'] and all(rel(vs[k], vsc[k]) < 1e-3 for k in ('D_over_G', 'D_over_h', 'dstar_over_h')))
# 2. The Monte-Carlo fabrication tolerance on a planned opencif design
T = tempfile.mkdtemp(prefix='opencif_dc_test_')
d = SP.load(os.path.join(ROOT, 'examples', 'matrix', 'm15_linear_shear_nomeander_filter_teardrop.yaml'))
S.clear_caches()
spec, opts, sc, g, env_bad, notes = PL.prepare(d)
kind, n_legs, spec, sc = PL.fit_layout(spec, sc, opts, notes, lambda m_: None)
fluid, sc_f, info, _, _, _, _ = BU.plan_and_draw(spec, kind, n_legs, opts['leg_spacing'], opts['failsafe'],
                                                 lambda m_: None, 'm15', 0.0)
f0 = DCR.fab_tolerance(sc_f, info, n=3, sigma_um=0.0)
F_design_net = float(DCR.NW.solve(sc_f.Rc[:-1], sc_f.Rs, sc_f.Rg, sc_f.Rc_out, sc_f.Rs_out)['F'])
check('fabrication tolerance at sigma 0: every draw the design network\'s F, yield 1',
      f0['F_std'] == 0.0 and rel(f0['F_mean'], F_design_net) < 1e-6 and f0['yield_F'] == 1.0, 'F %s' % f0['F_mean'])
fa, fb = DCR.fab_tolerance(sc_f, info), DCR.fab_tolerance(sc_f, info)
fw = DCR.fab_tolerance(sc_f, info, sigma_um=1.0)
check('fabrication tolerance (24 draws, sigma 0.25 um, seed 0): deterministic; a wider spread spreads F',
      fa == fb and fw['F_std'] > fa['F_std'] > 0,
      'F std %.2e (0.25 um), %.2e (1 um); yield %.2f' % (fa['F_std'], fw['F_std'], fa['yield_F']))
# 3. The reports on the drawn device
dxf = os.path.join(T, 'm15.dxf')
from opencif.geometry import dxf as io_dxf                      # noqa: E402
io_dxf.write_dxf(fluid, dxf)
m = M.measure_device(fluid, sc_f, info)
rep = DCR.reports(sc_f, info, m, sc_f.spec, float(M.network(sc_f, m)['F']))
keys = {'B3_pressure_drop', 'B6_F_real', 'B7_transfer_bias', 'B8_sensitivity_map', 'B9_wall_porosity',
        'V_operating_regime', 'V_steric', 'V_fab_tolerance', 'role'}
b3 = rep['B3_pressure_drop']
check('reports on m15 (linear, no head): every report; the drawn network\'s pressure drop within 1 % of the trunk '
      'estimate; eta F; phi = G/Lambda',
      set(rep) == keys and b3['trunk_kPa'] > 0 and abs(b3['drawn_device_kPa'] / b3['trunk_kPa'] - 1) < 0.01
      and abs(rep['B6_F_real']['eta_F'] - 0.903 * rep['B6_F_real']['F_drawn']) < 1e-5
      and abs(rep['B9_wall_porosity']['phi'] - spec.G / (spec.G + spec.g_p)) < 1e-4,
      'dP trunk %s, device %s kPa; transfer bias %s' % (b3['trunk_kPa'], b3['drawn_device_kPa'],
                                                         rep['B7_transfer_bias']['D_true_over_D_predicted']))
# 4. Determinism: the same photomask from two separate processes
prog = ('import sys; sys.path.insert(0, %r)\n'
        'from opencif import pipeline as PL, spec as SP, build as BU\n'
        'from opencif.physics import schedule as S\n'
        'from opencif.geometry import dxf\n'
        'S.clear_caches()\n'
        'sp, op, sc, g, e, n = PL.prepare(SP.load(%r))\n'
        'k, nl, sp, sc = PL.fit_layout(sp, sc, op, n, lambda m: None)\n'
        'fl = BU.plan_and_draw(sp, k, nl, op["leg_spacing"], op["failsafe"], lambda m: None, "m15", 0.0)[0]\n'
        'dxf.write_dxf(fl, sys.argv[1])\n') % (ROOT, os.path.join(ROOT, 'examples', 'matrix',
                                                                'm15_linear_shear_nomeander_filter_teardrop.yaml'))
hs = []
for i in range(2):
    p_ = os.path.join(T, 'det_%d.dxf' % i)
    subprocess.run([sys.executable, '-I', '-B', '-c', prog, p_], check=True)
    hs.append(hashlib.sha256(open(p_, 'rb').read()).hexdigest())
hs.append(hashlib.sha256(open(dxf, 'rb').read()).hexdigest())
check('determinism (R1): two processes and this one draw the same photomask, byte for byte', len(set(hs)) == 1,
      hs[0][:16])
# 5. Invariance (M9) of the schedule under changes of the flow and of the viscosity
d2 = SP.load(os.path.join(ROOT, 'examples', 'decoupling', 'dc_D4_F0p5.yaml'))
base = S.plan(SP.to_engine(d2)[0])
same = True
for flow, mu in ((300.0, 1e-3), (30.0, 2.5e-3), (3.0, 1.2e-3)):
    d3 = copy.deepcopy(d2)
    d3.setdefault('operation', {})['flow_uL_min'] = flow
    d3['fluid'] = dict(density_kg_m3=1025.0, viscosity_Pa_s=mu, rheology='newtonian')
    S.clear_caches()
    sc3 = S.plan(SP.to_engine(SP.load(d3))[0])
    for a in ('f', 'P', 'S', 'wc', 'ws', 'Rs', 'x'):
        same &= np.array_equal(getattr(sc3, a), getattr(base, a), equal_nan=True)
    same &= sc3.N == base.N and sc3.NM == base.NM
check('invariance (M9): the schedule (f, P, S, widths, closure) at 3, 30 and 300 uL/min and 1-2.5 mPa s is the same',
      same, 'N %d' % base.N)
# In M9's pitch column, a longer or shorter post (Lambda = G + g_p) changes the side segments' lengths and so their
# widths, and never N, F, or the fractions.
same_p = True
for gp in (20.0, 40.0):
    d3 = copy.deepcopy(d2)
    d3.setdefault('process', {})['post_length_um'] = gp
    S.clear_caches()
    sc3 = S.plan(SP.to_engine(SP.load(d3))[0])
    for a in ('f', 'P', 'S', 'wc'):
        same_p &= np.array_equal(getattr(sc3, a), getattr(base, a), equal_nan=True)
    same_p &= sc3.N == base.N and sc3.F == base.F and sc3.spec.Lam_L == float(sc3.spec.G) + gp
check('invariance (M9, the pitch): post length 20 and 40 um instead of 27: the same N, F, f, P, S and centre widths',
      same_p, 'N %d, F %.9f' % (base.N, base.F))

# 6. Gifford 2014, Fig. 2 (the design code's M19). Each plan asks for F = 2/3, the paper's endpoint w_s = w_c, and for
# the cutoff that f_gap implies at w_c = 100 um, so that the planner's constant f(i) equals f_gap. The process floor and
# the aspect-ratio limit lie far below the first side width, and the gap's entrance term is zero, because Gifford's
# eqn (4) contains no gap resistance.
import math                                                      # noqa: E402
import opencif                                                   # noqa: E402


def fig2_plan(f, h, wc=100.0):
    """Return the plan of the Fig. 2 device with the constant fraction f, the depth h, and the center width wc (in um),
    planned without a head in evaluate mode."""
    D = S.KAPPA_G * S.x_of(f) * wc
    d = opencif.build_spec(D, 2.0 / 3.0, process=dict(depth_um=h, gap_um=20.5, post_length_um=25.6,
                                                      post_width_um=50.6, min_feature_um=0.2,
                                                      aspect_ratio_max=5000.0, gap_entrance_c=0.0),
                           width_in_um=wc, meander=None, mode='evaluate', name='fig2')
    S.clear_caches()
    return S.plan(SP.to_engine(SP.load(d))[0])


def i_T(f):
    """Return Gifford's i_T = ln 3 / (-ln(1 - 2 f)), the number of rows at which F reaches 2/3."""
    return math.log(3.0) / -math.log(1.0 - 2.0 * f)


F_AB = [k * 9.6e-5 for k in (6, 5, 4, 3, 2, 1)]                  # the caption's six values, 5.76e-4 down to 9.6e-5
ok_n = ok_f = ok_sig = ok_end = True
worst_sig = 0.0
for f in F_AB:
    pl = fig2_plan(f, 100.0)
    i = np.arange(1, pl.N + 1, dtype=float)
    sig = pl.S / pl.P
    worst_sig = max(worst_sig, float(np.max(np.abs(sig - ((1 - 2 * f) ** (-i) - 1) / 2))))
    ok_n &= pl.N == math.ceil(i_T(f))
    ok_f &= float(np.max(np.abs(pl.f / f - 1))) < 1e-12
    ok_end &= abs(pl.ws[pl.N - 2] - 100.0) < 0.1
check('M19, panel A: N = ceil(i_T) and f(i) = f_gap at every row, for the six f_gap = k x 9.6e-5', ok_n and ok_f,
      'N %s' % [math.ceil(i_T(f)) for f in F_AB])
check('M19, panel B: Q_s(i)/Q_c(i) = [(1 - 2 f_gap)^(-i) - 1]/2 at every row (within 1e-12)', ok_sig and
      worst_sig < 1e-12, '%.1e' % worst_sig)
check('M19, the endpoint: the last side segment is within 0.1 um of w_c = 100 um', ok_end)
CD = [fig2_plan(1.76e-4, h) for h in (150.0, 100.0, 50.0, 25.0, 10.0, 5.0)]
n_cd = {pl.N for pl in CD}
same_sig = all(np.array_equal(pl.S / pl.P, CD[0].S / CD[0].P) for pl in CD)
fan = CD[0].ws[999] - CD[-1].ws[999]
check('M19, panels C and D: at f_gap = 1.76e-4 the six depths give one N and one flow-ratio curve, whereas the widths '
      'fan out', len(n_cd) == 1 and same_sig and fan > 20.0,
      'N %s; w_s(1000) %.1f um at 150 um deep and %.1f um at 5 um' % (sorted(n_cd), CD[0].ws[999], CD[-1].ws[999]))
# DIGITIZED holds the last rows of the six curves of panel A, read from the published figure (Gifford 2014, Fig. 2A,
# digitized at 3544 x 2778 px from the embedded image of p. 4499 by extrapolating each curve from w_s = 97 and 99 um to
# 100 um).
# The caption and the panel's label print the largest value as 5.76e-3, which would end at row 95; the figure's first
# curve ends near row 950, the i_T of 5.76e-4, so the printed exponent is a misprint. The design code's own list,
# CIF_v2.gifford_fig2's (5.76, 4.83, 3.90, 2.97, 2.04, 0.96) x 1e-4, ends 5% early on the fifth curve.
DIGITIZED = [951, 1143, 1426, 1898, 2840, 5627]
dev = [abs(i_T(f) / n - 1) for f, n in zip(F_AB, DIGITIZED)]
cif_list = [5.76e-4, 4.83e-4, 3.90e-4, 2.97e-4, 2.04e-4, 9.6e-5]
dev_cif = [abs(i_T(f) / n - 1) for f, n in zip(cif_list, DIGITIZED)]
check('M19 against the figure: the six k x 9.6e-5 end within 2% of the digitized curves, CIF_v2\'s list does not',
      max(dev) < 0.02 and max(dev_cif) > 0.04,
      'largest deviation %.1f%% (CIF_v2\'s list %.1f%%)' % (100 * max(dev), 100 * max(dev_cif)))
print('ALL OK' if not fails else 'FAILED: %s' % fails)
sys.exit(1 if fails else 0)
